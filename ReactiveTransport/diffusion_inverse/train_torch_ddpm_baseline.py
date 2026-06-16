"""Train a conditional DDPM baseline for bidirectional pore evolution.

This model learns p(target image | source image, progress/RTM conditions). It is
intended as a local, reproducible diffusion baseline while the external BBDM
repository remains a handoff dependency.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import random
from pathlib import Path

import numpy as np
from PIL import Image

try:
    import torch
    from torch import nn
    from torch.utils.data import DataLoader, Dataset
except ImportError as exc:  # pragma: no cover
    raise SystemExit("PyTorch is required. Use C:\\Users\\imgw\\Documents\\Codex\\NMR-agent\\.venv\\Scripts\\python.exe") from exc


COND_FIELDS = [
    "source_progress",
    "target_progress",
    "signed_progress_delta",
    "source_time_s",
    "Da",
    "Pe",
    "source_porosity",
    "source_k_k0",
    "source_tortuosity",
]


def read_rows(path: Path, limit: int = 0, seed: int = 20260607) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if limit and len(rows) > limit:
        rng = random.Random(seed)
        rows = rng.sample(rows, limit)
    return rows


def condition_vector(row: dict) -> np.ndarray:
    values = []
    for field in COND_FIELDS:
        value = float(row.get(field, 0.0) or 0.0)
        if field in {"source_time_s", "Da", "Pe", "source_k_k0"}:
            value = math.log10(max(value, 1e-12))
        values.append(value)
    values.append(1.0 if row.get("layout") == "hex" else 0.0)
    values.append(1.0 if row.get("direction") == "forward" else -1.0)
    return np.asarray(values, dtype=np.float32)


def compute_condition_stats(rows: list[dict]) -> tuple[np.ndarray, np.ndarray]:
    matrix = np.stack([condition_vector(row) for row in rows])
    mean = matrix.mean(axis=0).astype(np.float32)
    std = matrix.std(axis=0).astype(np.float32)
    std[std == 0] = 1.0
    return mean, std


class PairDataset(Dataset):
    def __init__(self, rows: list[dict], mean: np.ndarray, std: np.ndarray, image_size: int):
        self.rows = rows
        self.mean = mean
        self.std = std
        self.image_size = image_size

    def __len__(self) -> int:
        return len(self.rows)

    def __getitem__(self, index: int):
        row = self.rows[index]
        source = load_image(Path(row["source_image"]), self.image_size)
        target = load_image(Path(row["target_image"]), self.image_size)
        cond = (condition_vector(row) - self.mean) / self.std
        return source, target, torch.tensor(cond, dtype=torch.float32), row["pair_id"]


def load_image(path: Path, image_size: int) -> torch.Tensor:
    image = Image.open(path).convert("RGB")
    if image.size != (image_size, image_size):
        image = image.resize((image_size, image_size), Image.Resampling.BILINEAR)
    arr = np.asarray(image, dtype=np.float32) / 255.0
    arr = arr * 2.0 - 1.0
    return torch.tensor(arr.transpose(2, 0, 1), dtype=torch.float32)


def save_image_tensor(chw: torch.Tensor, output: Path) -> None:
    arr = chw.detach().cpu().clamp(-1, 1).numpy()
    arr = ((arr + 1.0) * 127.5).clip(0, 255).astype(np.uint8).transpose(1, 2, 0)
    output.parent.mkdir(parents=True, exist_ok=True)
    Image.fromarray(arr, mode="RGB").save(output)


class SinusoidalTimeEmbedding(nn.Module):
    def __init__(self, dim: int):
        super().__init__()
        self.dim = dim

    def forward(self, t: torch.Tensor) -> torch.Tensor:
        half = self.dim // 2
        freqs = torch.exp(
            torch.arange(half, device=t.device, dtype=torch.float32)
            * -(math.log(10000.0) / max(half - 1, 1))
        )
        args = t.float()[:, None] * freqs[None, :]
        emb = torch.cat([torch.sin(args), torch.cos(args)], dim=1)
        if self.dim % 2:
            emb = torch.cat([emb, torch.zeros_like(emb[:, :1])], dim=1)
        return emb


class ConditionalDenoiser(nn.Module):
    def __init__(self, cond_dim: int, time_dim: int = 64):
        super().__init__()
        self.time_embed = nn.Sequential(
            SinusoidalTimeEmbedding(time_dim),
            nn.Linear(time_dim, time_dim),
            nn.SiLU(),
            nn.Linear(time_dim, time_dim),
        )
        in_channels = 6 + cond_dim + time_dim
        self.enc1 = block(in_channels, 48)
        self.enc2 = block(48, 96)
        self.enc3 = block(96, 160)
        self.pool = nn.MaxPool2d(2)
        self.up2 = nn.ConvTranspose2d(160, 96, 2, stride=2)
        self.dec2 = block(192, 96)
        self.up1 = nn.ConvTranspose2d(96, 48, 2, stride=2)
        self.dec1 = block(96, 48)
        self.out = nn.Conv2d(48, 3, 1)

    def forward(self, noisy_target: torch.Tensor, source: torch.Tensor, cond: torch.Tensor, t: torch.Tensor) -> torch.Tensor:
        b, _, h, w = source.shape
        t_emb = self.time_embed(t)
        cond_map = cond[:, :, None, None].expand(b, cond.shape[1], h, w)
        time_map = t_emb[:, :, None, None].expand(b, t_emb.shape[1], h, w)
        x = torch.cat([noisy_target, source, cond_map, time_map], dim=1)
        e1 = self.enc1(x)
        e2 = self.enc2(self.pool(e1))
        e3 = self.enc3(self.pool(e2))
        d2 = self.up2(e3)
        d2 = self.dec2(torch.cat([d2, e2], dim=1))
        d1 = self.up1(d2)
        d1 = self.dec1(torch.cat([d1, e1], dim=1))
        return self.out(d1)


def block(in_channels: int, out_channels: int) -> nn.Sequential:
    return nn.Sequential(
        nn.Conv2d(in_channels, out_channels, 3, padding=1),
        nn.GroupNorm(8, out_channels),
        nn.SiLU(),
        nn.Conv2d(out_channels, out_channels, 3, padding=1),
        nn.GroupNorm(8, out_channels),
        nn.SiLU(),
    )


class DiffusionSchedule:
    def __init__(self, timesteps: int, device: torch.device):
        self.timesteps = timesteps
        betas = torch.linspace(1e-4, 0.02, timesteps, device=device)
        alphas = 1.0 - betas
        alpha_bars = torch.cumprod(alphas, dim=0)
        self.betas = betas
        self.alphas = alphas
        self.alpha_bars = alpha_bars
        self.sqrt_alpha_bars = torch.sqrt(alpha_bars)
        self.sqrt_one_minus_alpha_bars = torch.sqrt(1.0 - alpha_bars)

    def add_noise(self, x0: torch.Tensor, t: torch.Tensor, noise: torch.Tensor) -> torch.Tensor:
        a = self.sqrt_alpha_bars[t][:, None, None, None]
        b = self.sqrt_one_minus_alpha_bars[t][:, None, None, None]
        return a * x0 + b * noise


def train_model(args: argparse.Namespace) -> dict:
    train_rows = read_rows(args.train_csv, args.max_train_pairs, args.seed)
    mean, std = compute_condition_stats(train_rows)
    dataset = PairDataset(train_rows, mean, std, args.image_size)
    loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=True, num_workers=0, pin_memory=True)
    device = torch.device(args.device if args.device else ("cuda" if torch.cuda.is_available() else "cpu"))
    model = ConditionalDenoiser(cond_dim=len(mean)).to(device)
    schedule = DiffusionSchedule(args.timesteps, device)
    optim = torch.optim.AdamW(model.parameters(), lr=args.lr)
    loss_fn = nn.MSELoss()
    losses = []
    step = 0
    model.train()
    for _epoch in range(args.epochs):
        for source, target, cond, _pair_ids in loader:
            source = source.to(device, non_blocking=True)
            target = target.to(device, non_blocking=True)
            cond = cond.to(device, non_blocking=True)
            t = torch.randint(0, args.timesteps, (target.shape[0],), device=device)
            noise = torch.randn_like(target)
            noisy = schedule.add_noise(target, t, noise)
            pred = model(noisy, source, cond, t)
            loss = loss_fn(pred, noise)
            optim.zero_grad(set_to_none=True)
            loss.backward()
            optim.step()
            losses.append(float(loss.item()))
            step += 1
            if args.max_steps and step >= args.max_steps:
                break
        if args.max_steps and step >= args.max_steps:
            break
    args.checkpoint.parent.mkdir(parents=True, exist_ok=True)
    torch.save(
        {
            "model_state": model.state_dict(),
            "cond_mean": mean,
            "cond_std": std,
            "image_size": args.image_size,
            "timesteps": args.timesteps,
            "cond_fields": COND_FIELDS,
            "train_csv": str(args.train_csv),
        },
        args.checkpoint,
    )
    summary = {
        "checkpoint": str(args.checkpoint),
        "train_rows": len(train_rows),
        "steps": step,
        "timesteps": args.timesteps,
        "mean_train_mse": float(np.mean(losses)) if losses else float("nan"),
        "final_train_mse": losses[-1] if losses else float("nan"),
        "device": str(device),
    }
    args.checkpoint.with_suffix(".train_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


@torch.no_grad()
def sample_batch(
    model: ConditionalDenoiser,
    source: torch.Tensor,
    cond: torch.Tensor,
    schedule: DiffusionSchedule,
    sample_steps: int,
    eta: float,
) -> torch.Tensor:
    device = source.device
    x = torch.randn_like(source)
    steps = torch.linspace(schedule.timesteps - 1, 0, sample_steps, device=device).long()
    for idx, t_value in enumerate(steps):
        t = torch.full((source.shape[0],), int(t_value.item()), device=device, dtype=torch.long)
        eps = model(x, source, cond, t)
        alpha_bar = schedule.alpha_bars[t][:, None, None, None]
        pred_x0 = (x - torch.sqrt(1 - alpha_bar) * eps) / torch.sqrt(alpha_bar)
        pred_x0 = pred_x0.clamp(-1, 1)
        if idx == len(steps) - 1:
            x = pred_x0
            continue
        t_next_value = int(steps[idx + 1].item())
        alpha_next = schedule.alpha_bars[torch.full_like(t, t_next_value)][:, None, None, None]
        if eta == 0:
            x = torch.sqrt(alpha_next) * pred_x0 + torch.sqrt(1 - alpha_next) * eps
        else:
            noise = torch.randn_like(x)
            sigma = eta * torch.sqrt((1 - alpha_next) / (1 - alpha_bar)) * torch.sqrt(1 - alpha_bar / alpha_next)
            direction = torch.sqrt(torch.clamp(1 - alpha_next - sigma**2, min=0.0)) * eps
            x = torch.sqrt(alpha_next) * pred_x0 + direction + sigma * noise
    return x.clamp(-1, 1)


def predict(args: argparse.Namespace) -> dict:
    device = torch.device(args.device if args.device else ("cuda" if torch.cuda.is_available() else "cpu"))
    checkpoint = torch.load(args.checkpoint, map_location=device, weights_only=False)
    mean = np.asarray(checkpoint["cond_mean"], dtype=np.float32)
    std = np.asarray(checkpoint["cond_std"], dtype=np.float32)
    image_size = int(checkpoint["image_size"])
    timesteps = int(checkpoint["timesteps"])
    rows = read_rows(args.test_csv, args.max_test_pairs, args.seed)
    dataset = PairDataset(rows, mean, std, image_size)
    loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=False, num_workers=0)
    model = ConditionalDenoiser(cond_dim=len(mean)).to(device)
    model.load_state_dict(checkpoint["model_state"])
    model.eval()
    schedule = DiffusionSchedule(timesteps, device)
    args.generated_dir.mkdir(parents=True, exist_ok=True)
    pred_rows = []
    for source, _target, cond, pair_ids in loader:
        source = source.to(device)
        cond = cond.to(device)
        samples = sample_batch(model, source, cond, schedule, args.sample_steps, args.eta)
        for image, pair_id in zip(samples, pair_ids):
            output = args.generated_dir / f"{pair_id}.png"
            save_image_tensor(image, output)
            pred_rows.append({"pair_id": pair_id, "generated_image": str(output)})
    args.predictions_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.predictions_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["pair_id", "generated_image"])
        writer.writeheader()
        writer.writerows(pred_rows)
    return {"predictions": len(pred_rows), "generated_dir": str(args.generated_dir)}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--train-csv", type=Path, required=True)
    parser.add_argument("--test-csv", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--generated-dir", type=Path, required=True)
    parser.add_argument("--predictions-csv", type=Path, required=True)
    parser.add_argument("--image-size", type=int, default=64)
    parser.add_argument("--epochs", type=int, default=1)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--lr", type=float, default=2e-4)
    parser.add_argument("--timesteps", type=int, default=100)
    parser.add_argument("--sample-steps", type=int, default=25)
    parser.add_argument("--eta", type=float, default=0.0)
    parser.add_argument("--max-train-pairs", type=int, default=0)
    parser.add_argument("--max-test-pairs", type=int, default=0)
    parser.add_argument("--max-steps", type=int, default=0)
    parser.add_argument("--seed", type=int, default=20260607)
    parser.add_argument("--device", default="")
    parser.add_argument("--predict-only", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    random.seed(args.seed)
    np.random.seed(args.seed)
    torch.manual_seed(args.seed)
    if not args.predict_only:
        summary = train_model(args)
        print(f"Trained DDPM steps: {summary['steps']} final MSE: {summary['final_train_mse']}")
    pred = predict(args)
    print(f"Generated predictions: {pred['predictions']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
