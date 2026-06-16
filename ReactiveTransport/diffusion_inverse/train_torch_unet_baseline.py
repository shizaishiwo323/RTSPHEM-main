"""Train a lightweight conditional U-Net baseline for interface evolution.

This is a GPU-capable neural baseline for the full RTM manifest. It is not a
DDPM/BBDM implementation; its purpose is to produce genuine model-generated
images through the same test/evaluation/figure pipeline while the external BBDM
code is being prepared.
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
    raise SystemExit("PyTorch is required. Use C:\\Users\\imgw\\.conda\\envs\\ml\\python.exe") from exc


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
    mean = matrix.mean(axis=0)
    std = matrix.std(axis=0)
    std[std == 0] = 1.0
    return mean.astype(np.float32), std.astype(np.float32)


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
        return source, torch.tensor(cond, dtype=torch.float32), target, row["pair_id"]


def load_image(path: Path, image_size: int) -> torch.Tensor:
    image = Image.open(path).convert("RGB")
    if image.size != (image_size, image_size):
        image = image.resize((image_size, image_size), Image.Resampling.BILINEAR)
    arr = np.asarray(image, dtype=np.float32) / 255.0
    return torch.tensor(arr.transpose(2, 0, 1), dtype=torch.float32)


class ConditionalUNet(nn.Module):
    def __init__(self, cond_dim: int):
        super().__init__()
        channels = 3 + cond_dim
        self.enc1 = block(channels, 32)
        self.enc2 = block(32, 64)
        self.enc3 = block(64, 128)
        self.pool = nn.MaxPool2d(2)
        self.up2 = nn.ConvTranspose2d(128, 64, 2, stride=2)
        self.dec2 = block(128, 64)
        self.up1 = nn.ConvTranspose2d(64, 32, 2, stride=2)
        self.dec1 = block(64, 32)
        self.out = nn.Conv2d(32, 3, 1)

    def forward(self, image: torch.Tensor, cond: torch.Tensor) -> torch.Tensor:
        b, _, h, w = image.shape
        cond_map = cond[:, :, None, None].expand(b, cond.shape[1], h, w)
        x = torch.cat([image, cond_map], dim=1)
        e1 = self.enc1(x)
        e2 = self.enc2(self.pool(e1))
        e3 = self.enc3(self.pool(e2))
        d2 = self.up2(e3)
        d2 = self.dec2(torch.cat([d2, e2], dim=1))
        d1 = self.up1(d2)
        d1 = self.dec1(torch.cat([d1, e1], dim=1))
        return torch.sigmoid(self.out(d1))


def block(in_channels: int, out_channels: int) -> nn.Sequential:
    return nn.Sequential(
        nn.Conv2d(in_channels, out_channels, 3, padding=1),
        nn.GroupNorm(4, out_channels),
        nn.SiLU(),
        nn.Conv2d(out_channels, out_channels, 3, padding=1),
        nn.GroupNorm(4, out_channels),
        nn.SiLU(),
    )


def train_model(args: argparse.Namespace) -> dict:
    train_rows = read_rows(args.train_csv, args.max_train_pairs, args.seed)
    mean, std = compute_condition_stats(train_rows)
    dataset = PairDataset(train_rows, mean, std, args.image_size)
    loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=True, num_workers=0, pin_memory=True)
    device = torch.device(args.device if args.device else ("cuda" if torch.cuda.is_available() else "cpu"))
    model = ConditionalUNet(cond_dim=len(mean)).to(device)
    optim = torch.optim.AdamW(model.parameters(), lr=args.lr)
    loss_fn = nn.L1Loss()
    step = 0
    losses = []
    model.train()
    for epoch in range(args.epochs):
        for source, cond, target, _ in loader:
            source = source.to(device, non_blocking=True)
            cond = cond.to(device, non_blocking=True)
            target = target.to(device, non_blocking=True)
            pred = model(source, cond)
            loss = loss_fn(pred, target)
            optim.zero_grad(set_to_none=True)
            loss.backward()
            optim.step()
            step += 1
            losses.append(float(loss.item()))
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
            "cond_fields": COND_FIELDS,
            "train_csv": str(args.train_csv),
        },
        args.checkpoint,
    )
    summary = {
        "checkpoint": str(args.checkpoint),
        "train_rows": len(train_rows),
        "steps": step,
        "epochs_requested": args.epochs,
        "mean_train_l1": float(np.mean(losses)) if losses else float("nan"),
        "final_train_l1": losses[-1] if losses else float("nan"),
        "device": str(device),
    }
    args.checkpoint.with_suffix(".train_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def predict(args: argparse.Namespace) -> dict:
    device = torch.device(args.device if args.device else ("cuda" if torch.cuda.is_available() else "cpu"))
    checkpoint = torch.load(args.checkpoint, map_location=device, weights_only=False)
    mean = np.asarray(checkpoint["cond_mean"], dtype=np.float32)
    std = np.asarray(checkpoint["cond_std"], dtype=np.float32)
    image_size = int(checkpoint["image_size"])
    rows = read_rows(args.test_csv, args.max_test_pairs, args.seed)
    dataset = PairDataset(rows, mean, std, image_size)
    loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=False, num_workers=0)
    model = ConditionalUNet(cond_dim=len(mean)).to(device)
    model.load_state_dict(checkpoint["model_state"])
    model.eval()
    args.generated_dir.mkdir(parents=True, exist_ok=True)
    prediction_rows = []
    with torch.no_grad():
        for source, cond, _, pair_ids in loader:
            pred = model(source.to(device), cond.to(device)).cpu().numpy()
            for image, pair_id in zip(pred, pair_ids):
                output = args.generated_dir / f"{pair_id}.png"
                save_tensor_image(image, output)
                prediction_rows.append({"pair_id": pair_id, "generated_image": str(output)})
    args.predictions_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.predictions_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["pair_id", "generated_image"])
        writer.writeheader()
        writer.writerows(prediction_rows)
    return {"predictions": len(prediction_rows), "generated_dir": str(args.generated_dir)}


def save_tensor_image(chw: np.ndarray, output: Path) -> None:
    arr = np.clip(chw.transpose(1, 2, 0) * 255.0, 0, 255).astype(np.uint8)
    Image.fromarray(arr, mode="RGB").save(output)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--train-csv", type=Path, required=True)
    parser.add_argument("--test-csv", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--generated-dir", type=Path, required=True)
    parser.add_argument("--predictions-csv", type=Path, required=True)
    parser.add_argument("--image-size", type=int, default=128)
    parser.add_argument("--epochs", type=int, default=1)
    parser.add_argument("--batch-size", type=int, default=16)
    parser.add_argument("--lr", type=float, default=1e-3)
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
        train_summary = train_model(args)
        print(f"Trained U-Net steps: {train_summary['steps']} final L1: {train_summary['final_train_l1']}")
    pred_summary = predict(args)
    print(f"Generated predictions: {pred_summary['predictions']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
