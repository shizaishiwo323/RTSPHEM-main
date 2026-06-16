"""Train RTSPHEM pore evolution with OpenAI improved-diffusion.

This V2 runner keeps openai/improved-diffusion as an external dependency and
uses its UNetModel plus GaussianDiffusion training/sampling code. The local
wrapper only adapts RTSPHEM paired interface images into a conditional
image-to-image problem:

    target_t ~ q(target_0, t)
    eps_theta([target_t, source, condition_maps], t)
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import random
import subprocess
import sys
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


def add_improved_diffusion_to_path(root: Path) -> str:
    if not (root / "improved_diffusion" / "unet.py").is_file():
        raise FileNotFoundError(f"Not an improved-diffusion checkout: {root}")
    root_str = str(root.resolve())
    if root_str not in sys.path:
        sys.path.insert(0, root_str)
    try:
        return subprocess.check_output(["git", "-C", root_str, "rev-parse", "HEAD"], text=True).strip()
    except Exception:
        return "unknown"


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


class RTSPHEMConditionedUNet(nn.Module):
    def __init__(
        self,
        *,
        cond_dim: int,
        image_size: int,
        num_channels: int,
        num_res_blocks: int,
        attention_resolutions: str,
        dropout: float,
        num_heads: int,
        use_scale_shift_norm: bool,
        use_checkpoint: bool,
    ):
        super().__init__()
        from improved_diffusion.unet import UNetModel

        if image_size == 64:
            channel_mult = (1, 2, 3, 4)
        elif image_size == 32:
            channel_mult = (1, 2, 2, 2)
        else:
            raise ValueError(f"unsupported image size for V2 runner: {image_size}")
        attention_ds = tuple(image_size // int(res) for res in attention_resolutions.split(",") if res)
        self.cond_dim = cond_dim
        self.model = UNetModel(
            in_channels=6 + cond_dim,
            model_channels=num_channels,
            out_channels=3,
            num_res_blocks=num_res_blocks,
            attention_resolutions=attention_ds,
            dropout=dropout,
            channel_mult=channel_mult,
            num_classes=None,
            use_checkpoint=use_checkpoint,
            num_heads=num_heads,
            num_heads_upsample=-1,
            use_scale_shift_norm=use_scale_shift_norm,
        )

    def forward(self, x: torch.Tensor, timesteps: torch.Tensor, source: torch.Tensor, cond: torch.Tensor):
        b, _, h, w = x.shape
        cond_map = cond[:, :, None, None].expand(b, cond.shape[1], h, w)
        return self.model(torch.cat([x, source, cond_map], dim=1), timesteps)


def make_diffusion(args: argparse.Namespace, *, sample: bool = False):
    from improved_diffusion.script_util import create_gaussian_diffusion

    return create_gaussian_diffusion(
        steps=args.diffusion_steps,
        learn_sigma=False,
        sigma_small=False,
        noise_schedule=args.noise_schedule,
        use_kl=False,
        predict_xstart=False,
        rescale_timesteps=args.rescale_timesteps,
        rescale_learned_sigmas=False,
        timestep_respacing=(args.timestep_respacing if sample else ""),
    )


def make_model(args: argparse.Namespace, cond_dim: int) -> nn.Module:
    return RTSPHEMConditionedUNet(
        cond_dim=cond_dim,
        image_size=args.image_size,
        num_channels=args.num_channels,
        num_res_blocks=args.num_res_blocks,
        attention_resolutions=args.attention_resolutions,
        dropout=args.dropout,
        num_heads=args.num_heads,
        use_scale_shift_norm=args.use_scale_shift_norm,
        use_checkpoint=args.use_checkpoint,
    )


def train_model(args: argparse.Namespace, commit: str) -> dict:
    train_rows = read_rows(args.train_csv, args.max_train_pairs, args.seed)
    mean, std = compute_condition_stats(train_rows)
    dataset = PairDataset(train_rows, mean, std, args.image_size)
    loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=True, num_workers=0, pin_memory=True)
    device = torch.device(args.device if args.device else ("cuda" if torch.cuda.is_available() else "cpu"))
    model = make_model(args, cond_dim=len(mean)).to(device)
    diffusion = make_diffusion(args)
    opt = torch.optim.AdamW(model.parameters(), lr=args.lr)
    losses = []
    history_rows = []
    step = 0
    model.train()
    for _epoch in range(args.epochs):
        for source, target, cond, _pair_ids in loader:
            source = source.to(device, non_blocking=True)
            target = target.to(device, non_blocking=True)
            cond = cond.to(device, non_blocking=True)
            t = torch.randint(0, diffusion.num_timesteps, (target.shape[0],), device=device)
            terms = diffusion.training_losses(model, target, t, model_kwargs={"source": source, "cond": cond})
            loss = terms["loss"].mean()
            opt.zero_grad(set_to_none=True)
            loss.backward()
            opt.step()
            loss_value = float(loss.item())
            losses.append(loss_value)
            step += 1
            if args.loss_history_csv and (step == 1 or step % args.log_interval == 0):
                window = losses[-args.log_interval :]
                history_rows.append(
                    {
                        "step": step,
                        "loss": loss_value,
                        "loss_window_mean": float(np.mean(window)),
                        "lr": args.lr,
                    }
                )
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
            "improved_diffusion_root": str(args.improved_diffusion_root),
            "improved_diffusion_commit": commit,
            "model_args": {
                "num_channels": args.num_channels,
                "num_res_blocks": args.num_res_blocks,
                "attention_resolutions": args.attention_resolutions,
                "dropout": args.dropout,
                "num_heads": args.num_heads,
                "use_scale_shift_norm": args.use_scale_shift_norm,
                "use_checkpoint": args.use_checkpoint,
                "diffusion_steps": args.diffusion_steps,
                "noise_schedule": args.noise_schedule,
                "rescale_timesteps": args.rescale_timesteps,
            },
        },
        args.checkpoint,
    )
    summary = {
        "checkpoint": str(args.checkpoint),
        "train_rows": len(train_rows),
        "steps": step,
        "mean_train_loss": float(np.mean(losses)) if losses else float("nan"),
        "final_train_loss": losses[-1] if losses else float("nan"),
        "device": str(device),
        "improved_diffusion_root": str(args.improved_diffusion_root),
        "improved_diffusion_commit": commit,
    }
    args.checkpoint.with_suffix(".train_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    if args.loss_history_csv:
        write_csv(args.loss_history_csv, history_rows, ["step", "loss", "loss_window_mean", "lr"])
    return summary


@torch.no_grad()
def predict(args: argparse.Namespace, commit: str) -> dict:
    device = torch.device(args.device if args.device else ("cuda" if torch.cuda.is_available() else "cpu"))
    checkpoint = torch.load(args.checkpoint, map_location=device, weights_only=False)
    model_args = checkpoint.get("model_args", {})
    for key, value in model_args.items():
        if hasattr(args, key):
            setattr(args, key, value)
    mean = np.asarray(checkpoint["cond_mean"], dtype=np.float32)
    std = np.asarray(checkpoint["cond_std"], dtype=np.float32)
    image_size = int(checkpoint["image_size"])
    rows = read_rows(args.test_csv, args.max_test_pairs, args.seed)
    dataset = PairDataset(rows, mean, std, image_size)
    loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=False, num_workers=0)
    model = make_model(args, cond_dim=len(mean)).to(device)
    model.load_state_dict(checkpoint["model_state"])
    model.eval()
    diffusion = make_diffusion(args, sample=True)
    sample_fn = diffusion.ddim_sample_loop if args.use_ddim else diffusion.p_sample_loop
    args.generated_dir.mkdir(parents=True, exist_ok=True)
    pred_rows = []
    for source, _target, cond, pair_ids in loader:
        source = source.to(device)
        cond = cond.to(device)
        samples = sample_fn(
            model,
            (source.shape[0], 3, image_size, image_size),
            clip_denoised=True,
            model_kwargs={"source": source, "cond": cond},
            progress=False,
        )
        for image, pair_id in zip(samples, pair_ids):
            output = args.generated_dir / f"{pair_id}.png"
            save_image_tensor(image, output)
            pred_rows.append({"pair_id": pair_id, "generated_image": str(output)})
    args.predictions_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.predictions_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["pair_id", "generated_image"])
        writer.writeheader()
        writer.writerows(pred_rows)
    return {
        "predictions": len(pred_rows),
        "generated_dir": str(args.generated_dir),
        "improved_diffusion_commit": commit,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--improved-diffusion-root", type=Path, required=True)
    parser.add_argument("--train-csv", type=Path, required=True)
    parser.add_argument("--test-csv", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--generated-dir", type=Path, required=True)
    parser.add_argument("--predictions-csv", type=Path, required=True)
    parser.add_argument("--image-size", type=int, default=64)
    parser.add_argument("--num-channels", type=int, default=32)
    parser.add_argument("--num-res-blocks", type=int, default=2)
    parser.add_argument("--attention-resolutions", default="16,8")
    parser.add_argument("--dropout", type=float, default=0.0)
    parser.add_argument("--num-heads", type=int, default=1)
    parser.add_argument("--use-scale-shift-norm", action="store_true")
    parser.add_argument("--use-checkpoint", action="store_true")
    parser.add_argument("--diffusion-steps", type=int, default=100)
    parser.add_argument("--noise-schedule", default="linear")
    parser.add_argument("--rescale-timesteps", action="store_true")
    parser.add_argument("--timestep-respacing", default="ddim20")
    parser.add_argument("--use-ddim", action="store_true")
    parser.add_argument("--epochs", type=int, default=1)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--lr", type=float, default=2e-4)
    parser.add_argument("--loss-history-csv", type=Path, default=None)
    parser.add_argument("--log-interval", type=int, default=100)
    parser.add_argument("--max-train-pairs", type=int, default=0)
    parser.add_argument("--max-test-pairs", type=int, default=0)
    parser.add_argument("--max-steps", type=int, default=0)
    parser.add_argument("--seed", type=int, default=20260607)
    parser.add_argument("--device", default="")
    parser.add_argument("--predict-only", action="store_true")
    return parser.parse_args()


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    commit = add_improved_diffusion_to_path(args.improved_diffusion_root)
    random.seed(args.seed)
    np.random.seed(args.seed)
    torch.manual_seed(args.seed)
    if not args.predict_only:
        summary = train_model(args, commit)
        print(f"Trained improved-diffusion V2 steps: {summary['steps']} final loss: {summary['final_train_loss']}")
    pred = predict(args, commit)
    print(f"Generated predictions: {pred['predictions']}")
    print(f"improved-diffusion commit: {commit}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
