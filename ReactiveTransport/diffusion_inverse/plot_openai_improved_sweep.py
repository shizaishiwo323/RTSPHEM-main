"""Plot OpenAI improved-diffusion hyperparameter sweep results."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def collect_runs(sweep_root: Path) -> list[dict]:
    rows = []
    for run_dir in sorted(path for path in sweep_root.iterdir() if path.is_dir()):
        config_path = run_dir / "run_config.json"
        summary_path = run_dir / "checkpoint.train_summary.json"
        metrics_path = run_dir / "metrics.summary.json"
        if not (config_path.is_file() and summary_path.is_file() and metrics_path.is_file()):
            continue
        config = read_json(config_path)
        summary = read_json(summary_path)
        metrics = read_json(metrics_path)
        rows.append(
            {
                "run_id": run_dir.name,
                "max_steps": int(config["max_steps"]),
                "lr": float(config["lr"]),
                "noise_schedule": config["noise_schedule"],
                "num_channels": int(config["num_channels"]),
                "diffusion_steps": int(config["diffusion_steps"]),
                "sample_respacing": config["timestep_respacing"],
                "final_train_loss": float(summary["final_train_loss"]),
                "mean_train_loss": float(summary["mean_train_loss"]),
                "mean_pixel_mae": float(metrics["mean_pixel_mae"]),
                "mean_porosity_abs_error": float(metrics["mean_porosity_abs_error"]),
            }
        )
    return rows


def plot_loss_histories(sweep_root: Path, output: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    for run_dir in sorted(path for path in sweep_root.iterdir() if path.is_dir()):
        history_path = run_dir / "loss_history.csv"
        config_path = run_dir / "run_config.json"
        if not (history_path.is_file() and config_path.is_file()):
            continue
        config = read_json(config_path)
        history = read_csv(history_path)
        if not history:
            continue
        steps = [int(row["step"]) for row in history]
        loss = [float(row["loss_window_mean"]) for row in history]
        label = run_dir.name
        if config.get("sweep_axis") == "steps":
            label = f"steps={config['max_steps']}"
        elif config.get("sweep_axis") == "lr":
            label = f"lr={config['lr']}"
        elif config.get("sweep_axis") == "schedule":
            label = f"{config['noise_schedule']}"
        elif config.get("sweep_axis") == "channels":
            label = f"ch={config['num_channels']}"
        ax.plot(steps, loss, linewidth=1.4, label=label)
    ax.set_xlabel("training step")
    ax.set_ylabel("window mean training loss")
    ax.set_title("Improved-diffusion training loss curves")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_step_tradeoff(rows: list[dict], output: Path) -> None:
    step_rows = [
        row
        for row in rows
        if row["run_id"].startswith("steps_")
        or row["run_id"] in {"channels_64"}
        or row["run_id"].startswith("channels64_steps")
    ]
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.8))
    for channels, color in [(32, "#2b5c8a"), (64, "#b24a3b")]:
        group = sorted([row for row in step_rows if row["num_channels"] == channels], key=lambda r: r["max_steps"])
        if not group:
            continue
        x = [row["max_steps"] for row in group]
        train = [row["final_train_loss"] for row in group]
        mae = [row["mean_pixel_mae"] for row in group]
        axes[0].plot(x, train, marker="o", color=color, label=f"{channels} ch")
        axes[1].plot(x, mae, marker="s", color=color, label=f"{channels} ch")
    axes[0].set_xlabel("training steps")
    axes[0].set_ylabel("final train loss")
    axes[0].set_title("Loss vs steps")
    axes[1].set_xlabel("training steps")
    axes[1].set_ylabel("validation RGB MAE")
    axes[1].set_title("Validation error vs steps")
    for ax in axes:
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_hyperparams(rows: list[dict], output: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(11.0, 3.6))
    groups = [
        ("lr", [row for row in rows if row["run_id"].startswith("lr_") or row["run_id"] == "steps_1000"]),
        ("noise_schedule", [row for row in rows if row["run_id"].startswith("schedule_") or row["run_id"] == "steps_1000"]),
        ("num_channels", [row for row in rows if row["run_id"].startswith("channels_") or row["run_id"] == "steps_1000"]),
    ]
    for ax, (field, group_rows) in zip(axes, groups):
        labels = [str(row[field]) for row in group_rows]
        values = [row["mean_pixel_mae"] for row in group_rows]
        ax.bar(labels, values, color="#4f7d6a")
        ax.set_title(field)
        ax.set_ylabel("validation RGB MAE")
        ax.tick_params(axis="x", rotation=25)
        ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220)
    plt.close(fig)


def write_report(rows: list[dict], output: Path) -> dict:
    best = min(rows, key=lambda row: row["mean_pixel_mae"])
    lines = [
        "# OpenAI Improved-Diffusion Hyperparameter Sweep",
        "",
        f"Best validation RGB MAE: `{best['mean_pixel_mae']:.6f}`",
        f"Best run: `{best['run_id']}`",
        "",
        "| run_id | steps | lr | schedule | channels | train_loss | val_rgb_mae | val_porosity_err |",
        "|---|---:|---:|---|---:|---:|---:|---:|",
    ]
    for row in sorted(rows, key=lambda item: item["mean_pixel_mae"]):
        lines.append(
            f"| {row['run_id']} | {row['max_steps']} | {row['lr']:.1e} | {row['noise_schedule']} | "
            f"{row['num_channels']} | {row['final_train_loss']:.6f} | {row['mean_pixel_mae']:.6f} | "
            f"{row['mean_porosity_abs_error']:.6f} |"
        )
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return best


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sweep-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    rows = collect_runs(args.sweep_root)
    if not rows:
        raise SystemExit(f"No complete runs found in {args.sweep_root}")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.output_dir / "sweep_results.csv", rows, list(rows[0].keys()))
    plot_loss_histories(args.sweep_root, args.output_dir / "loss_curves.png")
    plot_step_tradeoff(rows, args.output_dir / "step_loss_validation_tradeoff.png")
    plot_hyperparams(rows, args.output_dir / "hyperparameter_comparison.png")
    best = write_report(rows, args.output_dir / "sweep_report.md")
    (args.output_dir / "best_run.json").write_text(json.dumps(best, indent=2), encoding="utf-8")
    print(f"Best run: {best['run_id']}  val RGB MAE: {best['mean_pixel_mae']:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
