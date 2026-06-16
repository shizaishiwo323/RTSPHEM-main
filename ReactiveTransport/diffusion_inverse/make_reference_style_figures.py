"""Create validation figures inspired by the pore-editing reference paper.

Figures include real/generated image panels, metric distributions, and physical
trajectory summaries from RTM labels. They can be used with retrieval, BBDM, or
any model that writes generated images as `<pair_id>.png`.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from PIL import Image  # noqa: E402


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def make_figures(
    pairs_csv: Path,
    generated_dir: Path,
    metrics_csv: Path,
    output_dir: Path,
    *,
    max_examples: int = 12,
) -> dict:
    output_dir.mkdir(parents=True, exist_ok=True)
    pairs = read_csv(pairs_csv)
    metrics = read_csv(metrics_csv)
    metric_by_pair = {row["pair_id"]: row for row in metrics}
    examples = [row for row in pairs if (generated_dir / f"{row['pair_id']}.png").is_file()][:max_examples]
    if examples:
        plot_examples(examples, generated_dir, output_dir / "real_vs_generated_examples.png")
    if metrics:
        plot_metric_summary(metrics, output_dir / "metric_summary.png")
    plot_physical_summary(pairs, output_dir / "physical_trajectory_summary.png")
    summary = {
        "pairs": len(pairs),
        "metrics": len(metrics),
        "examples": len(examples),
        "mean_pixel_mae": mean_float(row.get("pixel_mae") for row in metrics),
        "mean_porosity_abs_error": mean_float(row.get("porosity_abs_error") for row in metrics),
        "figures": [
            str(output_dir / "real_vs_generated_examples.png"),
            str(output_dir / "metric_summary.png"),
            str(output_dir / "physical_trajectory_summary.png"),
        ],
    }
    (output_dir / "figure_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def plot_examples(examples: list[dict], generated_dir: Path, output_png: Path) -> None:
    cols = len(examples)
    fig, axes = plt.subplots(3, cols, figsize=(max(6, cols * 1.6), 4.8), squeeze=False)
    for col, row in enumerate(examples):
        source = Image.open(row["source_image"]).convert("RGB")
        true = Image.open(row["target_image"]).convert("RGB")
        generated = Image.open(generated_dir / f"{row['pair_id']}.png").convert("RGB")
        for ax, image, title in [
            (axes[0][col], source, f"src {float(row['source_progress']):.2f}"),
            (axes[1][col], true, f"true {float(row['target_progress']):.2f}"),
            (axes[2][col], generated, row.get("direction", "")),
        ]:
            ax.imshow(image)
            ax.set_title(title, fontsize=8)
            ax.axis("off")
    fig.tight_layout()
    fig.savefig(output_png, dpi=220)
    plt.close(fig)


def plot_metric_summary(metrics: list[dict], output_png: Path) -> None:
    progress = np.asarray([float(row["target_progress"]) for row in metrics], dtype=float)
    pixel_mae = np.asarray([float(row["pixel_mae"]) for row in metrics], dtype=float)
    porosity_error = np.asarray([float(row["porosity_abs_error"]) for row in metrics], dtype=float)
    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.2))
    axes[0].scatter(progress, pixel_mae, s=10, alpha=0.55, c="#3b6ea8")
    axes[0].set_xlabel("target progress")
    axes[0].set_ylabel("pixel MAE")
    axes[0].set_title("Image error")
    axes[1].scatter(progress, porosity_error, s=10, alpha=0.55, c="#b54a3a")
    axes[1].set_xlabel("target progress")
    axes[1].set_ylabel("porosity abs. error")
    axes[1].set_title("Porosity consistency")
    fig.tight_layout()
    fig.savefig(output_png, dpi=220)
    plt.close(fig)


def plot_physical_summary(pairs: list[dict], output_png: Path) -> None:
    progress = np.asarray([float(row["target_progress"]) for row in pairs], dtype=float)
    porosity = np.asarray([float(row.get("target_porosity", "nan")) for row in pairs], dtype=float)
    kk0 = np.asarray([float(row.get("target_k_k0", "nan")) for row in pairs], dtype=float)
    tort = np.asarray([float(row.get("target_tortuosity", "nan")) for row in pairs], dtype=float)
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.1))
    for ax, values, ylabel, color in [
        (axes[0], porosity, "RTM porosity", "#367c48"),
        (axes[1], kk0, "k/k0", "#7652a3"),
        (axes[2], tort, "tortuosity", "#a76b27"),
    ]:
        mask = np.isfinite(values)
        ax.scatter(progress[mask], values[mask], s=8, alpha=0.4, c=color)
        ax.set_xlabel("target progress")
        ax.set_ylabel(ylabel)
    fig.tight_layout()
    fig.savefig(output_png, dpi=220)
    plt.close(fig)


def mean_float(values) -> float:
    vals = []
    for value in values:
        try:
            vals.append(float(value))
        except (TypeError, ValueError):
            continue
    return float(np.mean(vals)) if vals else float("nan")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pairs-csv", type=Path, required=True)
    parser.add_argument("--generated-dir", type=Path, required=True)
    parser.add_argument("--metrics-csv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--max-examples", type=int, default=12)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    summary = make_figures(
        args.pairs_csv,
        args.generated_dir,
        args.metrics_csv,
        args.output_dir,
        max_examples=args.max_examples,
    )
    print(f"Wrote figures: {args.output_dir}")
    print(f"Examples: {summary['examples']}  Mean pixel MAE: {summary['mean_pixel_mae']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
