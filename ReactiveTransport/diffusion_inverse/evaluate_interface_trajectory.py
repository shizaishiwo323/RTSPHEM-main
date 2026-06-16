"""Evaluate generated interface-image trajectories against RTM references.

This evaluator intentionally uses image-derived morphology metrics plus the
RTM manifest columns. It does not rerun RTM or COMSOL; permeability/k/k0 and
tortuosity comparisons are manifest-level checks when reference records exist.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import deque
from pathlib import Path
from typing import Iterable

import numpy as np
from PIL import Image

try:
    from scipy import ndimage  # type: ignore
except ImportError:  # pragma: no cover - fallback path
    ndimage = None


def pore_mask_from_image(path: Path, *, white_threshold: int = 245) -> np.ndarray:
    arr = np.asarray(Image.open(path).convert("RGB"))
    yellow = (arr[:, :, 0] >= 200) & (arr[:, :, 1] >= 180) & (arr[:, :, 2] <= 80)
    if np.any(yellow):
        return yellow
    return np.any(arr < white_threshold, axis=2)


def rgb_float(path: Path) -> np.ndarray:
    return np.asarray(Image.open(path).convert("RGB"), dtype=np.float32) / 255.0


def image_metrics(path: Path) -> dict[str, float | int]:
    mask = pore_mask_from_image(path)
    porosity = float(mask.mean())
    components, largest = connected_components(mask)
    perimeter = interface_length(mask)
    s2 = two_point_correlation(mask)
    return {
        "porosity": porosity,
        "connected_components": components,
        "connected_pore_fraction": float(largest / max(1, int(mask.sum()))),
        "interface_length_px": float(perimeter),
        "s2_x_lag1": float(s2["x_lag1"]),
        "s2_y_lag1": float(s2["y_lag1"]),
    }


def connected_components(mask: np.ndarray) -> tuple[int, int]:
    if ndimage is not None:
        structure = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]], dtype=bool)
        labels, count = ndimage.label(mask, structure=structure)
        if count == 0:
            return 0, 0
        sizes = np.bincount(labels.ravel())
        largest = int(sizes[1:].max()) if len(sizes) > 1 else 0
        return int(count), largest
    visited = np.zeros(mask.shape, dtype=bool)
    components = 0
    largest = 0
    height, width = mask.shape
    for y in range(height):
        for x in range(width):
            if not mask[y, x] or visited[y, x]:
                continue
            components += 1
            size = 0
            queue: deque[tuple[int, int]] = deque([(y, x)])
            visited[y, x] = True
            while queue:
                cy, cx = queue.popleft()
                size += 1
                for ny, nx in ((cy - 1, cx), (cy + 1, cx), (cy, cx - 1), (cy, cx + 1)):
                    if 0 <= ny < height and 0 <= nx < width and mask[ny, nx] and not visited[ny, nx]:
                        visited[ny, nx] = True
                        queue.append((ny, nx))
            largest = max(largest, size)
    return components, largest


def interface_length(mask: np.ndarray) -> int:
    vertical = np.count_nonzero(mask[1:, :] != mask[:-1, :])
    horizontal = np.count_nonzero(mask[:, 1:] != mask[:, :-1])
    return int(vertical + horizontal)


def two_point_correlation(mask: np.ndarray) -> dict[str, float]:
    data = mask.astype(np.float32)
    if data.shape[1] > 1:
        x_lag1 = float(np.mean(data[:, :-1] * data[:, 1:]))
    else:
        x_lag1 = float("nan")
    if data.shape[0] > 1:
        y_lag1 = float(np.mean(data[:-1, :] * data[1:, :]))
    else:
        y_lag1 = float("nan")
    return {"x_lag1": x_lag1, "y_lag1": y_lag1}


def compare_images(true_path: Path, generated_path: Path) -> dict:
    true_mask = pore_mask_from_image(true_path).astype(np.float32)
    generated_mask = pore_mask_from_image(generated_path).astype(np.float32)
    true_rgb = rgb_float(true_path)
    generated_rgb = rgb_float(generated_path)
    if true_mask.shape != generated_mask.shape:
        generated_mask = np.asarray(
            Image.fromarray((generated_mask * 255).astype(np.uint8)).resize(
                (true_mask.shape[1], true_mask.shape[0]), Image.Resampling.NEAREST
            )
        ).astype(np.float32) / 255.0
    if true_rgb.shape != generated_rgb.shape:
        generated_rgb = (
            np.asarray(
                Image.fromarray((generated_rgb * 255).astype(np.uint8)).resize(
                    (true_rgb.shape[1], true_rgb.shape[0]), Image.Resampling.BILINEAR
                ),
                dtype=np.float32,
            )
            / 255.0
        )
    true_metrics = image_metrics(true_path)
    generated_metrics = image_metrics(generated_path)
    return {
        "true": true_metrics,
        "generated": generated_metrics,
        "pixel_mae": float(np.mean(np.abs(true_rgb - generated_rgb))),
        "pore_mask_mae": float(np.mean(np.abs(true_mask - generated_mask))),
        "porosity_abs_error": abs(float(true_metrics["porosity"]) - float(generated_metrics["porosity"])),
        "component_count_abs_error": abs(
            int(true_metrics["connected_components"]) - int(generated_metrics["connected_components"])
        ),
        "interface_length_abs_error": abs(
            float(true_metrics["interface_length_px"]) - float(generated_metrics["interface_length_px"])
        ),
    }


def count_monotonicity_violations(values: Iterable[float], direction: str) -> int:
    sequence = list(values)
    violations = 0
    for prev, current in zip(sequence, sequence[1:]):
        if direction == "forward" and current < prev:
            violations += 1
        elif direction == "backward" and current > prev:
            violations += 1
    return violations


def evaluate_pairs(pairs_csv: Path, generated_dir: Path, output_csv: Path) -> dict:
    rows: list[dict] = []
    with pairs_csv.open(newline="", encoding="utf-8") as handle:
        for pair in csv.DictReader(handle):
            generated = generated_dir / f"{pair['pair_id']}.png"
            if not generated.is_file():
                continue
            result = compare_images(Path(pair["target_image"]), generated)
            rows.append(
                {
                    "pair_id": pair["pair_id"],
                    "direction": pair["direction"],
                    "target_progress": pair["target_progress"],
                    "pixel_mae": result["pixel_mae"],
                    "pore_mask_mae": result["pore_mask_mae"],
                    "porosity_abs_error": result["porosity_abs_error"],
                    "component_count_abs_error": result["component_count_abs_error"],
                    "interface_length_abs_error": result["interface_length_abs_error"],
                    "true_porosity_image": result["true"]["porosity"],
                    "generated_porosity_image": result["generated"]["porosity"],
                    "target_rtm_porosity": pair.get("target_porosity", ""),
                    "target_permeability_mD": pair.get("target_permeability_mD", ""),
                    "target_k_k0": pair.get("target_k_k0", ""),
                    "target_tortuosity": pair.get("target_tortuosity", ""),
                }
            )
    fieldnames = [
        "pair_id",
        "direction",
        "target_progress",
        "pixel_mae",
        "pore_mask_mae",
        "porosity_abs_error",
        "component_count_abs_error",
        "interface_length_abs_error",
        "true_porosity_image",
        "generated_porosity_image",
        "target_rtm_porosity",
        "target_permeability_mD",
        "target_k_k0",
        "target_tortuosity",
    ]
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    summary = {
        "evaluated_pairs": len(rows),
        "mean_pixel_mae": _mean(row["pixel_mae"] for row in rows),
        "mean_porosity_abs_error": _mean(row["porosity_abs_error"] for row in rows),
    }
    output_csv.with_suffix(".summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def _mean(values: Iterable[float]) -> float:
    vals = [float(value) for value in values]
    if not vals:
        return math.nan
    return float(sum(vals) / len(vals))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pairs-csv", type=Path, required=True)
    parser.add_argument("--generated-dir", type=Path, required=True)
    parser.add_argument("--output-csv", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    summary = evaluate_pairs(args.pairs_csv, args.generated_dir, args.output_csv)
    print(f"Evaluated pairs: {summary['evaluated_pairs']}")
    print(f"Mean pixel MAE: {summary['mean_pixel_mae']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
