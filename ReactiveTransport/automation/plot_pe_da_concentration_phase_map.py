"""Build a Pe-Da phase map from RTM concentration snapshots.

Inputs
------
- A batch directory containing ``exp_XXX`` folders.
- Each experiment folder must include ``run_metadata.json``,
  ``global_evolution_log.csv``, and ``individual_plots/concentration``.

Outputs
-------
- A PNG/PDF phase-map figure with one cropped concentration snapshot per
  Pe-Da point.
- A CSV manifest recording the selected timestep for each experiment.

Selection rule
--------------
For each experiment, choose the timestep whose porosity is closest to a target
porosity, by default 0.80.
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image


@dataclass(frozen=True)
class Snapshot:
    exp: str
    da: float
    pe: float
    timestep: int
    porosity_initial: float
    porosity_final: float
    porosity_target: float
    porosity_selected: float
    image_path: Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot Pe-Da concentration snapshot phase map."
    )
    parser.add_argument(
        "--batch-dir",
        type=Path,
        required=True,
        help="Batch output directory containing exp_XXX folders.",
    )
    parser.add_argument(
        "--first-n",
        type=int,
        default=70,
        help="Number of experiments to include after sorting by folder name.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory. Defaults to <batch-dir>/figures/pe_da_concentration_phase_map.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Figure export DPI.",
    )
    parser.add_argument(
        "--thumbnail-pad-frac",
        type=float,
        default=0.08,
        help="Fractional whitespace around each thumbnail cell.",
    )
    parser.add_argument(
        "--target-porosity",
        type=float,
        default=0.80,
        help="Absolute porosity target used to select each snapshot.",
    )
    parser.add_argument(
        "--no-crop",
        action="store_true",
        help="Use each full original concentration PNG without cropping.",
    )
    return parser.parse_args()


def read_porosity_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    return [
        row
        for row in rows
        if row.get("timestep") and row.get("porosity") and row["porosity"].strip()
    ]


def choose_target_porosity_snapshot(exp_dir: Path, target_porosity: float) -> Snapshot:
    metadata_path = exp_dir / "run_metadata.json"
    evolution_path = exp_dir / "global_evolution_log.csv"

    with metadata_path.open("r", encoding="utf-8") as handle:
        metadata = json.load(handle)

    rows = read_porosity_rows(evolution_path)
    if not rows:
        raise ValueError(f"No porosity rows found in {evolution_path}")

    porosities = np.array([float(row["porosity"]) for row in rows], dtype=float)
    timesteps = np.array([int(float(row["timestep"])) for row in rows], dtype=int)
    phi_initial = float(porosities[0])
    phi_final = float(porosities[-1])
    phi_target = float(target_porosity)
    selected_index = int(np.argmin(np.abs(porosities - phi_target)))
    timestep = int(timesteps[selected_index])

    image_path = (
        exp_dir
        / "individual_plots"
        / "concentration"
        / f"concentration_{timestep:04d}.png"
    )
    if not image_path.exists():
        raise FileNotFoundError(f"Missing selected concentration image: {image_path}")

    params = metadata["parameters"]
    return Snapshot(
        exp=exp_dir.name,
        da=float(params["Da"]),
        pe=float(params["Pe"]),
        timestep=timestep,
        porosity_initial=phi_initial,
        porosity_final=phi_final,
        porosity_target=phi_target,
        porosity_selected=float(porosities[selected_index]),
        image_path=image_path,
    )


def find_plot_area(image: Image.Image) -> tuple[int, int, int, int]:
    """Crop to the main colored concentration panel, excluding labels/colorbar."""
    rgb = np.asarray(image.convert("RGB"))
    height, width, _ = rgb.shape

    # The concentration panel is the large saturated block. Text and axes are
    # dark but not saturated, while the colorbar is narrow and rightmost.
    maxc = rgb.max(axis=2).astype(float)
    minc = rgb.min(axis=2).astype(float)
    saturation = maxc - minc
    colored = saturation > 40

    x_counts = colored.sum(axis=0)
    y_counts = colored.sum(axis=1)

    x_threshold = max(10, int(0.03 * height))
    y_threshold = max(10, int(0.08 * width))
    xs = np.where(x_counts > x_threshold)[0]
    ys = np.where(y_counts > y_threshold)[0]
    if xs.size == 0 or ys.size == 0:
        return (0, 0, width, height)

    # Drop the narrow colorbar by keeping the wide colored run that starts first.
    x_runs = contiguous_runs(xs)
    main_run = max(x_runs, key=lambda run: (run[1] - run[0], -run[0]))
    x0, x1 = main_run

    y_runs = contiguous_runs(ys)
    main_y_run = max(y_runs, key=lambda run: run[1] - run[0])
    y0, y1 = main_y_run

    pad_x = max(2, int(0.005 * width))
    pad_y = max(2, int(0.005 * height))
    return (
        max(0, x0 - pad_x),
        max(0, y0 - pad_y),
        min(width, x1 + 1 + pad_x),
        min(height, y1 + 1 + pad_y),
    )


def contiguous_runs(values: np.ndarray) -> list[tuple[int, int]]:
    if values.size == 0:
        return []
    breaks = np.where(np.diff(values) > 1)[0]
    starts = np.r_[values[0], values[breaks + 1]]
    ends = np.r_[values[breaks], values[-1]]
    return [(int(start), int(end)) for start, end in zip(starts, ends)]


def infer_common_plot_area(paths: Iterable[Path]) -> tuple[int, int, int, int]:
    """Infer one shared plot crop so all thumbnails keep the same geometry."""
    boxes: list[tuple[int, int, int, int]] = []
    image_size: tuple[int, int] | None = None

    for path in paths:
        with Image.open(path) as image:
            if image_size is None:
                image_size = image.size
            if image.size != image_size:
                continue
            boxes.append(find_plot_area(image))

    if not boxes:
        raise ValueError("Could not infer a common plot crop from snapshot images.")

    coords = np.array(boxes, dtype=float)
    x0, y0, x1, y1 = np.median(coords, axis=0).round().astype(int)
    if x1 <= x0 or y1 <= y0:
        raise ValueError(f"Invalid common crop box: {(x0, y0, x1, y1)}")
    return (x0, y0, x1, y1)


def full_image_box(path: Path) -> tuple[int, int, int, int]:
    with Image.open(path) as image:
        width, height = image.size
    return (0, 0, width, height)


def load_cropped_image(
    path: Path,
    crop_box: tuple[int, int, int, int],
) -> np.ndarray:
    with Image.open(path) as image:
        cropped = image.crop(crop_box).convert("RGB")
    return np.asarray(cropped)


def format_tick(value: float) -> str:
    if value >= 1:
        return f"{value:g}"
    return f"{value:.2g}"


def format_porosity_label(value: float) -> str:
    return f"p{value:.2f}".replace(".", "p")


def write_manifest(path: Path, snapshots: Iterable[Snapshot]) -> None:
    fields = [
        "exp",
        "Da",
        "Pe",
        "selected_timestep",
        "porosity_initial",
        "porosity_final",
        "porosity_target",
        "porosity_selected",
        "image_path",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for item in snapshots:
            writer.writerow(
                {
                    "exp": item.exp,
                    "Da": f"{item.da:.12g}",
                    "Pe": f"{item.pe:.12g}",
                    "selected_timestep": item.timestep,
                    "porosity_initial": f"{item.porosity_initial:.12g}",
                    "porosity_final": f"{item.porosity_final:.12g}",
                    "porosity_target": f"{item.porosity_target:.12g}",
                    "porosity_selected": f"{item.porosity_selected:.12g}",
                    "image_path": str(item.image_path),
                }
            )


def plot_phase_map(
    snapshots: list[Snapshot],
    output_dir: Path,
    dpi: int,
    thumbnail_pad_frac: float,
    target_porosity: float,
    no_crop: bool,
) -> Path:
    da_values = sorted({item.da for item in snapshots})
    pe_values = sorted({item.pe for item in snapshots})
    da_index = {value: idx for idx, value in enumerate(da_values)}
    pe_index = {value: idx for idx, value in enumerate(pe_values)}

    n_cols = len(da_values)
    n_rows = len(pe_values)
    fig_width = max(12.0, 1.55 * n_cols + 2.6)
    fig_height = max(8.0, 1.55 * n_rows + 2.4)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height), constrained_layout=False)
    ax.set_xlim(-0.5, n_cols - 0.5)
    ax.set_ylim(-0.5, n_rows - 0.5)
    ax.set_aspect("equal")
    ax.set_xlabel("Da", fontsize=16, labelpad=12)
    ax.set_ylabel("Pe", fontsize=16, labelpad=12)
    ax.set_title(
        f"Concentration snapshots near porosity {target_porosity:.2f}",
        fontsize=18,
        pad=20,
        weight="bold",
    )
    ax.set_xticks(range(n_cols), [format_tick(value) for value in da_values], fontsize=11)
    ax.set_yticks(range(n_rows), [format_tick(value) for value in pe_values], fontsize=11)
    ax.grid(True, color="#d8d8d8", linewidth=0.6, zorder=0)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    pad = max(0.0, min(0.35, thumbnail_pad_frac))
    half_w = 0.5 - pad
    half_h = 0.5 - pad
    if no_crop:
        common_crop_box = full_image_box(snapshots[0].image_path)
    else:
        common_crop_box = infer_common_plot_area(item.image_path for item in snapshots)

    for item in snapshots:
        col = da_index[item.da]
        row = pe_index[item.pe]
        image = load_cropped_image(item.image_path, common_crop_box)
        img_h, img_w = image.shape[:2]
        image_aspect = img_w / img_h
        box_w = 2 * half_w
        box_h = 2 * half_h
        box_aspect = box_w / box_h
        if image_aspect >= box_aspect:
            draw_w = box_w
            draw_h = draw_w / image_aspect
        else:
            draw_h = box_h
            draw_w = draw_h * image_aspect
        ax.imshow(
            image,
            extent=(
                col - 0.5 * draw_w,
                col + 0.5 * draw_w,
                row - 0.5 * draw_h,
                row + 0.5 * draw_h,
            ),
            aspect="equal",
            interpolation="lanczos",
            zorder=2,
        )
        ax.text(
            col,
            row - 0.5 * draw_h - 0.045,
            f"{item.exp}  t={item.timestep:04d}",
            ha="center",
            va="top",
            fontsize=5.8,
            color="#333333",
            clip_on=False,
        )

    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.09, top=0.93)
    porosity_label = format_porosity_label(target_porosity)
    crop_label = "full_originals" if no_crop else "phase_map"
    output_png = output_dir / f"pe_da_concentration_{porosity_label}_{crop_label}.png"
    output_pdf = output_dir / f"pe_da_concentration_{porosity_label}_{crop_label}.pdf"
    fig.savefig(output_png, dpi=dpi, facecolor="white")
    fig.savefig(output_pdf, facecolor="white")
    plt.close(fig)
    return output_png


def main() -> None:
    args = parse_args()
    batch_dir = args.batch_dir.resolve()
    output_dir = (
        args.output_dir.resolve()
        if args.output_dir
        else batch_dir / "figures" / "pe_da_concentration_phase_map"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    exp_dirs = sorted(
        path
        for path in batch_dir.iterdir()
        if path.is_dir() and path.name.startswith("exp_")
    )[: args.first_n]
    if not exp_dirs:
        raise ValueError(f"No exp_XXX directories found in {batch_dir}")

    snapshots = [
        choose_target_porosity_snapshot(exp_dir, args.target_porosity)
        for exp_dir in exp_dirs
    ]
    snapshots.sort(key=lambda item: (item.pe, item.da, item.exp))

    manifest_path = output_dir / (
        f"selected_{format_porosity_label(args.target_porosity)}_snapshots.csv"
    )
    write_manifest(manifest_path, snapshots)
    output_png = plot_phase_map(
        snapshots,
        output_dir=output_dir,
        dpi=args.dpi,
        thumbnail_pad_frac=args.thumbnail_pad_frac,
        target_porosity=args.target_porosity,
        no_crop=args.no_crop,
    )
    print(f"Wrote {output_png}")
    print(f"Wrote {manifest_path}")
    print(
        "Included "
        f"{len(snapshots)} snapshots, "
        f"{len({item.da for item in snapshots})} Da values, "
        f"{len({item.pe for item in snapshots})} Pe values."
    )


if __name__ == "__main__":
    main()
