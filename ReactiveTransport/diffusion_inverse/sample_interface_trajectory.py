"""Prepare target-progress queries for pore-evolution trajectory sampling.

When generated images already exist from an external BBDM run, this script can
assemble a trajectory grid. Without a model command it writes a query manifest
that records source state, target progress, nearest RTM reference state, and
estimated physical time from the RTM trajectory.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

from PIL import Image, ImageDraw

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from common import crop_resize_image, find_nearest_state, read_json, StateRecord, write_csv  # noqa: E402

try:  # noqa: SIM105 - optional script-local dependency
    from fit_time_estimator import load_model, predict_time_s  # type: ignore
except ImportError:  # pragma: no cover
    load_model = None
    predict_time_s = None


def load_states(states_csv: Path, exp_id: str) -> list[StateRecord]:
    states: list[StateRecord] = []
    with states_csv.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row["exp_id"] != exp_id:
                continue
            states.append(
                StateRecord(
                    exp_id=row["exp_id"],
                    timestep=int(row["timestep"]),
                    time_s=float(row["time_s"]),
                    log_time_s=float(row["log_time_s"]),
                    porosity=float(row["porosity"]),
                    permeability_mD=float(row["permeability_mD"]),
                    k_k0=float(row["k_k0"]),
                    tortuosity=float(row["tortuosity"]),
                    surface_area_cm2=float(row["surface_area_cm2"]),
                    progress=float(row["progress"]),
                    image_path=row["image_path"],
                    cropped_image_path=row["cropped_image_path"],
                    Da=float(row["Da"]),
                    Pe=float(row["Pe"]),
                    layout=row["layout"],
                )
            )
    return sorted(states, key=lambda state: state.progress)


def build_query_manifest(
    *,
    dataset_root: Path,
    exp_id: str,
    source_timestep: int,
    target_progresses: list[float],
    output_dir: Path,
    image_size: int = 256,
    time_estimator_path: Path | None = None,
) -> Path:
    states = load_states(dataset_root / "manifests" / "states.csv", exp_id)
    source = next((state for state in states if state.timestep == source_timestep), None)
    if source is None:
        raise ValueError(f"Source timestep {source_timestep} for {exp_id} was not found in states.csv")
    time_model = None
    if time_estimator_path is not None:
        if load_model is None or predict_time_s is None:
            raise RuntimeError("fit_time_estimator.py could not be imported")
        time_model = load_model(time_estimator_path)

    source_query_image = output_dir / "source" / f"{exp_id}_t{source_timestep:04d}.png"
    crop_resize_image(Path(source.image_path), source_query_image, image_size=image_size)
    rows = []
    for target_progress in target_progresses:
        reference = find_nearest_state(states, target_progress)
        predicted_time_s = reference.time_s
        time_source = "nearest_rtm_state"
        if time_model is not None:
            predicted_time_s = predict_time_s(
                time_model,
                {
                    "source_progress": source.progress,
                    "target_progress": target_progress,
                    "signed_progress_delta": target_progress - source.progress,
                    "source_time_s": source.time_s,
                    "Da": source.Da,
                    "Pe": source.Pe,
                    "layout": source.layout,
                },
            )
            time_source = "time_estimator"
        rows.append(
            {
                "query_id": f"{exp_id}_t{source_timestep:04d}_p{target_progress:.4f}".replace(".", "p"),
                "exp_id": exp_id,
                "source_timestep": source_timestep,
                "source_image": str(source_query_image),
                "source_progress": source.progress,
                "target_progress": target_progress,
                "signed_progress_delta": target_progress - source.progress,
                "direction": "forward" if target_progress > source.progress else "backward",
                "nearest_reference_timestep": reference.timestep,
                "nearest_reference_image": reference.cropped_image_path,
                "estimated_target_time_s": predicted_time_s,
                "estimated_target_log_time_s": math.log10(max(predicted_time_s, 1e-12)),
                "time_estimate_source": time_source,
                "Da": source.Da,
                "Pe": source.Pe,
                "layout": source.layout,
            }
        )
    manifest_path = output_dir / "query_manifest.csv"
    write_csv(manifest_path, rows, list(rows[0].keys()) if rows else [])
    return manifest_path


def write_trajectory_grid(query_manifest: Path, generated_dir: Path, output_png: Path) -> None:
    rows = list(csv.DictReader(query_manifest.open(newline="", encoding="utf-8")))
    if not rows:
        raise ValueError("Query manifest is empty")
    thumbs = []
    labels = []
    for row in rows:
        generated = generated_dir / f"{row['query_id']}.png"
        image_path = generated if generated.is_file() else Path(row["nearest_reference_image"])
        thumbs.append(Image.open(image_path).convert("RGB").resize((160, 160)))
        labels.append(f"p={float(row['target_progress']):.2f} t~{float(row['estimated_target_time_s']):.3g}s")
    width = 160 * len(thumbs)
    canvas = Image.new("RGB", (width, 190), "white")
    draw = ImageDraw.Draw(canvas)
    for idx, image in enumerate(thumbs):
        x = idx * 160
        canvas.paste(image, (x, 0))
        draw.text((x + 4, 166), labels[idx], fill=(0, 0, 0))
    output_png.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(output_png)


def parse_progresses(text: str) -> list[float]:
    values = [float(item.strip()) for item in text.split(",") if item.strip()]
    for value in values:
        if not 0.0 <= value <= 1.0:
            raise ValueError(f"target progress must be in [0, 1], got {value}")
    return values


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-root", type=Path, required=True)
    parser.add_argument("--exp-id", required=True)
    parser.add_argument("--source-timestep", type=int, required=True)
    parser.add_argument("--target-progresses", required=True, help="Comma-separated values in [0, 1].")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--image-size", type=int, default=256)
    parser.add_argument("--generated-dir", type=Path, default=None)
    parser.add_argument("--time-estimator", type=Path, default=None)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    _ = read_json(args.dataset_root / "manifests" / "summary.json")
    manifest = build_query_manifest(
        dataset_root=args.dataset_root,
        exp_id=args.exp_id,
        source_timestep=args.source_timestep,
        target_progresses=parse_progresses(args.target_progresses),
        output_dir=args.output_dir,
        image_size=args.image_size,
        time_estimator_path=args.time_estimator,
    )
    print(f"Wrote query manifest: {manifest}")
    if args.generated_dir:
        grid = args.output_dir / "trajectory_grid.png"
        write_trajectory_grid(manifest, args.generated_dir, grid)
        print(f"Wrote trajectory grid: {grid}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
