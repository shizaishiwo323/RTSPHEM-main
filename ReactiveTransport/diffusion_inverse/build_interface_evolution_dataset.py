"""Build a bidirectional interface-image dataset for diffusion pore evolution.

The dataset includes all timesteps whose global porosity is at or below the
configured cutoff. Each valid state is cropped to the non-white interface
content and resized. Pairs are then formed within the same RTM experiment so a
model can learn both backward history reconstruction and forward dissolution.
"""

from __future__ import annotations

import argparse
import math
import random
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from common import (  # noqa: E402
    StateRecord,
    copy_or_overwrite,
    crop_resize_image,
    ensure_dir,
    read_evolution_table,
    read_json,
    split_exp_ids,
    state_to_row,
    timestep_from_path,
    write_csv,
    write_json,
)


STATE_FIELDS = [
    "exp_id",
    "timestep",
    "time_s",
    "log_time_s",
    "porosity",
    "permeability_mD",
    "k_k0",
    "tortuosity",
    "surface_area_cm2",
    "progress",
    "image_path",
    "cropped_image_path",
    "Da",
    "Pe",
    "layout",
]

PAIR_FIELDS = [
    "pair_id",
    "split",
    "direction",
    "exp_id",
    "source_timestep",
    "target_timestep",
    "source_image",
    "target_image",
    "bbdm_source_image",
    "bbdm_target_image",
    "source_progress",
    "target_progress",
    "signed_progress_delta",
    "source_time_s",
    "target_time_s",
    "delta_time_s",
    "source_log_time_s",
    "target_log_time_s",
    "source_porosity",
    "target_porosity",
    "source_permeability_mD",
    "target_permeability_mD",
    "source_k_k0",
    "target_k_k0",
    "source_tortuosity",
    "target_tortuosity",
    "Da",
    "Pe",
    "layout",
]


def load_exp_states(
    exp_dir: Path,
    output_root: Path,
    *,
    porosity_max: float,
    image_size: int,
    white_threshold: int,
) -> tuple[list[StateRecord], dict[str, int]]:
    exp_id = exp_dir.name
    metadata = read_json(exp_dir / "run_metadata.json")
    params = metadata.get("parameters", {})
    da = float(params.get("Da", math.nan))
    pe = float(params.get("Pe", math.nan))
    layout = str(metadata.get("layoutType", params.get("layoutType", "")))

    rows = read_evolution_table(exp_dir)
    row_by_step = {int(row["TimeStep"]): row for row in rows}
    image_dir = exp_dir / "interface_images"
    image_paths = sorted(image_dir.glob("timestep_*.png"))
    valid_rows = [
        row
        for row in rows
        if math.isfinite(float(row["Porosity"])) and float(row["Porosity"]) <= porosity_max
    ]
    skipped_high = len(rows) - len(valid_rows)
    if len(valid_rows) < 2:
        return [], {"skipped_high_porosity": skipped_high, "missing_images": 0}

    phi_initial = float(valid_rows[0]["Porosity"])
    phi_cutoff = min(porosity_max, max(float(row["Porosity"]) for row in valid_rows))
    denom = phi_cutoff - phi_initial
    if denom <= 0:
        return [], {"skipped_high_porosity": skipped_high, "missing_images": 0}

    states: list[StateRecord] = []
    missing_images = 0
    image_by_step = {timestep_from_path(path): path for path in image_paths}
    for row in valid_rows:
        timestep = int(row["TimeStep"])
        source_image = image_by_step.get(timestep)
        if source_image is None or not source_image.is_file():
            missing_images += 1
            continue
        cropped_path = output_root / "images" / "states" / f"{exp_id}_t{timestep:04d}.png"
        crop_resize_image(
            source_image,
            cropped_path,
            image_size=image_size,
            white_threshold=white_threshold,
        )
        time_s = float(row["Time_s"])
        progress = (float(row["Porosity"]) - phi_initial) / denom
        states.append(
            StateRecord(
                exp_id=exp_id,
                timestep=timestep,
                time_s=time_s,
                log_time_s=math.log10(max(time_s, 1e-12)),
                porosity=float(row["Porosity"]),
                permeability_mD=float(row["Permeability_mD"]),
                k_k0=float(row["k_k0"]),
                tortuosity=float(row["Tortuosity"]),
                surface_area_cm2=float(row["SurfaceArea_cm2"]),
                progress=max(0.0, min(1.0, progress)),
                image_path=str(source_image),
                cropped_image_path=str(cropped_path),
                Da=da,
                Pe=pe,
                layout=layout,
            )
        )
    return states, {"skipped_high_porosity": skipped_high, "missing_images": missing_images}


def make_pairs(
    states_by_exp: dict[str, list[StateRecord]],
    split_map: dict[str, str],
    output_root: Path,
    *,
    seed: int,
    max_pairs_per_exp: int,
    materialize_bbdm: bool,
) -> list[dict]:
    rng = random.Random(seed)
    pairs: list[dict] = []
    for exp_id, states in sorted(states_by_exp.items()):
        candidates = [(src, dst) for src in states for dst in states if src.timestep != dst.timestep]
        if max_pairs_per_exp > 0 and len(candidates) > max_pairs_per_exp:
            candidates = rng.sample(candidates, max_pairs_per_exp)
        for index, (source, target) in enumerate(candidates, start=1):
            pair_id = f"{exp_id}_p{index:06d}"
            split = split_map[exp_id]
            source_pair_path = output_root / "images" / "A" / split / f"{pair_id}.png"
            target_pair_path = output_root / "images" / "B" / split / f"{pair_id}.png"
            if materialize_bbdm:
                copy_or_overwrite(Path(source.cropped_image_path), source_pair_path)
                copy_or_overwrite(Path(target.cropped_image_path), target_pair_path)
            signed_progress_delta = target.progress - source.progress
            pairs.append(
                {
                    "pair_id": pair_id,
                    "split": split,
                    "direction": "forward" if signed_progress_delta > 0 else "backward",
                    "exp_id": exp_id,
                    "source_timestep": source.timestep,
                    "target_timestep": target.timestep,
                    "source_image": source.cropped_image_path,
                    "target_image": target.cropped_image_path,
                    "bbdm_source_image": str(source_pair_path) if materialize_bbdm else "",
                    "bbdm_target_image": str(target_pair_path) if materialize_bbdm else "",
                    "source_progress": source.progress,
                    "target_progress": target.progress,
                    "signed_progress_delta": signed_progress_delta,
                    "source_time_s": source.time_s,
                    "target_time_s": target.time_s,
                    "delta_time_s": target.time_s - source.time_s,
                    "source_log_time_s": source.log_time_s,
                    "target_log_time_s": target.log_time_s,
                    "source_porosity": source.porosity,
                    "target_porosity": target.porosity,
                    "source_permeability_mD": source.permeability_mD,
                    "target_permeability_mD": target.permeability_mD,
                    "source_k_k0": source.k_k0,
                    "target_k_k0": target.k_k0,
                    "source_tortuosity": source.tortuosity,
                    "target_tortuosity": target.tortuosity,
                    "Da": source.Da,
                    "Pe": source.Pe,
                    "layout": source.layout,
                }
            )
    return pairs


def build_dataset(
    *,
    batch_root: Path,
    output_root: Path,
    porosity_max: float = 0.80,
    image_size: int = 256,
    white_threshold: int = 245,
    seed: int = 20260607,
    max_pairs_per_exp: int = 2000,
    materialize_bbdm: bool = True,
) -> dict:
    ensure_dir(output_root)
    states_by_exp: dict[str, list[StateRecord]] = {}
    skipped_high_porosity = 0
    missing_images = 0
    for exp_dir in sorted(path for path in batch_root.iterdir() if path.is_dir() and path.name.startswith("exp_")):
        try:
            states, counts = load_exp_states(
                exp_dir,
                output_root,
                porosity_max=porosity_max,
                image_size=image_size,
                white_threshold=white_threshold,
            )
        except (FileNotFoundError, RuntimeError, ValueError):
            continue
        skipped_high_porosity += counts["skipped_high_porosity"]
        missing_images += counts["missing_images"]
        if len(states) >= 2:
            states_by_exp[exp_dir.name] = states

    split_map = split_exp_ids(list(states_by_exp), seed)
    states = [state for exp_states in states_by_exp.values() for state in exp_states]
    pairs = make_pairs(
        states_by_exp,
        split_map,
        output_root,
        seed=seed,
        max_pairs_per_exp=max_pairs_per_exp,
        materialize_bbdm=materialize_bbdm,
    )

    write_csv(output_root / "manifests" / "states.csv", [state_to_row(state) for state in states], STATE_FIELDS)
    write_csv(output_root / "manifests" / "pairs.csv", pairs, PAIR_FIELDS)
    for split in ("train", "val", "test"):
        write_csv(
            output_root / "splits" / f"{split}.csv",
            [row for row in pairs if row["split"] == split],
            PAIR_FIELDS,
        )

    summary = {
        "batch_root": str(batch_root),
        "output_root": str(output_root),
        "porosity_max": porosity_max,
        "image_size": image_size,
        "experiments": len(states_by_exp),
        "states": len(states),
        "pairs": len(pairs),
        "skipped_high_porosity": skipped_high_porosity,
        "missing_images": missing_images,
        "max_pairs_per_exp": max_pairs_per_exp,
        "materialize_bbdm": materialize_bbdm,
    }
    write_json(output_root / "manifests" / "summary.json", summary)
    write_qc_report(output_root, summary, pairs)
    return summary


def write_qc_report(output_root: Path, summary: dict, pairs: list[dict]) -> None:
    by_split = {split: sum(1 for row in pairs if row["split"] == split) for split in ("train", "val", "test")}
    forward = sum(1 for row in pairs if row["direction"] == "forward")
    backward = sum(1 for row in pairs if row["direction"] == "backward")
    lines = [
        "# Interface Evolution Dataset QC",
        "",
        f"- experiments: {summary['experiments']}",
        f"- states: {summary['states']}",
        f"- pairs: {summary['pairs']}",
        f"- porosity_max: {summary['porosity_max']}",
        f"- skipped_high_porosity: {summary['skipped_high_porosity']}",
        f"- missing_images: {summary['missing_images']}",
        f"- forward_pairs: {forward}",
        f"- backward_pairs: {backward}",
        f"- split_train_pairs: {by_split['train']}",
        f"- split_val_pairs: {by_split['val']}",
        f"- split_test_pairs: {by_split['test']}",
        "",
    ]
    qc_path = output_root / "qc" / "dataset_qc.md"
    ensure_dir(qc_path.parent)
    qc_path.write_text("\n".join(lines), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--porosity-max", type=float, default=0.80)
    parser.add_argument("--image-size", type=int, default=256)
    parser.add_argument("--white-threshold", type=int, default=245)
    parser.add_argument("--seed", type=int, default=20260607)
    parser.add_argument("--max-pairs-per-exp", type=int, default=2000, help="0 means all directed pairs.")
    parser.add_argument("--no-materialize-bbdm", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    summary = build_dataset(
        batch_root=args.batch_root,
        output_root=args.output_root,
        porosity_max=args.porosity_max,
        image_size=args.image_size,
        white_threshold=args.white_threshold,
        seed=args.seed,
        max_pairs_per_exp=args.max_pairs_per_exp,
        materialize_bbdm=not args.no_materialize_bbdm,
    )
    print(f"Wrote dataset: {summary['output_root']}")
    print(f"States: {summary['states']}  Pairs: {summary['pairs']}  Experiments: {summary['experiments']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
