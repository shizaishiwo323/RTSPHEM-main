"""Build selected forward/backward test-case comparisons.

Inputs:
    A prepared diffusion_inverse dataset root with manifests/states.csv.

Outputs:
    selected_pairs.csv with three geometry cases x three Pe/Da settings x
    forward/backward trajectories, plus optional real-vs-generated comparison
    figures and metrics when generated images are available.

This script does not run RTM, COMSOL, or model inference. Use the emitted
selected_pairs.csv as --test-csv for a trained generator, then rerun this script
with --generated-dir to assemble figures and metrics.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterable

import numpy as np
from PIL import Image, ImageDraw, ImageFilter

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from common import find_nearest_state, StateRecord, write_csv  # noqa: E402
from evaluate_interface_trajectory import compare_images  # noqa: E402


PARAMETER_CASES = [
    ("pe1_da0p01", 1.0, 0.01),
    ("pe0p1_da1", 0.1, 1.0),
    ("pe10_da1", 10.0, 1.0),
]

GEOMETRY_CASE_NAMES = ["hex", "random_realization_1", "random_realization_2"]
FORWARD_STAGES = [0.25, 0.50, 0.75, 1.00]
BACKWARD_STAGES = [0.75, 0.50, 0.25, 0.00]
PORE_RED = (255, 0, 0)
SOLID_YELLOW = (255, 255, 0)

PAIR_FIELDNAMES = [
    "pair_id",
    "split",
    "direction",
    "case_id",
    "geometry_case",
    "geometry_index",
    "parameter_case",
    "requested_target_progress",
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


def load_states(states_csv: Path) -> dict[str, list[StateRecord]]:
    by_exp: dict[str, list[StateRecord]] = defaultdict(list)
    with states_csv.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            state = StateRecord(
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
            by_exp[state.exp_id].append(state)
    for exp_id in by_exp:
        by_exp[exp_id] = sorted(by_exp[exp_id], key=lambda state: state.progress)
    return dict(sorted(by_exp.items(), key=lambda item: exp_number(item[0])))


def exp_number(exp_id: str) -> int:
    try:
        return int(exp_id.split("_")[-1])
    except ValueError:
        return 10**9


def build_selected_pairs(dataset_root: Path) -> tuple[list[dict], list[dict], list[str]]:
    by_exp = load_states(dataset_root / "manifests" / "states.csv")
    pairs: list[dict] = []
    case_rows: list[dict] = []
    notes: list[str] = []

    for parameter_label, pe, da in PARAMETER_CASES:
        matching = [
            (exp_id, states)
            for exp_id, states in by_exp.items()
            if states and math.isclose(states[0].Pe, pe) and math.isclose(states[0].Da, da)
        ]
        matching = sorted(matching, key=lambda item: exp_number(item[0]))
        if len(matching) < 3:
            notes.append(
                f"{parameter_label} has {len(matching)} matching experiments; expected three geometry cases."
            )
            continue

        for geometry_index, (exp_id, states) in enumerate(matching[:3], start=1):
            geometry_case = GEOMETRY_CASE_NAMES[geometry_index - 1]
            layout = states[0].layout
            case_id = f"{geometry_case}_{parameter_label}"
            case_rows.append(
                {
                    "case_id": case_id,
                    "geometry_case": geometry_case,
                    "geometry_index": geometry_index,
                    "layout": layout,
                    "parameter_case": parameter_label,
                    "Pe": pe,
                    "Da": da,
                    "exp_id": exp_id,
                    "valid_states": len(states),
                    "progress_min": states[0].progress,
                    "progress_max": states[-1].progress,
                    "porosity_min": states[0].porosity,
                    "porosity_max": states[-1].porosity,
                }
            )
            pairs.extend(
                build_direction_pairs(
                    case_id=case_id,
                    geometry_case=geometry_case,
                    geometry_index=geometry_index,
                    parameter_label=parameter_label,
                    source=find_nearest_state(states, 0.0),
                    stages=FORWARD_STAGES,
                    states=states,
                    direction="forward",
                )
            )
            pairs.extend(
                build_direction_pairs(
                    case_id=case_id,
                    geometry_case=geometry_case,
                    geometry_index=geometry_index,
                    parameter_label=parameter_label,
                    source=find_nearest_state(states, 1.0),
                    stages=BACKWARD_STAGES,
                    states=states,
                    direction="backward",
                )
            )
    return pairs, case_rows, notes


def build_direction_pairs(
    *,
    case_id: str,
    geometry_case: str,
    geometry_index: int,
    parameter_label: str,
    source: StateRecord,
    stages: list[float],
    states: list[StateRecord],
    direction: str,
) -> list[dict]:
    rows = []
    for stage_index, requested_progress in enumerate(stages, start=1):
        target = find_nearest_state(states, requested_progress)
        pair_id = f"{case_id}_{direction}_s{stage_index:02d}_p{progress_token(requested_progress)}"
        rows.append(
            {
                "pair_id": pair_id,
                "split": "selected_test_case",
                "direction": direction,
                "case_id": case_id,
                "geometry_case": geometry_case,
                "geometry_index": geometry_index,
                "parameter_case": parameter_label,
                "requested_target_progress": requested_progress,
                "exp_id": source.exp_id,
                "source_timestep": source.timestep,
                "target_timestep": target.timestep,
                "source_image": source.cropped_image_path,
                "target_image": target.cropped_image_path,
                "bbdm_source_image": "",
                "bbdm_target_image": "",
                "source_progress": source.progress,
                "target_progress": target.progress,
                "signed_progress_delta": target.progress - source.progress,
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
    return rows


def progress_token(value: float) -> str:
    return f"{value:.2f}".replace(".", "p")


def write_case_table(output_dir: Path, case_rows: list[dict]) -> Path:
    path = output_dir / "selected_cases.csv"
    write_csv(
        path,
        case_rows,
        [
            "case_id",
            "geometry_case",
            "geometry_index",
            "layout",
            "parameter_case",
            "Pe",
            "Da",
            "exp_id",
            "valid_states",
            "progress_min",
            "progress_max",
            "porosity_min",
            "porosity_max",
        ],
    )
    return path


def write_manifest_summary(output_dir: Path, pairs: list[dict], case_rows: list[dict], notes: list[str]) -> None:
    summary = {
        "cases": len(case_rows),
        "pairs": len(pairs),
        "directions": sorted({row["direction"] for row in pairs}),
        "stages_per_case_direction": {
            "forward": len(FORWARD_STAGES),
            "backward": len(BACKWARD_STAGES),
        },
        "parameter_cases": PARAMETER_CASES,
        "geometry_case_names": GEOMETRY_CASE_NAMES,
        "notes": notes,
    }
    (output_dir / "selection_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


def load_pairs(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_metrics(pairs: list[dict], generated_dir: Path, output_csv: Path) -> list[dict]:
    rows = []
    for row in pairs:
        generated = generated_dir / f"{row['pair_id']}.png"
        if not generated.is_file():
            continue
        result = compare_images(Path(row["target_image"]), generated)
        rows.append(
            {
                "pair_id": row["pair_id"],
                "case_id": row["case_id"],
                "geometry_case": row["geometry_case"],
                "parameter_case": row["parameter_case"],
                "direction": row["direction"],
                "requested_target_progress": row["requested_target_progress"],
                "target_progress": row["target_progress"],
                "source_timestep": row["source_timestep"],
                "target_timestep": row["target_timestep"],
                "pixel_mae": result["pixel_mae"],
                "pore_mask_mae": result["pore_mask_mae"],
                "porosity_abs_error": result["porosity_abs_error"],
                "component_count_abs_error": result["component_count_abs_error"],
                "interface_length_abs_error": result["interface_length_abs_error"],
                "true_porosity_image": result["true"]["porosity"],
                "generated_porosity_image": result["generated"]["porosity"],
                "target_rtm_porosity": row["target_porosity"],
                "target_permeability_mD": row["target_permeability_mD"],
                "target_k_k0": row["target_k_k0"],
                "target_tortuosity": row["target_tortuosity"],
            }
        )
    fieldnames = [
        "pair_id",
        "case_id",
        "geometry_case",
        "parameter_case",
        "direction",
        "requested_target_progress",
        "target_progress",
        "source_timestep",
        "target_timestep",
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
    write_csv(output_csv, rows, fieldnames)
    summary = {
        "evaluated_pairs": len(rows),
        "mean_pixel_mae": mean(row["pixel_mae"] for row in rows),
        "mean_porosity_abs_error": mean(row["porosity_abs_error"] for row in rows),
        "by_direction": summarize_by(rows, "direction"),
        "by_geometry_case": summarize_by(rows, "geometry_case"),
        "by_parameter_case": summarize_by(rows, "parameter_case"),
    }
    output_csv.with_suffix(".summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return rows


def mean(values: Iterable[float | str]) -> float:
    vals = [float(value) for value in values]
    return float(sum(vals) / len(vals)) if vals else math.nan


def summarize_by(rows: list[dict], key: str) -> dict:
    grouped: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        grouped[str(row[key])].append(row)
    return {
        name: {
            "pairs": len(items),
            "mean_pixel_mae": mean(item["pixel_mae"] for item in items),
            "mean_porosity_abs_error": mean(item["porosity_abs_error"] for item in items),
        }
        for name, items in sorted(grouped.items())
    }


def create_figures(pairs: list[dict], generated_dir: Path, figure_dir: Path) -> list[dict]:
    figure_rows = []
    grouped: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for row in pairs:
        grouped[(row["case_id"], row["direction"])].append(row)
    for (case_id, direction), rows in sorted(grouped.items()):
        rows = sorted(rows, key=lambda row: int(row["pair_id"].split("_s")[-1].split("_")[0]))
        missing = [row["pair_id"] for row in rows if not (generated_dir / f"{row['pair_id']}.png").is_file()]
        if missing:
            continue
        output_png = figure_dir / f"{case_id}_{direction}_real_vs_generated.png"
        draw_case_direction_grid(rows, generated_dir, output_png, binarized=False)
        figure_rows.append(
            {
                "case_id": case_id,
                "direction": direction,
                "figure": str(output_png),
                "stages": len(rows),
                "postprocess": "raw",
            }
        )
        binary_output_png = figure_dir.parent / "figures_binarized" / f"{case_id}_{direction}_real_vs_generated_binarized.png"
        draw_case_direction_grid(rows, generated_dir, binary_output_png, binarized=True)
        figure_rows.append(
            {
                "case_id": case_id,
                "direction": direction,
                "figure": str(binary_output_png),
                "stages": len(rows),
                "postprocess": "binarized",
            }
        )
    return figure_rows


def draw_case_direction_grid(rows: list[dict], generated_dir: Path, output_png: Path, *, binarized: bool) -> None:
    thumb = 150
    col_w = 205
    row_h = 220
    label_w = 100
    header_h = 110
    cols = len(rows) + 1
    width = label_w + cols * col_w
    height = header_h + 2 * row_h + 34
    canvas = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(canvas)
    first = rows[0]
    title = (
        f"{first['case_id']} | {first['direction']} | "
        f"Pe={float(first['Pe']):g}, Da={float(first['Da']):g}, layout={first['layout']}"
    )
    draw.text((14, 12), title, fill=(0, 0, 0))
    subtitle = "Top row: real RTM target. Bottom row: generated target."
    if binarized:
        subtitle += " Binarized: solid=yellow, pore=red."
    draw.text((14, 40), subtitle, fill=(40, 40, 40))

    source = open_thumb(Path(first["source_image"]), thumb, binarized=binarized)
    sx = label_w + 20
    canvas.paste(source, (sx, header_h + 62))
    draw_label_lines(draw, (sx, header_h + 2), source_label_lines(first), fill=(0, 0, 0))
    draw.text((18, header_h + 96), "real", fill=(0, 0, 0))
    draw.text((18, header_h + row_h + 86), "generated", fill=(0, 0, 0))

    for idx, row in enumerate(rows, start=1):
        x = label_w + idx * col_w + 20
        real = open_thumb(Path(row["target_image"]), thumb, binarized=binarized)
        generated_path = generated_dir / f"{row['pair_id']}.png"
        generated_full = Image.open(generated_path).convert("RGB")
        generated = open_thumb(generated_path, thumb, binarized=binarized)
        canvas.paste(real, (x, header_h + 62))
        canvas.paste(generated, (x, header_h + row_h + 46))
        draw_label_lines(draw, (x, header_h + 2), target_label_lines(row, phase="real"), fill=(0, 0, 0))
        draw_label_lines(
            draw,
            (x, header_h + row_h + 10),
            target_label_lines(row, phase="generated", generated_image=generated_full),
            fill=(0, 0, 0),
        )
    output_png.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(output_png)


def open_thumb(path: Path, size: int, *, binarized: bool = False) -> Image.Image:
    image = Image.open(path).convert("RGB")
    if binarized:
        image = binarize_phase_image(image, median_size=3)
    return image.resize((size, size), Image.Resampling.NEAREST if binarized else Image.Resampling.BILINEAR)


def binarize_phase_image(image: Image.Image, *, green_threshold: int = 128, median_size: int = 0) -> Image.Image:
    """Map a red/yellow RTM-style image to pure pore-red and solid-yellow phases."""
    rgb = image.convert("RGB")
    if median_size > 1:
        if median_size % 2 == 0:
            raise ValueError("median_size must be odd")
        rgb = rgb.filter(ImageFilter.MedianFilter(size=median_size))
    arr = np.asarray(rgb)
    solid = arr[:, :, 1] >= green_threshold
    out = np.empty_like(arr)
    out[:, :] = PORE_RED
    out[solid] = SOLID_YELLOW
    return Image.fromarray(out, mode="RGB")


def phase_porosity_from_image(image: Image.Image, *, green_threshold: int = 128, median_size: int = 0) -> float:
    """Return the red pore-phase fraction after the same phase cleanup used for display."""
    binary = binarize_phase_image(image, green_threshold=green_threshold, median_size=median_size)
    arr = np.asarray(binary)
    pore = np.all(arr == PORE_RED, axis=2)
    return float(pore.mean())


def source_label_lines(row: dict) -> list[str]:
    lines = ["source", f"p={float(row['source_progress']):.3f}"]
    if row.get("source_porosity", "") != "":
        lines.append(f"phi={float(row['source_porosity']):.3f}")
    if row.get("source_time_s", "") != "":
        lines.append(f"t={format_time_s(row['source_time_s'])}")
    return lines


def target_label_lines(row: dict, *, phase: str, generated_image: Image.Image | None = None) -> list[str]:
    if phase == "real":
        return [
            f"target p={float(row['target_progress']):.3f}",
            f"phi={float(row['target_porosity']):.3f}",
            f"t={format_time_s(row['target_time_s'])}",
        ]
    if phase == "generated":
        if generated_image is None:
            raise ValueError("generated_image is required for generated target labels")
        generated_porosity = phase_porosity_from_image(generated_image, median_size=3)
        return [
            f"phi={generated_porosity:.3f}",
            f"t={format_time_s(row['target_time_s'])}",
        ]
    raise ValueError(f"unknown label phase: {phase}")


def format_time_s(value: str | float) -> str:
    return f"{float(value):.3g}s"


def draw_label_lines(draw: ImageDraw.ImageDraw, xy: tuple[int, int], lines: list[str], *, fill: tuple[int, int, int]) -> None:
    x, y = xy
    for idx, line in enumerate(lines):
        draw.text((x, y + idx * 15), line, fill=fill)


def write_report(
    *,
    output_dir: Path,
    dataset_root: Path,
    pairs: list[dict],
    case_rows: list[dict],
    notes: list[str],
    metrics_rows: list[dict],
    figure_rows: list[dict],
    generated_dir: Path | None,
) -> Path:
    metrics_summary = {}
    metrics_summary_path = output_dir / "selected_metrics.summary.json"
    if metrics_summary_path.is_file():
        metrics_summary = json.loads(metrics_summary_path.read_text(encoding="utf-8"))
    lines = [
        "# Selected Forward/Backward Interface Evolution Test Cases",
        "",
        f"Dataset root: `{dataset_root}`",
        f"Selected cases: {len(case_rows)}",
        f"Selected source-target pairs: {len(pairs)}",
        f"Generated image directory: `{generated_dir}`" if generated_dir else "Generated image directory: not provided",
        "",
        "## Selection Logic",
        "",
        "- Parameter cases: Pe=1 Da=0.01; Pe=0.1 Da=1; Pe=10 Da=1.",
        "- Geometry cases are assigned from the three matching experiments for each parameter case in exp-id order.",
        "- The manifest only has two layout labels (`hex`, `random`), so the second and third cases are treated as separate random realizations.",
        "- Forward source is the nearest progress 0.0 state; target stages are 0.25, 0.50, 0.75, 1.00.",
        "- Backward source is the nearest progress 1.0 state; target stages are 0.75, 0.50, 0.25, 0.00.",
        "",
        "## Case Table",
        "",
        "| case_id | exp_id | geometry_case | layout | Pe | Da | states | porosity range |",
        "|---|---:|---|---|---:|---:|---:|---|",
    ]
    for row in case_rows:
        lines.append(
            "| {case_id} | {exp_id} | {geometry_case} | {layout} | {Pe:g} | {Da:g} | {valid_states} | {porosity_min:.4f}-{porosity_max:.4f} |".format(
                **row
            )
        )
    if notes:
        lines.extend(["", "## Notes", ""])
        lines.extend(f"- {note}" for note in notes)
    if metrics_rows:
        lines.extend(
            [
                "",
                "## Metric Summary",
                "",
                f"- Evaluated generated pairs: {metrics_summary.get('evaluated_pairs', len(metrics_rows))}",
                f"- Mean RGB pixel MAE: {metrics_summary.get('mean_pixel_mae', math.nan):.6f}",
                f"- Mean image-derived porosity absolute error: {metrics_summary.get('mean_porosity_abs_error', math.nan):.6f}",
                "",
                "### By Direction",
                "",
                "| direction | pairs | mean RGB MAE | mean porosity abs error |",
                "|---|---:|---:|---:|",
            ]
        )
        for direction, summary in metrics_summary.get("by_direction", {}).items():
            lines.append(
                f"| {direction} | {summary['pairs']} | {summary['mean_pixel_mae']:.6f} | {summary['mean_porosity_abs_error']:.6f} |"
            )
    if figure_rows:
        lines.extend(
            [
                "",
                "## Comparison Figures",
                "",
                "Raw figures preserve the generator colors. Binarized figures map solid to pure yellow and pore to pure red for easier phase comparison.",
                "",
            ]
        )
        for row in figure_rows:
            fig_path = Path(row["figure"])
            rel = fig_path.relative_to(output_dir).as_posix()
            lines.append(f"![{row['case_id']} {row['direction']} {row['postprocess']}]({rel})")
            lines.append("")
    report = output_dir / "selected_test_case_report.md"
    report.write_text("\n".join(lines), encoding="utf-8")
    return report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--generated-dir", type=Path, default=None)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    pairs, case_rows, notes = build_selected_pairs(args.dataset_root)
    selected_pairs = args.output_dir / "selected_pairs.csv"
    write_csv(selected_pairs, pairs, PAIR_FIELDNAMES)
    write_case_table(args.output_dir, case_rows)
    write_manifest_summary(args.output_dir, pairs, case_rows, notes)

    metrics_rows: list[dict] = []
    figure_rows: list[dict] = []
    if args.generated_dir is not None:
        metrics_rows = write_metrics(pairs, args.generated_dir, args.output_dir / "selected_metrics.csv")
        figure_rows = create_figures(pairs, args.generated_dir, args.output_dir / "figures")
        write_csv(
            args.output_dir / "comparison_figures.csv",
            figure_rows,
            ["case_id", "direction", "figure", "stages", "postprocess"],
        )

    report = write_report(
        output_dir=args.output_dir,
        dataset_root=args.dataset_root,
        pairs=pairs,
        case_rows=case_rows,
        notes=notes,
        metrics_rows=metrics_rows,
        figure_rows=figure_rows,
        generated_dir=args.generated_dir,
    )
    print(f"Wrote selected pairs: {selected_pairs}")
    print(f"Wrote report: {report}")
    if args.generated_dir is not None:
        print(f"Wrote comparison figures: {args.output_dir / 'figures'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
