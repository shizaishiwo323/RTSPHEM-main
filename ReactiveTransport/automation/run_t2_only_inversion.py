"""Standalone T2 inversion runner for existing NMR/COMSOL decay workbooks.

This script does not run RTM or COMSOL. It only scans an existing result folder
for NMR decay Excel files, then runs T2 inversion into a separate run-ID folder.

Typical input folder:
    C:\\Users\\imgw\\Documents\\Codex\\RTSPHEM-main\\outputs\\rtm_tests\\visual_test_matlab

Expected input files:
    <target_folder>\\comsol_results\\T2*.xlsx

Outputs:
    <target_folder>\\inversion_results\\<run_id>\\
        *_T2.png
        *_T2.mat
        *_T2_lcurve.png      # only when REGULARIZATION_MODE = "lcurve"
        inversion_summary.csv
        run_config.json

You can either edit the USER CONFIG section below and run this file directly,
or pass command-line arguments to override the defaults.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

import pandas as pd


# =============================================================================
# USER CONFIG - edit these values when running manually
# =============================================================================

TARGET_FOLDER = r"C:\Users\imgw\Documents\Codex\RTSPHEM-main\outputs\rtm_runs\rtm_20260502_220315_980_hex"

# T2 axis settings, in milliseconds.
T2_MIN_MS = 1.0e-2
T2_MAX_MS = 1.0e4
NUM_T2_BINS = 200

# Regularization mode:
#   "fixed"  : use FIXED_REGULARIZATION
#   "lcurve" : search the best factor in [LCURVE_ALPHA_MIN, LCURVE_ALPHA_MAX]
REGULARIZATION_MODE = "fixed"
FIXED_REGULARIZATION = 0.5

# L-curve search settings. Used only when REGULARIZATION_MODE = "lcurve".
LCURVE_ALPHA_MIN = 1e-6
LCURVE_ALPHA_MAX = 1e4
LCURVE_ALPHA_COUNT = 120

# COMSOL workbook time unit conversion. If the workbook time column is seconds,
# keep 1000. If it is already milliseconds, set this to 1.
TIME_TO_MS_SCALE = 1000.0

# Calibration:
#   True  : first workbook is recalculated using the first porosity value from
#           global_evolution.xlsx, and subsequent workbooks reuse that factor.
#   False : use MANUAL_CALIBRATION_FACTOR; if None, use T2_process default.
USE_FIRST_POROSITY_CALIBRATION = True
MANUAL_CALIBRATION_FACTOR: float | None = None

# Only files matching this pattern are treated as NMR decay workbooks.
EXCEL_NAME_PATTERN = "T2*.xlsx"

# Leave empty to auto-generate a unique ID from date/time and parameters.
RUN_ID = ""


# =============================================================================
# Implementation
# =============================================================================


@dataclass
class T2OnlyConfig:
    target_folder: Path
    t2_min_ms: float = T2_MIN_MS
    t2_max_ms: float = T2_MAX_MS
    num_t2_bins: int = NUM_T2_BINS
    regularization_mode: str = REGULARIZATION_MODE
    fixed_regularization: float = FIXED_REGULARIZATION
    lcurve_alpha_min: float = LCURVE_ALPHA_MIN
    lcurve_alpha_max: float = LCURVE_ALPHA_MAX
    lcurve_alpha_count: int = LCURVE_ALPHA_COUNT
    time_to_ms_scale: float = TIME_TO_MS_SCALE
    use_first_porosity_calibration: bool = USE_FIRST_POROSITY_CALIBRATION
    manual_calibration_factor: float | None = MANUAL_CALIBRATION_FACTOR
    excel_name_pattern: str = EXCEL_NAME_PATTERN
    run_id: str = RUN_ID


def add_t2_process_to_path() -> Path:
    automation_dir = Path(__file__).resolve().parent
    reactive_transport_dir = automation_dir.parent
    package_root = reactive_transport_dir / "T2_process"
    sys.path.insert(0, str(package_root))
    return package_root


add_t2_process_to_path()

from run_t2_process_inversion import _build_matlab_compatible_outputs  # noqa: E402


def sanitize_id(text: str) -> str:
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", text)
    text = re.sub(r"_+", "_", text).strip("_")
    return text or "run"


def make_run_id(config: T2OnlyConfig) -> str:
    if config.run_id:
        return sanitize_id(config.run_id)

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    if config.regularization_mode.lower() == "fixed":
        mode_text = f"fixed_reg_{config.fixed_regularization:g}"
    else:
        mode_text = (
            f"lcurve_{config.lcurve_alpha_min:g}_{config.lcurve_alpha_max:g}_"
            f"n{config.lcurve_alpha_count}"
        )
    return sanitize_id(
        f"t2only_{stamp}_{mode_text}_T2_{config.t2_min_ms:g}_{config.t2_max_ms:g}_bins_{config.num_t2_bins}"
    )


def find_nmr_decay_excels(target_folder: Path, pattern: str) -> list[Path]:
    comsol_dir = target_folder / "comsol_results"
    if comsol_dir.is_dir():
        hits = sorted(comsol_dir.glob(pattern))
    else:
        hits = []

    if not hits:
        hits = sorted(target_folder.rglob(pattern))

    filtered: list[Path] = []
    for path in hits:
        name = path.name.lower()
        if any(skip in name for skip in ("global_evolution", "tortuosity", "summary")):
            continue
        if "_t2." in name or name.endswith("_t2.xlsx"):
            continue
        filtered.append(path)
    return filtered


def read_first_porosity(target_folder: Path) -> float | None:
    global_file = target_folder / "global_evolution.xlsx"
    if not global_file.is_file():
        return None

    try:
        data = pd.read_excel(global_file)
    except Exception as exc:
        print(f"[warning] Could not read {global_file}: {exc}")
        return None

    porosity_columns = [col for col in data.columns if "porosity" in str(col).lower()]
    if not porosity_columns or data.empty:
        return None

    value = data.loc[data.index[0], porosity_columns[0]]
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None

    return value if math.isfinite(value) else None


def run_one_inversion(
    *,
    excel_file: Path,
    output_dir: Path,
    config: T2OnlyConfig,
    calibration_factor: float | None,
) -> dict[str, Any]:
    return _build_matlab_compatible_outputs(
        input_excel=excel_file,
        output_dir=output_dir,
        time_to_ms_scale=float(config.time_to_ms_scale),
        regularization=float(config.fixed_regularization),
        regularization_mode=config.regularization_mode.lower(),
        alpha_min=float(config.lcurve_alpha_min),
        alpha_max=float(config.lcurve_alpha_max),
        alpha_count=int(config.lcurve_alpha_count),
        calibration_factor=calibration_factor,
        num_bins=int(config.num_t2_bins),
        t2_min_ms=float(config.t2_min_ms),
        t2_max_ms=float(config.t2_max_ms),
    )


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)


def write_summary_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fieldnames = [
        "input_excel",
        "success",
        "total_water",
        "raw_spectrum_sum",
        "calibration_factor",
        "regularization",
        "regularization_mode",
        "output_png",
        "output_mat",
        "output_lcurve_png",
        "error",
    ]
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run_batch(config: T2OnlyConfig) -> Path:
    target_folder = config.target_folder.resolve()
    if not target_folder.is_dir():
        raise FileNotFoundError(f"target_folder does not exist: {target_folder}")

    mode = config.regularization_mode.lower().strip()
    if mode not in {"fixed", "lcurve"}:
        raise ValueError("regularization_mode must be 'fixed' or 'lcurve'.")
    config.regularization_mode = mode

    if not (0 < config.t2_min_ms < config.t2_max_ms):
        raise ValueError("Require 0 < t2_min_ms < t2_max_ms.")
    if config.num_t2_bins < 3:
        raise ValueError("num_t2_bins must be at least 3.")
    if config.lcurve_alpha_count < 3:
        raise ValueError("lcurve_alpha_count must be at least 3.")

    excel_files = find_nmr_decay_excels(target_folder, config.excel_name_pattern)
    if not excel_files:
        raise FileNotFoundError(
            "No NMR decay Excel files were found.\n"
            f"  target_folder: {target_folder}\n"
            f"  preferred location: {target_folder / 'comsol_results'}\n"
            f"  pattern: {config.excel_name_pattern}\n"
            "Please run COMSOL NMR simulation first, or point target_folder to a folder that already contains T2*.xlsx."
        )

    run_id = make_run_id(config)
    output_dir = target_folder / "inversion_results" / run_id
    output_dir.mkdir(parents=True, exist_ok=True)

    config_payload = asdict(config)
    config_payload["target_folder"] = str(target_folder)
    config_payload["run_id"] = run_id
    config_payload["created_at"] = datetime.now().isoformat(timespec="seconds")
    config_payload["input_excels"] = [str(path) for path in excel_files]
    write_json(output_dir / "run_config.json", config_payload)

    print(f"T2-only inversion run ID: {run_id}")
    print(f"Input folder : {target_folder}")
    print(f"Output folder: {output_dir}")
    print(f"Found {len(excel_files)} NMR decay workbook(s).")
    print(f"Mode: {config.regularization_mode}")
    if config.regularization_mode == "fixed":
        print(f"Fixed regularization: {config.fixed_regularization:g}")
    else:
        print(
            "L-curve alpha search: "
            f"[{config.lcurve_alpha_min:g}, {config.lcurve_alpha_max:g}], "
            f"count={config.lcurve_alpha_count}"
        )
    print("")

    first_porosity = read_first_porosity(target_folder)
    calibration_factor = config.manual_calibration_factor
    summary_rows: list[dict[str, Any]] = []

    for index, excel_file in enumerate(excel_files, start=1):
        print(f"========== [{index}/{len(excel_files)}] {excel_file.name} ==========")
        try:
            result = run_one_inversion(
                excel_file=excel_file,
                output_dir=output_dir,
                config=config,
                calibration_factor=calibration_factor,
            )

            if (
                index == 1
                and config.use_first_porosity_calibration
                and first_porosity is not None
                and result.get("raw_spectrum_sum")
                and float(result["raw_spectrum_sum"]) > 0
            ):
                calibration_factor = float(first_porosity) / float(result["raw_spectrum_sum"])
                print(
                    "[calibration] first RTM porosity "
                    f"{first_porosity:.6g} -> calibration_factor={calibration_factor:.12g}"
                )
                result = run_one_inversion(
                    excel_file=excel_file,
                    output_dir=output_dir,
                    config=config,
                    calibration_factor=calibration_factor,
                )

            summary_rows.append(
                {
                    "input_excel": str(excel_file),
                    "success": bool(result.get("success", False)),
                    "total_water": result.get("total_water"),
                    "raw_spectrum_sum": result.get("raw_spectrum_sum"),
                    "calibration_factor": result.get("calibration_factor"),
                    "regularization": result.get("regularization"),
                    "regularization_mode": result.get("regularization_mode", config.regularization_mode),
                    "output_png": result.get("output_png"),
                    "output_mat": result.get("output_mat"),
                    "output_lcurve_png": result.get("output_lcurve_png"),
                    "error": "",
                }
            )
        except Exception as exc:
            print(f"[error] {excel_file}: {exc}")
            summary_rows.append(
                {
                    "input_excel": str(excel_file),
                    "success": False,
                    "total_water": "",
                    "raw_spectrum_sum": "",
                    "calibration_factor": calibration_factor,
                    "regularization": "",
                    "regularization_mode": config.regularization_mode,
                    "output_png": "",
                    "output_mat": "",
                    "output_lcurve_png": "",
                    "error": str(exc),
                }
            )

    write_summary_csv(output_dir / "inversion_summary.csv", summary_rows)
    print("")
    print("Done.")
    print(f"Summary: {output_dir / 'inversion_summary.csv'}")
    return output_dir


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run standalone T2 inversion for existing NMR decay Excel files.")
    parser.add_argument("--target-folder", type=Path, default=Path(TARGET_FOLDER))
    parser.add_argument("--t2-min-ms", type=float, default=T2_MIN_MS)
    parser.add_argument("--t2-max-ms", type=float, default=T2_MAX_MS)
    parser.add_argument("--num-t2-bins", type=int, default=NUM_T2_BINS)
    parser.add_argument("--regularization-mode", choices=["fixed", "lcurve"], default=REGULARIZATION_MODE)
    parser.add_argument("--fixed-regularization", type=float, default=FIXED_REGULARIZATION)
    parser.add_argument("--lcurve-alpha-min", type=float, default=LCURVE_ALPHA_MIN)
    parser.add_argument("--lcurve-alpha-max", type=float, default=LCURVE_ALPHA_MAX)
    parser.add_argument("--lcurve-alpha-count", type=int, default=LCURVE_ALPHA_COUNT)
    parser.add_argument("--time-to-ms-scale", type=float, default=TIME_TO_MS_SCALE)
    parser.add_argument("--excel-name-pattern", default=EXCEL_NAME_PATTERN)
    parser.add_argument("--run-id", default=RUN_ID)
    parser.add_argument("--manual-calibration-factor", type=float, default=MANUAL_CALIBRATION_FACTOR)
    parser.add_argument(
        "--no-first-porosity-calibration",
        action="store_true",
        help="Disable calibration by the first porosity value in global_evolution.xlsx.",
    )
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    config = T2OnlyConfig(
        target_folder=args.target_folder,
        t2_min_ms=args.t2_min_ms,
        t2_max_ms=args.t2_max_ms,
        num_t2_bins=args.num_t2_bins,
        regularization_mode=args.regularization_mode,
        fixed_regularization=args.fixed_regularization,
        lcurve_alpha_min=args.lcurve_alpha_min,
        lcurve_alpha_max=args.lcurve_alpha_max,
        lcurve_alpha_count=args.lcurve_alpha_count,
        time_to_ms_scale=args.time_to_ms_scale,
        use_first_porosity_calibration=not args.no_first_porosity_calibration,
        manual_calibration_factor=args.manual_calibration_factor,
        excel_name_pattern=args.excel_name_pattern,
        run_id=args.run_id,
    )
    run_batch(config)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
