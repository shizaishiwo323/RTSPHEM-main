"""Run NMR T2 inversion through the standardized T2_process package.

Inputs:
    --input-excel: COMSOL exported spin-echo workbook.
    --output-dir: Folder where automation-compatible outputs are written.
    --calibration-factor: Optional porosity calibration factor.

Outputs:
    <input_stem>_T2.mat with legacy MATLAB automation fields.
    <input_stem>_T2.png using the T2_process plotting style.

This wrapper can run either fixed-regularization NNLS or L-curve
regularization selection.
"""

from __future__ import annotations

import argparse
import json
import math
import shutil
import sys
import tempfile
import traceback
from pathlib import Path
from typing import Any

import numpy as np
from scipy.io import savemat
from scipy.signal import find_peaks


DEFAULT_REGULARIZATION = 0.5
DEFAULT_CALIBRATION_FACTOR = 0.00001 * 2.15544


def _add_t2_process_to_path() -> Path:
    reactive_transport_dir = Path(__file__).resolve().parents[1]
    package_root = reactive_transport_dir / "T2_process"
    sys.path.insert(0, str(package_root))
    return package_root


_add_t2_process_to_path()

from nmr_t2 import NnlsConfig, PlotConfig, invert_single_signal_nnls  # noqa: E402
from nmr_t2.config import LCurveConfig  # noqa: E402
from nmr_t2.io_utils import load_decay_table_multi_column, sort_and_filter_signal, trim_signal_from_global_peak  # noqa: E402
from nmr_t2.lcurve import invert_single_signal_lcurve  # noqa: E402
from nmr_t2.plotting import plot_decay_and_spectrum_pair, plot_lcurve_result  # noqa: E402


def _json_number(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _copy_input_if_needed(input_excel: Path) -> Path:
    """Copy a long Windows input path to temp for libraries with path limits."""

    if sys.platform != "win32" or len(str(input_excel)) <= 220:
        return input_excel

    temp_path = Path(tempfile.gettempdir()) / input_excel.name
    shutil.copy2(input_excel, temp_path)
    print(f"        [long path] Copied input workbook to temp: {temp_path}")
    return temp_path


def _build_matlab_compatible_outputs(
    *,
    input_excel: Path,
    output_dir: Path,
    time_to_ms_scale: float,
    regularization: float,
    regularization_mode: str,
    alpha_min: float,
    alpha_max: float,
    alpha_count: int,
    calibration_factor: float | None,
    num_bins: int,
    t2_min_ms: float,
    t2_max_ms: float,
) -> dict[str, Any]:
    workbook_for_reading = _copy_input_if_needed(input_excel)
    output_dir.mkdir(parents=True, exist_ok=True)

    time_ms, signal_matrix, signal_names, valid_excel_columns = load_decay_table_multi_column(
        workbook_for_reading,
        time_to_ms_scale=float(time_to_ms_scale),
        signal_name_prefix="col",
    )

    regularization_mode = str(regularization_mode).lower().strip()
    if regularization_mode not in {"fixed", "lcurve"}:
        raise ValueError("regularization_mode must be either 'fixed' or 'lcurve'.")

    nnls_config = NnlsConfig(
        num_bins=int(num_bins),
        regularization=float(regularization),
        t2_min_ms=float(t2_min_ms),
        t2_max_ms=float(t2_max_ms),
        min_points_after_trim=10,
    )
    lcurve_config = LCurveConfig(
        num_bins=int(num_bins),
        t2_min_ms=float(t2_min_ms),
        t2_max_ms=float(t2_max_ms),
        alpha_min=float(alpha_min),
        alpha_max=float(alpha_max),
        alpha_count=int(alpha_count),
        min_points_after_trim=10,
    )

    spectra: list[np.ndarray] = []
    fit_times: list[np.ndarray] = []
    fit_amplitudes: list[np.ndarray] = []
    residuals: list[np.ndarray] = []
    trimmed_times: list[np.ndarray] = []
    trimmed_amplitudes: list[np.ndarray] = []
    selected_regularizations: list[float] = []
    lcurve_results = []
    t2_bins_ms: np.ndarray | None = None

    if regularization_mode == "fixed":
        print(f"        T2_process fixed NNLS regularization: {regularization:g}")
    else:
        print(
            "        T2_process L-curve search: "
            f"alpha=[{alpha_min:g}, {alpha_max:g}], count={int(alpha_count)}"
        )
    print(f"        Valid signal columns: {len(signal_names)}")

    for idx, signal_name in enumerate(signal_names):
        signal = signal_matrix[:, idx]
        try:
            trimmed = trim_signal_from_global_peak(
                signal_name,
                time_ms,
                signal,
                min_points_after_trim=int(nnls_config.min_points_after_trim),
            )
            inv_time_ms = trimmed.trimmed_time_ms
            inv_amplitude = trimmed.trimmed_amplitude
        except ValueError:
            inv_time_ms, inv_amplitude = sort_and_filter_signal(time_ms, signal, minimum_points=2)

        if regularization_mode == "fixed":
            inversion = invert_single_signal_nnls(
                inv_time_ms,
                inv_amplitude,
                signal_name=signal_name,
                config=nnls_config,
            )
            selected_regularizations.append(float(inversion.regularization))
        else:
            inversion = invert_single_signal_lcurve(
                inv_time_ms,
                inv_amplitude,
                signal_name=signal_name,
                config=lcurve_config,
            )
            selected_regularizations.append(float(inversion.best_regularization))
            lcurve_results.append(inversion)

        if t2_bins_ms is None:
            t2_bins_ms = inversion.t2_bins_ms

        spectra.append(np.asarray(inversion.spectrum, dtype=float))
        fit_times.append(np.asarray(inversion.fit_time_ms, dtype=float))
        fit_amplitudes.append(np.asarray(inversion.fit_amplitude, dtype=float))
        residuals.append(np.asarray(inversion.residual, dtype=float))
        trimmed_times.append(np.asarray(inv_time_ms, dtype=float))
        trimmed_amplitudes.append(np.asarray(inv_amplitude, dtype=float))

    if t2_bins_ms is None or not spectra:
        raise RuntimeError("No valid T2 inversion result was produced.")

    spectra_matrix = np.column_stack(spectra)
    combined_spectrum = np.sum(spectra_matrix, axis=1)
    raw_spectrum_sum = float(np.sum(combined_spectrum))

    if calibration_factor is None or not math.isfinite(float(calibration_factor)):
        calibration_factor = DEFAULT_CALIBRATION_FACTOR

    water_inc = combined_spectrum * float(calibration_factor)
    water_cum = np.cumsum(water_inc)
    total_water = float(water_cum[-1]) if water_cum.size else float("nan")

    if combined_spectrum.size and float(np.max(combined_spectrum)) > 0:
        locs, peak_properties = find_peaks(
            combined_spectrum,
            height=float(np.max(combined_spectrum)) * 0.05,
            prominence=float(np.max(combined_spectrum)) * 0.01,
        )
        peaks = np.asarray(peak_properties.get("peak_heights", combined_spectrum[locs]), dtype=float)
    else:
        locs = np.array([], dtype=int)
        peaks = np.array([], dtype=float)

    output_png = output_dir / f"{input_excel.stem}_T2.png"
    output_mat = output_dir / f"{input_excel.stem}_T2.mat"
    output_lcurve_png = output_dir / f"{input_excel.stem}_T2_lcurve.png"

    raw_combined = np.nansum(signal_matrix, axis=1)
    sorted_time_ms, sorted_raw_combined = sort_and_filter_signal(time_ms, raw_combined, minimum_points=2)

    plot_config = PlotConfig()
    plot_decay_and_spectrum_pair(
        signal_name=f"{input_excel.stem} combined",
        raw_time_ms=sorted_time_ms,
        raw_amplitude=sorted_raw_combined,
        t2_bins_ms=t2_bins_ms,
        spectrum=combined_spectrum,
        output_path=output_png,
        config=plot_config,
    )

    lcurve_metric_payload = {}
    if regularization_mode == "lcurve" and lcurve_results:
        # For multi-column workbooks, save the first L-curve plot and retain all
        # selected regularizations in the MAT payload.
        plot_lcurve_result(lcurve_results[0], output_path=output_lcurve_png, config=plot_config)
        first_lcurve = lcurve_results[0]
        lcurve_metric_payload = {
            "lcurve_alpha_values": first_lcurve.alpha_values.reshape(-1, 1),
            "lcurve_residual_norms": first_lcurve.residual_norms.reshape(-1, 1),
            "lcurve_roughness_norms": first_lcurve.roughness_norms.reshape(-1, 1),
            "lcurve_zeta_values": first_lcurve.zeta_values.reshape(-1, 1),
            "lcurve_eta_values": first_lcurve.eta_values.reshape(-1, 1),
            "lcurve_slope_reciprocal_values": first_lcurve.slope_reciprocal_values.reshape(-1, 1),
            "lcurve_best_index": np.array([[float(first_lcurve.best_index + 1)]], dtype=float),
            "lcurve_used_range_filter": np.array([[float(first_lcurve.used_range_filter)]], dtype=float),
        }

    mat_payload = {
        "T2_bins": (t2_bins_ms / 1000.0).reshape(-1, 1),
        "T2_bins_ms": t2_bins_ms.reshape(-1, 1),
        "T2_log": np.log10(t2_bins_ms / 1000.0).reshape(-1, 1),
        "T2_log_ms": np.log10(t2_bins_ms).reshape(-1, 1),
        "combined_spectrum": combined_spectrum.reshape(-1, 1),
        "spectra": spectra_matrix,
        "peaks": peaks.reshape(-1, 1),
        "locs": (locs.astype(float) + 1.0).reshape(-1, 1),
        "peak_T2_ms": t2_bins_ms[locs].reshape(-1, 1) if locs.size else np.empty((0, 1)),
        "water_cum": water_cum.reshape(-1, 1),
        "total_water": np.array([[total_water]], dtype=float),
        "t": (time_ms / 1000.0).reshape(-1, 1),
        "time_ms": time_ms.reshape(-1, 1),
        "y_matrix": signal_matrix,
        "signal_names": np.asarray(signal_names, dtype=object),
        "valid_excel_columns": np.asarray(valid_excel_columns, dtype=float).reshape(-1, 1),
        "calibration_factor": np.array([[float(calibration_factor)]], dtype=float),
        "raw_spectrum_sum": np.array([[raw_spectrum_sum]], dtype=float),
        "regularization": np.array([[float(np.nanmean(selected_regularizations))]], dtype=float),
        "selected_regularizations": np.asarray(selected_regularizations, dtype=float).reshape(-1, 1),
        "regularization_mode": regularization_mode,
        "alpha_min": np.array([[float(alpha_min)]], dtype=float),
        "alpha_max": np.array([[float(alpha_max)]], dtype=float),
        "alpha_count": np.array([[float(alpha_count)]], dtype=float),
        "t2_min_ms_config": np.array([[float(t2_min_ms)]], dtype=float),
        "t2_max_ms_config": np.array([[float(t2_max_ms)]], dtype=float),
        "num_bins": np.array([[float(num_bins)]], dtype=float),
        "inversion_backend": f"ReactiveTransport/T2_process/nmr_t2 {regularization_mode}",
        "fit_times_ms": np.asarray(fit_times, dtype=object),
        "fit_amplitudes": np.asarray(fit_amplitudes, dtype=object),
        "fit_residuals": np.asarray(residuals, dtype=object),
        "trimmed_times_ms": np.asarray(trimmed_times, dtype=object),
        "trimmed_amplitudes": np.asarray(trimmed_amplitudes, dtype=object),
    }
    mat_payload.update(lcurve_metric_payload)
    savemat(output_mat, mat_payload, do_compression=True)

    print(f"        [OK] Figure saved: {output_png}")
    if regularization_mode == "lcurve" and lcurve_results:
        print(f"        [OK] L-curve figure saved: {output_lcurve_png}")
    print(f"        [OK] Data saved: {output_mat}")
    print(f"        Total water content: {total_water:.6g}")

    return {
        "success": True,
        "total_water": total_water,
        "raw_spectrum_sum": raw_spectrum_sum,
        "calibration_factor": float(calibration_factor),
        "output_png": str(output_png),
        "output_mat": str(output_mat),
        "regularization": float(np.nanmean(selected_regularizations)),
        "regularization_mode": regularization_mode,
        "selected_regularizations": selected_regularizations,
        "output_lcurve_png": str(output_lcurve_png) if regularization_mode == "lcurve" and lcurve_results else None,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Run T2_process NNLS inversion.")
    parser.add_argument("--input-excel", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--calibration-factor", type=float, default=None)
    parser.add_argument("--regularization", type=float, default=DEFAULT_REGULARIZATION)
    parser.add_argument("--regularization-mode", choices=["fixed", "lcurve"], default="fixed")
    parser.add_argument("--alpha-min", type=float, default=1e-6)
    parser.add_argument("--alpha-max", type=float, default=1e4)
    parser.add_argument("--alpha-count", type=int, default=120)
    parser.add_argument("--time-to-ms-scale", type=float, default=1000.0)
    parser.add_argument("--num-bins", type=int, default=200)
    parser.add_argument("--t2-min-ms", type=float, default=1.0)
    parser.add_argument("--t2-max-ms", type=float, default=1e4)
    args = parser.parse_args()

    try:
        result = _build_matlab_compatible_outputs(
            input_excel=args.input_excel,
            output_dir=args.output_dir,
            time_to_ms_scale=args.time_to_ms_scale,
            regularization=args.regularization,
            regularization_mode=args.regularization_mode,
            alpha_min=args.alpha_min,
            alpha_max=args.alpha_max,
            alpha_count=args.alpha_count,
            calibration_factor=args.calibration_factor,
            num_bins=args.num_bins,
            t2_min_ms=args.t2_min_ms,
            t2_max_ms=args.t2_max_ms,
        )
        print("RESULT_JSON=" + json.dumps(result, ensure_ascii=False, allow_nan=False))
        return 0
    except Exception as exc:
        traceback.print_exc()
        result = {
            "success": False,
            "error": str(exc),
            "total_water": _json_number(float("nan")),
            "raw_spectrum_sum": _json_number(float("nan")),
            "calibration_factor": _json_number(args.calibration_factor),
            "regularization": float(args.regularization),
            "regularization_mode": args.regularization_mode,
        }
        print("RESULT_JSON=" + json.dumps(result, ensure_ascii=False, allow_nan=False))
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
