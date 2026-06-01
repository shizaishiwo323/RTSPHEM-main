"""PNG-based NMR simulation bridge for RTSPHEM RTM runs.

Inputs:
  - one RTSPHEM experiment directory containing run_metadata.json
  - one interface PNG such as interface_images/timestep_0001.png
  - a JSON configuration written from NMRSimulationConfig.m

Outputs:
  - PNG-derived decay CSV/PNG files under png_nmr_results/<method>/decay
  - fixed-alpha T2 inversion CSV/MAT/PNG under png_nmr_results/<method>/inversion
  - one per-step manifest JSON

This module intentionally imports the mature NMR scripts from
C:/Users/imgw/Documents/Codex/NMR模拟/advanced_tools instead of copying their
solver code into RTSPHEM.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from dataclasses import asdict
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg", force=True)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image
from scipy.io import savemat


METHOD_TO_SOLVER = {
    "png_pixel_cpu": "pixel",
    "png_pixel_gpu": "pixel_gpu",
    "png_mesh": "triangular",
}

DEFAULT_OUTSIDE_BOUNDARY = {
    "left": "gas",
    "right": "gas",
    "top": "solid",
    "bottom": "solid",
}


def method_to_solver(method: str) -> str:
    key = str(method).strip().lower()
    if key not in METHOD_TO_SOLVER:
        raise ValueError(f"Unsupported PNG NMR method: {method}")
    return METHOD_TO_SOLVER[key]


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    if isinstance(value, np.ndarray):
        return json_safe(value.tolist())
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, (np.floating, float)):
        value = float(value)
        return value if math.isfinite(value) else None
    if isinstance(value, (np.integer, int)):
        return int(value)
    return value


def ensure_on_path(path: str | Path) -> None:
    path = str(Path(path))
    if path and path not in sys.path:
        sys.path.insert(0, path)


def metadata_lengths_um(metadata: dict[str, Any]) -> tuple[float, float]:
    params = metadata.get("parameters", {})
    length_x_cm = float(params["lengthXAxis_cm"])
    length_y_cm = float(params["lengthYAxis_cm"])
    return length_x_cm * 10000.0, length_y_cm * 10000.0


def normalize_boundary_kind(value: Any, *, default: str) -> str:
    value = str(scalar_value(value) or default).strip().lower()
    aliases = {
        "solid": "solid",
        "solid_liquid": "solid",
        "solid-liquid": "solid",
        "water_solid": "solid",
        "water-solid": "solid",
        "gas": "gas",
        "gas_liquid": "gas",
        "gas-liquid": "gas",
        "water_gas": "gas",
        "water-gas": "gas",
    }
    if value not in aliases:
        raise ValueError(f"Unsupported outside boundary kind: {value!r}; use 'solid' or 'gas'.")
    return aliases[value]


def parse_outside_boundary_config(config: dict[str, Any]) -> dict[str, str]:
    mode = str(scalar_value(config.get("outside_boundary_mode", "left_right_gas_top_bottom_solid"))).strip().lower()
    if mode in {"default", "left_right_gas_top_bottom_solid", "lr_gas_tb_solid", "rtsphem_default"}:
        boundaries = dict(DEFAULT_OUTSIDE_BOUNDARY)
    elif mode in {"all_solid", "solid", "solid_liquid"}:
        boundaries = {"left": "solid", "right": "solid", "top": "solid", "bottom": "solid"}
    elif mode in {"all_gas", "gas", "gas_liquid"}:
        boundaries = {"left": "gas", "right": "gas", "top": "gas", "bottom": "gas"}
    elif mode in {"custom", "per_side", "manual"}:
        boundaries = dict(DEFAULT_OUTSIDE_BOUNDARY)
    else:
        raise ValueError(
            f"Unsupported outside_boundary_mode: {mode!r}; use default, all_solid, all_gas, or custom."
        )

    for side in ("left", "right", "top", "bottom"):
        key = f"outside_boundary_{side}"
        if key in config and scalar_value(config.get(key)) is not None:
            boundaries[side] = normalize_boundary_kind(config[key], default=boundaries[side])
    return boundaries


def boundary_kind_for_pixel_side(
    row: int,
    col: int,
    neighbor_row: int,
    neighbor_col: int,
    bbox: tuple[int, int, int, int],
    boundaries: dict[str, str],
) -> str:
    r_min, r_max, c_min, c_max = bbox
    if neighbor_col < c_min or col == c_min:
        return boundaries["left"]
    if neighbor_col > c_max or col == c_max:
        return boundaries["right"]
    if neighbor_row < r_min or row == r_min:
        return boundaries["top"]
    if neighbor_row > r_max or row == r_max:
        return boundaries["bottom"]
    return boundaries["left"]


def boundary_kind_for_mesh_side(row: int, col: int, bbox: tuple[int, int, int, int], boundaries: dict[str, str]) -> str:
    r_min, r_max, c_min, c_max = bbox
    if col <= c_min + 2:
        return boundaries["left"]
    if col >= c_max - 2:
        return boundaries["right"]
    if row <= r_min + 2:
        return boundaries["top"]
    if row >= r_max - 2:
        return boundaries["bottom"]
    return boundaries["left"]


def apply_outside_boundary_convention(config: dict[str, Any]) -> dict[str, str]:
    import advanced_tools.png_phase_nmr_decay as png_decay

    boundaries = parse_outside_boundary_config(config)

    def configured_pixel_boundary(row, col, neighbor_row, neighbor_col, bbox):
        return boundary_kind_for_pixel_side(row, col, neighbor_row, neighbor_col, bbox, boundaries)

    def configured_mesh_boundary(labels, midpoint, params):
        row, col = png_decay.xy_to_row_col(labels, float(midpoint[0]), float(midpoint[1]), params)
        r0, r1 = max(0, row - 2), min(labels.shape[0], row + 3)
        c0, c1 = max(0, col - 2), min(labels.shape[1], col + 3)
        local = labels[r0:r1, c0:c1]
        if np.any(local == png_decay.SOLID):
            return "solid"
        return boundary_kind_for_mesh_side(row, col, png_decay.sample_bbox(labels), boundaries)

    png_decay.boundary_kind_for_outside_neighbor = configured_pixel_boundary
    png_decay.classify_boundary_edge = configured_mesh_boundary
    return boundaries


def nonwhite_bbox(rgb: np.ndarray) -> dict[str, int]:
    nonwhite = np.any(rgb[:, :, :3] < 250, axis=2)
    yy, xx = np.nonzero(nonwhite)
    if yy.size == 0 or xx.size == 0:
        raise ValueError("No non-white sample bbox was found in the interface PNG.")
    return {
        "row_min": int(yy.min()),
        "row_max": int(yy.max()),
        "col_min": int(xx.min()),
        "col_max": int(xx.max()),
        "bbox_height_px": int(yy.max() - yy.min() + 1),
        "bbox_width_px": int(xx.max() - xx.min() + 1),
    }


def compute_pixel_size_from_png(
    *,
    exp_dir: Path,
    metadata: dict[str, Any],
    use_nonwhite_bbox: bool,
    png_path: Path | None = None,
) -> tuple[float, float, dict[str, Any]]:
    if png_path is None:
        png_path = next((exp_dir / "interface_images").glob("timestep_*.png"))
    rgb = np.asarray(Image.open(png_path).convert("RGB"))
    height, width = rgb.shape[:2]
    length_x_um, length_y_um = metadata_lengths_um(metadata)

    info: dict[str, Any] = {
        "raw_height_px": int(height),
        "raw_width_px": int(width),
        "length_x_um": float(length_x_um),
        "length_y_um": float(length_y_um),
    }
    if use_nonwhite_bbox:
        bbox = nonwhite_bbox(rgb)
        info.update(bbox)
        info["geometry_size_mode"] = "nonwhite_bbox"
        pixel_size_x_um = length_x_um / float(bbox["bbox_width_px"])
        pixel_size_y_um = length_y_um / float(bbox["bbox_height_px"])
    else:
        info["geometry_size_mode"] = "full_png"
        pixel_size_x_um = length_x_um / float(width)
        pixel_size_y_um = length_y_um / float(height)

    info["pixel_size_x_um"] = float(pixel_size_x_um)
    info["pixel_size_y_um"] = float(pixel_size_y_um)
    return float(pixel_size_x_um), float(pixel_size_y_um), info


def bool_from_config(config: dict[str, Any], key: str, default: bool) -> bool:
    value = scalar_value(config.get(key, default))
    if value is None:
        return bool(default)
    if isinstance(value, str):
        return value.strip().lower() in {"1", "true", "yes", "on"}
    return bool(value)


def scalar_value(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        value = value.tolist()
    if isinstance(value, (list, tuple)):
        if len(value) == 0:
            return None
        if len(value) == 1:
            return scalar_value(value[0])
        raise ValueError(f"Expected a scalar or empty value, got {value!r}")
    return value


def optional_int(value: Any) -> int | None:
    value = scalar_value(value)
    if value is None:
        return None
    if isinstance(value, str) and value.strip().lower() in {"", "none", "null", "original", "full"}:
        return None
    if isinstance(value, float) and math.isnan(value):
        return None
    return int(value)


def config_float(config: dict[str, Any], key: str, default: float) -> float:
    value = scalar_value(config.get(key, default))
    if value is None:
        value = default
    return float(value)


def config_int(config: dict[str, Any], key: str, default: int) -> int:
    value = scalar_value(config.get(key, default))
    if value is None:
        value = default
    return int(value)


def build_params(config: dict[str, Any], pixel_size_x_um: float, pixel_size_y_um: float, solver: str):
    from advanced_tools.png_phase_nmr_decay import SimulationParams

    return SimulationParams(
        pixel_size_x_um=float(pixel_size_x_um),
        pixel_size_y_um=float(pixel_size_y_um),
        diffusion_um2_per_ms=config_float(config, "diffusion_um2_per_ms", 2.0),
        bulk_t2_ms=config_float(config, "bulk_t2_ms", 3000.0),
        rho_solid_um_per_ms=config_float(config, "rho_solid_um_per_ms", 0.05),
        rho_gas_um_per_ms=config_float(config, "rho_gas_um_per_ms", 0.0),
        dt_ms=config_float(config, "decay_dt_ms", 5.0),
        t_max_ms=config_float(config, "decay_t_max_ms", 5500.0),
        max_grid_size=optional_int(config.get("max_grid_size", 500)),
        solver=solver,
        mesh_bulk_size_um=config_float(config, "mesh_bulk_size_um", 16.0),
        mesh_boundary_size_um=config_float(config, "mesh_boundary_size_um", 5.0),
        mesh_max_points=config_int(config, "mesh_max_points", 25000),
    )


def normalize(values: np.ndarray) -> np.ndarray:
    first = float(values[0]) if values.size else 0.0
    return values / max(abs(first), 1e-30)


def trim_to_decay(time_ms: np.ndarray, signal: np.ndarray) -> tuple[np.ndarray, np.ndarray, int]:
    peak_idx = int(np.nanargmax(signal))
    return time_ms[peak_idx:] - time_ms[peak_idx], signal[peak_idx:], peak_idx


def inversion_spectrum_plot_values(spectrum: np.ndarray) -> tuple[np.ndarray, str]:
    return np.asarray(spectrum, dtype=float), "spectrum amplitude"


def alpha_token(alpha: float) -> str:
    return f"alpha{alpha:g}".replace(".", "p").replace("-", "m")


def save_inversion_plot(
    *,
    time_ms: np.ndarray,
    signal: np.ndarray,
    t2_ms: np.ndarray,
    spectrum: np.ndarray,
    output_path: Path,
    title: str,
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.8), constrained_layout=True)
    axes[0].plot(time_ms, normalize(signal), color="black", lw=1.8)
    axes[0].set_xlabel("elapsed time after decay peak (ms)")
    axes[0].set_ylabel("normalized signal")
    axes[0].grid(alpha=0.25)
    axes[0].set_title("T2 decay")

    spectrum_for_plot, spectrum_ylabel = inversion_spectrum_plot_values(spectrum)
    axes[1].plot(t2_ms, spectrum_for_plot, color="#d62728", lw=1.8)
    axes[1].set_xscale("log")
    axes[1].set_xlabel("T2 (ms)")
    axes[1].set_ylabel(spectrum_ylabel)
    axes[1].grid(alpha=0.25, which="both")
    axes[1].set_title("NNLS spectrum")
    fig.suptitle(title)
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def invert_decay(decay_csv: Path, inversion_dir: Path, config: dict[str, Any], timestep: int) -> dict[str, Any]:
    inversion_dir.mkdir(parents=True, exist_ok=True)
    alpha = config_float(config, "nnls_alpha", 0.1)
    token = alpha_token(alpha)
    spectrum_csv = inversion_dir / f"timestep_{timestep:04d}_t2_spectrum_{token}.csv"
    spectrum_mat = inversion_dir / f"timestep_{timestep:04d}_t2_spectrum_{token}.mat"
    figure_png = inversion_dir / f"timestep_{timestep:04d}_t2_inversion_{token}.png"
    if bool_from_config(config, "skip_existing", True) and spectrum_csv.exists() and spectrum_mat.exists():
        return {
            "success": True,
            "skipped_existing_inversion": True,
            "spectrum_csv": str(spectrum_csv.resolve()),
            "spectrum_mat": str(spectrum_mat.resolve()),
            "figure_png": str(figure_png.resolve()) if figure_png.exists() else "",
        }

    from nmr_t2.config import NnlsConfig
    from nmr_t2.nnls import invert_single_signal_nnls

    frame = pd.read_csv(decay_csv)
    time_ms_raw = frame["time_ms"].to_numpy(dtype=float)
    signal_raw = frame["signal"].to_numpy(dtype=float)
    valid = np.isfinite(time_ms_raw) & np.isfinite(signal_raw) & (time_ms_raw >= 0)
    time_ms, signal, peak_idx = trim_to_decay(time_ms_raw[valid], signal_raw[valid])

    t2_min = config_float(config, "t2_min_ms", 1.0)
    t2_max = config_float(config, "t2_max_ms", 100000.0)
    num_bins = config_int(config, "t2_num_bins", 240)
    nnls_config = NnlsConfig(
        num_bins=num_bins,
        regularization=alpha,
        t2_min_ms=t2_min,
        t2_max_ms=t2_max,
        min_points_after_trim=config_int(config, "min_points_after_trim", 10),
    )
    result = invert_single_signal_nnls(
        time_ms,
        signal,
        signal_name=f"timestep_{timestep:04d}",
        config=nnls_config,
    )

    spectrum_norm = result.spectrum / max(float(np.nanmax(result.spectrum)), 1e-30)
    pd.DataFrame(
        {
            "t2_ms": result.t2_bins_ms,
            "spectrum": result.spectrum,
            "spectrum_normalized": spectrum_norm,
        }
    ).to_csv(spectrum_csv, index=False)
    savemat(
        spectrum_mat,
        {
            "T2_bins_ms": result.t2_bins_ms.reshape(-1, 1),
            "T2_log_ms": np.log10(result.t2_bins_ms).reshape(-1, 1),
            "combined_spectrum": result.spectrum.reshape(-1, 1),
            "spectrum": result.spectrum.reshape(-1, 1),
            "spectrum_normalized": spectrum_norm.reshape(-1, 1),
            "time_ms": time_ms.reshape(-1, 1),
            "signal": signal.reshape(-1, 1),
            "regularization": np.array([[alpha]], dtype=float),
            "peak_t2_ms": np.array([[float(result.t2_bins_ms[int(np.nanargmax(result.spectrum))])]], dtype=float),
        },
        do_compression=True,
    )
    save_inversion_plot(
        time_ms=time_ms,
        signal=signal,
        t2_ms=result.t2_bins_ms,
        spectrum=result.spectrum,
        output_path=figure_png,
        title=f"timestep {timestep:04d}, alpha={alpha:g}",
    )

    peak_idx_t2 = int(np.nanargmax(result.spectrum)) if result.spectrum.size else 0
    return {
        "success": True,
        "trim_peak_index": int(peak_idx),
        "peak_t2_ms": float(result.t2_bins_ms[peak_idx_t2]) if result.t2_bins_ms.size else float("nan"),
        "residual_norm": float(result.residual_norm),
        "roughness_norm": float(result.roughness_norm),
        "spectrum_csv": str(spectrum_csv.resolve()),
        "spectrum_mat": str(spectrum_mat.resolve()),
        "figure_png": str(figure_png.resolve()),
    }


def run_decay(
    *,
    method: str,
    png_path: Path,
    decay_dir: Path,
    params: Any,
    config: dict[str, Any],
) -> dict[str, Any]:
    decay_dir.mkdir(parents=True, exist_ok=True)
    solver = method_to_solver(method)
    if solver == "pixel_gpu":
        import advanced_tools.run_four_way_nmr_comparison_gpu as gpu_workflow
        import advanced_tools.png_phase_nmr_decay as png_decay

        gpu_workflow.build_water_operator = png_decay.build_water_operator
        gpu_workflow.PIXEL_GPU_SOLVER = str(config.get("pixel_gpu_solver", "cg"))
        gpu_workflow.PIXEL_GPU_CG_TOL = config_float(config, "pixel_gpu_cg_tol", 1e-8)
        gpu_workflow.PIXEL_GPU_CG_MAXITER = config_int(config, "pixel_gpu_cg_maxiter", 2000)
        return gpu_workflow.simulate_png_gpu_pixel(png_path, decay_dir, params)

    from advanced_tools.png_phase_nmr_decay import simulate_png

    return simulate_png(png_path, decay_dir, params)


def run_one(args: argparse.Namespace) -> dict[str, Any]:
    exp_dir = args.exp_dir.resolve()
    png_path = args.png.resolve()
    config = load_json(args.config.resolve())
    method = str(config.get("method", "png_mesh")).strip().lower()
    solver = method_to_solver(method)

    advanced_tools_dir = Path(config.get("advanced_tools_dir", "")).resolve()
    if not advanced_tools_dir.exists():
        raise FileNotFoundError(f"advanced_tools_dir does not exist: {advanced_tools_dir}")
    ensure_on_path(advanced_tools_dir.parent)

    t2_process_path = Path(config.get("t2_process_path", "")).resolve()
    if not t2_process_path.exists():
        raise FileNotFoundError(f"t2_process_path does not exist: {t2_process_path}")
    ensure_on_path(t2_process_path)
    outside_boundaries = apply_outside_boundary_convention(config)

    metadata_path = exp_dir / "run_metadata.json"
    metadata = load_json(metadata_path)
    use_bbox = bool_from_config(config, "use_nonwhite_bbox_for_geometry_size", True)
    pixel_size_x_um, pixel_size_y_um, geometry_info = compute_pixel_size_from_png(
        exp_dir=exp_dir,
        metadata=metadata,
        use_nonwhite_bbox=use_bbox,
        png_path=png_path,
    )

    output_root = Path(config.get("output_root", exp_dir / "png_nmr_results"))
    if not output_root.is_absolute():
        output_root = exp_dir / output_root
    method_root = output_root / method
    decay_dir = method_root / "decay"
    inversion_dir = method_root / "inversion"
    manifest_dir = method_root / "manifests"
    manifest_dir.mkdir(parents=True, exist_ok=True)

    timestep = int(args.timestep)
    stem = png_path.stem
    decay_csv = decay_dir / f"{stem}_nmr_decay.csv"
    skip_existing = bool_from_config(config, "skip_existing", True)
    if skip_existing and decay_csv.exists():
        decay_summary = {
            "input_path": str(png_path),
            "curve_csv": str(decay_csv.resolve()),
            "skipped_existing_decay": True,
        }
        decay_wall_time_s = 0.0
    else:
        params = build_params(config, pixel_size_x_um, pixel_size_y_um, solver)
        start = time.perf_counter()
        decay_summary = run_decay(
            method=method,
            png_path=png_path,
            decay_dir=decay_dir,
            params=params,
            config=config,
        )
        decay_wall_time_s = time.perf_counter() - start
        decay_summary["params_requested"] = asdict(params)

    inversion_summary = {"success": False, "skipped": True}
    if bool_from_config(config, "run_inversion", True):
        inversion_summary = invert_decay(decay_csv, inversion_dir, config, timestep)

    manifest = {
        "success": bool(inversion_summary.get("success", False) or not bool_from_config(config, "run_inversion", True)),
        "method": method,
        "solver": solver,
        "timestep": timestep,
        "time_s": float(args.time_s) if args.time_s is not None else None,
        "porosity": float(args.porosity) if args.porosity is not None else None,
        "exp_dir": str(exp_dir),
        "png_path": str(png_path),
        "geometry_info": geometry_info,
        "outside_boundary": outside_boundaries,
        "decay_wall_time_s": float(decay_wall_time_s),
        "decay": decay_summary,
        "inversion": inversion_summary,
    }
    manifest_path = manifest_dir / f"timestep_{timestep:04d}_manifest.json"
    manifest = json_safe(manifest)
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
    manifest["manifest_path"] = str(manifest_path.resolve())
    return manifest


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exp-dir", required=True, type=Path)
    parser.add_argument("--png", required=True, type=Path)
    parser.add_argument("--timestep", required=True, type=int)
    parser.add_argument("--time-s", type=float, default=None)
    parser.add_argument("--porosity", type=float, default=None)
    parser.add_argument("--config", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    try:
        result = run_one(parse_args(argv))
    except Exception as exc:
        result = {"success": False, "error": f"{type(exc).__name__}: {exc}"}
        print("RESULT_JSON=" + json.dumps(json_safe(result), ensure_ascii=False, allow_nan=False), flush=True)
        return 1
    print("RESULT_JSON=" + json.dumps(json_safe(result), ensure_ascii=False, allow_nan=False), flush=True)
    return 0 if result.get("success") else 1


if __name__ == "__main__":
    raise SystemExit(main())
