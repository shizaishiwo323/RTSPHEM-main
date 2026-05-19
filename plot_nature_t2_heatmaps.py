#!/usr/bin/env python3
"""Draw Nature-style T2 spectrum sequence heatmaps for each scenario."""

from __future__ import annotations

import argparse
import json
import os
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/matplotlib_codex_cache")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.ticker import LogFormatterMathtext, LogLocator, ScalarFormatter
import numpy as np
import pandas as pd
import scipy.io as sio


def apply_publication_style(font_size: int = 8) -> None:
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
    plt.rcParams["svg.fonttype"] = "none"
    plt.rcParams["font.size"] = font_size
    plt.rcParams["axes.spines.right"] = False
    plt.rcParams["axes.spines.top"] = False
    plt.rcParams["axes.linewidth"] = 0.8
    plt.rcParams["legend.frameon"] = False
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42


def edges_from_centers(values: np.ndarray, *, log: bool = False) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if values.size == 1:
        if log:
            return np.array([values[0] / np.sqrt(10.0), values[0] * np.sqrt(10.0)])
        return np.array([values[0] - 0.5, values[0] + 0.5])

    if log:
        log_values = np.log10(values)
        mids = (log_values[:-1] + log_values[1:]) / 2
        first = log_values[0] - (mids[0] - log_values[0])
        last = log_values[-1] + (log_values[-1] - mids[-1])
        return 10 ** np.r_[first, mids, last]

    mids = (values[:-1] + values[1:]) / 2
    first = values[0] - (mids[0] - values[0])
    last = values[-1] + (values[-1] - mids[-1])
    return np.r_[first, mids, last]


def load_sample_times(samples_csv: Path) -> dict[str, float]:
    samples = pd.read_csv(samples_csv)
    required = {"sample_id", "time_s"}
    missing = required.difference(samples.columns)
    if missing:
        missing_text = ", ".join(sorted(missing))
        raise ValueError(f"{samples_csv} is missing required columns: {missing_text}")
    return dict(zip(samples["sample_id"].astype(str), samples["time_s"].astype(float)))


def load_scenario(scenedir: Path, time_by_sample: dict[str, float]) -> dict[str, object]:
    spectrum_path = scenedir / "scenario_spectrum_sequence.xlsx"
    config_path = scenedir / "scenario_config.json"

    spectra = pd.read_excel(spectrum_path, sheet_name="spectra_by_sample")

    columns = [col for col in spectra.columns if col != "t2_ms"]
    col_info: list[tuple[float, str]] = []
    for col in columns:
        sample_id = str(col).split("__")[0]
        if sample_id in time_by_sample:
            col_info.append((float(time_by_sample[sample_id]), col))
    if not col_info:
        raise ValueError(f"No sample columns in {spectrum_path} matched {len(time_by_sample)} sample times.")
    col_info.sort(key=lambda item: item[0])

    ordered_cols = [col for _, col in col_info]
    times_s = np.array([time_s for time_s, _ in col_info], dtype=float)
    t2_ms = spectra["t2_ms"].to_numpy(dtype=float)
    amplitude = spectra[ordered_cols].to_numpy(dtype=float)
    amplitude = np.nan_to_num(amplitude, nan=0.0, posinf=0.0, neginf=0.0)
    amplitude[amplitude < 0] = 0.0

    meta: dict[str, object] = {}
    if config_path.exists():
        meta = json.loads(config_path.read_text(encoding="utf-8"))

    return {
        "path": scenedir,
        "t2_ms": t2_ms,
        "times_s": times_s,
        "amplitude": amplitude,
        "meta": meta,
    }


def _read_exp_metadata(expdir: Path) -> dict[str, object]:
    metadata_path = expdir / "run_metadata.json"
    if not metadata_path.exists():
        return {"scenario_id": expdir.name}

    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    params = metadata.get("parameters") or {}
    if isinstance(params, dict):
        metadata.setdefault("Da", params.get("Da"))
        metadata.setdefault("Pe", params.get("Pe"))
    metadata.setdefault("geometry", metadata.get("layoutType"))
    metadata.setdefault("scenario_id", metadata.get("run_id") or expdir.name)
    return metadata


def _timestep_from_inversion_path(path: Path) -> int | None:
    match = re.search(r"T2_t(\d+)_T2\.mat$", path.name)
    if not match:
        return None
    return int(match.group(1))


def load_exp_scenario(expdir: Path, *, porosity_max: float) -> dict[str, object]:
    evolution_path = expdir / "global_evolution.xlsx"
    inversion_dir = expdir / "surrogate_inversion_results"
    if not evolution_path.exists():
        raise FileNotFoundError(f"{expdir} is missing global_evolution.xlsx")
    if not inversion_dir.exists():
        raise FileNotFoundError(f"{expdir} is missing surrogate_inversion_results")

    evolution = pd.read_excel(evolution_path)
    required = {"TimeStep", "Time_s", "Porosity"}
    missing = required.difference(evolution.columns)
    if missing:
        missing_text = ", ".join(sorted(missing))
        raise ValueError(f"{evolution_path} is missing required columns: {missing_text}")

    evolution = evolution[["TimeStep", "Time_s", "Porosity"]].copy()
    evolution["TimeStep"] = evolution["TimeStep"].astype(int)
    by_timestep = evolution.set_index("TimeStep")

    records: list[tuple[int, float, np.ndarray, np.ndarray]] = []
    skipped_full_porosity = 0
    skipped_missing_metadata = 0
    for mat_path in sorted(inversion_dir.glob("T2_t*_T2.mat")):
        timestep = _timestep_from_inversion_path(mat_path)
        if timestep is None or timestep not in by_timestep.index:
            skipped_missing_metadata += 1
            continue

        row = by_timestep.loc[timestep]
        porosity = float(row["Porosity"])
        if porosity >= porosity_max:
            skipped_full_porosity += 1
            continue

        data = sio.loadmat(mat_path)
        t2_ms = np.ravel(data.get("T2_bins_ms", data.get("T2_bins"))).astype(float)
        spectrum = np.ravel(data["combined_spectrum"]).astype(float)
        spectrum = np.nan_to_num(spectrum, nan=0.0, posinf=0.0, neginf=0.0)
        spectrum[spectrum < 0] = 0.0
        records.append((timestep, float(row["Time_s"]), t2_ms, spectrum))

    if not records:
        raise ValueError(f"No usable T2 spectra remain in {expdir} after porosity filtering.")

    records.sort(key=lambda item: item[0])
    t2_ms = records[0][2]
    for timestep, _, current_t2_ms, spectrum in records:
        if current_t2_ms.shape != t2_ms.shape or not np.allclose(current_t2_ms, t2_ms):
            raise ValueError(f"T2 bins differ at {expdir} timestep {timestep}.")
        if spectrum.shape != t2_ms.shape:
            raise ValueError(f"Spectrum length does not match T2 bins at {expdir} timestep {timestep}.")

    meta = _read_exp_metadata(expdir)
    meta["loaded_timestep_count"] = len(records)
    meta["skipped_full_porosity_count"] = skipped_full_porosity
    meta["skipped_missing_metadata_count"] = skipped_missing_metadata

    return {
        "path": expdir,
        "t2_ms": t2_ms,
        "times_s": np.array([time_s for _, time_s, _, _ in records], dtype=float),
        "amplitude": np.column_stack([spectrum for _, _, _, spectrum in records]),
        "meta": meta,
    }


def select_expdirs(run_dir: Path, *, exp_start: int, exp_count: int | None) -> list[Path]:
    expdirs = sorted(path for path in run_dir.glob("exp_*") if path.is_dir())
    selected = []
    for expdir in expdirs:
        match = re.fullmatch(r"exp_(\d+)", expdir.name)
        if not match:
            continue
        exp_number = int(match.group(1))
        if exp_number < exp_start:
            continue
        if exp_count is not None and exp_number >= exp_start + exp_count:
            continue
        selected.append(expdir)
    return selected


def scenario_title(scenedir: Path, meta: dict[str, object]) -> str:
    scenario_id = str(meta.get("scenario_id") or scenedir.name.split("_Da_")[0])
    da = meta.get("Da")
    pe = meta.get("Pe")
    geometry = meta.get("geometry") or meta.get("layoutType")
    if da is None or pe is None or geometry is None:
        return scenario_id
    return f"{scenario_id}  |  Da={float(da):g}, Pe={float(pe):g}, {geometry}"


def plot_scenario(
    scenario: dict[str, object],
    outdir: Path,
    *,
    vmax: float,
    cmap: LinearSegmentedColormap,
    xscale: str,
    y_min_ms: float,
    y_max_ms: float,
) -> list[Path]:
    scenedir = scenario["path"]
    assert isinstance(scenedir, Path)
    t2_ms = scenario["t2_ms"]
    times_s = scenario["times_s"]
    amplitude = scenario["amplitude"]
    meta = scenario["meta"]
    assert isinstance(meta, dict)
    assert isinstance(t2_ms, np.ndarray)
    assert isinstance(times_s, np.ndarray)
    assert isinstance(amplitude, np.ndarray)

    scenario_id = str(meta.get("scenario_id") or scenedir.name.split("_Da_")[0])
    x_edges = edges_from_centers(times_s, log=(xscale == "log"))
    y_edges = edges_from_centers(t2_ms, log=True)

    fig, ax = plt.subplots(figsize=(3.75, 2.55), constrained_layout=True)
    mesh = ax.pcolormesh(
        x_edges,
        y_edges,
        amplitude,
        shading="auto",
        cmap=cmap,
        norm=Normalize(vmin=0.0, vmax=vmax),
        rasterized=True,
    )

    ax.set_xscale(xscale)
    ax.set_yscale("log")
    ax.set_xlim(float(x_edges.min()), float(x_edges.max()))
    ax.set_ylim(y_min_ms, y_max_ms)
    ax.set_xlabel("Time, t (s)", labelpad=3)
    ax.set_ylabel("T2 relaxation time (ms)", labelpad=3)
    ax.set_title(scenario_title(scenedir, meta), loc="left", pad=5, fontsize=8.5)
    ax.set_frame_on(False)
    ax.tick_params(axis="both", which="both", length=0, pad=2)
    if xscale == "log":
        ax.xaxis.set_major_locator(LogLocator(base=10))
        ax.xaxis.set_major_formatter(LogFormatterMathtext(base=10))
        ax.xaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10) * 0.1))
        ax.xaxis.set_minor_formatter(plt.NullFormatter())
    else:
        ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.yaxis.set_major_locator(LogLocator(base=10))
    ax.yaxis.set_major_formatter(LogFormatterMathtext(base=10))
    ax.yaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10) * 0.1))
    ax.yaxis.set_minor_formatter(plt.NullFormatter())

    cbar = fig.colorbar(mesh, ax=ax, fraction=0.046, pad=0.025)
    cbar.set_label("T2 amplitude (a.u.)", labelpad=4)
    cbar.ax.tick_params(length=0, pad=2)
    cbar.outline.set_linewidth(0.6)

    for spine in cbar.ax.spines.values():
        spine.set_linewidth(0.6)

    outputs = [
        outdir / f"{scenario_id}_t2_time_s_heatmap_nature.svg",
        outdir / f"{scenario_id}_t2_time_s_heatmap_nature.pdf",
        outdir / f"{scenario_id}_t2_time_s_heatmap_nature.png",
    ]
    fig.savefig(outputs[0], bbox_inches="tight")
    fig.savefig(outputs[1], bbox_inches="tight")
    fig.savefig(outputs[2], dpi=300, bbox_inches="tight")
    plt.close(fig)
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", type=Path, help="Directory containing scenario_* outputs or exp_* batch outputs.")
    parser.add_argument(
        "--outdir",
        type=Path,
        default=None,
        help="Output directory. Defaults to run_dir/figures/nature_t2_time_s_heatmaps.",
    )
    parser.add_argument(
        "--samples-csv",
        type=Path,
        default=Path("Data/metadata/samples.csv"),
        help="Sample metadata table containing sample_id and time_s.",
    )
    parser.add_argument(
        "--xscale",
        choices=["log", "linear"],
        default="linear",
        help="Scale for the physical-time x-axis. Use log when time_s spans orders of magnitude and early-time detail matters.",
    )
    parser.add_argument(
        "--y-min-ms",
        type=float,
        default=1.0,
        help="Lower T2 axis limit in ms.",
    )
    parser.add_argument(
        "--y-max-ms",
        type=float,
        default=10000.0,
        help="Upper T2 axis limit in ms.",
    )
    parser.add_argument(
        "--vmax-percentile",
        type=float,
        default=99.5,
        help="Global percentile used as the shared colour scale maximum.",
    )
    parser.add_argument(
        "--exp-start",
        type=int,
        default=1,
        help="First exp number to include when run_dir contains exp_* batch outputs.",
    )
    parser.add_argument(
        "--exp-count",
        type=int,
        default=None,
        help="Number of exp_* directories to include. For example, 70 includes exp_001 through exp_070 when --exp-start is 1.",
    )
    parser.add_argument(
        "--porosity-max",
        type=float,
        default=0.999999,
        help="Exclude exp_* timesteps whose simulated porosity is greater than or equal to this value.",
    )
    args = parser.parse_args()

    apply_publication_style()
    run_dir = args.run_dir
    default_dirname = "nature_t2_time_s_heatmaps" if args.xscale == "linear" else "nature_t2_time_s_heatmaps_log"
    outdir = args.outdir or run_dir / "figures" / default_dirname
    outdir.mkdir(parents=True, exist_ok=True)

    scenedirs = sorted(path for path in run_dir.glob("scenario_*") if path.is_dir())
    if scenedirs:
        time_by_sample = load_sample_times(args.samples_csv)
        scenarios = [load_scenario(scenedir, time_by_sample) for scenedir in scenedirs]
        mode = "scenario"
    else:
        expdirs = select_expdirs(run_dir, exp_start=args.exp_start, exp_count=args.exp_count)
        if not expdirs:
            raise ValueError(f"No scenario_* or selected exp_* directories found in {run_dir}")
        scenarios = [load_exp_scenario(expdir, porosity_max=args.porosity_max) for expdir in expdirs]
        mode = "exp"

    all_values = np.concatenate([scenario["amplitude"].ravel() for scenario in scenarios])
    vmax = float(np.nanpercentile(all_values, args.vmax_percentile))
    if not np.isfinite(vmax) or vmax <= 0:
        vmax = float(np.nanmax(all_values)) if all_values.size else 1.0

    cmap = LinearSegmentedColormap.from_list(
        "nmr_t2_amplitude",
        ["#F7F7F7", "#DDEEEE", "#77D7D1", "#33B5A5", "#3775BA", "#0F4D92", "#272727"],
        N=256,
    )
    cmap.set_bad(color="white")

    written: list[Path] = []
    for scenario in scenarios:
        written.extend(
            plot_scenario(
                scenario,
                outdir,
                vmax=vmax,
                cmap=cmap,
                xscale=args.xscale,
                y_min_ms=args.y_min_ms,
                y_max_ms=args.y_max_ms,
            )
        )

    print(f"Output directory: {outdir.resolve()}")
    print(f"Input mode: {mode}")
    if mode == "scenario":
        print(f"Sample metadata: {args.samples_csv.resolve()}")
    else:
        skipped_full = sum(int(scenario["meta"].get("skipped_full_porosity_count", 0)) for scenario in scenarios)
        loaded_steps = sum(int(scenario["meta"].get("loaded_timestep_count", 0)) for scenario in scenarios)
        print(f"EXP selection: start={args.exp_start}, count={args.exp_count or 'all'}")
        print(f"Porosity filter: kept timesteps with Porosity < {args.porosity_max:g}")
        print(f"Loaded spectra: {loaded_steps}; skipped 100%-porosity spectra: {skipped_full}")
    print(f"X axis: time_s ({args.xscale} scale)")
    print(f"Y axis: T2 from {args.y_min_ms:g} to {args.y_max_ms:g} ms")
    print(f"Shared colour scale vmax ({args.vmax_percentile:g}th percentile): {vmax:.6g}")
    for path in written:
        print(path.resolve())


if __name__ == "__main__":
    main()
