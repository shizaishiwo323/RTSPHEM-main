#!/usr/bin/env python3
"""Draw a Da-Pe phase map of T2 time-sequence heatmaps from exp_* batch outputs.

Inputs
------
- A batch directory containing exp_* folders.
- Each exp_* folder must contain run_metadata.json, global_evolution.xlsx, and
  surrogate_inversion_results/T2_t*_T2.mat files.

Outputs
-------
- A Nature-style Da-Pe grid figure saved as SVG, PDF, and PNG.
- Timesteps with simulated porosity >= --porosity-max are excluded before plotting.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.ticker import LogFormatterMathtext, LogLocator
import numpy as np
import pandas as pd

from plot_nature_t2_heatmaps import (
    apply_publication_style,
    edges_from_centers,
    load_exp_scenario,
    select_expdirs,
)


PANEL_BOX_ASPECT = 2.55 / 3.75


def format_axis_value(value: float) -> str:
    return f"{value:g}"


def metadata_value(scenario: dict[str, object], key: str) -> float:
    meta = scenario["meta"]
    assert isinstance(meta, dict)
    return float(meta[key])


def plot_phase_map(
    scenarios: list[dict[str, object]],
    outdir: Path,
    *,
    vmax_percentile: float,
    y_min_ms: float,
    y_max_ms: float,
    wspace: float,
    hspace: float,
) -> list[Path]:
    da_values = sorted({metadata_value(scenario, "Da") for scenario in scenarios})
    pe_values = sorted({metadata_value(scenario, "Pe") for scenario in scenarios})
    by_coord = {
        (metadata_value(scenario, "Da"), metadata_value(scenario, "Pe")): scenario
        for scenario in scenarios
    }

    all_values = np.concatenate([scenario["amplitude"].ravel() for scenario in scenarios])
    vmax = float(np.nanpercentile(all_values, vmax_percentile))
    if not np.isfinite(vmax) or vmax <= 0:
        vmax = float(np.nanmax(all_values)) if all_values.size else 1.0

    cmap = LinearSegmentedColormap.from_list(
        "nmr_t2_amplitude",
        ["#F7F7F7", "#DDEEEE", "#77D7D1", "#33B5A5", "#3775BA", "#0F4D92", "#272727"],
        N=256,
    )
    cmap.set_bad(color="white")
    norm = Normalize(vmin=0.0, vmax=vmax)

    ncols = len(da_values)
    nrows = len(pe_values)
    cell_width_in = 1.05
    cell_height_in = cell_width_in * PANEL_BOX_ASPECT
    fig_width = 1.35 + ncols * cell_width_in
    fig_height = 1.05 + nrows * cell_height_in
    fig = plt.figure(figsize=(fig_width, fig_height))
    gs = fig.add_gridspec(
        nrows=nrows,
        ncols=ncols,
        left=0.075,
        right=0.875,
        bottom=0.115,
        top=0.895,
        wspace=wspace,
        hspace=hspace,
    )

    axes = np.empty((nrows, ncols), dtype=object)
    last_mesh = None
    row_pe_values = list(reversed(pe_values))
    for row, pe in enumerate(row_pe_values):
        for col, da in enumerate(da_values):
            ax = fig.add_subplot(gs[row, col])
            axes[row, col] = ax
            ax.set_box_aspect(PANEL_BOX_ASPECT)
            scenario = by_coord.get((da, pe))
            if scenario is None:
                ax.set_axis_off()
                continue

            t2_ms = scenario["t2_ms"]
            times_s = scenario["times_s"]
            amplitude = scenario["amplitude"]
            assert isinstance(t2_ms, np.ndarray)
            assert isinstance(times_s, np.ndarray)
            assert isinstance(amplitude, np.ndarray)

            x_edges = edges_from_centers(times_s, log=False)
            y_edges = edges_from_centers(t2_ms, log=True)
            last_mesh = ax.pcolormesh(
                x_edges,
                y_edges,
                amplitude,
                shading="auto",
                cmap=cmap,
                norm=norm,
                rasterized=True,
            )
            ax.set_yscale("log")
            ax.set_xlim(float(x_edges.min()), float(x_edges.max()))
            ax.set_ylim(y_min_ms, y_max_ms)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.tick_params(length=0)
            for spine in ax.spines.values():
                spine.set_linewidth(0.35)
                spine.set_color("#B8B8B8")

            if row == 0:
                ax.set_title(format_axis_value(da), fontsize=6.4, pad=4)
            if col == 0:
                ax.set_ylabel(format_axis_value(pe), rotation=0, ha="right", va="center", labelpad=13, fontsize=6.4)

    fig.text(0.5, 0.047, "Damkohler number, Da", ha="center", va="center", fontsize=8)
    fig.text(0.026, 0.505, "Peclet number, Pe", ha="center", va="center", rotation=90, fontsize=8)
    fig.text(0.5, 0.945, "T2 relaxation phase map", ha="center", va="center", fontsize=9)

    if last_mesh is not None:
        cbar_ax = fig.add_axes([0.905, 0.22, 0.014, 0.56])
        cbar = fig.colorbar(last_mesh, cax=cbar_ax)
        cbar.set_label("T2 amplitude (a.u.)", labelpad=4, fontsize=7)
        cbar.ax.tick_params(length=0, pad=2, labelsize=6.2)
        cbar.outline.set_linewidth(0.5)

    outputs = [
        outdir / "da_pe_t2_phase_map_nature.svg",
        outdir / "da_pe_t2_phase_map_nature.pdf",
        outdir / "da_pe_t2_phase_map_nature.png",
    ]
    fig.savefig(outputs[0])
    fig.savefig(outputs[1])
    fig.savefig(outputs[2], dpi=450)
    plt.close(fig)
    return outputs


def write_manifest(scenarios: list[dict[str, object]], outdir: Path) -> Path:
    rows = []
    for scenario in scenarios:
        meta = scenario["meta"]
        path = scenario["path"]
        assert isinstance(meta, dict)
        assert isinstance(path, Path)
        rows.append(
            {
                "exp": path.name,
                "Da": float(meta["Da"]),
                "Pe": float(meta["Pe"]),
                "geometry": meta.get("geometry") or meta.get("layoutType"),
                "loaded_timestep_count": int(meta.get("loaded_timestep_count", 0)),
                "skipped_full_porosity_count": int(meta.get("skipped_full_porosity_count", 0)),
            }
        )
    manifest = outdir / "da_pe_t2_phase_map_manifest.csv"
    pd.DataFrame(rows).sort_values(["Da", "Pe"]).to_csv(manifest, index=False)
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", type=Path, help="Directory containing exp_* batch outputs.")
    parser.add_argument(
        "--outdir",
        type=Path,
        default=None,
        help="Output directory. Defaults to run_dir/figures/da_pe_t2_phase_map.",
    )
    parser.add_argument("--exp-start", type=int, default=1, help="First exp number to include.")
    parser.add_argument("--exp-count", type=int, default=70, help="Number of exp_* directories to include.")
    parser.add_argument(
        "--porosity-max",
        type=float,
        default=0.999999,
        help="Exclude timesteps whose simulated porosity is greater than or equal to this value.",
    )
    parser.add_argument("--y-min-ms", type=float, default=1.0, help="Lower T2 axis limit in ms.")
    parser.add_argument("--y-max-ms", type=float, default=10000.0, help="Upper T2 axis limit in ms.")
    parser.add_argument(
        "--vmax-percentile",
        type=float,
        default=99.5,
        help="Global percentile used as the shared colour scale maximum.",
    )
    parser.add_argument("--wspace", type=float, default=0.10, help="Horizontal spacing between panels.")
    parser.add_argument("--hspace", type=float, default=0.12, help="Vertical spacing between panels.")
    args = parser.parse_args()

    apply_publication_style(font_size=7)
    outdir = args.outdir or args.run_dir / "figures" / "da_pe_t2_phase_map"
    outdir.mkdir(parents=True, exist_ok=True)

    expdirs = select_expdirs(args.run_dir, exp_start=args.exp_start, exp_count=args.exp_count)
    if not expdirs:
        raise ValueError(f"No selected exp_* directories found in {args.run_dir}")

    scenarios = [load_exp_scenario(expdir, porosity_max=args.porosity_max) for expdir in expdirs]
    outputs = plot_phase_map(
        scenarios,
        outdir,
        vmax_percentile=args.vmax_percentile,
        y_min_ms=args.y_min_ms,
        y_max_ms=args.y_max_ms,
        wspace=args.wspace,
        hspace=args.hspace,
    )
    manifest = write_manifest(scenarios, outdir)

    print(f"Output directory: {outdir.resolve()}")
    print(f"EXP selection: start={args.exp_start}, count={args.exp_count}")
    print(f"Porosity filter: kept timesteps with Porosity < {args.porosity_max:g}")
    print(f"Grid: {len({metadata_value(s, 'Da') for s in scenarios})} Da values x {len({metadata_value(s, 'Pe') for s in scenarios})} Pe values")
    print(f"Panel spacing: wspace={args.wspace:g}, hspace={args.hspace:g}")
    print(f"Manifest: {manifest.resolve()}")
    for output in outputs:
        print(output.resolve())


if __name__ == "__main__":
    main()
