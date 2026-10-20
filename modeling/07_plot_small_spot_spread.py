"""
07_plot_small_spot_spread.py - Plot the spread of voltage and calcium from a small light spot

Panels, in the style of Figure 4:
- Calcium and voltage maps on the reconstructed cell for control and 25% retained
  inhibition at each contrast (calcium maps share one scale, voltage maps another)
- Spatial attenuation of voltage (dashed) and calcium (solid) from the spot center
  at each contrast, with the Figure 4B focal-synapse curves for comparison
- Voltage traces at five distances from the spot center (same y-axis limits)

Prerequisites:
    - output/small_spot_<diameter>um/ from 06_run_small_spot_spread.py
    - output/fig4/nominal25/figure_data.json from 01_run_fig4_simulations.py

Usage:
    python 07_plot_small_spot_spread.py [diameter_um]     (default 50)

Output:
    - Panels in figures/small_spot_<diameter>um/

Author: Pratyush Ramakrishna
Date: 2026
"""

import os
import sys
import json
import importlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, Normalize

# =============================================================================
# CONFIGURATION
# =============================================================================

DIAMETER = float(sys.argv[1]) if len(sys.argv) > 1 else 50.0
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(SCRIPT_DIR, 'output', f'small_spot_{DIAMETER:g}um')
FIG4_DIR = os.path.join(SCRIPT_DIR, 'output', 'fig4', 'nominal25')
FIGURES_DIR = os.path.join(SCRIPT_DIR, 'figures', f'small_spot_{DIAMETER:g}um')

CONDITIONS = [("ctrl", "Control", "#333333"),
              ("reduced", "25% retained", "#d55e00")]
START_MS = 100.0
PLOT_FROM = -20.0     # ms relative to onset

fig4_plots = importlib.import_module('02_plot_fig4_panels')   # save_panel, nice_max
fig4_run = importlib.import_module('01_run_fig4_simulations')  # spatial_metrics
fig4_plots.FIGURES_DIR = FIGURES_DIR

plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                     "axes.spines.right": False, "svg.fonttype": "none",
                     "pdf.fonttype": 42})

# =============================================================================
# FUNCTIONS
# =============================================================================

def map_panel(values, name, title, edges, index, limits, palette, vmax, center_xy, label):
    """Plot a map on the morphology (one line per segment edge) with the spot outline."""
    values = np.clip(values, 0, None)
    fig, ax = plt.subplots(figsize=(3.5, 3.0), layout="constrained")
    collection = LineCollection(edges[:, :, :2], array=values[index], cmap=palette,
                                norm=Normalize(0, vmax), linewidths=.65,
                                capstyle="round", joinstyle="round")
    ax.add_collection(collection)
    ax.add_patch(plt.Circle(center_xy, DIAMETER / 2, fill=False, ec="#079ba8", lw=.8, ls="--"))
    x, y = limits[1] - 27, limits[2] + 7
    ax.plot([x, x + 20], [y, y], color="#333333", lw=1)
    ax.text(x + 10, y - 2.4, "20 µm", ha="center", va="top")
    ax.set(xlim=limits[:2], ylim=limits[2:], aspect="equal", title=title)
    ax.set_axis_off()
    bar = fig.colorbar(collection, ax=ax, orientation="horizontal", ticks=[0, vmax / 2, vmax])
    bar.set_label(label)
    fig4_plots.save_panel(fig, name)


def plot_attenuation(maps, contrast, focal_curves):
    """Voltage (dashed) and calcium (solid) relative to the spot center against path distance."""
    distance, length = maps['path_distance_from_source_um'], maps['cable_lengths']
    source = int(maps['source_segment'])
    fig, ax = plt.subplots(figsize=(3.5, 3.0), layout="constrained")
    for key, label, color in CONDITIONS:
        dv, ca = maps[f'{key}_{contrast}_dv'], maps[f'{key}_{contrast}_ca']
        for y, style, kind in [(ca, "-", "Ca"), (dv, "--", "V")]:
            if y[source] <= 1e-6:
                continue
            m = fig4_run.spatial_metrics(y / y[source], distance, length)
            ax.plot(m['distance_um'], m['normalized_response'], style, color=color, label=f"{label}: {kind}")
        f = focal_curves[key]
        ax.plot(f['distance_um'], f['transfer_v'], ":", color=color, lw=.9, alpha=.8,
                label=f"{label}: V, focal synapse (Fig. 4B)")
    ax.axhline(1 / np.e, color=".7", lw=.7, linestyle=":")
    ax.set(xlim=(0, 105), ylim=(0, 1.1), xlabel="Path distance (µm)", ylabel="Remote / center response",
           title=f"{DIAMETER:g} µm spot, {contrast}% contrast")
    ax.legend(frameon=False, fontsize=6)
    fig4_plots.save_panel(fig, f"spatial_attenuation_{contrast}")


def plot_traces(traces, contrast, ylim):
    """Voltage traces at the five distances from the spot center, control and 25% retained."""
    t = traces['t_ms'] - START_MS
    show = t >= PLOT_FROM
    distances = traces['site_distance_um']
    fig, axes = plt.subplots(1, len(distances), figsize=(3.0 * len(distances), 3.2), sharey=True,
                             layout="constrained")
    for col, distance in enumerate(distances):
        ax = axes[col]
        for key, label, color in CONDITIONS:
            ax.plot(t[show], traces[f'{key}_{contrast}'][col, show], color=color,
                    label="Control" if key == "ctrl" else "25% retained inhibition")
        ax.set(xlabel="Time from onset (ms)", ylim=ylim)
        ax.set_title("spot center (0 µm)" if col == 0 else f"{distance:.0f} µm from spot center", fontsize=9)
    axes[0].set_ylabel("Voltage (mV)")
    axes[0].legend(frameon=False, fontsize=8, loc="upper right")
    fig.suptitle(f"Voltage spread from a {DIAMETER:g} µm OFF spot at {contrast}% contrast", fontsize=10)
    for extension in ("png", "svg"):
        fig.savefig(os.path.join(FIGURES_DIR, f"voltage_traces_{contrast}.{extension}"), dpi=200,
                    bbox_inches="tight")
    plt.close(fig)
    print(f"      Saved: voltage_traces_{contrast}")


# =============================================================================
# MAIN SCRIPT
# =============================================================================

def main():
    print("=" * 60)
    print(f"VG3 Amacrine Cell Model - {DIAMETER:g} um Spot Spread Panels")
    print("=" * 60)
    os.makedirs(FIGURES_DIR, exist_ok=True)

    print("\n[1/4] Loading results")
    with np.load(os.path.join(RESULTS_DIR, 'small_spot_maps.npz')) as archive:
        maps = {key: archive[key] for key in archive.files}
    with np.load(os.path.join(RESULTS_DIR, 'small_spot_traces.npz')) as archive:
        traces = {key: archive[key] for key in archive.files}
    with open(os.path.join(RESULTS_DIR, 'small_spot_summary.json')) as f:
        summary = json.load(f)
    with open(os.path.join(FIG4_DIR, 'figure_data.json')) as f:
        focal_curves = json.load(f)['spatial']
    contrasts = summary['contrasts_pct']

    edges, index = maps["geometry_edges"], maps["edge_segment_index"]
    order = np.argsort(edges[:, :, 2].mean(1), kind="stable")
    edges, index = edges[order], index[order]
    xymin, xymax = edges[:, :, :2].min((0, 1)), edges[:, :, :2].max((0, 1))
    limits = (xymin[0]-4, xymax[0]+4, xymin[1]-9, xymax[1]+6)
    palette = ListedColormap(plt.get_cmap("magma")(np.linspace(0, .82, 256)))

    print("\n[2/4] Plotting calcium and voltage maps")
    ca_max = fig4_plots.nice_max(max(maps[f'{k}_{c}_ca'].max() for k, _, _ in CONDITIONS for c in contrasts))
    dv_max = fig4_plots.nice_max(max(maps[f'{k}_{c}_dv'].max() for k, _, _ in CONDITIONS for c in contrasts))
    print(f"      Calcium maps 0 to {ca_max:g}; voltage maps 0 to {dv_max:g} mV")
    for key, label, _ in CONDITIONS:
        for c in contrasts:
            map_panel(maps[f'{key}_{c}_ca'], f"calcium_{key}_{c}", f"{label}, {c}% contrast",
                      edges, index, limits, palette, ca_max, maps['spot_center_xy'], "Ca proxy (a.u.)")
            map_panel(maps[f'{key}_{c}_dv'], f"voltage_{key}_{c}", f"{label}, {c}% contrast",
                      edges, index, limits, palette, dv_max, maps['spot_center_xy'], "Peak ΔV (mV)")

    print("\n[3/4] Plotting spatial attenuation")
    for c in contrasts:
        plot_attenuation(maps, c, focal_curves)

    print("\n[4/4] Plotting voltage traces")
    t = traces['t_ms'] - START_MS
    all_v = np.concatenate([traces[f'{k}_{c}'][:, t >= PLOT_FROM].ravel()
                            for k, _, _ in CONDITIONS for c in contrasts])
    ylim = (all_v.min() - 1, all_v.max() + 1)
    for c in contrasts:
        plot_traces(traces, c, ylim)

    print("\n" + "=" * 60)
    print("PLOTTING COMPLETE")
    print("=" * 60)
    print(f"\n  Figures: {FIGURES_DIR}")


if __name__ == '__main__':
    main()
