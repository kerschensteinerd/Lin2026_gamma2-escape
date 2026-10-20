"""
02_plot_fig4_panels.py - Plot the Figure 4 model panels

This script plots the model results as separate panels:
- Calcium maps on the reconstructed cell for the focal input and for the spot at 15%, 30% and 70% contrast
- Spatial attenuation of calcium and voltage away from the focal input
- Contrast-response curves (peak voltage, peak calcium, normalized peak calcium)

Control is gray and 25% retained inhibition is orange. The spot calcium maps
share one color scale (0 to 1, activation units); the two focal maps share
their own scale. Each panel is saved as SVG, PDF and PNG. The excitatory and
inhibitory contrast curves are plotted by 00_plot_contrast_gains.py.

Prerequisites:
    - output/fig4/nominal25/ from 01_run_fig4_simulations.py

Usage:
    python 02_plot_fig4_panels.py

Output:
    - Panels in figures/fig4_panels/

Author: Pratyush Ramakrishna
Date: 2026
"""

import os
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, Normalize

# =============================================================================
# CONFIGURATION
# =============================================================================

# Directory paths
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(SCRIPT_DIR, 'output', 'fig4', 'nominal25')
FIGURES_DIR = os.path.join(SCRIPT_DIR, 'figures', 'fig4_panels')

# Conditions: key in the data files, label, line color
CONDITIONS = [("ctrl", "Control", "#333333"),
              ("reduced", "25% retained", "#d55e00")]

MAP_CONTRASTS = (15, 30, 70)   # % contrast
CALCIUM_MAX = 1.0              # top of the shared spot calcium color scale (activation units)

plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                     "axes.spines.right": False, "svg.fonttype": "none",
                     "pdf.fonttype": 42})

# =============================================================================
# FUNCTIONS
# =============================================================================

def save_panel(fig, name):
    """Save a figure as SVG, PDF and PNG in the figures folder, then close it."""
    for extension in ("svg", "pdf", "png"):
        path = os.path.join(FIGURES_DIR, f"{name}.{extension}")
        fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"      Saved: {name}")


def nice_max(value):
    """Round a positive value up to 1, 2 or 5 times a power of ten."""
    exponent = np.floor(np.log10(value))
    for m in (1, 2, 5, 10):
        if value <= m * 10 ** exponent:
            return float(m * 10 ** exponent)


def map_panel(values, name, title, edges, index, limits, palette, vmax, source_xy=None, center_xy=None):
    """
    Plot a calcium map on the morphology, one line per segment edge.

    Parameters
    ----------
    values : ndarray
        Calcium value of each segment
    name : str
        File name without extension
    title : str
        Panel title
    edges : ndarray
        Segment edges, shape (n_edges, 2, 3), sorted by depth
    index : ndarray
        Segment index of each edge
    limits : tuple
        (xmin, xmax, ymin, ymax) of the plot (um)
    palette : Colormap
        Color map for the calcium values
    vmax : float
        Top of the color scale
    source_xy : ndarray, optional
        Position of the focal input, drawn as a cyan circle
    center_xy : ndarray, optional
        Center of the spot stimulus, drawn as a gray cross
    """
    values = np.clip(values, 0, None)      # hyperpolarized segments have no calcium increase
    assert np.isfinite(values).all()
    assert values.max() <= vmax + 1e-5, "Calcium values exceed the shared scale."

    fig, ax = plt.subplots(figsize=(3.5, 3.0), layout="constrained")
    collection = LineCollection(edges[:, :, :2], array=values[index], cmap=palette,
                                norm=Normalize(0, vmax), linewidths=.65,
                                capstyle="round", joinstyle="round")
    ax.add_collection(collection)

    if source_xy is not None:
        ax.plot(source_xy[0], source_xy[1], "o", ms=4, mfc="white", mec="#079ba8", mew=1)
    else:
        ax.plot(*center_xy, "+", color="#777777", ms=5)

    # 20 um scale bar
    x, y = limits[1] - 27, limits[2] + 7
    ax.plot([x, x + 20], [y, y], color="#333333", lw=1)
    ax.text(x + 10, y - 2.4, "20 µm", ha="center", va="top")
    ax.set(xlim=limits[:2], ylim=limits[2:], aspect="equal", title=title)
    ax.set_axis_off()

    bar = fig.colorbar(collection, ax=ax, orientation="horizontal", ticks=[0, vmax / 2, vmax])
    bar.set_label("Ca proxy (a.u.)")
    save_panel(fig, name)


def plot_spatial_attenuation(curves):
    """Plot calcium (solid) and voltage (dashed) relative to the source against path distance."""
    fig, ax = plt.subplots(figsize=(3.5, 3.0), layout="constrained")
    for key, label, color in CONDITIONS:
        curve = curves["spatial"][key]
        ax.plot(curve["distance_um"], curve["transfer_ca"], color=color,
                label=f"{label}: Ca")
        ax.plot(curve["distance_um"], curve["transfer_v"], color=color,
                linestyle="--", label=f"{label}: V")
    ax.axhline(1/np.e, color=".7", lw=.7, linestyle=":")
    ax.set(xlim=(0, 105), ylim=(0, 1), xlabel="Path distance (µm)",
           ylabel="Remote / source response")
    ax.legend(frameon=False, fontsize=7)
    save_panel(fig, "spatial_attenuation")


def plot_contrast_response(curves, field, name, ylabel, imaging=None):
    """Plot one contrast-response curve per condition, with dots at the mapped contrasts."""
    fig, ax = plt.subplots(figsize=(3.5, 3.0), layout="constrained")
    for i, (key, label, color) in enumerate(CONDITIONS):
        curve = curves["contrast"][key]
        contrast, values = np.asarray(curve["contrast_pct"]), np.asarray(curve[field])
        if field == "normalized_peak_ca":
            # Normalized between the 15% and 70% responses: only that range is defined
            keep = (contrast >= 15) & (contrast <= 70)
            contrast, values = contrast[keep], values[keep]
        ax.plot(contrast, values, color=color, label=label)
        selected = np.isin(contrast, MAP_CONTRASTS)
        ax.plot(contrast[selected], values[selected], "o", color=color, ms=3.3)
        if imaging is not None:
            ax.plot([30], [imaging[i]], "D", ms=5, mfc="white", mec=color, mew=1.2,
                    label=f"{label}: imaging")
    ax.set(xlim=(0, 73), xlabel="Contrast (%)", ylabel=ylabel)
    if field == "mean_peak_ca_au":
        ax.set(ylim=(0, CALCIUM_MAX), yticks=[0, CALCIUM_MAX / 2, CALCIUM_MAX])
    elif field == "mean_peak_dv_mV":
        ax.set_ylim(bottom=0)
    elif field == "normalized_peak_ca":
        ax.set(xlim=(15, 70), xticks=[15, 30, 50, 70], ylim=(-0.03, 1.03), yticks=[0, 0.5, 1])
    ax.legend(frameon=False, fontsize=7 if imaging is not None else 9)
    save_panel(fig, name)


# =============================================================================
# MAIN SCRIPT
# =============================================================================

def main():
    print("=" * 60)
    print("VG3 Amacrine Cell Model - Figure 4 Panels")
    print("=" * 60)
    os.makedirs(FIGURES_DIR, exist_ok=True)

    # Load results
    print("\n[1/4] Loading results")
    with np.load(os.path.join(RESULTS_DIR, "anatomy_maps.npz")) as archive:
        anatomy = {key: archive[key] for key in archive.files}
    with np.load(os.path.join(RESULTS_DIR, "finite_local_maps.npz")) as archive:
        focal = {key: archive[key] for key in archive.files}
    assert np.array_equal(anatomy["geometry_edges"], focal["geometry_edges"])
    with open(os.path.join(RESULTS_DIR, "figure_data.json")) as f:
        curves = json.load(f)

    # Draw the deepest edges first and set the plot limits
    edges, index = anatomy["geometry_edges"], anatomy["edge_segment_index"]
    order = np.argsort(edges[:, :, 2].mean(1), kind="stable")
    edges, index = edges[order], index[order]
    xymin, xymax = edges[:, :, :2].min((0, 1)), edges[:, :, :2].max((0, 1))
    limits = (xymin[0]-4, xymax[0]+4, xymin[1]-9, xymax[1]+6)
    palette = ListedColormap(plt.get_cmap("magma")(np.linspace(0, .82, 256)))

    # Calcium maps
    print("\n[2/4] Plotting calcium maps")
    source = int(focal["local_source_segment"])
    source_xy = anatomy["xyz_segments"][source, :2]
    focal_max = nice_max(max(focal[f"local_{key}_ca"].max() for key, _, _ in CONDITIONS))
    print(f"      Focal maps share a color scale of 0 to {focal_max:g}")
    for key, label, _ in CONDITIONS:
        map_panel(focal[f"local_{key}_ca"], f"focal_calcium_{key}", label,
                  edges, index, limits, palette, focal_max, source_xy=source_xy)
        for contrast in MAP_CONTRASTS:
            map_panel(anatomy[f"spot_{key}_{contrast}_ca"],
                      f"spot_calcium_{key}_{contrast}", f"{label}, {contrast}% contrast",
                      edges, index, limits, palette, CALCIUM_MAX, center_xy=anatomy["stimulus_center_xy"])

    # Spatial attenuation
    print("\n[3/4] Plotting spatial attenuation")
    plot_spatial_attenuation(curves)

    # Contrast responses
    print("\n[4/4] Plotting contrast responses")
    for field, name, ylabel, imaging in [
        ("mean_peak_dv_mV", "voltage_contrast", "Peak ΔV (mV)", None),
        ("mean_peak_ca_au", "calcium_contrast", "Ca proxy (a.u.)", None),
        ("normalized_peak_ca", "normalized_calcium_contrast", "Normalized peak calcium",
         curves["imaging_OFF_normalized30"]),
    ]:
        plot_contrast_response(curves, field, name, ylabel, imaging)

    print("\n" + "=" * 60)
    print("PLOTTING COMPLETE")
    print("=" * 60)
    print(f"\n  Figures: {FIGURES_DIR}")
    print("\nNext step: Run 03_run_calcium_sweep.py")


if __name__ == '__main__':
    main()
