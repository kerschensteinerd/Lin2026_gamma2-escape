"""
05_plot_example_voltage_traces.py - Reproduce the Figure 4B spatial-attenuation curve as voltage traces

Figure 4B plots the peak voltage and calcium at each point along the dendrites,
normalized to their value at a single focal excitatory synapse, against path
distance from that synapse. This script uses the same stimulus (one excitatory
synapse at the same source segment, producing a 15 mV peak depolarization at the
source in control, together with the phasic inhibition of the 100 um spot at 30%
contrast) and records the full time course of the voltage at five points taken
directly from the Figure 4B curve's distance axis, from the source out to about
100 um. Every panel uses the same y-axis limits, so the drop in peak voltage with
distance is easy to see. Control inhibition is compared with 25% retained
inhibition.

Every input is an Exp2Syn conductance synapse; no voltage clamps are used. The
peak reached at each site is checked against the saved Figure 4 focal maps
(local_ctrl_finite_dv and local_reduced_finite_dv).

Prerequisites:
    - output/fig4/calibration.json and output/fig4/nominal25/ from 01_run_fig4_simulations.py

Usage:
    python 05_plot_example_voltage_traces.py

Output:
    - figures/example_traces/voltage_trace_signal_spread.png/.svg

Author: Pratyush Ramakrishna
Date: 2026
"""

import os
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from base_functions import (
    load_json, build_model, path_distance, spot_inhibition_weights, focal_response,
    START_MS, SCRIPT_DIR,
)

# =============================================================================
# CONFIGURATION
# =============================================================================

N_SITES = 5                   # points taken from the Figure 4B spatial-attenuation curve
RETAINED_FRACTIONS = [(1.0, "#333333", "Control"),
                      (0.25, "#d55e00", "25% retained inhibition")]
TSTOP = 600.0                 # ms, long enough to show the decay of the 300 ms synaptic conductances
PLOT_FROM = -20.0             # ms relative to stimulus onset

FIG4_DIR = os.path.join(SCRIPT_DIR, "output", "fig4")
NOMINAL_DIR = os.path.join(FIG4_DIR, "nominal25")
FIGURES_DIR = os.path.join(SCRIPT_DIR, "figures", "example_traces")

plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})

# =============================================================================
# FUNCTIONS
# =============================================================================

def attenuation_curve_targets(n_sites):
    """
    Evenly spaced distances taken from the saved Figure 4B spatial-attenuation curve.

    Parameters
    ----------
    n_sites : int
        Number of distances to take

    Returns
    -------
    ndarray
        Target path distances from the focal point (um)
    """
    with open(os.path.join(NOMINAL_DIR, "figure_data.json")) as f:
        curve = json.load(f)["spatial"]["ctrl"]["distance_um"]
    index = np.linspace(0, len(curve) - 1, n_sites).round().astype(int)
    return np.asarray(curve)[index]


def pick_distance_sites(distance_from_focal_point, targets):
    """
    Find the segment nearest each target path distance from the focal point.

    Parameters
    ----------
    distance_from_focal_point : ndarray
        Path distance of every segment from the stimulation focal point (um)
    targets : iterable of float
        Distances to match (um)

    Returns
    -------
    sites : list of int
        Segment index nearest each target distance
    distances : list of float
        The actual path distance of each chosen segment (um)
    """
    sites, distances = [], []
    for target in targets:
        site = int(np.argmin(np.abs(distance_from_focal_point - target)))
        sites.append(site)
        distances.append(float(distance_from_focal_point[site]))
    return sites, distances


# =============================================================================
# MAIN SCRIPT
# =============================================================================

def main():
    print("=" * 60)
    print("VG3 Amacrine Cell Model - Spatial Attenuation as Voltage Traces")
    print("=" * 60)
    os.makedirs(FIGURES_DIR, exist_ok=True)

    print("\n[1/4] Building the cell model")
    calibration = load_json(os.path.join(FIG4_DIR, "calibration.json"))
    focal_point = int(calibration["focal"]["source_segment"])
    source_g = float(calibration["focal"]["conductance_uS"])
    model, fits = build_model()
    inh_nodes, inh_weights = spot_inhibition_weights(model, fits, calibration["focal"]["inhibition_contrast"])

    focal_distance = path_distance(model, focal_point)
    targets = attenuation_curve_targets(N_SITES)
    sites, distances = pick_distance_sites(focal_distance, targets)
    print(f"      Focal synapse: segment {focal_point}, peak conductance {source_g * 1000:.3f} nS")
    for site, distance in zip(sites, distances):
        print(f"      Recording site: segment {site}, {distance:.1f} um from the focal point")

    print("\n[2/4] Simulating the focal synapse at each retained fraction")
    traces = {}
    for fraction, color, label in RETAINED_FRACTIONS:
        r = focal_response(model, inh_nodes, inh_weights, fraction, focal_point, source_g,
                           record=sites, trace_sites=sites, tstop=TSTOP)
        traces[fraction] = (r["t"] - START_MS, r["traces"], r["peak_v"] - r["base_v"])
        print(f"      {label}: peak dV at each site (mV) = "
              + ", ".join(f"{p:.2f}" for p in traces[fraction][2]))

    print("\n[3/4] Checking the peak voltages against the saved Figure 4B focal maps")
    focal_maps = np.load(os.path.join(NOMINAL_DIR, "finite_local_maps.npz"))
    archived = {1.0: focal_maps["local_ctrl_finite_dv"], 0.25: focal_maps["local_reduced_finite_dv"]}
    for fraction, color, label in RETAINED_FRACTIONS:
        error = np.max(np.abs(traces[fraction][2] - archived[fraction][sites]))
        print(f"      {label}: max difference from the archived Figure 4B map = {error:.2e} mV")

    print("\n[4/4] Plotting the traces")
    show = traces[1.0][0] >= PLOT_FROM
    ymin = min(v[:, show].min() for t, v, _ in traces.values()) - 1
    ymax = max(v[:, show].max() for t, v, _ in traces.values()) + 1
    fig, axes = plt.subplots(1, N_SITES, figsize=(3.0 * N_SITES, 3.2), sharey=True, layout="constrained")

    for col, distance in enumerate(distances):
        ax = axes[col]
        for fraction, color, label in RETAINED_FRACTIONS:
            t, v, _ = traces[fraction]
            ax.plot(t[show], v[col, show], color=color, label=label)
        ax.set(xlabel="Time from onset (ms)", ylim=(ymin, ymax))
        label = "focal point (0 µm)" if distance == distances[0] else f"{distance:.0f} µm from focal point"
        ax.set_title(label, fontsize=9)

    axes[0].set_ylabel("Voltage (mV)")
    axes[0].legend(frameon=False, fontsize=8, loc="upper right")
    fig.suptitle("Voltage spread from the Figure 4B focal synapse (with phasic 30%-contrast spot inhibition)",
                 fontsize=10)

    for extension in ("png", "svg"):
        fig.savefig(os.path.join(FIGURES_DIR, f"voltage_trace_signal_spread.{extension}"),
                    dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("      Saved: voltage_trace_signal_spread")

    print("\n" + "=" * 60)
    print("PLOTTING COMPLETE")
    print("=" * 60)
    print(f"\n  Figures: {FIGURES_DIR}")


if __name__ == '__main__':
    main()
