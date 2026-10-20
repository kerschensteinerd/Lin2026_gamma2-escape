"""
00_plot_contrast_gains.py - Plot the excitatory and inhibitory contrast-response curves

This script plots the synaptic drive of the model against stimulus contrast for a
100 um OFF spot centered on the soma:
- Excitation: bipolar-cell input, Naka-Rushton contrast gain
- Inhibition: phasic amacrine-cell input, Naka-Rushton contrast gain with a lower
  half-maximal contrast, for control and for 25% retained inhibition

Each curve is the summed peak synaptic conductance (sum of the Exp2Syn weights;
all synapses fire together). At 100% contrast the somatic voltage clamp of the
model reads Kim et al. (2015) Fig. 1F (1.0 nS excitation, 1.2 nS inhibition).

Prerequisites:
    - Files in data/ and src/, and a native build of mod/khh_vg3.mod (see README.md)

Usage:
    python 00_plot_contrast_gains.py

Output:
    - figures/fig4_panels/contrast_gains.svg/.pdf/.png

Author: Pratyush Ramakrishna
Date: 2026
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from base_functions import (build_model, gaussian_spot_weights, exc_eff_weight, inh_eff_weight,
                            G_EXC, G_INH, SPOT_DIAMETER, SPOT_SIGMA, POLARITY, SCRIPT_DIR)

# =============================================================================
# CONFIGURATION
# =============================================================================

FIGURES_DIR = os.path.join(SCRIPT_DIR, 'figures', 'fig4_panels')
NOMINAL_FRACTION = 0.25

EXCITATION_COLOR = "#079ba8"
CONDITIONS = [("ctrl", "Control", "#333333", 1.0),
              ("reduced", "25% retained", "#d55e00", NOMINAL_FRACTION)]

plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                     "axes.spines.right": False, "svg.fonttype": "none",
                     "pdf.fonttype": 42})

# =============================================================================
# FUNCTIONS
# =============================================================================

def naka_rushton(c, c50, n):
    """Contrast gain rescaled so that 100% contrast gives 1."""
    c = np.asarray(c, dtype=float)
    return c ** n / (c50 ** n + c ** n) * (c50 ** n + 1.0)


def save_panel(fig, name):
    """Save a figure as SVG, PDF and PNG in the figures folder, then close it."""
    for extension in ("svg", "pdf", "png"):
        fig.savefig(os.path.join(FIGURES_DIR, f"{name}.{extension}"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"      Saved: {name}")


# =============================================================================
# MAIN SCRIPT
# =============================================================================

def main():
    print("=" * 60)
    print("VG3 Amacrine Cell Model - Excitatory and Inhibitory Contrast Curves")
    print("=" * 60)
    os.makedirs(FIGURES_DIR, exist_ok=True)

    print("\n[1/2] Summing the synaptic weights of the 100 um OFF spot")
    model, fits = build_model()
    p, coords = model['p'], model['node_coords']
    exc_spot = gaussian_spot_weights(coords[model['off_nodes'], :2], model['soma_xy'], SPOT_DIAMETER, SPOT_SIGMA)
    inh_spot = gaussian_spot_weights(coords[model['inh_nodes'], :2], model['soma_xy'], SPOT_DIAMETER, SPOT_SIGMA)
    g_exc = exc_eff_weight(fits, POLARITY, SPOT_DIAMETER) * G_EXC * exc_spot[exc_spot > 1e-6].sum() * 1000
    g_inh = inh_eff_weight(fits, POLARITY, SPOT_DIAMETER) * G_INH * inh_spot[inh_spot > 1e-6].sum() * 1000
    c50_e, n_e, c50_i, n_i = p['contrast_c50'], p['contrast_n'], p['inh_c50'], p['inh_cn']
    print(f"      100% contrast: excitation {g_exc:.2f} nS, inhibition {g_inh:.2f} nS")
    print(f"      Excitation half-max {c50_e:.0%} (n {n_e:g}); inhibition half-max {c50_i:.0%} (n {n_i:g})")

    print("\n[2/2] Plotting")
    contrast = np.linspace(0, 1, 401)
    fig, ax = plt.subplots(figsize=(3.5, 3.0), layout="constrained")
    ax.plot(contrast * 100, g_exc * naka_rushton(contrast, c50_e, n_e), color=EXCITATION_COLOR,
            label="Excitation")
    for key, label, color, fraction in CONDITIONS:
        ax.plot(contrast * 100, fraction * g_inh * naka_rushton(contrast, c50_i, n_i), color=color,
                linestyle="-" if key == "ctrl" else "--", label=f"Inhibition, {label.lower()}")
    for c in (15, 30, 70):
        ax.axvline(c, color=".7", lw=.7, linestyle=":")
    ax.set(xlim=(0, 100), ylim=(0, None), xticks=[0, 15, 30, 70, 100],
           xlabel="Contrast (%)", ylabel="Synaptic conductance (nS)")
    ax.legend(frameon=False, fontsize=7)
    save_panel(fig, "contrast_gains")

    print("\n" + "=" * 60)
    print("PLOTTING COMPLETE")
    print("=" * 60)
    print(f"\n  Figures: {FIGURES_DIR}")
    print("\nNext step: Run 01_run_fig4_simulations.py")


if __name__ == '__main__':
    main()
