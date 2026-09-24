"""
plotting.py - shared plotting helpers: overlay traces, 3D skeleton spread,
tuning curves.
"""
import os, sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401, registers 3d projection

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
from run_functions import save_fig

BLACK = "#000000"   # Ctrl, full inhibition
ORANGE = "#EB752E"  # -75% inh. removed


def plot_overlay_trace(t_ctrl, v_off_ctrl, v_on_ctrl, t_ko, v_off_ko, v_on_ko,
                        rest_v, title, out_noext, xlabel="t (ms)", ylabel="V (mV)"):
    """Ctrl and -75% inh. as separate subplots (shared y-axis); OFF solid /
    ON dashed within each, black=Ctrl / orange=-75% inh."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharey=True)
    panels = [(axes[0], t_ctrl, v_off_ctrl, v_on_ctrl, BLACK, "Ctrl"),
              (axes[1], t_ko, v_off_ko, v_on_ko, ORANGE, "-75% inh.")]
    for ax, t, v_off, v_on, color, label in panels:
        ax.plot(t, v_off, color=color, ls="-", lw=1.6, label="OFF")
        ax.plot(t, v_on, color=color, ls="--", lw=1.6, label="ON")
        ax.axhline(rest_v, color="grey", ls=":", lw=0.8)
        ax.set_xlabel(xlabel)
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
        ax.legend(fontsize=8, frameon=False)
        ax.set_title(label, fontsize=11)
    axes[0].set_ylabel(ylabel)
    fig.suptitle(title, fontsize=12)
    fig.tight_layout()
    save_fig(fig, out_noext)
    plt.close(fig)


def render_3d_skeleton(node_coords, soma_idx, depol_by_condition, condition_labels,
                        title, out_noext, overlay_xy=None, vmin=None, vmax=None,
                        cmap="inferno", formats=("png", "pdf")):
    """(n_conditions x 3) grid: 3D scatter | XY | XZ, shared color scale.

    depol_by_condition: dict condition_key -> 1D array of peak depolarization,
    computed live by the caller as all_v.max(axis=1) - rest (no saved-file
    round trip, so there's no column-index to get wrong).
    """
    nc = node_coords
    keys = list(depol_by_condition.keys())
    if vmin is None:
        vmin = min(d.min() for d in depol_by_condition.values())
    if vmax is None:
        vmax = max(d.max() for d in depol_by_condition.values())
    n = len(keys)
    fig = plt.figure(figsize=(15, 5 * n))
    for r, key in enumerate(keys):
        depol = depol_by_condition[key]
        label = condition_labels.get(key, key)

        ax3d = fig.add_subplot(n, 3, r * 3 + 1, projection="3d")
        ax3d.scatter(nc[:, 0], nc[:, 1], nc[:, 2], c=depol, cmap=cmap, s=5, vmin=vmin, vmax=vmax)
        ax3d.scatter(*nc[soma_idx], c="lime", s=80, marker="*", edgecolor="k")
        ax3d.set_title(f"{label} — 3D")

        axxy = fig.add_subplot(n, 3, r * 3 + 2)
        axxy.scatter(nc[:, 0], nc[:, 1], c=depol, cmap=cmap, s=7, vmin=vmin, vmax=vmax)
        axxy.scatter(*nc[soma_idx, :2], c="lime", s=120, marker="*", edgecolor="k")
        if overlay_xy is not None:
            overlay_xy(axxy)
        axxy.set_aspect("equal"); axxy.set_title("XY (top)")
        axxy.set_xlabel("x (µm)"); axxy.set_ylabel("y (µm)")

        axxz = fig.add_subplot(n, 3, r * 3 + 3)
        p = axxz.scatter(nc[:, 0], nc[:, 2], c=depol, cmap=cmap, s=7, vmin=vmin, vmax=vmax)
        axxz.scatter(nc[soma_idx, 0], nc[soma_idx, 2], c="lime", s=120, marker="*", edgecolor="k")
        axxz.set_title("XZ (IPL depth)"); axxz.set_xlabel("x (µm)"); axxz.set_ylabel("z (µm)")
        fig.colorbar(p, ax=axxz, label="peak depol (mV)")

    fig.suptitle(title, fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    save_fig(fig, out_noext, formats=formats)
    plt.close(fig)


def plot_tuning_curve(ax, x, y_ctrl, y_ko, xlabel, ylabel, title=None):
    """Draw one black(Ctrl)/orange(-75% inh.) curve onto a provided Axes."""
    ax.plot(x, y_ctrl, color=BLACK, ls="-", marker="o", lw=2, label="Ctrl")
    ax.plot(x, y_ko, color=ORANGE, ls="-", marker="o", lw=2, label="-75% inh.")
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.legend(fontsize=9, frameon=False)
    if title:
        ax.set_title(title, fontsize=11)
    return ax
