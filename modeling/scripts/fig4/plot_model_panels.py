#!/usr/bin/env python3
"""Plot the model results as separate panel assets."""
import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, Normalize
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
CONDITIONS = [("ctrl", "Control", "#333333"),
              ("reduced", "25% retained", "#d55e00")]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path,
                        default=ROOT / "outputs/fig4_v9_retained_inhibition_20260916/nominal25")
    parser.add_argument("--outdir", type=Path, default=ROOT / "outputs/model_panels")
    args = parser.parse_args()
    data, out = args.results.resolve(), args.outdir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    with np.load(data / "anatomy_maps.npz") as archive:
        anatomy = {key: archive[key] for key in archive.files}
    with np.load(data / "finite_local_maps.npz") as archive:
        focal = {key: archive[key] for key in archive.files}
    assert np.array_equal(anatomy["geometry_edges"], focal["geometry_edges"])
    curves = json.loads((data / "figure_data.json").read_text())
    edges, index = anatomy["geometry_edges"], anatomy["edge_segment_index"]
    order = np.argsort(edges[:, :, 2].mean(1), kind="stable")
    edges, index = edges[order], index[order]
    xymin, xymax = edges[:, :, :2].min((0, 1)), edges[:, :, :2].max((0, 1))
    limits = (xymin[0]-4, xymax[0]+4, xymin[1]-9, xymax[1]+6)
    palette = ListedColormap(plt.get_cmap("magma")(np.linspace(0, .82, 256)))
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                         "axes.spines.right": False, "svg.fonttype": "none",
                         "pdf.fonttype": 42})
    outputs = []

    def save(fig, name):
        for extension in ("svg", "pdf", "png"):
            path = out / f"{name}.{extension}"
            fig.savefig(path, dpi=300, bbox_inches="tight")
            outputs.append(path.name)
        plt.close(fig)

    def map_panel(values, name, title, focal_input=False):
        assert np.isfinite(values).all() and values.min() > -1e-5
        assert values.max() <= 50 + 1e-5, "Calcium values exceed the shared scale."
        fig, ax = plt.subplots(figsize=(3.5, 3.0), layout="constrained")
        collection = LineCollection(edges[:, :, :2], array=values[index], cmap=palette,
                                    norm=Normalize(0, 50), linewidths=.65,
                                    capstyle="round", joinstyle="round")
        ax.add_collection(collection)
        if focal_input:
            source = int(focal["local_source_segment"])
            sx, sy = anatomy["xyz_segments"][source, :2]
            ax.plot(sx, sy, "o", ms=4, mfc="white", mec="#079ba8", mew=1)
        else:
            ax.plot(*anatomy["stimulus_center_xy"], "+", color="#777777", ms=5)
        x, y = limits[1]-27, limits[2]+7
        ax.plot([x, x+20], [y, y], color="#333333", lw=1)
        ax.text(x+10, y-2.4, "20 µm", ha="center", va="top")
        ax.set(xlim=limits[:2], ylim=limits[2:], aspect="equal", title=title)
        ax.set_axis_off()
        bar = fig.colorbar(collection, ax=ax, orientation="horizontal", ticks=[0, 25, 50])
        bar.set_label("Ca proxy (a.u.)")
        save(fig, name)

    for key, label, _ in CONDITIONS:
        map_panel(focal[f"local_{key}_ca"], f"focal_calcium_{key}", label,
                  focal_input=True)
        for contrast in (15, 30, 70):
            map_panel(anatomy[f"spot_{key}_{contrast}_ca"],
                      f"spot_calcium_{key}_{contrast}", f"{label}, {contrast}% contrast")

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
    save(fig, "spatial_attenuation")

    for field, name, ylabel in [
        ("mean_peak_dv_mV", "voltage_contrast", "Peak ΔV (mV)"),
        ("mean_peak_ca_au", "calcium_contrast", "Ca proxy (a.u.)"),
        ("auc_shape", "normalized_auc_contrast", "Normalized calcium AUC"),
    ]:
        fig, ax = plt.subplots(figsize=(3.5, 3.0), layout="constrained")
        for key, label, color in CONDITIONS:
            curve = curves["contrast"][key]
            contrast, values = np.asarray(curve["contrast_pct"]), np.asarray(curve[field])
            line, = ax.plot(contrast, values, color=color, label=label)
            assert np.array_equal(line.get_ydata(), values)
            selected = np.isin(contrast, [15, 30, 70])
            ax.plot(contrast[selected], values[selected], "o", color=color, ms=3.3)
        ax.set(xlim=(0, 73), xlabel="Contrast (%)", ylabel=ylabel)
        if field == "mean_peak_ca_au":
            ax.set(ylim=(0, 50), yticks=[0, 25, 50])
        elif field == "mean_peak_dv_mV":
            ax.set_ylim(bottom=0)
        ax.legend(frameon=False)
        save(fig, name)

    inputs = [data / name for name in
              ("anatomy_maps.npz", "finite_local_maps.npz", "figure_data.json")]
    manifest = {
        "inputs": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
        "outputs": outputs, "panels": len(outputs)//3,
        "calcium_scale": [0, 50],
        "geometry": "Depth-ordered pt3d XY projection; equal aspect; no smoothing",
        "normalized_auc": "Model curves only, using the saved 15%/70% endpoint normalization",
    }
    (out / "plot_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"panels": manifest["panels"], "outdir": str(out)}))


if __name__ == "__main__":
    main()
