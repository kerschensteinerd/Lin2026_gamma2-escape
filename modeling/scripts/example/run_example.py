#!/usr/bin/env python
"""
run_example.py - reproduce the reference panel for one stimulus (spot or bar):
3D voltage-spread skeleton, OFF/ON trace (Ctrl vs -75% inh. subplots),
contrast+calcium tuning curve. Compares full inhibition (Ctrl) vs
75%-inhibition-removed.

Prerequisites:
    ./setup.sh (installs NEURON + dependencies)

Usage:
    python scripts/example/run_example.py [--stim spot|bar] [--outdir DIR]

Output:
    outputs/example/trace_overlay.{png,svg}
    outputs/example/skeleton_3d.{png,pdf}
    outputs/example/contrast_tuning.{png,svg}
    outputs/example/contrast_tuning_data.txt
"""
import os, sys, argparse, time
import numpy as np
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
SRC = os.path.join(ROOT, "src")
sys.path.insert(0, SRC)
from run_functions import (build, load_fits, load_named_config, run_spot, run_bar,
                            off_field_corner, region_voltages, simulate_vgcc_calcium,
                            save_fig)
from plotting import plot_overlay_trace, render_3d_skeleton, plot_tuning_curve

CONFIG = "sep_g22"
SPOT_DIAMETER_UM = 100.0
BAR_W, BAR_H, BAR_DEFF = 50.0, 80.0, 100.0
POLARITY = "OFF"
T_STOP_MS = 2000.0    # long enough for the trace to decay back to baseline
CONTRASTS = [round(0.1 * i, 1) for i in range(11)]   # 0.0 ... 1.0
KO = {"ctrl": 1.0, "ko75": 0.25}                       # condition -> keep fraction
KO_LABEL = {"ctrl": "Ctrl", "ko75": "-75% inh."}
KO_SEED = {"ctrl": 0, "ko75": 75}


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--stim", choices=["spot", "bar"], default="spot")
    ap.add_argument("--outdir", default=os.path.join(ROOT, "outputs", "example"))
    return ap.parse_args()


class Progress:
    """Prints '[i/N] label ... done (Xs, ~Ys left)' after each simulation run."""
    def __init__(self, total):
        self.total = total
        self.done = 0
        self.t0 = time.time()

    def tick(self, label):
        start = time.time()
        def finish():
            self.done += 1
            elapsed = time.time() - self.t0
            per_run = elapsed / self.done
            remaining = per_run * (self.total - self.done)
            print(f"  [{self.done}/{self.total}] {label} done "
                  f"({time.time() - start:.1f}s, ~{remaining:.0f}s left)")
        return finish


def run_condition(model, fits, G_exc, G_inh, stim, keep, seed, contrast=1.0):
    """Uniform (t, all_v) wrapper over run_spot / run_bar."""
    if stim == "spot":
        _, _, t, all_v = run_spot(model, fits, SPOT_DIAMETER_UM, POLARITY, record_v=True,
                                   exc_scale=G_exc, inh_scale=G_inh,
                                   inh_keep_fraction=keep, inh_seed=seed, contrast=contrast)
        return t, all_v
    corner = off_field_corner(model)
    t, all_v, _, _ = run_bar(model, fits, corner, bar_w=BAR_W, bar_h=BAR_H, D_eff=BAR_DEFF,
                              polarity=POLARITY, record_v=True,
                              exc_scale=G_exc, inh_scale=G_inh,
                              inh_keep_fraction=keep, inh_seed=seed, contrast=contrast)
    return t, all_v


def stim_overlay(model, stim):
    """callable(ax) drawing the stimulus outline on an XY-plane axes."""
    if stim == "spot":
        sx, sy = model["soma_xy"]; R = SPOT_DIAMETER_UM / 2.0
        th = np.linspace(0, 2 * np.pi, 200)
        return lambda ax: ax.plot(sx + R * np.cos(th), sy + R * np.sin(th),
                                   color="cyan", ls="--", lw=1.5)
    corner = off_field_corner(model)
    return lambda ax: ax.add_patch(plt.Rectangle((corner[0], corner[1]), BAR_W, BAR_H,
                                                  fill=False, edgecolor="cyan", lw=1.5))


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    cfg, desc = load_named_config(CONFIG)
    print(f"Config: {desc}")
    model = build(g_pas=cfg["g_pas"])
    model["p"]["inh_rev"] = cfg.get("inh_rev", -60.0)
    model["p"]["t_stop"] = T_STOP_MS
    fits = load_fits()
    G_exc, G_inh, rest = cfg["G_exc"], cfg["G_inh"], model["rest_v"]

    total_runs = len(KO) + len(KO) * len(CONTRASTS)   # 2 trace runs + contrast sweep
    print(f"Running {total_runs} simulations at t_stop={T_STOP_MS:.0f} ms "
          f"(each run takes longer with a longer t_stop)...")
    prog = Progress(total_runs)

    # === run Ctrl / -75% inh, record everything ===
    traces, depol = {}, {}
    for key, keep in KO.items():
        finish = prog.tick(f"{args.stim} trace, {KO_LABEL[key]}")
        t, all_v = run_condition(model, fits, G_exc, G_inh, args.stim, keep, KO_SEED[key])
        finish()
        traces[key] = (t, region_voltages(model, all_v))
        depol[key] = all_v.max(axis=1) - rest

    # === OFF/ON trace, Ctrl and -75% inh. as separate subplots (shared y-axis) ===
    t_c, V_c = traces["ctrl"]; t_k, V_k = traces["ko75"]
    plot_overlay_trace(t_c, V_c["OFF"], V_c["ON"], t_k, V_k["OFF"], V_k["ON"], rest_v=rest,
                        title=f"{args.stim.upper()} stimulus, {POLARITY} polarity [{CONFIG}]",
                        out_noext=os.path.join(args.outdir, "trace_overlay"))

    # === 3D skeleton, shared color scale ===
    render_3d_skeleton(model["node_coords"], model["soma_idx"], depol, KO_LABEL,
                        title=f"{args.stim.upper()} stimulus — peak depolarization spread",
                        out_noext=os.path.join(args.outdir, "skeleton_3d"),
                        overlay_xy=stim_overlay(model, args.stim))

    # === contrast + calcium tuning curve (OFF layer) ===
    off_layer = model["off_layer"]
    dv_curve, ca_curve = {k: [] for k in KO}, {k: [] for k in KO}
    rows = ["contrast\tcondition\tOFFlayer_peak_dV_mV\tOFFlayer_calcium_norm\n"]
    for key, keep in KO.items():
        raw_ca = []
        for C in CONTRASTS:
            finish = prog.tick(f"contrast={C:.1f}, {KO_LABEL[key]}")
            _, all_v = run_condition(model, fits, G_exc, G_inh, args.stim, keep, KO_SEED[key], contrast=C)
            finish()
            dv = float(all_v[off_layer].max(axis=1).max() - rest) if off_layer.any() else 0.0
            dv_curve[key].append(dv)
            raw_ca.append(simulate_vgcc_calcium(dv, rest_v=rest))
        lo, hi = min(raw_ca), max(raw_ca)
        ca_curve[key] = [(c - lo) / (hi - lo) if hi > lo else 0.0 for c in raw_ca]
        for C, dv, ca in zip(CONTRASTS, dv_curve[key], ca_curve[key]):
            rows.append(f"{C:.1f}\t{key}\t{dv:.4f}\t{ca:.4f}\n")
    with open(os.path.join(args.outdir, "contrast_tuning_data.txt"), "w") as f:
        f.writelines(rows)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    plot_tuning_curve(axes[0], CONTRASTS, dv_curve["ctrl"], dv_curve["ko75"],
                       xlabel="contrast", ylabel="OFF-layer peak ΔV (mV)", title="Voltage response")
    plot_tuning_curve(axes[1], CONTRASTS, ca_curve["ctrl"], ca_curve["ko75"],
                       xlabel="contrast", ylabel="Calcium influx (norm.)", title="Calcium response")
    fig.suptitle(f"{args.stim.upper()} {POLARITY} — contrast tuning [{CONFIG}]")
    fig.tight_layout()
    save_fig(fig, os.path.join(args.outdir, "contrast_tuning"))
    plt.close(fig)

    print(f"Done. Outputs written to {args.outdir}")


if __name__ == "__main__":
    main()
