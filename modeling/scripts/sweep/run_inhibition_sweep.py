#!/usr/bin/env python
"""
run_inhibition_sweep.py - sweep the fraction of inhibitory synapses removed
(0/25/50/75/100%) at a fixed spot stimulus and gain, and plot max OFF-layer
depolarization vs inhibition removed.

Prerequisites:
    ./setup.sh (installs NEURON + dependencies)

Usage:
    python scripts/sweep/run_inhibition_sweep.py [--outdir DIR]

Output:
    outputs/sweep/inhibition_sweep.{png,svg}
    outputs/sweep/inhibition_sweep_data.txt
"""
import os, sys, argparse, time
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
SRC = os.path.join(ROOT, "src")
sys.path.insert(0, SRC)
from run_functions import build, load_fits, load_named_config, run_spot, save_fig
from plotting import BLACK

CONFIG = "sep_g22"
SPOT_DIAMETER_UM = 100.0
POLARITY = "OFF"
CONTRAST = 1.0
PCT_REMOVED = [0, 25, 50, 75, 100]   # % of inhibitory synapses removed


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", default=os.path.join(ROOT, "outputs", "sweep"))
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
            remaining = (elapsed / self.done) * (self.total - self.done)
            print(f"  [{self.done}/{self.total}] {label} done "
                  f"({time.time() - start:.1f}s, ~{remaining:.0f}s left)")
        return finish


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    cfg, desc = load_named_config(CONFIG)
    print(f"Config: {desc}")
    model = build(g_pas=cfg["g_pas"])
    model["p"]["inh_rev"] = cfg.get("inh_rev", -60.0)
    fits = load_fits()
    G_exc, G_inh, rest = cfg["G_exc"], cfg["G_inh"], model["rest_v"]
    off_layer = model["off_layer"]

    print(f"Running {len(PCT_REMOVED)} simulations...")
    prog = Progress(len(PCT_REMOVED))

    dv = []
    rows = ["pct_inhibition_removed\tOFFlayer_peak_dV_mV\n"]
    for pct in PCT_REMOVED:
        finish = prog.tick(f"{pct}% inh. removed")
        keep = 1.0 - pct / 100.0
        _, _, t, all_v = run_spot(model, fits, SPOT_DIAMETER_UM, POLARITY, record_v=True,
                                   exc_scale=G_exc, inh_scale=G_inh,
                                   inh_keep_fraction=keep, inh_seed=pct, contrast=CONTRAST)
        finish()
        d = float(all_v[off_layer].max(axis=1).max() - rest) if off_layer.any() else 0.0
        dv.append(d)
        rows.append(f"{pct}\t{d:.4f}\n")
    with open(os.path.join(args.outdir, "inhibition_sweep_data.txt"), "w") as f:
        f.writelines(rows)

    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.plot(PCT_REMOVED, dv, color=BLACK, marker="o", lw=2)
    ax.set_xlabel("inhibition removed (%)")
    ax.set_ylabel("OFF-layer peak ΔV (mV)")
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.set_title(f"{POLARITY} spot {SPOT_DIAMETER_UM:.0f}µm, contrast={CONTRAST} [{CONFIG}]")
    fig.tight_layout()
    save_fig(fig, os.path.join(args.outdir, "inhibition_sweep"))
    plt.close(fig)
    print(f"Done. Outputs written to {args.outdir}")


if __name__ == "__main__":
    main()
