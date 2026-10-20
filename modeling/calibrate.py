#!/usr/bin/env python
"""
calibrate.py -- check the model against Kim 2015 and Grimes 2011 (and, with
--solve, re-derive the Kv density).

Validated (these are what the model is fitted to):
  1. Ctrl OFF voltage size tuning at the soma, 100% contrast, spot diameters
     50-600 um                                       vs Kim 2015 Fig. 1C (OFF)
  2. Kv: leak-subtracted clamp at the true soma, hold -90, 200 ms step to 0 mV
                                                     vs Grimes 2011 Fig. 3E: ~993 pA
Reported somatic measurements and their reference values:
  3. somatic input resistance                       vs Grimes 2011: 1075 +/- 477 MOhm
  4. somatic v-clamp EPSC (hold -60) / IPSC (hold 0), 100 um OFF spot, 100%,
     Kv off (Kim's Cs/TEA internal)                 vs Kim 2015 Fig. 1F: 1.0 / 1.2 nS
The clamp is placed on model["soma_idx"] (max native diameter, high-z INL
cluster); the SWC type-1 'soma' label is a thin terminal and must not be used.

Usage: python calibrate.py [--solve]
"""
import os, sys, math, argparse
from concurrent.futures import ProcessPoolExecutor

ROOT = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(ROOT, "src")
sys.path.insert(0, SRC)
sys.path.insert(0, ROOT)
from base_functions import G_EXC, G_INH, INH_REV, SPOT_DIAMETER as SPOT_DIAMETER_UM

D, START = SPOT_DIAMETER_UM, 100.0
KIM_EXC_NS, KIM_INH_NS = 1.0, 1.2
GRIMES_RIN, GRIMES_RIN_SD, GRIMES_KV_PA = 1075.0, 477.0, 993.0
# Kim 2015 Fig. 1C, OFF, somatic dV (mV) vs spot diameter (radius 25-300 um), read from the figure
KIM_1C_DIAMS = [50, 100, 150, 200, 300, 400, 500, 600]
KIM_1C_DV = [8.7, 11.6, 9.1, 5.8, 0.5, -3.0, -4.1, -5.9]


def _setup(active_kv, gkbar=None):
    sys.path.insert(0, SRC)
    from run_functions import build
    ov = {"active_kv": active_kv}
    if gkbar is not None:
        ov["gkbar"] = gkbar
    m = build(overrides=ov)
    m["p"].update(inh_rev=INH_REV, start_time=START, t_stop=300.0)
    return m


def secant(measure, x_ref, target, it=8, tol=0.01):
    x0 = math.log(x_ref); y0 = math.log(measure(x_ref))
    x1 = x0 + math.log(target) - y0; y1 = math.log(measure(math.exp(x1)))
    for _ in range(it):
        if abs(math.exp(y1) / target - 1) < tol:
            break
        sl = (y1 - y0) / (x1 - x0)
        x0, y0, x1 = x1, y1, x1 + (math.log(target) - y1) / max(sl, 0.05)
        y1 = math.log(measure(math.exp(x1)))
    return math.exp(x1), math.exp(y1)


def task_size(diameter):
    """Signed peak somatic dV for a Ctrl OFF spot at 100% contrast."""
    from run_functions import load_fits, run_spot, region_voltages
    m = _setup(active_kv=True)
    _, _, t, all_v = run_spot(m, load_fits(), float(diameter), "OFF", record_v=True,
                              exc_scale=G_EXC, inh_scale=G_INH, contrast=1.0)
    s = region_voltages(m, all_v)["soma"]
    d = s[t >= START] - s[(t >= START - 10) & (t < START)].mean()
    return diameter, float(d.max() if abs(d.max()) >= abs(d.min()) else d.min())


def task_rin(_):
    from neuron import h
    m = _setup(active_kv=True)
    s = m["segments"][m["soma_idx"]]
    h.finitialize(-45.0); imp = h.Impedance(); imp.loc(s.x, sec=s.sec); imp.compute(0, 1)
    return {"rin_MOhm": imp.input(s.x, sec=s.sec)}


def task_synapses(_):
    import numpy as np
    from run_functions import load_fits, exc_eff_weight, inh_eff_weight
    from vg3_base_functions import run_vclamp_spot
    m = _setup(active_kv=False)
    fits = load_fits(); soma = m["segments"][m["soma_idx"]]
    p = dict(m["p"]); p.update(d50_psi=1e9, d50_inh=1e-9)    # same disabled global scales as run_spot

    def clamp(gexc, ginh, vhold):
        t, i, *_ = run_vclamp_spot(m["cell"], m["segments"], m["node_coords"], m["off_nodes"], m["inh_nodes"],
                                   p, m["soma_xy"], D, v_hold=vhold,
                                   exc_weight=exc_eff_weight(fits, "OFF", D) * gexc,
                                   inh_weight=inh_eff_weight(fits, "OFF", D) * ginh,
                                   n_events=1, start_time=START, t_stop=250.0, dt=p["dt"], clamp_seg=soma)
        i = np.asarray(i) * 1000.0
        d = i[t >= START] - i[(t >= START - 10) & (t < START)].mean()
        return abs(d.min() if vhold < -30 else d.max()) / 60.0

    return {"epsc_nS": clamp(G_EXC, 0.0, -60.0), "ipsc_nS": clamp(0.0, G_INH, 0.0)}


def task_kv(solve):
    import numpy as np
    from neuron import h
    m = _setup(active_kv=True)
    secs = list(m["cell"].all); soma = m["segments"][m["soma_idx"]]
    gk0 = float(m["p"]["gkbar"])

    def step(gk):
        for s in secs:
            s.gkbar_khh_vg3 = gk
        c = h.SEClamp(soma.x, sec=soma.sec)
        c.dur1, c.amp1, c.dur2, c.amp2, c.dur3, c.amp3, c.rs = 100.0, -90.0, 200.0, 0.0, 50.0, -90.0, 0.001
        iv = h.Vector().record(c._ref_i); tv = h.Vector().record(h._ref_t)
        h.CVode().active(0); h.dt = m["p"]["dt"]; h.finitialize(-90.0)
        while h.t < 350.0:
            h.fadvance()
        t, i = np.array(tv), np.array(iv) * 1000.0
        del c
        return float(i[(t >= 280.0) & (t < 300.0)].mean())

    leak = step(0.0)
    kv = lambda gk: step(gk) - leak
    out = {"kv_0mV_pA": kv(gk0), "gkbar": gk0}
    if solve:
        out["gkbar_solved"], _ = secant(kv, gk0, GRIMES_KV_PA)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--solve", action="store_true", help="re-derive the Kv density (gkbar)")
    a = ap.parse_args()
    with ProcessPoolExecutor(max_workers=6) as ex:
        sizes = ex.map(task_size, KIM_1C_DIAMS)
        fr = ex.submit(task_rin, None); fs = ex.submit(task_synapses, None); fk = ex.submit(task_kv, a.solve)
        size = dict(sizes)
        r = {**fr.result(), **fs.result(), **fk.result()}

    print("1. Ctrl OFF size tuning at the soma, 100% contrast (model / Kim 2015 Fig. 1C, mV):")
    sse = 0.0
    for d, k in zip(KIM_1C_DIAMS, KIM_1C_DV):
        sse += (size[d] - k) ** 2
        print(f"     {d:4d} um: {size[d]:+6.2f} / {k:+5.1f}")
    print(f"     SSE {sse:.1f} mV^2")
    ok_kv = abs(r["kv_0mV_pA"] / GRIMES_KV_PA - 1) < 0.05
    print(f"2. Kv leak-subtracted current at 0 mV: {r['kv_0mV_pA']:.0f} pA (Grimes ~{GRIMES_KV_PA:.0f}) "
          f"{'OK' if ok_kv else 'CHECK'}; gkbar {r['gkbar']:.5f}")
    z = (r["rin_MOhm"] - GRIMES_RIN) / GRIMES_RIN_SD
    print(f"3. somatic input resistance: {r['rin_MOhm']:.0f} MOhm (Grimes {GRIMES_RIN:.0f} +/- {GRIMES_RIN_SD:.0f}; "
          f"{z:+.1f} SD) -- known deviation")
    print(f"4. somatic clamp EPSC / IPSC: {r['epsc_nS']:.2f} / {r['ipsc_nS']:.2f} nS (Kim {KIM_EXC_NS} / {KIM_INH_NS}; "
          f"x{r['epsc_nS'] / KIM_EXC_NS:.1f} / x{r['ipsc_nS'] / KIM_INH_NS:.1f}) -- known deviation")
    if a.solve:
        print(f"\nre-derived gkbar for ~{GRIMES_KV_PA:.0f} pA: {r['gkbar_solved']:.5f} (current {r['gkbar']})")


if __name__ == "__main__":
    main()
