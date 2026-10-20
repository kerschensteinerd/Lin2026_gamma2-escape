"""
run_functions.py - build the cid3 cell from fitted parameters and run a spot or
bar stimulus on it. Used by the model simulation and calibration scripts.
"""
import os, sys, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))          # .../src
ROOT = os.path.dirname(HERE)                                 # repo root
DATA = os.path.join(ROOT, "data")
FITTED = os.path.join(DATA, "fitted")
if HERE not in sys.path:
    sys.path.insert(0, HERE)
from vg3_base_functions import (
    load_locked_params, setup_cell, load_and_classify_synapses,
    find_true_soma_segment, run_polarity_spot_g, hill, define_layer_masks,
    contrast_gain,
)

OFF_JSON = os.path.join(FITTED, "fitted_off_epsc_params.json")
ON_JSON = os.path.join(FITTED, "fitted_on_epsc_params.json")
INH_JSON = os.path.join(FITTED, "fitted_inh_params.json")
LOCKED = os.path.join(DATA, "locked_params.json")
SWC = os.path.join(DATA, "cid3.swc")
SYN = os.path.join(DATA, "sm2nrn_cid3.mat")
ONB = os.path.join(DATA, "vg3_onbip_list.txt")
OFB = os.path.join(DATA, "vg3_offbip_list.txt")

EXC_HILL_OFF_D50 = 1e9    # disable internal descending exc Hill (-> 1)
INH_HILL_OFF_D50 = 1e-9   # disable internal ascending inh Hill (-> 1)
N_EVENTS = 1              # single stimulus event, not a pulse train
GCAMP_VTH = -45.0         # mV, VGCC threshold for the simple GCaMP proxy


def load_fits():
    """Load the OFF EPSC, ON EPSC, and inhibition fits."""
    with open(OFF_JSON) as f: off = json.load(f)
    with open(ON_JSON) as f: on = json.load(f)
    with open(INH_JSON) as f: inh = json.load(f)
    return {"off": off, "on": on, "inh": inh}


def build(g_pas=None, layer_radius=4.0, overrides=None):
    """Build the cid3 cell once. g_pas overrides the locked passive leak.
    overrides: optional dict merged into the locked params before building
    (e.g. {"active_kv": True, "gkbar": 0.0012}) -- for validation/tuning runs
    that need something beyond g_pas changed without touching locked_params.json.
    """
    p = load_locked_params(LOCKED)
    if g_pas is not None:
        p["g_pas"] = float(g_pas)
    if overrides:
        p.update(overrides)
    p["n_events"] = N_EVENTS
    cell, node_coords, segments = setup_cell(SWC, p)
    syn = load_and_classify_synapses(SYN, ONB, OFB, p["cell_id"], node_coords)
    soma_idx, soma_xyz, _ = find_true_soma_segment(segments, node_coords, cell=cell)
    on_nodes = np.asarray(syn["on_nodes"], dtype=int)
    off_nodes = np.asarray(syn["off_nodes"], dtype=int)
    inh_nodes = np.asarray(syn["inh_nodes"], dtype=int)
    off_m, on_m = define_layer_masks(node_coords, on_nodes, off_nodes,
                                     radius_um=layer_radius)
    return {
        "p": p, "cell": cell, "node_coords": node_coords, "segments": segments,
        "soma_idx": soma_idx, "soma_xy": soma_xyz[:2], "rest_v": float(p["rest_v"]),
        "on_nodes": on_nodes, "off_nodes": off_nodes, "inh_nodes": inh_nodes,
        "on_layer": on_m & ~off_m, "off_layer": off_m & ~on_m, "g_pas": p["g_pas"],
    }


# === fitted size-tuning gains ===
def exc_size_gain(fits, pol, D):
    if pol.upper() == "OFF":
        f = fits["off"]
        return (hill(D, f["d50_exc_off"], f["n_exc_off"], ascending=True) *
                hill(D, f["d50_supp_off"], f["n_supp_off"], ascending=False))
    f = fits["on"]
    return (hill(D, f["d50_exc_on"], f["n_exc_on"], ascending=True) *
            hill(D, f["d50_supp_on"], f["n_supp_on"], ascending=False))


def exc_eff_weight(fits, pol, D):
    if pol.upper() == "OFF":
        return fits["off"]["exc_weight"] * exc_size_gain(fits, "OFF", D)
    return (fits["on"]["exc_weight"] * fits["on"].get("on_scale", 2.0) *
            exc_size_gain(fits, "ON", D))


def inh_eff_weight(fits, pol, D):
    inh = fits["inh"]
    A = inh["A_off"] if pol.upper() == "OFF" else inh["A_on"]
    d50 = inh["d50_inh_off"] if pol.upper() == "OFF" else inh["d50_inh_on"]
    n = inh["n_inh_off"] if pol.upper() == "OFF" else inh["n_inh_on"]
    return inh["inh_weight"] * A * hill(D, d50, n, ascending=True)


def subsample_inh(inh_nodes, keep_fraction, seed=0):
    """Randomly keep keep_fraction of the inhibitory synapses (seeded).
    1.0 = all, 0.25 = 75% removed, 0.0 = all removed."""
    inh_nodes = np.asarray(inh_nodes, dtype=int)
    if keep_fraction >= 1.0 or len(inh_nodes) == 0:
        return inh_nodes
    if keep_fraction <= 0.0:
        return np.array([], dtype=int)
    rng = np.random.default_rng(seed)
    k = max(1, int(round(len(inh_nodes) * float(keep_fraction))))
    return np.sort(rng.choice(inh_nodes, size=k, replace=False))


def inh_contrast_gain(contrast, p):
    """Contrast gain for PHASIC inhibition. Same as excitation unless p has
    'inh_c50': then a saturating Naka-Rushton c^n/(c50^n + c^n), rescaled so
    100% contrast still gives 1 (the 100%-contrast fits are unchanged)."""
    if "inh_c50" not in p:
        return contrast_gain(contrast, p)
    c = max(float(contrast), 0.0)
    c50 = float(p["inh_c50"]); n = float(p.get("inh_cn", 2.0))
    return (c ** n / (c50 ** n + c ** n)) * (c50 ** n + 1.0)


def run_spot(model, fits, D, polarity, center_xy=None, record_v=False,
             exc_scale=1.0, inh_scale=1.0, inh_keep_fraction=1.0, inh_seed=0,
             contrast=1.0):
    """Run a combined excitation+inhibition spot stimulus with fitted gains.

    exc_scale/inh_scale multiply the fitted weights. contrast (0-1) scales both
    via contrast_gain. inh_keep_fraction<1 randomly removes inhibitory synapses.
    Only the chosen polarity's bipolar cells fire; inhibition is active for both.
    Returns (peak_g_exc_nS, peak_g_inh_nS[, t, all_v]).
    """
    p = dict(model["p"])
    p["d50_psi"] = EXC_HILL_OFF_D50
    p["d50_inh_on"] = INH_HILL_OFF_D50
    p["d50_inh_off"] = INH_HILL_OFF_D50
    center = model["soma_xy"] if center_xy is None else center_xy
    cg = contrast_gain(contrast, p)
    exc_eff = exc_eff_weight(fits, polarity, float(D)) * float(exc_scale) * cg
    inh_eff = inh_eff_weight(fits, polarity, float(D)) * float(inh_scale) * inh_contrast_gain(contrast, p)
    inh_use = subsample_inh(model["inh_nodes"], inh_keep_fraction, inh_seed)
    return run_polarity_spot_g(
        model["cell"], model["segments"], model["node_coords"],
        model["on_nodes"], model["off_nodes"], inh_use, p,
        center, float(D), polarity=polarity.upper(),
        exc_weight=exc_eff, inh_weight=inh_eff, record_v=record_v)


def region_voltages(model, all_v):
    """Per-region voltage traces: soma, ON-dendrite region, OFF-dendrite region."""
    soma = all_v[model["soma_idx"]]
    on = all_v[model["on_layer"]].mean(axis=0) if model["on_layer"].any() else soma
    off = all_v[model["off_layer"]].mean(axis=0) if model["off_layer"].any() else soma
    return {"soma": soma, "ON": on, "OFF": off}


def gcamp_signal(v, vth=GCAMP_VTH):
    """Half-wave rectification: max(0, V - vth). No kinetics."""
    return np.maximum(0.0, np.asarray(v) - vth)


def peak_dv(v, t, rest, t_on):
    """Signed peak deviation from rest after stimulus onset."""
    m = np.asarray(t) >= t_on
    d = np.asarray(v)[m] - rest
    return float(d[int(np.argmax(np.abs(d)))]) if d.size else 0.0


def simulate_vgcc_calcium(delta_v, rest_v=-60.0, v_half=-40.0, k=4.0):
    """Boltzmann transform: peak depolarization (mV from rest) -> VGCC open-
    probability / Ca2+ influx proxy in [0, 1]."""
    absolute_v = rest_v + delta_v
    return 1.0 / (1.0 + np.exp(-(absolute_v - v_half) / k))


def run_combined(model, fits, D, center_xy=None, record_v=True,
                 exc_scale=1.0, inh_scale=1.0, inh_keep_fraction=1.0, inh_seed=0,
                 contrast=1.0):
    """Activate BOTH ON and OFF bipolar classes at once, plus shared inhibition
    (additive recruitment over both pathways' Hills). Returns (t, all_v)."""
    from vg3_base_functions import gaussian_spot_weights
    from neuron import h
    p = dict(model["p"]); p["n_events"] = N_EVENTS
    nc = model["node_coords"]; segs = model["segments"]
    center = model["soma_xy"] if center_xy is None else center_xy
    sigma = p.get("spot_sigma_um", 25.0)
    cg = contrast_gain(contrast, p)

    inh_use = subsample_inh(model["inh_nodes"], inh_keep_fraction, inh_seed)
    on_w = fits["on"]["exc_weight"] * fits["on"].get("on_scale", 2.0) * exc_size_gain(fits, "ON", D) * exc_scale * cg
    off_w = fits["off"]["exc_weight"] * exc_size_gain(fits, "OFF", D) * exc_scale * cg
    inh = fits["inh"]
    inh_drive = (inh["A_on"] * hill(D, inh["d50_inh_on"], inh["n_inh_on"], ascending=True)
                 + inh["A_off"] * hill(D, inh["d50_inh_off"], inh["n_inh_off"], ascending=True))
    inh_w = inh["inh_weight"] * inh_drive * inh_scale * cg

    syns = []
    for nodes, w, e, t1, t2 in [
            (model["on_nodes"], on_w, p["exc_rev"], p["exc_tau1"], p["exc_tau2"]),
            (model["off_nodes"], off_w, p["exc_rev"], p["exc_tau1"], p["exc_tau2"]),
            (inh_use, inh_w, p["inh_rev"], p["inh_tau1"], p["inh_tau2"])]:
        xy = nc[np.asarray(nodes, dtype=int), :2]
        sw = gaussian_spot_weights(xy, center, float(D), sigma)
        for idx, gw in zip(nodes, sw):
            if gw <= 1e-6:
                continue
            seg = segs[int(idx)]
            s = h.Exp2Syn(seg.x, sec=seg.sec)
            s.e = e; s.tau1 = t1; s.tau2 = t2
            syns.append((s, float(w) * float(gw)))

    stim = h.NetStim(); stim.number = p["n_events"]
    stim.start = p.get("start_time", 50.0)
    stim.interval = p.get("event_interval_ms", 10.0); stim.noise = 0
    ncons = []
    for s, w in syns:
        c = h.NetCon(stim, s); c.weight[0] = w; c.delay = 0; ncons.append(c)

    t_vec = h.Vector().record(h._ref_t)
    seg_vecs = [h.Vector().record(seg._ref_v) for seg in segs]
    h.CVode().active(0); h.dt = p["dt"]; h.finitialize(p["rest_v"])
    while h.t < p["t_stop"]:
        h.fadvance()
    t = np.array(t_vec); all_v = np.array([np.array(v) for v in seg_vecs])
    del ncons, syns, stim, t_vec, seg_vecs
    return t, all_v


def save_fig(fig, path_noext, dpi=130, formats=("png", "svg")):
    """Save a figure in the given formats. Use ('png','pdf') for editable PDF."""
    for ext in formats:
        if ext == "pdf":
            fig.savefig(f"{path_noext}.pdf", format="pdf")
        else:
            fig.savefig(f"{path_noext}.{ext}", dpi=dpi)


def off_field_corner(model, margin=5.0):
    """A corner of the OFF-synapse footprint, for placing a bar so its signal
    must spread across the cell to reach the soma."""
    xy = model["node_coords"][model["off_nodes"], :2]
    return np.array([xy[:, 0].min() - margin, xy[:, 1].min() - margin])


def run_bar(model, fits, corner_xy, bar_w=50.0, bar_h=80.0, D_eff=100.0,
            polarity="OFF", record_v=True, exc_scale=1.0, inh_scale=1.0,
            inh_keep_fraction=1.0, inh_seed=0, contrast=1.0):
    """A bar (bar_w x bar_h um) activating the chosen polarity's bipolar
    synapses inside the rectangle at corner_xy, plus matching-polarity
    inhibition over the same rectangle. Returns (t, all_v, exc_in, inh_in)."""
    from neuron import h
    pol = polarity.upper()
    p = dict(model["p"]); p["n_events"] = N_EVENTS
    nc = model["node_coords"]; segs = model["segments"]
    x0, y0 = float(corner_xy[0]), float(corner_xy[1])
    cg = contrast_gain(contrast, p)

    def in_bar(nodes):
        xy = nc[np.asarray(nodes, dtype=int), :2]
        m = ((xy[:, 0] >= x0) & (xy[:, 0] <= x0 + bar_w) &
             (xy[:, 1] >= y0) & (xy[:, 1] <= y0 + bar_h))
        return np.asarray(nodes, dtype=int)[m]

    exc_nodes = model["off_nodes"] if pol == "OFF" else model["on_nodes"]
    exc_in = in_bar(exc_nodes)
    inh_in = in_bar(subsample_inh(model["inh_nodes"], inh_keep_fraction, inh_seed))
    exc_w = exc_eff_weight(fits, pol, D_eff) * exc_scale * cg
    inh = fits["inh"]
    A = inh["A_off"] if pol == "OFF" else inh["A_on"]
    d50 = inh["d50_inh_off"] if pol == "OFF" else inh["d50_inh_on"]
    n = inh["n_inh_off"] if pol == "OFF" else inh["n_inh_on"]
    inh_w = inh["inh_weight"] * A * hill(D_eff, d50, n, ascending=True) * inh_scale * cg

    syns = []
    for nodes, w, e, t1, t2 in [
            (exc_in, exc_w, p["exc_rev"], p["exc_tau1"], p["exc_tau2"]),
            (inh_in, inh_w, p["inh_rev"], p["inh_tau1"], p["inh_tau2"])]:
        for idx in nodes:
            seg = segs[int(idx)]
            s = h.Exp2Syn(seg.x, sec=seg.sec)
            s.e = e; s.tau1 = t1; s.tau2 = t2
            syns.append((s, float(w)))

    stim = h.NetStim(); stim.number = p["n_events"]
    stim.start = p.get("start_time", 50.0)
    stim.interval = p.get("event_interval_ms", 10.0); stim.noise = 0
    ncons = []
    for s, w in syns:
        c = h.NetCon(stim, s); c.weight[0] = w; c.delay = 0; ncons.append(c)

    t_vec = h.Vector().record(h._ref_t)
    seg_vecs = [h.Vector().record(seg._ref_v) for seg in segs]
    h.CVode().active(0); h.dt = p["dt"]; h.finitialize(p["rest_v"])
    while h.t < p["t_stop"]:
        h.fadvance()
    t = np.array(t_vec); all_v = np.array([np.array(v) for v in seg_vecs])
    del ncons, syns, stim, t_vec, seg_vecs
    return t, all_v, exc_in, inh_in
