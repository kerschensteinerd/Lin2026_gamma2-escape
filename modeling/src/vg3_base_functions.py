"""
vg3_base_functions.py - morphology, synapse mapping, and spot-stimulus helpers
for the cid3 VG3 amacrine cell model. Used by run_functions.py.
"""
from neuron import h
import numpy as np
import json
import os
import scipy.io as sio
import scipy.ndimage as ndi
import matplotlib.pyplot as plt

# Load NEURON standard libraries
h.load_file('stdrun.hoc')
h.load_file('stdlib.hoc')
h.load_file('import3d.hoc')


# =============================================================================
# SWC + segment loaders (unchanged from the older pipeline)
# =============================================================================
class Cell:
    def __init__(self, name='neuron'):
        self.name = name
        self.soma, self.apic, self.dend, self.axon, self.all = [], [], [], [], []


def load_swc(filename, xshift=0, yshift=0, zshift=0):
    cell = Cell(name="".join(filename.split('.')[:-1]))
    morph = h.Import3d_SWC_read()
    morph.input(filename)
    i3d = h.Import3d_GUI(morph, 0)
    swc_secs = i3d.swc.sections
    swc_secs = [swc_secs.object(i) for i in range(int(swc_secs.count()))]
    name_form = {1: 'soma[%d]', 2: 'axon[%d]', 3: 'dend[%d]', 4: 'apic[%d]'}
    sec_list = {1: cell.soma, 2: cell.axon, 3: cell.dend, 4: cell.apic}
    real_secs = {}
    for swc_sec in swc_secs:
        cell_part = int(swc_sec.type)
        if swc_sec.is_subsidiary:
            continue
        if cell_part not in name_form:
            cell_part = 3
        sec = h.Section(cell=cell)
        if swc_sec.parentsec is not None:
            sec.connect(real_secs[swc_sec.parentsec.hname()](swc_sec.parentx))
        if swc_sec.first == 1:
            h.pt3dstyle(1, swc_sec.raw.getval(0, 0), swc_sec.raw.getval(1, 0),
                        swc_sec.raw.getval(2, 0), sec=sec)
        j = swc_sec.first
        xx, yy, zz = [swc_sec.raw.getrow(i).c(j) for i in range(3)]
        dd = swc_sec.d.c(j)
        if dd.size() == 1:
            x, y, z, d = [dim.x[0] for dim in [xx, yy, zz, dd]]
            d = max(d, 0.1)
            for xp in [x - d/2., x, x + d/2.]:
                h.pt3dadd(xp + xshift, y + yshift, z + zshift, d, sec=sec)
        else:
            for x, y, z, d in zip(xx, yy, zz, dd):
                d = max(d, 0.1)
                h.pt3dadd(x + xshift, y + yshift, z + zshift, d, sec=sec)
        sec_list[cell_part].append(sec)
        real_secs[swc_sec.hname()] = sec
        sec.nseg = max(1, int(sec.L / 1))
    cell.all = cell.soma + cell.apic + cell.dend + cell.axon
    return cell, real_secs


def extract_node_coords(cell):
    node_coords, segments = [], []
    for sec in cell.all:
        n3d = int(h.n3d(sec=sec))
        x3d = np.array([h.x3d(i, sec=sec) for i in range(n3d)])
        y3d = np.array([h.y3d(i, sec=sec) for i in range(n3d)])
        z3d = np.array([h.z3d(i, sec=sec) for i in range(n3d)])
        lengths = np.sqrt(np.diff(x3d)**2 + np.diff(y3d)**2 + np.diff(z3d)**2)
        cum = np.insert(np.cumsum(lengths), 0, 0.0)
        total = cum[-1] if len(cum) > 0 else 0
        for seg in sec:
            tl = seg.x * sec.L
            if tl <= total and total > 0:
                idx = max(0, np.searchsorted(cum, tl) - 1)
                sl = lengths[idx] if idx < len(lengths) else 0
                frac = (tl - cum[idx]) / sl if sl > 0 else 0.0
                xp = x3d[idx] + frac * (x3d[idx+1] - x3d[idx]) if idx+1 < len(x3d) else x3d[idx]
                yp = y3d[idx] + frac * (y3d[idx+1] - y3d[idx]) if idx+1 < len(y3d) else y3d[idx]
                zp = z3d[idx] + frac * (z3d[idx+1] - z3d[idx]) if idx+1 < len(z3d) else z3d[idx]
                node_coords.append([xp, yp, zp])
                segments.append(seg)
    return np.array(node_coords), segments


def map_synapses_to_nodes(syn_pos, node_coords):
    n_syn = syn_pos.shape[0]
    nearest = np.zeros(n_syn, dtype=int)
    available = np.ones(node_coords.shape[0], dtype=bool)
    for s in range(n_syn):
        d = np.sqrt(np.sum((node_coords - syn_pos[s])**2, axis=1))
        d[~available] = np.inf
        m = np.argmin(d)
        nearest[s] = m
        available[m] = False
    return nearest


# =============================================================================
# Locked-params helpers
# =============================================================================
DEFAULT_LOCKED_PARAMS = {
    "cell_id": 3,
    "g_pas": 1e-4,
    "ra": 100,
    "cm": 1,
    "rest_v": -45,          # VGluT3 RMP per Hsiang/Friedrichsen (Kim 2015 reported ~-38 but model uses -45)
    "temp": 33,
    "exc_rev": 0,
    "inh_rev": -60,         # chloride reversal at RMP → no IPSC at V_hold=-60
    "exc_weight": 0.001,    # uS (will be tuned to Kim 2015 EPSC targets)
    "inh_weight": 1.12e-4,  # uS (will be tuned to Kim 2015 IPSC targets)
    "on_scale": 0.3,
    "exc_tau1": 1.0,
    "exc_tau2": 1.0,
    "inh_tau1": 1.0,
    "inh_tau2": 3.0,
    "gaba_fraction": 0.7,
    # Voltage-clamp protocol
    "v_hold_epsc": -60.0,
    "v_hold_ipsc": 0.0,
    "spot_sigma_um": 25.0,  # Gaussian-blur SD (bipolar RF size)
    "start_time": 50.0,
    "t_stop": 300.0,
    "dt": 0.1,
    "n_events": 10,         # short release train per spot
    "event_interval_ms": 10.0
}


def load_locked_params(path):
    if not os.path.exists(path):
        return dict(DEFAULT_LOCKED_PARAMS)
    with open(path, 'r') as f:
        loaded = json.load(f)
    merged = dict(DEFAULT_LOCKED_PARAMS)
    merged.update(loaded)
    return merged


def save_locked_params(path, params):
    with open(path, 'w') as f:
        json.dump(params, f, indent=2)


# =============================================================================
# Cell + synapse loaders
# =============================================================================
def setup_cell(swc_file, p):
    """Build the VG3 cell. Passive only — graded, non-spiking (no Na channels).
    Effective leak = g_pas + g_tonic_shunt (uniform background Cl- conductance)
    with reversal at rest_v (-60). Stage 10's tonic-shunt lever uses this.
    """
    h('forall delete_section()')
    cell, _ = load_swc(swc_file)
    g_eff = float(p["g_pas"]) + float(p.get("g_tonic_shunt", 0.0))

    # Identify the true soma BEFORE any diameter flattening, by native max
    # diameter (the cid3 SWC type-1 marker is misplaced at a thin terminal).
    node_coords, segments = extract_node_coords(cell)
    native_diams = np.array([seg.diam for seg in segments])
    cell.soma_idx = int(np.argmax(native_diams))
    cell.soma_xyz = node_coords[cell.soma_idx].copy()

    # Optional uniform cable diameter (Friedrichsen modeling choice: setDiam=0.5).
    # This is a biophysical parameter, NOT a morphology change — branch points,
    # segment lengths, XYZ coords and synapse locations are all unchanged; only
    # the (EM-tracing-unreliable) cable thickness is set to a literature value.
    # Shortens the passive length constant uniformly (λ ∝ √diam).
    set_diam = p.get("set_diam", None)
    if set_diam is not None:
        sd = float(set_diam)
        for sec in cell.all:
            for i in range(int(h.n3d(sec=sec))):
                h.pt3dchange(i, sd, sec=sec)

    for sec in cell.all:
        sec.insert('pas')
        sec.g_pas = g_eff
        sec.e_pas = p["rest_v"]
        sec.Ra = p["ra"]
        sec.cm = p["cm"]
    h.celsius = p["temp"]
    node_coords, segments = extract_node_coords(cell)
    return cell, node_coords, segments


def find_true_soma_segment(segments, node_coords, cell=None):
    """Find the segment representing the actual cell body.

    The cid3 SWC has a misplaced type=1 marker pointing to a tiny terminal
    varicosity at (175.8, 145.3, 13.6) instead of the real soma at
    (135.7, 145.7, ~38.5) with diameter 6.6 µm. So `segments[0]` is wrong.

    Strategy: find the segment with the largest diameter — soma is the
    thickest part of the cell. Returns (soma_seg_index, soma_xyz, soma_diam).

    If `cell` is given and has a `soma_idx` attribute (set by setup_cell before
    any diameter flattening), use that — needed when set_diam has made all
    diameters uniform so argmax(diam) would be meaningless.
    """
    import numpy as _np
    if cell is not None and hasattr(cell, "soma_idx"):
        idx = int(cell.soma_idx)
        return idx, node_coords[idx], float(segments[idx].diam)
    diams = _np.array([seg.diam for seg in segments])
    idx_max = int(_np.argmax(diams))
    return idx_max, node_coords[idx_max], float(diams[idx_max])


def load_and_classify_synapses(mat_file, on_bip_file, off_bip_file, cell_id, node_coords):
    mat_data = sio.loadmat(mat_file, squeeze_me=True, struct_as_record=False)
    syn_data = mat_data['sm2nrn'].syn
    syn_pos = syn_data.pos
    syn_pre = syn_data.pre
    syn_type = syn_data.synType

    with open(on_bip_file, 'r') as f:
        on_bip_list = [int(x) for x in f.read().split()]
    with open(off_bip_file, 'r') as f:
        off_bip_list = [int(x) for x in f.read().split()]

    on_syn_idx, off_syn_idx, inh_syn_idx = [], [], []
    for i in range(len(syn_pre)):
        pre_id = int(syn_pre[i])
        stype = int(syn_type[i])
        if pre_id in on_bip_list:
            on_syn_idx.append(i)
        elif pre_id in off_bip_list:
            off_syn_idx.append(i)
        if stype == 1 and pre_id != cell_id:
            inh_syn_idx.append(i)

    mapped_nodes = map_synapses_to_nodes(syn_pos, node_coords)
    return {
        "syn_pos": syn_pos,
        "mapped_nodes": mapped_nodes,
        "on_nodes":  mapped_nodes[on_syn_idx],
        "off_nodes": mapped_nodes[off_syn_idx],
        "inh_nodes": mapped_nodes[inh_syn_idx],
        "on_syn_idx":  np.array(on_syn_idx),
        "off_syn_idx": np.array(off_syn_idx),
        "inh_syn_idx": np.array(inh_syn_idx),
    }


# =============================================================================
# Gaussian-blurred spot stimulus
# =============================================================================
def gaussian_spot_weights(synapse_xy, center_xy, spot_diameter_um,
                          sigma_um=25.0, pixel_um=2.0, pad_um=None, normalize=False):
    """Per-synapse activation weight in [0, 1] for a circular spot of given
    diameter, blurred by a Gaussian kernel of SD sigma_um.

    Computes the 2D convolution numerically on a regular grid, then bilinearly
    samples at each synapse XY.

    If normalize=True, min-max normalizes the blurred grid to [0, 1] before
    sampling, so that center synapses read ~1.0 at every spot diameter.
    If normalize=False (default), leaves the grid un-normalized so small spots
    give lower peak intensity than large spots.
    """
    if pad_um is None:
        pad_um = 4 * sigma_um
    syn_xy = np.asarray(synapse_xy)
    if len(syn_xy) == 0:
        return np.zeros(0)
    R = spot_diameter_um / 2.0

    xmin = min(syn_xy[:, 0].min(), center_xy[0]) - pad_um - R
    xmax = max(syn_xy[:, 0].max(), center_xy[0]) + pad_um + R
    ymin = min(syn_xy[:, 1].min(), center_xy[1]) - pad_um - R
    ymax = max(syn_xy[:, 1].max(), center_xy[1]) + pad_um + R

    nx = int(np.ceil((xmax - xmin) / pixel_um)) + 1
    ny = int(np.ceil((ymax - ymin) / pixel_um)) + 1
    xs = xmin + np.arange(nx) * pixel_um
    ys = ymin + np.arange(ny) * pixel_um
    X, Y = np.meshgrid(xs, ys)
    disk = (((X - center_xy[0])**2 + (Y - center_xy[1])**2) <= R**2).astype(float)
    blurred = ndi.gaussian_filter(disk, sigma=sigma_um / pixel_um, mode='constant', cval=0.0)
    
    # Optionally normalize the blurred grid to [0, 1] (matches th2 convention)
    if normalize:
        blurred_min = blurred.min()
        blurred_max = blurred.max()
        if blurred_max > blurred_min:
            blurred = (blurred - blurred_min) / (blurred_max - blurred_min)
        else:
            blurred = np.zeros_like(blurred)
    # else: Leave un-normalized — small spots SHOULD give smaller drive than large spots.

    fx = (syn_xy[:, 0] - xmin) / pixel_um
    fy = (syn_xy[:, 1] - ymin) / pixel_um
    ix = np.clip(fx.astype(int), 0, nx - 2)
    iy = np.clip(fy.astype(int), 0, ny - 2)
    dx = fx - ix; dy = fy - iy
    w = ((1 - dx) * (1 - dy) * blurred[iy, ix]
         + dx       * (1 - dy) * blurred[iy, ix + 1]
         + (1 - dx) * dy       * blurred[iy + 1, ix]
         + dx       * dy       * blurred[iy + 1, ix + 1])
    return np.clip(w, 0.0, None)


# =============================================================================
# Network surround model — global scaling factors applied per spot diameter.
#
# These capture two mechanisms that are NOT in the cell's anatomy:
#   1) Presynaptic inhibition of bipolar terminals: wide-field amacrines GABA-
#      ergically suppress BC glutamate release when the surround is lit. As D
#      grows, every BC fires less → S_exc(D) DECREASES.
#   2) Wide-field amacrine recruitment: the amacrines that release inhibition
#      onto VG3 have larger RFs than VG3's dendritic field. Small spots don't
#      drive them; large spots do → S_inh(D) INCREASES.
#
# Both are applied as GLOBAL multipliers on every per-synapse NetCon weight
# (on top of the per-synapse Gaussian spatial weight from gaussian_spot_weights).
# Total free parameters: 4 (D50_psi, n_psi, D50_inh, n_inh).
# =============================================================================
def surround_scale_exc(D, D50=150.0, n=4.0):
    """Descending Hill — presynaptic inhibition of BC terminals."""
    D = max(float(D), 1e-6)
    return 1.0 / (1.0 + (D / D50) ** n)


def surround_scale_inh(D, D50=250.0, n=2.0):
    """Ascending Hill — wide-field amacrine surround recruitment (legacy: pre-polarity model)."""
    D = max(float(D), 1e-6)
    return 1.0 / (1.0 + (D50 / D) ** n)


def surround_scale_inh_polarity(D, polarity, p):
    """Polarity-specific AC Hill scale.
        polarity = 'ON'  → uses d50_inh_on, n_inh_on
        polarity = 'OFF' → uses d50_inh_off, n_inh_off
    Falls back to legacy d50_inh / n_inh keys if polarity-specific not present.
    """
    D = max(float(D), 1e-6)
    if polarity.upper() == "ON":
        D50 = float(p.get("d50_inh_on", p.get("d50_inh", 250.0)))
        n   = float(p.get("n_inh_on",   p.get("n_inh",   2.0)))
    elif polarity.upper() == "OFF":
        D50 = float(p.get("d50_inh_off", p.get("d50_inh", 250.0)))
        n   = float(p.get("n_inh_off",   p.get("n_inh",   2.0)))
    else:
        raise ValueError(f"polarity must be 'ON' or 'OFF'; got {polarity!r}")
    return 1.0 / (1.0 + (D50 / D) ** n)


# =============================================================================
# Presynaptic EPSC size-tuning gain: TWO Hill curves + a contrast hook.
#
# The bipolar -> VG3 EPSC spot-size tuning is modeled as a CONTRAST-INDEPENDENT
# product of two Hill curves, applied as a single global per-spot multiplier on
# top of the anatomical recruitment R0(D) (the un-normalized Gaussian drive
# summed over the recruited synapses):
#
#     gain(D, C) = contrast_gain(C) * H_exc(D) * H_supp(D)
#     peak_g_exc(D, C) = exc_weight * gain(D, C) * R0(D)
#
#   H_exc(D)  ASCENDING Hill  — center excitation recruitment of BC release.
#             ~0 for tiny spots -> 1 for large spots. Sharpens the rise and
#             removes the small-spot overshoot left by raw recruitment alone.
#                 H_exc(D)  = 1 / (1 + (d50_exc  / D)^n_exc)
#   H_supp(D) DESCENDING Hill — surround suppression / presynaptic inhibition of
#             BC terminals. 1 for small spots -> 0 for large spots. The fall.
#                 H_supp(D) = 1 / (1 + (D / d50_supp)^n_supp)
#
#   contrast_gain(C) maps stimulus contrast (0..1) to a DRIVE multiplier
#   (equivalent to scaling the stimulus disk amplitude). The SAME H_exc / H_supp
#   curves are reused at every contrast — only this orthogonal factor changes.
#   Default is linear (returns C; 100% contrast = 1.0 -> 1.0); drop in a
#   Naka-Rushton C^q/(C50^q + C^q) later by setting contrast_c50 / contrast_n,
#   WITHOUT touching the size-tuning Hills.
#
# Parameter keys are polarity-tagged so ON and OFF carry independent curves:
#     d50_exc_off,  n_exc_off,  d50_supp_off,  n_supp_off    (and the *_on set)
# This is the generic definition Part 2 fits at 100% contrast and that the
# planned contrast sweep (20/50/100%) will reuse unchanged.
# =============================================================================
def hill(x, x50, n, ascending=False):
    """Generic Hill function returning a value in [0, 1].

    ascending=False (default): DESCENDING, 1 -> 0 as x grows: 1/(1+(x/x50)^n).
    ascending=True:            ASCENDING,  0 -> 1 as x grows: 1/(1+(x50/x)^n).
    """
    x = max(float(x), 1e-9)
    x50 = max(float(x50), 1e-9)
    if ascending:
        return 1.0 / (1.0 + (x50 / x) ** float(n))
    return 1.0 / (1.0 + (x / x50) ** float(n))


def contrast_gain(contrast, p=None):
    """Stimulus contrast (0..1) -> drive multiplier. Generic, contrast-axis hook.

    Default: linear (returns `contrast`), so 100% contrast (1.0) -> 1.0 and the
    Part-2 fit is unaffected. If `p` carries 'contrast_c50' and 'contrast_n', a
    saturating Naka-Rushton contrast-response is used instead. The size-tuning
    Hills are independent of this factor.
    """
    c = max(float(contrast), 0.0)
    if p is not None and "contrast_c50" in p and "contrast_n" in p:
        c50 = max(float(p["contrast_c50"]), 1e-9)
        q = float(p["contrast_n"])
        return c ** q / (c50 ** q + c ** q)
    return c


def presynaptic_exc_gain(D, p, polarity="OFF", contrast=1.0):
    """Per-spot excitatory gain = contrast_gain(C) * H_exc(D) * H_supp(D).

    Reads d50_exc_{pol}, n_exc_{pol}, d50_supp_{pol}, n_supp_{pol} from `p`
    (pol = 'on'/'off'). Returns the global multiplier applied to every recruited
    OFF/ON synapse's weight at spot diameter D and contrast C.
    """
    pol = polarity.lower()
    d50_exc  = float(p[f"d50_exc_{pol}"]);  n_exc  = float(p[f"n_exc_{pol}"])
    d50_supp = float(p[f"d50_supp_{pol}"]); n_supp = float(p[f"n_supp_{pol}"])
    h_exc  = hill(D, d50_exc,  n_exc,  ascending=True)
    h_supp = hill(D, d50_supp, n_supp, ascending=False)
    return contrast_gain(contrast, p) * h_exc * h_supp


# =============================================================================
# Voltage-clamp spot simulation
# =============================================================================
def run_vclamp_spot(cell, segments, node_coords,
                    exc_nodes, inh_nodes, p,
                    spot_center_xy, spot_diameter_um,
                    sigma_um=None, v_hold=None,
                    exc_weight=None, inh_weight=None,
                    on_scale=1.0,
                    exc_tau1=None, exc_tau2=None,
                    inh_tau1=None, inh_tau2=None,
                    n_events=None, interval=None,
                    start_time=None, t_stop=None, dt=None,
                    rs=0.001, record_all_traces=False):
    """Voltage-clamp the soma, present a Gaussian-blurred spot stimulus.

    Each synapse fires once per event (n_events total), with NetCon weight =
    base_weight * gaussian_spot_intensity_at_synapse_XY.

    Returns: t (ms), clamp_i_nA, soma_v_mV, all_v (or None),
             exc_w (per-syn intensity), inh_w (per-syn intensity).

    SEClamp.i sign convention: positive = electrode pushes current INTO cell.
    EPSC (depolarizing drive at V_hold=-60): clamp pulls cell down → i < 0.
    IPSC (hyperpolarizing drive at V_hold=0): clamp holds cell up → i > 0.
    """
    if sigma_um   is None: sigma_um   = p.get("spot_sigma_um", 25.0)
    if v_hold     is None: v_hold     = p.get("v_hold_epsc", -60.0)
    if exc_weight is None: exc_weight = p["exc_weight"]
    if inh_weight is None: inh_weight = p["inh_weight"]
    if exc_tau1   is None: exc_tau1   = p["exc_tau1"]
    if exc_tau2   is None: exc_tau2   = p["exc_tau2"]
    if inh_tau1   is None: inh_tau1   = p["inh_tau1"]
    if inh_tau2   is None: inh_tau2   = p["inh_tau2"]
    if n_events   is None: n_events   = p.get("n_events", 1)
    if interval   is None: interval   = p.get("event_interval_ms", 20.0)
    if start_time is None: start_time = p["start_time"]
    if t_stop     is None: t_stop     = p["t_stop"]
    if dt         is None: dt         = p["dt"]

    exc_nodes = np.asarray(exc_nodes, dtype=int)
    inh_nodes = np.asarray(inh_nodes, dtype=int)
    exc_xy = node_coords[exc_nodes, :2] if len(exc_nodes) else np.zeros((0, 2))
    inh_xy = node_coords[inh_nodes, :2] if len(inh_nodes) else np.zeros((0, 2))
    exc_w = gaussian_spot_weights(exc_xy, spot_center_xy, spot_diameter_um, sigma_um)
    inh_w = gaussian_spot_weights(inh_xy, spot_center_xy, spot_diameter_um, sigma_um)

    # Global network-surround scales (descending Hill for exc, ascending for inh)
    exc_global = surround_scale_exc(spot_diameter_um,
                                    D50=p.get("d50_psi", 150.0),
                                    n=p.get("n_psi", 4.0))
    inh_global = surround_scale_inh(spot_diameter_um,
                                    D50=p.get("d50_inh", 250.0),
                                    n=p.get("n_inh", 2.0))

    instantiated = []
    for idx, w in zip(exc_nodes, exc_w):
        if w <= 1e-6:
            continue
        seg = segments[int(idx)]
        syn = h.Exp2Syn(seg.x, sec=seg.sec)
        syn.e = p["exc_rev"]; syn.tau1 = exc_tau1; syn.tau2 = exc_tau2
        instantiated.append((syn, float(exc_weight) * float(w) * exc_global))
    for idx, w in zip(inh_nodes, inh_w):
        if w <= 1e-6:
            continue
        seg = segments[int(idx)]
        syn = h.Exp2Syn(seg.x, sec=seg.sec)
        syn.e = p["inh_rev"]; syn.tau1 = inh_tau1; syn.tau2 = inh_tau2
        instantiated.append((syn, float(inh_weight) * float(on_scale) * float(w) * inh_global))

    stim = h.NetStim()
    stim.number = n_events
    stim.start = start_time
    stim.interval = interval
    stim.noise = 0

    netcons = []
    for syn, w in instantiated:
        nc = h.NetCon(stim, syn)
        nc.weight[0] = w
        nc.delay = 0
        netcons.append(nc)

    soma_sec = cell.soma[0] if len(cell.soma) > 0 else cell.dend[0]
    seclamp = None; i_vec = None
    if v_hold is not None:
        seclamp = h.SEClamp(0.5, sec=soma_sec)
        seclamp.dur1 = t_stop
        seclamp.amp1 = float(v_hold)
        seclamp.rs   = float(rs)
        i_vec = h.Vector().record(seclamp._ref_i)        # nA

    t_vec = h.Vector().record(h._ref_t)
    v_vec = h.Vector().record(soma_sec(0.5)._ref_v)      # mV

    seg_vecs = []
    if record_all_traces:
        seg_vecs = [h.Vector().record(seg._ref_v) for seg in segments]

    h.CVode().active(0)
    h.dt = dt
    h.finitialize(p["rest_v"])
    while h.t < t_stop:
        h.fadvance()

    t_arr = np.array(t_vec)
    i_arr = np.array(i_vec) if i_vec is not None else np.zeros_like(t_arr)
    v_arr = np.array(v_vec)
    all_v = np.array([np.array(v) for v in seg_vecs]) if record_all_traces else None

    del netcons, instantiated, stim, t_vec, v_vec, seg_vecs
    if seclamp is not None:
        del seclamp, i_vec
    return t_arr, i_arr, v_arr, all_v, exc_w, inh_w, exc_global, inh_global


# =============================================================================
# Depth / layer helpers (kept for layer-trace plots)
# =============================================================================
def ipl_depth_pct(z, ipl_thickness=40.0):
    return (ipl_thickness - z) / ipl_thickness * 100.0


def define_layer_masks(node_coords, on_nodes, off_nodes, radius_um=15.0):
    """Layer masks by physical proximity to ON / OFF bipolar synapse mappings."""
    try:
        from scipy.spatial import cKDTree as _Tree
    except ImportError:
        from scipy.spatial import KDTree as _Tree
    on_idx  = np.asarray(on_nodes,  dtype=int)
    off_idx = np.asarray(off_nodes, dtype=int)
    n = len(node_coords)
    if len(on_idx) == 0 or len(off_idx) == 0:
        return np.zeros(n, bool), np.zeros(n, bool)
    on_dist  = _Tree(node_coords[on_idx]).query(node_coords)[0]
    off_dist = _Tree(node_coords[off_idx]).query(node_coords)[0]
    return (off_dist <= radius_um), (on_dist <= radius_um)


def plot_layer_traces(t_on, all_v_on, t_off, all_v_off,
                      node_coords, rest_v, off_mask, on_mask,
                      title, out_path,
                      t_max_plot=None, offset_ms=1000.0):
    """3-panel: Soma | OFF-dendrite region | ON-dendrite region.
    OFF trace x-displaced by offset_ms so amplitudes are visible side-by-side.
    """
    soma_on  = all_v_on[0]
    soma_off = all_v_off[0]
    if t_max_plot is None:
        t_max_plot = float(max(t_on[-1], t_off[-1]))
    on_keep  = t_on  <= t_max_plot
    off_keep = t_off <= t_max_plot
    t_on_c   = t_on[on_keep]
    t_off_c  = t_off[off_keep]
    t_off_s  = t_off_c + offset_ms
    v_on_data  = all_v_on[:, on_keep]
    v_off_data = all_v_off[:, off_keep]

    def avg(mask, v):
        return v[mask].mean(axis=0) if mask.any() else v[0]

    panels = [
        (soma_on[on_keep],         soma_off[off_keep],        "Soma"),
        (avg(off_mask, v_on_data), avg(off_mask, v_off_data),
         f"OFF-dendrite region (n={int(off_mask.sum())})"),
        (avg(on_mask, v_on_data),  avg(on_mask, v_off_data),
         f"ON-dendrite region (n={int(on_mask.sum())})"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.5), sharey=True)
    for ax, (v_on_tr, v_off_tr, lbl) in zip(axes, panels):
        ax.plot(t_on_c,  v_on_tr,  "r", lw=1.5, label="ON stim")
        ax.plot(t_off_s, v_off_tr, "b", lw=1.5, label="OFF stim")
        ax.axhline(rest_v, color="grey", ls=":", lw=0.8)
        ax.axvspan(t_max_plot, offset_ms, alpha=0.07, color="grey")
        ax.set_title(lbl)
        ax.set_xlabel(f"t (ms) — ON: 0–{t_max_plot:.0f}  |  "
                      f"OFF: {offset_ms:.0f}–{offset_ms+t_max_plot:.0f}")
    axes[0].set_ylabel("V (mV)"); axes[0].legend(fontsize=8)
    fig.suptitle(title, fontsize=11)
    fig.tight_layout(); fig.savefig(out_path); plt.close(fig)


# =============================================================================
# Stage 10 — Local-bar protocol helpers
# =============================================================================
def box_subset_indices(node_coords, mapped_nodes, center_xy, box_size_um,
                       sublamina=None, ipl_thickness=40.0):
    """Return subset of mapped_nodes inside an XY box (full width = box_size_um)
    centered on center_xy. If sublamina ∈ {'ON','OFF'}, also filter by IPL%.

    Convention: ON sublamina = IPL% > 50 (raw z < ipl_thickness/2).
                OFF sublamina = IPL% < 50 (raw z > ipl_thickness/2).
    """
    nodes = np.asarray(mapped_nodes, dtype=int)
    if len(nodes) == 0:
        return nodes
    half = box_size_um / 2.0
    xy = node_coords[nodes, :2]
    in_box = (np.abs(xy[:, 0] - center_xy[0]) <= half) & \
             (np.abs(xy[:, 1] - center_xy[1]) <= half)
    if sublamina is not None:
        z = node_coords[nodes, 2]
        pct = (ipl_thickness - z) / ipl_thickness * 100.0
        if sublamina.upper() == "ON":
            in_box &= pct > 50.0
        elif sublamina.upper() == "OFF":
            in_box &= pct < 50.0
    return nodes[in_box]


def select_target_varicosities(node_coords, on_nodes, off_nodes,
                                n_per_layer=2, min_dist_from_soma=40.0,
                                min_separation_um=25.0, ipl_thickness=40.0):
    """Pick distal BC-synapse positions as target varicosities, mixing ON/OFF.

    Sublamina is classified by Z-DEPTH at the synapse position (NOT by
    presynaptic BC type), matching the user's spec: ON = IPL% > 50, OFF = IPL% < 50.

    Strategy: pool all BC synapse nodes, split by z into z-ON and z-OFF pools,
    then within each pool pick the most-distal ones (sorted by Euclidean distance
    from soma) with a min_separation_um spacing rule.
    """
    soma = node_coords[0]
    all_bc = np.unique(np.concatenate([np.asarray(on_nodes, dtype=int),
                                        np.asarray(off_nodes, dtype=int)]))
    pct = (ipl_thickness - node_coords[all_bc, 2]) / ipl_thickness * 100.0
    pools = {"ON": all_bc[pct > 50.0], "OFF": all_bc[pct < 50.0]}

    out = []
    for sublamina, nodes in pools.items():
        if len(nodes) == 0:
            continue
        pts = node_coords[nodes]
        d_from_soma = np.linalg.norm(pts - soma, axis=1)
        order = np.argsort(-d_from_soma)
        chosen = []
        for k in order:
            if d_from_soma[k] < min_dist_from_soma:
                break
            cand_xy = pts[k, :2]
            too_close = False
            for c in chosen:
                if np.linalg.norm(node_coords[c['idx'], :2] - cand_xy) < min_separation_um:
                    too_close = True
                    break
            if too_close:
                continue
            chosen.append({
                "idx": int(nodes[k]),
                "xy": [float(pts[k, 0]), float(pts[k, 1])],
                "z":  float(pts[k, 2]),
                "sublamina": sublamina
            })
            if len(chosen) >= n_per_layer:
                break
        out.extend(chosen)
    return out


def find_measurement_nodes(target, node_coords, lateral_um,
                           ipl_thickness=40.0, xy_tol=10.0, r_tol=5.0):
    """For one target dict, find the segment indices for V_cross and V_lateral(r).

    Returns dict with keys 'center', 'cross', 'r{lateral_um[0]}', 'r{lateral_um[1]}'.
    Values that can't be found (no nearby segment satisfying criteria) are None.
    """
    cx, cy = target["xy"]
    cz = target["z"]
    pct_target = (ipl_thickness - cz) / ipl_thickness * 100.0
    sublamina = target["sublamina"]
    n_seg = len(node_coords)

    pct = (ipl_thickness - node_coords[:, 2]) / ipl_thickness * 100.0
    same_layer = (pct > 50.0) if sublamina == "ON" else (pct < 50.0)
    cross_layer = ~same_layer

    xy_d = np.linalg.norm(node_coords[:, :2] - np.array([cx, cy]), axis=1)

    # V_cross: closest segment with same XY (within xy_tol) but in opposite layer.
    cross_cands = np.where(cross_layer & (xy_d <= xy_tol))[0]
    cross_idx = int(cross_cands[np.argmin(xy_d[cross_cands])]) if len(cross_cands) else None

    result = {"center": int(target["idx"]), "cross": cross_idx}

    if not hasattr(lateral_um, "__iter__"):
        lateral_um = [lateral_um]
    for r in lateral_um:
        cands = np.where(same_layer & (np.abs(xy_d - r) <= r_tol))[0]
        # Exclude the target node itself
        cands = cands[cands != target["idx"]]
        if len(cands) == 0:
            result[f"r{int(r)}"] = None
        else:
            # Choose the candidate with smallest |xy_d - r|
            result[f"r{int(r)}"] = int(cands[np.argmin(np.abs(xy_d[cands] - r))])
    return result


def run_local_bar(cell, segments, node_coords,
                  exc_nodes_in_box, inh_nodes_in_box, p,
                  exc_weight=None, inh_weight=None,
                  exc_tau1=None, exc_tau2=None,
                  inh_tau1=None, inh_tau2=None,
                  n_events=None, interval=None,
                  start_time=None, t_stop=None, dt=None,
                  v_hold=None, rs=0.001):
    """Local-bar protocol: activate ALL synapses in the box at full weight,
    no Gaussian footprint, no global Hill scaling.

    v_hold:
        None  → no voltage clamp; returns (t, all_v).
        float → SEClamp at soma at this potential (mV); returns (t, all_v, i_clamp_nA).
                Sign: SEClamp.i positive = electrode pushes current INTO cell.
                So EPSCs at V_hold=-60 read negative, IPSCs at V_hold=0 read positive.
    """
    if exc_weight is None: exc_weight = p["exc_weight"]
    if inh_weight is None: inh_weight = p["inh_weight"]
    if exc_tau1   is None: exc_tau1   = p["exc_tau1"]
    if exc_tau2   is None: exc_tau2   = p["exc_tau2"]
    if inh_tau1   is None: inh_tau1   = p["inh_tau1"]
    if inh_tau2   is None: inh_tau2   = p["inh_tau2"]
    if n_events   is None: n_events   = p.get("n_events", 10)
    if interval   is None: interval   = p.get("event_interval_ms", 10.0)
    if start_time is None: start_time = p["start_time"]
    if t_stop     is None: t_stop     = p["t_stop"]
    if dt         is None: dt         = p["dt"]

    instantiated = []
    for idx in exc_nodes_in_box:
        seg = segments[int(idx)]
        syn = h.Exp2Syn(seg.x, sec=seg.sec)
        syn.e = p["exc_rev"]; syn.tau1 = exc_tau1; syn.tau2 = exc_tau2
        instantiated.append((syn, float(exc_weight)))
    for idx in inh_nodes_in_box:
        seg = segments[int(idx)]
        syn = h.Exp2Syn(seg.x, sec=seg.sec)
        syn.e = p["inh_rev"]; syn.tau1 = inh_tau1; syn.tau2 = inh_tau2
        instantiated.append((syn, float(inh_weight)))

    stim = h.NetStim()
    stim.number = n_events; stim.start = start_time
    stim.interval = interval; stim.noise = 0
    netcons = []
    for syn, w in instantiated:
        nc = h.NetCon(stim, syn); nc.weight[0] = w; nc.delay = 0
        netcons.append(nc)

    seclamp = None
    i_vec = None
    if v_hold is not None:
        soma_sec = cell.soma[0] if len(cell.soma) > 0 else cell.dend[0]
        seclamp = h.SEClamp(0.5, sec=soma_sec)
        seclamp.dur1 = t_stop; seclamp.amp1 = float(v_hold); seclamp.rs = float(rs)
        i_vec = h.Vector().record(seclamp._ref_i)

    t_vec = h.Vector().record(h._ref_t)
    seg_vecs = [h.Vector().record(seg._ref_v) for seg in segments]

    h.CVode().active(0)
    h.dt = dt
    h.finitialize(p["rest_v"])
    while h.t < t_stop:
        h.fadvance()

    t_arr = np.array(t_vec)
    all_v = np.array([np.array(v) for v in seg_vecs])
    if v_hold is None:
        del netcons, instantiated, stim, t_vec, seg_vecs
        return t_arr, all_v
    else:
        i_nA = np.array(i_vec)
        del netcons, instantiated, stim, t_vec, seg_vecs, seclamp, i_vec
        return t_arr, all_v, i_nA


# =============================================================================
# GCaMP6f post-hoc model.
# Half-wave rectifier above V_th (VGCC opening), then two single-exp lowpass
# stages: τ_Ca (cytosolic clearance) and τ_GCaMP (indicator binding kinetics).
# Standard model used in retina papers e.g. Borghuis 2013, Wei & Feller 2015.
# =============================================================================
def lowpass_exp(x, dt_ms, tau_ms):
    """Forward-Euler single-pole lowpass."""
    if tau_ms <= 0:
        return np.asarray(x, dtype=float).copy()
    alpha = dt_ms / (tau_ms + dt_ms)
    y = np.zeros_like(x, dtype=float)
    if len(x) == 0:
        return y
    y[0] = x[0]
    for i in range(1, len(x)):
        y[i] = y[i-1] + alpha * (x[i] - y[i-1])
    return y


def compute_gcamp_trace(t, V, v_th=-45.0, tau_ca=50.0, tau_gcamp=200.0):
    """Single-segment V → (Ca, GCaMP6f) traces in arbitrary units.

    t: ms (1-D); V: mV (1-D, same length).
    Returns (ca, gcamp) arrays in arb. units; baseline ≈ 0 when V ≤ v_th.
    """
    if len(t) < 2:
        return np.zeros_like(V), np.zeros_like(V)
    dt = float(np.median(np.diff(t)))
    ca_drive = np.maximum(np.asarray(V) - v_th, 0.0)
    ca = lowpass_exp(ca_drive, dt, tau_ca)
    gcamp = lowpass_exp(ca, dt, tau_gcamp)
    return ca, gcamp


def compute_gcamp_batch(t, V_all, v_th=-45.0, tau_ca=50.0, tau_gcamp=200.0):
    """Vectorised over segments. V_all shape (n_seg, n_t). Returns (ca, gcamp) same shape."""
    if V_all.shape[1] < 2:
        return np.zeros_like(V_all), np.zeros_like(V_all)
    dt = float(np.median(np.diff(t)))
    drive = np.maximum(V_all - v_th, 0.0)
    a_ca = dt / (tau_ca + dt)
    a_g  = dt / (tau_gcamp + dt)
    ca = np.zeros_like(V_all, dtype=float)
    g  = np.zeros_like(V_all, dtype=float)
    ca[:, 0] = drive[:, 0]
    g[:, 0]  = ca[:, 0]
    for j in range(1, V_all.shape[1]):
        ca[:, j] = ca[:, j-1] + a_ca * (drive[:, j] - ca[:, j-1])
        g[:, j]  = g[:, j-1]  + a_g  * (ca[:, j]    - g[:, j-1])
    return ca, g


# =============================================================================
# Polarity-specific spot stimulation + peak synaptic conductance measurement
# (Stage 13 — polarity-aware Kim 2015 fit).
#
# Records sum(Exp2Syn.g) over all active synapses → peak total g_exc and g_inh
# in nS. This bypasses the SEClamp settling artifact and matches the user's
# "fit nS directly, not pA" instruction.
# =============================================================================
def run_polarity_spot_g(cell, segments, node_coords,
                         on_nodes, off_nodes, inh_nodes, p,
                         spot_center_xy, spot_diameter_um, polarity,
                         exc_weight=None, inh_weight=None,
                         exc_tau1=None, exc_tau2=None,
                         inh_tau1=None, inh_tau2=None,
                         n_events=None, interval=None,
                         start_time=None, t_stop=None, dt=None,
                         sigma_um=None, record_v=False):
    """Polarity-specific spot. Only the chosen BC subset (ON or OFF) is activated.
    All AC synapses are activated, scaled by surround_scale_inh_polarity(polarity).
    Records sum-of-Exp2Syn-g across activated synapses → peak total conductance.

    Returns:
        peak_g_exc_nS, peak_g_inh_nS, optional (t, all_v) if record_v=True.
    """
    if sigma_um is None: sigma_um = p.get("spot_sigma_um", 25.0)
    if exc_weight is None: exc_weight = p["exc_weight"]
    if inh_weight is None: inh_weight = p["inh_weight"]
    if exc_tau1   is None: exc_tau1   = p["exc_tau1"]
    if exc_tau2   is None: exc_tau2   = p["exc_tau2"]
    if inh_tau1   is None: inh_tau1   = p["inh_tau1"]
    if inh_tau2   is None: inh_tau2   = p["inh_tau2"]
    if n_events   is None: n_events   = p.get("n_events", 10)
    if interval   is None: interval   = p.get("event_interval_ms", 10.0)
    if start_time is None: start_time = p["start_time"]
    if t_stop     is None: t_stop     = p["t_stop"]
    if dt         is None: dt         = p["dt"]

    # Pick the BC subset matching polarity
    if polarity.upper() == "ON":
        exc_subset = np.asarray(on_nodes, dtype=int)
    elif polarity.upper() == "OFF":
        exc_subset = np.asarray(off_nodes, dtype=int)
    else:
        raise ValueError(f"polarity must be 'ON' or 'OFF'; got {polarity!r}")
    inh_nodes = np.asarray(inh_nodes, dtype=int)

    exc_xy = node_coords[exc_subset, :2] if len(exc_subset) else np.zeros((0, 2))
    inh_xy = node_coords[inh_nodes,  :2] if len(inh_nodes)  else np.zeros((0, 2))
    exc_w = gaussian_spot_weights(exc_xy, spot_center_xy, spot_diameter_um, sigma_um)
    inh_w = gaussian_spot_weights(inh_xy, spot_center_xy, spot_diameter_um, sigma_um)

    # Global scales
    exc_global = surround_scale_exc(spot_diameter_um,
                                     D50=p.get("d50_psi", 150.0),
                                     n=p.get("n_psi", 4.0))
    inh_global = surround_scale_inh_polarity(spot_diameter_um, polarity, p)

    # Instantiate synapses, record per-syn g
    exc_syns_g = []   # list of h.Vector recording syn.g (µS)
    inh_syns_g = []
    instantiated = []
    for idx, w in zip(exc_subset, exc_w):
        if w <= 1e-6:
            continue
        seg = segments[int(idx)]
        syn = h.Exp2Syn(seg.x, sec=seg.sec)
        syn.e = p["exc_rev"]; syn.tau1 = exc_tau1; syn.tau2 = exc_tau2
        instantiated.append((syn, float(exc_weight) * float(w) * exc_global, "exc"))
        exc_syns_g.append(h.Vector().record(syn._ref_g))
    for idx, w in zip(inh_nodes, inh_w):
        if w <= 1e-6:
            continue
        seg = segments[int(idx)]
        syn = h.Exp2Syn(seg.x, sec=seg.sec)
        syn.e = p["inh_rev"]; syn.tau1 = inh_tau1; syn.tau2 = inh_tau2
        instantiated.append((syn, float(inh_weight) * float(w) * inh_global, "inh"))
        inh_syns_g.append(h.Vector().record(syn._ref_g))

    stim = h.NetStim(); stim.number = n_events; stim.start = start_time
    stim.interval = interval; stim.noise = 0
    netcons = []
    for syn, w, _ in instantiated:
        nc = h.NetCon(stim, syn); nc.weight[0] = w; nc.delay = 0
        netcons.append(nc)

    t_vec = h.Vector().record(h._ref_t)
    seg_vecs = [h.Vector().record(seg._ref_v) for seg in segments] if record_v else []

    h.CVode().active(0); h.dt = dt
    h.finitialize(p["rest_v"])
    while h.t < t_stop:
        h.fadvance()

    t_arr = np.array(t_vec)
    # Convert each g vector to µS array, sum across synapses, find peak
    if exc_syns_g:
        g_exc_sum = np.sum([np.array(g) for g in exc_syns_g], axis=0)  # µS
        peak_g_exc_uS = float(np.max(g_exc_sum))
    else:
        peak_g_exc_uS = 0.0
    if inh_syns_g:
        g_inh_sum = np.sum([np.array(g) for g in inh_syns_g], axis=0)
        peak_g_inh_uS = float(np.max(g_inh_sum))
    else:
        peak_g_inh_uS = 0.0

    peak_g_exc_nS = peak_g_exc_uS * 1000.0
    peak_g_inh_nS = peak_g_inh_uS * 1000.0

    all_v = np.array([np.array(v) for v in seg_vecs]) if record_v else None
    del netcons, instantiated, stim, t_vec, seg_vecs, exc_syns_g, inh_syns_g
    if record_v:
        return peak_g_exc_nS, peak_g_inh_nS, t_arr, all_v
    return peak_g_exc_nS, peak_g_inh_nS
