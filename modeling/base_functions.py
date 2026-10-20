"""
base_functions.py - Core functions for the VG3 amacrine cell (cid3) model

This module contains functions for:
- Building the cell model (passive leak, sustained Kv, Exp2Syn excitatory and
  inhibitory synapses at the EM contacts, sigmoid contrast gains) from src/ and data/
- Running the 100 um OFF spot and returning per-segment peak voltages
- Running a single focal excitatory synapse together with the phasic spot inhibition
- Exporting the morphology geometry used for the calcium maps
- Small helpers for writing JSON and TSV files

Every input is an Exp2Syn conductance synapse; no voltage clamps are used to
deliver excitation or inhibition. The model code itself lives in src/.
Parameter values are in data/locked_params.json and src/vg3_base_functions.py,
except the synaptic scales G_EXC / G_INH below.

Author: Pratyush Ramakrishna
Date: 2026
"""

from __future__ import division
import os
import sys
import csv
import json
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SRC_DIR = os.path.join(SCRIPT_DIR, 'src')
DATA_DIR = os.path.join(SCRIPT_DIR, 'data')
if SRC_DIR not in sys.path:
    sys.path.insert(0, SRC_DIR)

from neuron import h                                               # noqa: E402
from run_functions import (build, load_fits, run_spot, exc_eff_weight,  # noqa: E402
                           inh_eff_weight, inh_contrast_gain)
from vg3_base_functions import gaussian_spot_weights, contrast_gain    # noqa: E402

h.load_file('stdrun.hoc')

# =============================================================================
# FROZEN MODEL CONSTANTS
# =============================================================================

G_EXC, G_INH = 400.0, 600.0    # scales on the fitted weights; voltage fit to Kim 2015 Fig. 1C
INH_REV = -60.0                # mV, Cl- reversal (Kim 2015)
EXC_REV = 0.0                  # mV
REST = -45.0                   # mV
SPOT_DIAMETER = 100.0          # um
SPOT_SIGMA = 25.0              # um, Gaussian blur
POLARITY = 'OFF'
START_MS = 100.0               # spot onset (ms)
T_STOP_MS = 200.0              # peak is reached well inside 100 ms of onset
FOCAL_INH_CONTRAST = 0.30      # focal-spread panels use the spot inhibition at this contrast

# =============================================================================
# FILE HELPERS
# =============================================================================

def load_json(path):
    """Read a JSON file and return its contents."""
    with open(path, 'r') as f:
        return json.load(f)


def write_json(path, value):
    """
    Write a JSON file. NumPy values are converted and NaN becomes null.

    Parameters
    ----------
    path : str
        Output file
    value : dict or list
        Data to write
    """
    def clean(x):
        if isinstance(x, dict):
            return {str(k): clean(v) for k, v in x.items()}
        if isinstance(x, (list, tuple)):
            return [clean(v) for v in x]
        if isinstance(x, np.ndarray):
            return clean(x.tolist())
        if isinstance(x, np.integer):
            return int(x)
        if isinstance(x, np.bool_):
            return bool(x)
        if isinstance(x, (float, np.floating)):
            return float(x) if np.isfinite(x) else None
        return x

    with open(path, 'w') as f:
        f.write(json.dumps(clean(value), indent=2, allow_nan=False) + '\n')


def write_rows(path, rows):
    """Write a list of dictionaries as a tab-separated table."""
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


# =============================================================================
# MODEL
# =============================================================================

def build_model():
    """
    Build the frozen physiological cid3 model.

    Returns
    -------
    model : dict
        Model from run_functions.build, with the stimulus settings applied to
        model['p'] (inh_rev, start_time, t_stop)
    fits : dict
        Fitted OFF/ON excitation and inhibition size-tuning parameters
    """
    model = build()
    model['p'].update(inh_rev=INH_REV, start_time=START_MS, t_stop=T_STOP_MS)
    return model, load_fits()


def path_distance(model, source):
    """
    Path distance along the dendrites from one segment to every segment.

    Parameters
    ----------
    model : dict
        Model from build_model
    source : int
        Index of the source segment

    Returns
    -------
    ndarray
        Distance to each segment (um)
    """
    seg = model['segments'][int(source)]
    h.distance(0, seg.x, sec=seg.sec)
    return np.array([h.distance(s.x, sec=s.sec) for s in model['segments']])


# =============================================================================
# SIMULATION
# =============================================================================

def spot_response(model, fits, contrast, retained_fraction, record='off', trace_dt=1.0):
    """
    Simulate the 100 um OFF spot and return the peak response of each segment.

    Excitation and phasic inhibition are Exp2Syn conductances scaled by their
    contrast gains. Reduced inhibition scales every inhibitory contact by
    retained_fraction (no synapses are removed).

    Parameters
    ----------
    model, fits : from build_model
    contrast : float
        Stimulus contrast, between 0 and 1
    retained_fraction : float
        Fraction of the control inhibitory conductance kept at every contact
    record : str
        'off' to return the OFF contacts only, 'all' for every segment
    trace_dt : float
        Sampling interval of the returned mean OFF-contact trace (ms)

    Returns
    -------
    dict
        'peak_v' absolute peak voltage, 'base_v' pre-stimulus voltage (per
        returned segment), 'trace_t' and 'trace_off_mean' (mean OFF-contact
        voltage), and the summed peak conductances 'g_exc_nS', 'g_inh_nS'
    """
    g_exc, g_inh, t, all_v = run_spot(model, fits, SPOT_DIAMETER, POLARITY, record_v=True,
                                      exc_scale=G_EXC, inh_scale=G_INH * retained_fraction,
                                      contrast=contrast)
    off = np.asarray(model['off_nodes'], dtype=int)
    pre = (t >= START_MS - 10) & (t < START_MS)
    post = t >= START_MS
    rows = slice(None) if record == 'all' else off
    v = all_v[rows]
    step = max(1, int(round(trace_dt / (t[1] - t[0]))))
    out = dict(peak_v=v[:, post].max(axis=1), base_v=v[:, pre].mean(axis=1),
               trace_t=t[::step].copy(), trace_off_mean=all_v[off][:, ::step].mean(axis=0),
               g_exc_nS=float(g_exc), g_inh_nS=float(g_inh))
    del all_v
    return out


def spot_inhibition_weights(model, fits, contrast=FOCAL_INH_CONTRAST):
    """
    Synaptic weight of each inhibitory contact for the 100 um OFF spot at a given contrast.

    Same weights as the spot simulation: fitted size gain x G_INH x inhibitory
    contrast gain x Gaussian spot intensity. Contacts with negligible spot
    intensity are dropped, as in run_polarity_spot_g.

    Returns
    -------
    nodes : ndarray
        Segment index of each contact
    weights : ndarray
        Peak Exp2Syn conductance of each contact in control (uS)
    """
    nodes = np.asarray(model['inh_nodes'], dtype=int)
    xy = model['node_coords'][nodes, :2]
    spot = gaussian_spot_weights(xy, model['soma_xy'], SPOT_DIAMETER, SPOT_SIGMA)
    w = inh_eff_weight(fits, POLARITY, SPOT_DIAMETER) * G_INH * inh_contrast_gain(contrast, model['p']) * spot
    keep = spot > 1e-6
    return nodes[keep], w[keep]


def focal_response(model, inh_nodes, inh_weights, retained_fraction, focal_segment, focal_g,
                   record='all', trace_sites=None, tstop=T_STOP_MS, trace_dt=1.0):
    """
    A single focal excitatory synapse together with the phasic spot inhibition.

    Every input is an Exp2Syn conductance synapse; there are no clamps. The
    inhibitory synapses sit at the EM inhibitory contacts with the spot weights
    (scaled by retained_fraction for reduced inhibition), reversal and kinetics
    from model['p'] (inh_rev, inh_tau1, inh_tau2). The focal input is one
    excitatory synapse at focal_segment (exc_rev, exc_tau1, exc_tau2). All fire
    together at the stimulus onset (START_MS).

    Parameters
    ----------
    model : dict
        Model from build_model
    inh_nodes, inh_weights : ndarray
        From spot_inhibition_weights
    retained_fraction : float
        Fraction of the inhibitory synaptic weight kept at every contact
    focal_segment : int
        Segment index of the focal excitatory synapse
    focal_g : float
        Peak conductance of the focal synapse (uS)
    record : str or iterable of int
        'all' for every segment, or the segment indices whose peaks are returned
    trace_sites : iterable of int, optional
        Segments whose voltage traces are returned (sampled every trace_dt ms)
    tstop : float
        Simulation length (ms)

    Returns
    -------
    dict
        'peak_v' and 'base_v' (absolute peak after onset and pre-stimulus mean)
        for the recorded segments, and 't', 'traces' for trace_sites
    """
    p = model['p']
    segments = model['segments']
    objects = []
    stim = h.NetStim()
    stim.number, stim.start, stim.noise = 1, START_MS, 0
    objects.append(stim)

    def add_synapse(index, weight, rev, tau1, tau2):
        syn = h.Exp2Syn(segments[int(index)])
        syn.e, syn.tau1, syn.tau2 = rev, tau1, tau2
        nc = h.NetCon(stim, syn)
        nc.weight[0], nc.delay = float(weight), 0
        objects.extend([syn, nc])

    for index, weight in zip(inh_nodes, inh_weights):
        if weight * retained_fraction > 0:
            add_synapse(index, weight * retained_fraction, p['inh_rev'], p['inh_tau1'], p['inh_tau2'])
    add_synapse(focal_segment, focal_g, p['exc_rev'], p['exc_tau1'], p['exc_tau2'])

    rows = np.arange(len(segments)) if isinstance(record, str) else np.asarray(record, dtype=int)
    t_vec = h.Vector().record(h._ref_t)
    v_vecs = [h.Vector().record(segments[int(i)]._ref_v) for i in rows]
    h.CVode().active(0)
    h.dt = p['dt']
    h.finitialize(REST)
    h.continuerun(tstop)

    t = np.asarray(t_vec)
    v = np.array([np.asarray(x) for x in v_vecs])
    pre, post = (t >= START_MS - 10) & (t < START_MS), t >= START_MS
    out = dict(peak_v=v[:, post].max(axis=1), base_v=v[:, pre].mean(axis=1))
    if trace_sites is not None:
        step = max(1, int(round(trace_dt / (t[1] - t[0]))))
        where = [int(np.flatnonzero(rows == s)[0]) for s in trace_sites]
        out.update(t=t[::step].copy(), traces=v[where][:, ::step])
    del objects, v_vecs, v
    return out


# =============================================================================
# GEOMETRY EXPORT
# =============================================================================

def morphology_geometry(model):
    """
    Export the section polylines, split at segment borders, for plotting.

    Every edge lies on a consecutive 3D point pair of a section, and edges are
    cut at the segment boundaries so each edge belongs to exactly one segment.

    Parameters
    ----------
    model : dict
        Model from build_model

    Returns
    -------
    dict
        Arrays describing the geometry: edges, edge-to-segment index, segment
        coordinates, cable length of each segment, soma position, stimulus
        center and the OFF and inhibitory segment indices
    """
    segments = model['segments']
    sections = list(model['cell'].all)
    section_index = {sec.name(): i for i, sec in enumerate(sections)}
    segment_lookup = {
        (seg.sec.name(), min(seg.sec.nseg - 1, int(seg.x * seg.sec.nseg))): i
        for i, seg in enumerate(segments)
    }

    edges, edge_nodes, edge_sections = [], [], []
    for section_i, sec in enumerate(sections):
        n3d = int(h.n3d(sec=sec))
        xyz = np.array([
            [h.x3d(i, sec=sec), h.y3d(i, sec=sec), h.z3d(i, sec=sec)]
            for i in range(n3d)
        ])
        if len(xyz) < 2 or sec.L <= 0:
            continue
        arc = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(xyz, axis=0), axis=1))]
        boundaries = np.linspace(0.0, arc[-1], sec.nseg + 1)

        for j in range(n3d - 1):
            start, end = arc[j], arc[j + 1]
            if end <= start:
                continue
            inside = boundaries[(boundaries > start) & (boundaries < end)]
            cuts = np.r_[start, inside, end]
            for left, right in zip(cuts[:-1], cuts[1:]):
                if right - left <= 1e-12:
                    continue
                endpoints = xyz[j] + ((np.array([left, right]) - start) / (end - start))[:, None] * (xyz[j + 1] - xyz[j])
                local_segment = min(sec.nseg - 1, int(((left + right) / 2) / arc[-1] * sec.nseg))
                key = (sec.name(), local_segment)
                if key not in segment_lookup:
                    continue
                edges.append(endpoints)
                edge_nodes.append(segment_lookup[key])
                edge_sections.append(section_i)

    return {
        'xyz_segments': np.asarray(model['node_coords'], dtype=float),
        'geometry_edges': np.asarray(edges, dtype=float),
        'edge_segment_index': np.asarray(edge_nodes, dtype=np.int32),
        'edge_section_index': np.asarray(edge_sections, dtype=np.int32),
        'cable_lengths': np.array([seg.sec.L / seg.sec.nseg for seg in segments]),
        'segment_section_index': np.array([section_index[seg.sec.name()] for seg in segments], dtype=np.int32),
        'segment_x': np.array([seg.x for seg in segments]),
        'section_nseg': np.array([sec.nseg for sec in sections], dtype=np.int32),
        'soma_xyz': np.asarray(model['node_coords'][model['soma_idx']], dtype=float),
        'soma_segment_index': np.asarray(model['soma_idx'], dtype=np.int32),
        'stimulus_center_xy': np.asarray(model['soma_xy'], dtype=float),
        'off_nodes': np.asarray(model['off_nodes'], dtype=np.int32),
        'inh_nodes': np.asarray(model['inh_nodes'], dtype=np.int32),
    }
