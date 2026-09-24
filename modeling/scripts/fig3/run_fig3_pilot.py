#!/usr/bin/env python3
"""Shared model-loading, conductance simulation, and numerical export helpers."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import sys
import time

import numpy as np
from neuron import h

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'src'))
from run_functions import build, load_fits, load_named_config, exc_eff_weight, inh_eff_weight
from vg3_base_functions import gaussian_spot_weights
from calcium_readout import calcium_trace, influx_proxy

def write_rows(path, rows):
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)


def path_distance(model, source):
    h.distance(0, model['segments'][int(source)].x, sec=model['segments'][int(source)].sec)
    return np.array([h.distance(s.x, sec=s.sec) for s in model['segments']])


def simulate(model, exc_nodes, exc_weights, inh_weights, record_nodes,
             gaba_fraction=.7, gaba_remaining=1., tonic=False, e_inh=-60.,
             current_node=None, current_nA=.001, tstop=650., dt=.1, sample_dt=.5):
    """GABA and glycine share inherited kinetics/location; their weights stay explicit.

    Tonic protocol holds the fitted peak inhibitory load fixed at E_inh. It is
    an isolation experiment, not a simulation of the light stimulus time course.
    """
    p = model['p']; syns = []; ncons = []; clamps = []
    stim = h.NetStim(); stim.number = 1; stim.start = 50.; stim.noise = 0
    for idx, weight in zip(exc_nodes, exc_weights):
        if weight <= 0:
            continue
        s = h.Exp2Syn(model['segments'][int(idx)])
        s.e = p['exc_rev']; s.tau1 = p['exc_tau1']; s.tau2 = p['exc_tau2']
        nc = h.NetCon(stim, s); nc.weight[0] = float(weight); nc.delay = 0
        syns.append(s); ncons.append(nc)
    for fraction in [gaba_fraction * gaba_remaining, 1 - gaba_fraction]:
        for idx, total_weight in zip(model['inh_nodes'], inh_weights):
            weight = float(total_weight * fraction)
            if weight <= 1e-12:
                continue
            if tonic:
                s = h.SEClamp(model['segments'][int(idx)])
                s.dur1 = tstop + 1; s.amp1 = e_inh; s.rs = 1. / weight
                clamps.append(s)
            else:
                s = h.Exp2Syn(model['segments'][int(idx)])
                s.e = e_inh; s.tau1 = p['inh_tau1']; s.tau2 = p['inh_tau2']
                nc = h.NetCon(stim, s); nc.weight[0] = weight; nc.delay = 0
                syns.append(s); ncons.append(nc)
    if current_node is not None:
        s = h.IClamp(model['segments'][int(current_node)])
        s.delay = 50.; s.dur = 100.; s.amp = current_nA; clamps.append(s)
    tv = h.Vector().record(h._ref_t, sample_dt)
    vv = [h.Vector().record(model['segments'][int(idx)]._ref_v, sample_dt) for idx in record_nodes]
    h.CVode().active(0); h.dt = dt; h.finitialize(p['rest_v'])
    h.continuerun(tstop)
    return np.asarray(tv).copy(), np.array([np.asarray(v) for v in vv])
