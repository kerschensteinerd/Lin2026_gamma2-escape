"""Passive NEURON transient with a rigorously bounded constant-voltage tail.

Only the late voltage tail is approximated. Local calcium activation, the 2-ms
left-held clearance recurrence, and the 50--3050-ms trapezoidal AUC are preserved.
"""
from pathlib import Path
import sys
import numpy as np
from scipy.special import expit
from neuron import h

ROOT = Path(__file__).resolve().parents[3]
MODEL = ROOT
sys.path.insert(0, str(MODEL / "src"))
from calcium_readout import calcium_trace

REST = -45.0
VHALF = -37.47892660844856
SLOPE = 1.3984142481316566
TAU = 50.0
DT = 0.1
SAMPLE = 2.0
TEND = 3050.0


def simulate_spot(model, exc_weights, inh_weights, contrast=1.0,
                  excitation_scale=1.0, tolerance=1e-7, block_ms=10.0,
                  skip_baseline=True, return_trace=False):
    """Return local/mean metrics with a per-curve absolute AUC error bound.

    Weights are effective uS/contact, including the actual inhibitory condition.
    No hidden GABA/glycine split is applied. All conductances are held fixed.
    The remaining voltage infinity norm is bounded using the passive backward-
    Euler contraction and the difference over the last full block, measured at
    every modeled segment. The conservative bound also covers the lag between
    the last recorded 2-ms voltage and the current integration state.
    """
    segs = model["segments"]
    nodes = model["off_nodes"]
    if len(exc_weights) != len(nodes) or len(inh_weights) != len(model["inh_nodes"]):
        raise ValueError("Wrong contact-weight dimensions")
    if any(float(x) < 0 for x in np.r_[exc_weights, inh_weights]):
        raise ValueError("Nonnegative conductances required for contraction proof")
    if not all(float(sec.cm) == 1.0 for sec in model["cell"].all):
        raise ValueError("The declared tau bound requires cm=1 uF/cm2")
    leak = min(float(seg.g_pas) for seg in segs)
    if leak <= 0:
        raise ValueError("Positive baseline leak required")
    if not all(abs(float(seg.e_pas) - REST) < 1e-12 for seg in segs):
        raise ValueError("Pure shunt at shared rest required")
    clamps = []
    for index, weight in zip(nodes, np.asarray(exc_weights) * excitation_scale * contrast):
        if weight <= 0:
            continue
        c = h.SEClamp(segs[int(index)])
        c.dur1 = 50.0; c.amp1 = REST
        c.dur2 = 3000.0; c.amp2 = 0.0
        c.dur3 = 10.0; c.amp3 = REST; c.rs = 1.0 / float(weight)
        clamps.append(c)
    for index, weight in zip(model["inh_nodes"], inh_weights):
        if weight <= 0:
            continue
        c = h.SEClamp(segs[int(index)])
        c.dur1 = 3100.0; c.amp1 = REST; c.rs = 1.0 / float(weight)
        clamps.append(c)
    tv = h.Vector().record(h._ref_t, SAMPLE)
    vv = [h.Vector().record(segs[int(i)]._ref_v, SAMPLE) for i in nodes]
    h.CVode().active(0); h.secondorder = 0; h.dt = DT
    h.finitialize(REST)
    if skip_baseline:
        h.t = 50.0
        h.frecord_init()
    previous = np.array([seg.v for seg in segs])
    nstep = int(round(block_ms / DT))
    tau_max = 0.001 / leak
    q = (1 + DT / tau_max)**(-nstep)
    bound_v = float("inf")
    auc_bound = float("inf")
    scale = 1.0
    nblocks = 0
    while h.t < TEND:
        for _ in range(nstep):
            h.fadvance()
        nblocks += 1
        current = np.array([seg.v for seg in segs])
        block_change = float(np.max(np.abs(current - previous)))
        # Includes the complete latest block, hence conservatively bounds any
        # sample from that block as well as the current integration state.
        bound_v = block_change / (1.0 - q)
        previous = current
        sampled_t = np.asarray(tv)
        sampled_v = np.array([np.asarray(v) for v in vv])
        last_v = sampled_v[:, -1]
        remaining = max(0.0, TEND - float(sampled_t[-1]))
        auc_bound = TAU * remaining * bound_v / (4.0 * SLOPE)
        # This lower estimate is sufficient for a conservative relative error
        # criterion: monotonically rising calcium cannot exceed its final AUC.
        prefix_ca = calcium_trace(sampled_v, SAMPLE, tau_ms=TAU, mode="fixed", v_half=VHALF, slope=SLOPE, rest_v=REST)
        win = sampled_t >= 50 - 1e-6
        accumulated = float(np.trapezoid(prefix_ca[:, win], x=sampled_t[win], axis=1).mean())
        scale = max(1.0, accumulated)
        if auc_bound <= tolerance * scale and float(sampled_t[-1]) >= 60:
            break
    t = np.arange(50.0 if skip_baseline else 0.0, TEND + 0.1, SAMPLE)
    recorded_t = np.asarray(tv).copy()
    recorded_v = np.array([np.asarray(v).copy() for v in vv])
    keep = recorded_t <= TEND + 1e-6
    recorded_t, recorded_v = recorded_t[keep], recorded_v[:, keep]
    if not np.allclose(recorded_t, t[:len(recorded_t)], rtol=0, atol=1e-7):
        raise AssertionError("Recording does not align with historical sample times")
    voltage = np.empty((len(nodes), len(t)))
    voltage[:, :len(recorded_t)] = recorded_v
    voltage[:, len(recorded_t):] = recorded_v[:, -1, None]
    calcium = calcium_trace(voltage, SAMPLE, tau_ms=TAU, mode="fixed", v_half=VHALF, slope=SLOPE, rest_v=REST)
    win = (t >= 50 - 1e-6) & (t <= TEND + 1e-6)
    local_auc = np.trapezoid(calcium[:, win], x=t[win], axis=1)
    local_dv = voltage.max(axis=1) - REST
    local_ca = calcium.max(axis=1)
    result = dict(mean_ca_auc=float(local_auc.mean()), mean_peak_dv=float(local_dv.mean()),
                  mean_peak_ca=float(local_ca.mean()), max_peak_dv=float(local_dv.max()),
                  min_peak_dv=float(local_dv.min()), min_voltage=float(voltage.min()),
                  max_voltage=float(voltage.max()), baseline_mV=REST,
                  local_auc=local_auc, local_peak_dv=local_dv, local_peak_ca=local_ca,
                  recorded_t_ms=recorded_t, recorded_v_mV=recorded_v,
                  end_integration_ms=float(h.t), tail_start_ms=float(recorded_t[-1]),
                  tail_voltage_error_bound_mV=float(bound_v), tail_auc_error_bound=float(auc_bound),
                  requested_relative_tolerance=tolerance, bound_scale=scale,
                  tau_max_ms=tau_max, contraction_per_block=q, nblocks=nblocks,
                  activation_at_rest=float(expit((REST - VHALF) / SLOPE)),
                  max_activation=float(expit((float(voltage.max()) - VHALF) / SLOPE)))
    if return_trace:
        result.update(t_ms=t, v_mV=voltage, calcium=calcium)
    # Point processes must not survive into the next paired simulation.
    return result
