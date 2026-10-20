"""
01_run_fig4_simulations.py - Simulate the Figure 4 model with reduced inhibition

This script compares the control cell with cells that keep 20%, 25% and 30% of
the control inhibitory conductance at every inhibitory contact. The 25% case is
the nominal reduced (Gabrg2 cKO) condition. It:
- Simulates a 100 um OFF spot at 11 contrasts (0 to 70%) for control and each retained fraction
- Fits the calcium readout midpoint and slope to the OFF imaging
- Computes how far a single focal excitatory synapse spreads along the dendrites (peak
  voltage and calcium) together with the phasic spot inhibition at 30% contrast
- Converts voltage to the calcium proxy and saves the maps, curves and tables

Geometry, excitation, contrast gains and the readout are the same in every
condition. Only the inhibitory synaptic weight changes. Every input is an Exp2Syn
conductance synapse; no voltage clamps are used.

Prerequisites:
    - Files in data/ and src/, and a native build of mod/khh_vg3.mod (see README.md)

Usage:
    python 01_run_fig4_simulations.py

Output:
    - output/fig4/nominal25/ : maps, voltage traces and curves for the Figure 4 panels
    - output/fig4/ : calibration, OFF-contact peak voltages, sensitivity tables

Author: Pratyush Ramakrishna
Date: 2026
"""

from __future__ import division
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
import numpy as np

from base_functions import (
    load_json, write_json, write_rows, build_model, path_distance, morphology_geometry,
    spot_response, spot_inhibition_weights, focal_response,
    G_EXC, G_INH, INH_REV, REST, SPOT_DIAMETER, FOCAL_INH_CONTRAST, SCRIPT_DIR, DATA_DIR,
)
from calcium_functions import calcium_peak, normalize_15_70, fit_readout

# =============================================================================
# CONFIGURATION
# =============================================================================

# Retained fraction of the control inhibitory conductance (1.0 is control)
RETAINED_FRACTIONS = (0.20, 0.25, 0.30)
NOMINAL_FRACTION = 0.25
ALL_FRACTIONS = (1.0, *RETAINED_FRACTIONS)

# Contrasts (%)
MAP_CONTRASTS = (15, 30, 70)       # saved with full maps of every segment
ALL_CONTRASTS = (0, 5, 10, 15, 20, 25, 30, 40, 50, 60, 70)

N_WORKERS = 5                      # plus the main process for the focal-input search

# Directory paths
OUTPUT_DIR = os.path.join(SCRIPT_DIR, 'output', 'fig4')
NOMINAL_DIR = os.path.join(OUTPUT_DIR, 'nominal25')
FOCAL_FILE = os.path.join(DATA_DIR, 'published', 'focal_source.json')

# =============================================================================
# FUNCTIONS
# =============================================================================

_MODEL = None


def _init_worker():
    """Build the model once in each worker process."""
    global _MODEL
    _MODEL = build_model()


def _spot_job(args):
    """Run one spot simulation in a worker. Returns peak/baseline voltages and the mean OFF trace."""
    fraction, contrast, keep_all = args
    model, fits = _MODEL
    r = spot_response(model, fits, contrast / 100, fraction, record='all' if keep_all else 'off')
    return dict(fraction=fraction, contrast=contrast, keep_all=keep_all,
                peak_v=r['peak_v'].astype(np.float32), base_v=r['base_v'].astype(np.float32),
                trace_t=r['trace_t'], trace=r['trace_off_mean'].astype(np.float32),
                g_exc_nS=r['g_exc_nS'], g_inh_nS=r['g_inh_nS'])


def _focal_job(args):
    """Focal synapse plus phasic spot inhibition in a worker: peak and baseline voltage of every segment."""
    fraction, inh_nodes, inh_weights, source, focal_g = args
    model, _ = _MODEL
    r = focal_response(model, inh_nodes, inh_weights, fraction, source, focal_g)
    return fraction, r['base_v'], r['peak_v']


def spatial_metrics(y, distance, length):
    """
    Average a response in 5 um bins of path distance and find where it falls to 1/e.

    Bins are weighted by cable length. The crossing is the first point where the
    binned profile drops to 1/e of its value at the source, found by linear
    interpolation.

    Parameters
    ----------
    y : ndarray
        Response at each segment, normalized to the source
    distance : ndarray
        Path distance from the source to each segment (um)
    length : ndarray
        Cable length of each segment (um)

    Returns
    -------
    dict
        'path_1e_crossing_um' (None if there is no crossing), 'weighted_cable_um',
        'distance_um' and 'normalized_response' (the binned profile)
    """
    x, b = [0.0], [1.0]
    for low in np.arange(0, distance.max() + 5, 5):
        use = (distance >= low) & (distance < low + 5)
        if use.any():
            x.append(float(np.average(distance[use], weights=length[use])))
            b.append(float(np.average(y[use], weights=length[use])))

    xx, yy = np.asarray(x), np.asarray(b)
    crossing = np.flatnonzero((yy[:-1] > 1 / np.e) & (yy[1:] <= 1 / np.e))
    cross = None
    if len(crossing):
        i = crossing[0]
        cross = float(xx[i] + (1 / np.e - yy[i]) * (xx[i + 1] - xx[i]) / (yy[i + 1] - yy[i]))

    return dict(path_1e_crossing_um=cross, weighted_cable_um=float(np.sum(length * y)),
                distance_um=x, normalized_response=b)


def find_focal_conductance(model, inh_nodes, inh_weights, source, target_dv, g0):
    """
    Peak conductance of the focal synapse that gives a target_dv peak depolarization
    at the source in control, with the phasic spot inhibition present (secant in log g).

    Returns
    -------
    g : float
        Focal synaptic peak conductance (uS)
    """
    def source_dv(g):
        r = focal_response(model, inh_nodes, inh_weights, 1.0, source, g, record=[source])
        return float(r['peak_v'][0] - r['base_v'][0])

    x0, y0 = np.log(g0), source_dv(g0)
    x1 = x0 + 0.3
    y1 = source_dv(np.exp(x1))
    for _ in range(10):
        if abs(y1 - target_dv) < 0.01 or y1 == y0:
            break
        x0, y0, x1 = x1, y1, x1 + (target_dv - y1) * (x1 - x0) / (y1 - y0)
        y1 = source_dv(np.exp(x1))
        print(f"      focal g = {np.exp(x1) * 1000:.3f} nS -> source dV = {y1:.3f} mV")
    return float(np.exp(x1))


# =============================================================================
# MAIN SCRIPT
# =============================================================================

def main():
    print("=" * 60)
    print("VG3 Amacrine Cell Model - Figure 4 Simulations")
    print("=" * 60)
    started = time.time()

    if os.path.exists(os.path.join(OUTPUT_DIR, 'sensitivity_summary.json')):
        print(f"\nResults already exist in {OUTPUT_DIR}")
        print("  Delete the folder to run the simulations again")
        sys.exit(0)
    os.makedirs(NOMINAL_DIR, exist_ok=True)

    # Build the cell, the geometry and the focal source
    print("\n[1/6] Building the cell model")
    model, fits = build_model()
    geometry = morphology_geometry(model)
    coords = model['node_coords']
    off = np.asarray(model['off_nodes'], dtype=int)
    n_segments = len(coords)
    focal = load_json(FOCAL_FILE)
    source = int(np.argmin(np.linalg.norm(coords - np.asarray(focal['source_xyz_um']), axis=1)))
    nearest_off = float(np.min(np.linalg.norm(coords[off] - coords[source], axis=1)))
    print(f"      {n_segments} segments, {len(off)} OFF and {len(model['inh_nodes'])} inhibitory contacts")
    print(f"      Focal source: segment {source} at {np.round(coords[source], 2)} "
          f"({nearest_off:.2f} um from the nearest OFF contact)")

    # Spot simulations in worker processes
    print(f"\n[2/6] Simulating the spot at each contrast ({N_WORKERS} workers)")
    jobs = [(f, c, f in (1.0, NOMINAL_FRACTION) and c in MAP_CONTRASTS) for f in ALL_FRACTIONS for c in ALL_CONTRASTS]
    inh_nodes, inh_weights = spot_inhibition_weights(model, fits, FOCAL_INH_CONTRAST)
    print(f"      Phasic spot inhibition for the focal runs: {inh_weights.sum() * 1000:.3f} nS peak "
          f"({FOCAL_INH_CONTRAST:.0%} contrast spot, Exp2Syn, E_inh {INH_REV:g} mV)")
    with ProcessPoolExecutor(max_workers=N_WORKERS, initializer=_init_worker) as pool:
        spot_futures = [pool.submit(_spot_job, j) for j in jobs]

        # Focal conductance in the main process while the spots run
        print("\n[3/6] Focal synapse: peak conductance for a 15 mV source depolarization in control")
        # Start the conductance search from the 12.8 nS input estimate
        focal_g = find_focal_conductance(model, inh_nodes, inh_weights, source,
                                         focal['source_dv_target_mV'], focal['finite_sourceweight_uS'])
        print(f"      Focal synapse peak conductance: {focal_g * 1000:.3f} nS")
        focal_futures = [pool.submit(_focal_job, (f, inh_nodes, inh_weights, source, focal_g))
                         for f in ALL_FRACTIONS]

        spots = [fu.result() for fu in spot_futures]
        focal_runs = {f: (vb, vp) for f, vb, vp in (fu.result() for fu in focal_futures)}

    # Organize the OFF-contact peaks: (fraction, contrast, contact)
    peaks = np.zeros((len(ALL_FRACTIONS), len(ALL_CONTRASTS), len(off)))
    bases = np.zeros_like(peaks)
    full_maps, traces, conductance = {}, {}, {}
    for s in spots:
        i, j = ALL_FRACTIONS.index(s['fraction']), ALL_CONTRASTS.index(s['contrast'])
        pv, bv = s['peak_v'].astype(float), s['base_v'].astype(float)
        if s['keep_all']:
            full_maps[(s['fraction'], s['contrast'])] = (pv, bv)
            pv, bv = pv[off], bv[off]
        peaks[i, j], bases[i, j] = pv, bv
        traces[(s['fraction'], s['contrast'])] = (s['trace_t'], s['trace'])
        conductance[(s['fraction'], s['contrast'])] = (s['g_exc_nS'], s['g_inh_nS'])

    # Readout midpoint fitted to the OFF imaging (control and nominal 25%)
    print("\n[4/6] Fitting the calcium readout midpoint and slope")
    ic, ir = ALL_FRACTIONS.index(1.0), ALL_FRACTIONS.index(NOMINAL_FRACTION)
    targets = focal['imaging_OFF_normalized30_targets']
    v_half, slope, fit_error = fit_readout({'ctrl': peaks[ic], 'reduced': peaks[ir]},
                                           {'ctrl': bases[ic], 'reduced': bases[ir]},
                                           ALL_CONTRASTS, targets)
    calcium = dict(v_half_mV=v_half, slope_mV=slope, fit_error=fit_error, mode='peak_boltzmann')
    print(f"      Midpoint {v_half:.2f} mV, slope {slope:.1f} mV (fit error {fit_error:.4f})")

    # Contrast-response rows
    rows = []
    for i, f in enumerate(ALL_FRACTIONS):
        ca = calcium_peak(peaks[i], bases[i], v_half, slope).mean(axis=1)
        norm = normalize_15_70(ca, ALL_CONTRASTS)
        for j, c in enumerate(ALL_CONTRASTS):
            rows.append(dict(retained_effective_inhibition_fraction=f,
                             condition='ctrl' if f == 1.0 else 'reduced', polarity='OFF', contrast_pct=c,
                             n_model_sites=len(off),
                             mean_peak_dv_mV=float((peaks[i, j] - bases[i, j]).mean()),
                             mean_peak_ca_au=float(ca[j]), normalized_peak_ca=float(norm[j]),
                             mean_peak_absolute_v_mV=float(peaks[i, j].mean()),
                             total_g_exc_nS=conductance[(f, c)][0], total_g_inh_nS=conductance[(f, c)][1]))
            print(f"      Retained {f:.2f}, contrast {c:2d}%: mean peak dV = {rows[-1]['mean_peak_dv_mV']:.3f} mV, "
                  f"calcium = {ca[j]:.4f}, normalized = {norm[j]:.3f}")

    # Focal spread
    print("\n[5/6] Focal spread of calcium and voltage")
    distance = path_distance(model, source)
    length = geometry['cable_lengths']
    spatial_rows, spatial_arrays = [], {}
    for f in ALL_FRACTIONS:
        v_base, v_focal = focal_runs[f]          # baseline and peak voltage of every segment
        dv = v_focal - v_base
        ca = calcium_peak(v_focal, v_base, v_half, slope)
        row = dict(retained_effective_inhibition_fraction=f, source_conductance_uS=focal_g,
                   source_dv_mV=float(dv[source]), source_calcium_au=float(ca[source]),
                   source_baseline_mV=float(v_base[source]),
                   voltage=spatial_metrics(dv / dv[source], distance, length),
                   calcium=spatial_metrics(ca / ca[source], distance, length))
        spatial_rows.append(row)
        spatial_arrays[f] = dict(dv=dv, ca=ca, v_base=v_base, v_focal=v_focal,
                                 v_normalized=dv / dv[source], ca_normalized=ca / ca[source])
        print(f"      Retained {f:.2f}: source dV {dv[source]:.2f} mV, voltage 1/e = "
              f"{row['voltage']['path_1e_crossing_um']}, calcium 1/e = {row['calcium']['path_1e_crossing_um']}")

    # Save the nominal results used for the figure
    print("\n[6/6] Saving results")
    focal_maps = {**geometry, 'path_distance_from_source_um': distance,
                  'local_source_segment': np.asarray(source), 'source_segment': np.asarray(source),
                  'local_source_xyz': coords[source]}
    for condition, f in [('ctrl', 1.0), ('reduced', NOMINAL_FRACTION)]:
        a = spatial_arrays[f]
        focal_maps.update({f'local_{condition}_ca': a['ca'], f'local_{condition}_finite_dv': a['dv'],
                           f'local_{condition}_finite_v_normalized': a['v_normalized'],
                           f'local_{condition}_finite_ca_normalized': a['ca_normalized'],
                           f'local_{condition}_v_base': a['v_base'], f'local_{condition}_v_focal': a['v_focal']})
    np.savez_compressed(os.path.join(NOMINAL_DIR, 'finite_local_maps.npz'), **focal_maps)

    maps = {**geometry}
    for condition, f in [('ctrl', 1.0), ('reduced', NOMINAL_FRACTION)]:
        for c in MAP_CONTRASTS:
            pv, bv = full_maps[(f, c)]
            maps.update({f'spot_{condition}_{c}_dv': pv - bv, f'spot_{condition}_{c}_ca': calcium_peak(pv, bv, v_half, slope),
                         f'spot_{condition}_{c}_vpeak_mV': pv})
    np.savez_compressed(os.path.join(NOMINAL_DIR, 'anatomy_maps.npz'), **maps)

    t = traces[(1.0, 0)][0]
    trace_arrays = {f'OFF_{"ctrl" if f == 1.0 else "reduced"}_{c}': traces[(f, c)][1]
                    for f in (1.0, NOMINAL_FRACTION) for c in ALL_CONTRASTS}
    np.savez_compressed(os.path.join(NOMINAL_DIR, 'contrast_voltage_traces.npz'), t_ms=t, **trace_arrays)

    def curve(f):
        sub = [r for r in rows if r['retained_effective_inhibition_fraction'] == f]
        return {'contrast_pct': [r['contrast_pct'] for r in sub],
                'mean_peak_dv_mV': [r['mean_peak_dv_mV'] for r in sub],
                'mean_peak_ca_au': [r['mean_peak_ca_au'] for r in sub],
                'normalized_peak_ca': [r['normalized_peak_ca'] for r in sub]}

    def spatial(f):
        r = next(r for r in spatial_rows if r['retained_effective_inhibition_fraction'] == f)
        return dict(distance_um=r['voltage']['distance_um'], transfer_v=r['voltage']['normalized_response'],
                    transfer_ca=r['calcium']['normalized_response'])

    write_json(os.path.join(NOMINAL_DIR, 'figure_data.json'), dict(
        contrast={'ctrl': curve(1.0), 'reduced': curve(NOMINAL_FRACTION)},
        spatial={'ctrl': spatial(1.0), 'reduced': spatial(NOMINAL_FRACTION)},
        imaging_OFF_normalized30=targets, calcium=calcium,
        nominal_retained_effective_inhibition_fraction=NOMINAL_FRACTION))

    np.savez_compressed(os.path.join(OUTPUT_DIR, 'off_contact_peaks.npz'), peak_v=peaks, base_v=bases,
                        retained_fractions=np.asarray(ALL_FRACTIONS), contrast_pct=np.asarray(ALL_CONTRASTS),
                        off_nodes=off)
    write_json(os.path.join(OUTPUT_DIR, 'calibration.json'), dict(
        calcium=calcium, imaging_OFF_normalized30_targets=targets,
        synapses=dict(G_exc=G_EXC, G_inh=G_INH, inh_rev_mV=INH_REV, rest_mV=REST, spot_diameter_um=SPOT_DIAMETER),
        focal=dict(source_segment=source, source_xyz_um=coords[source], conductance_uS=focal_g,
                   inhibition='phasic Exp2Syn spot inhibition, no clamps',
                   inhibition_contrast=FOCAL_INH_CONTRAST,
                   inhibition_total_peak_nS=float(inh_weights.sum() * 1000))))
    write_json(os.path.join(OUTPUT_DIR, 'spatial_sensitivity.json'), spatial_rows)
    write_rows(os.path.join(OUTPUT_DIR, 'contrast_sensitivity.tsv'), rows)

    summary = []
    for f in RETAINED_FRACTIONS:
        sub = {r['contrast_pct']: r for r in rows if r['retained_effective_inhibition_fraction'] == f}
        sp = next(r for r in spatial_rows if r['retained_effective_inhibition_fraction'] == f)
        summary.append(dict(retained_effective_inhibition_fraction=f,
                            voltage_1e_um=sp['voltage']['path_1e_crossing_um'],
                            calcium_1e_um=sp['calcium']['path_1e_crossing_um'],
                            source_dv_mV=sp['source_dv_mV'],
                            normalized_peak_ca30=sub[30]['normalized_peak_ca'],
                            normalized_peak_ca30_minus_target=sub[30]['normalized_peak_ca'] - targets[1]))
    ctrl30 = next(r['normalized_peak_ca'] for r in rows
                  if r['retained_effective_inhibition_fraction'] == 1.0 and r['contrast_pct'] == 30)
    write_json(os.path.join(OUTPUT_DIR, 'sensitivity_summary.json'),
               dict(nominal_fraction=NOMINAL_FRACTION, calcium=calcium, control_normalized_peak_ca30=ctrl30,
                    imaging_OFF_normalized30_targets=targets, control_spatial=spatial_rows[0],
                    sensitivity=summary))

    print("\n" + "=" * 60)
    print("SIMULATIONS COMPLETE")
    print("=" * 60)
    print(f"\n  Elapsed: {time.time() - started:.0f} s")
    print(f"  Output dir: {OUTPUT_DIR}")
    print("\nNext step: Run 02_plot_fig4_panels.py to plot the Figure 4 model panels")


if __name__ == '__main__':
    main()
