"""
06_run_small_spot_spread.py - Spread of voltage and calcium from a small light spot

This is a more physiological version of the Figure 4A-B focal-spread test. Instead
of one strong excitatory synapse on top of the inhibition of a 100 um spot, the
cell is stimulated with a small OFF spot centered on the same dendritic site
(the published focal source). The default diameter is 50 um, the smallest spot in
Kim et al. 2015 Fig. 1C; at 25 um the fitted size tuning is extrapolated. The spot drives both the excitatory (bipolar) and
the inhibitory synapses under it, with the same Gaussian blur, fitted size
tuning and contrast gains as the 100 um spot of Figure 4. Nothing else is
active, so remote dendrites are not held under unopposed inhibition.

For each contrast and retained inhibition (control, 25%) it:
- Records the peak and pre-stimulus voltage of every segment
- Converts them to peak calcium with the readout fitted in 01_run_fig4_simulations.py
- Normalizes voltage and calcium to the spot-center segment and finds the 1/e
  path distance (same 5 um binning as Figure 4B)
- Records voltage traces at the five distances used by 05_plot_example_voltage_traces.py

Prerequisites:
    - output/fig4/calibration.json and output/fig4/nominal25/ from 01_run_fig4_simulations.py

Usage:
    python 06_run_small_spot_spread.py [diameter_um]     (default 50)

Output:
    - output/small_spot_<diameter>um/small_spot_maps.npz, small_spot_traces.npz, small_spot_summary.json

Author: Pratyush Ramakrishna
Date: 2026
"""

import os
import sys
import time
import json
import importlib
import numpy as np
from concurrent.futures import ProcessPoolExecutor

from base_functions import (load_json, write_json, build_model, path_distance, morphology_geometry,
                            run_spot, G_EXC, G_INH, POLARITY, START_MS, SCRIPT_DIR)
from calcium_functions import calcium_peak

# =============================================================================
# CONFIGURATION
# =============================================================================

SPOT_DIAMETER_SMALL = float(sys.argv[1]) if len(sys.argv) > 1 else 50.0   # um
CONTRASTS = (30, 70, 100)           # %
FRACTIONS = (1.0, 0.25)             # retained inhibition: control, nominal reduced
T_STOP = START_MS + 500.0           # ms, long enough for the traces
N_WORKERS = 6
N_TRACE_SITES = 5

FIG4_DIR = os.path.join(SCRIPT_DIR, 'output', 'fig4')
OUTPUT_DIR = os.path.join(SCRIPT_DIR, 'output', f'small_spot_{SPOT_DIAMETER_SMALL:g}um')

fig4_run = importlib.import_module('01_run_fig4_simulations')     # spatial_metrics
traces_05 = importlib.import_module('05_plot_example_voltage_traces')  # trace distances

# =============================================================================
# FUNCTIONS
# =============================================================================

_MODEL = None


def _init_worker():
    """Build the model once in each worker process."""
    global _MODEL
    _MODEL = build_model()
    _MODEL[0]['p']['t_stop'] = T_STOP


def _spot_job(args):
    """
    Run the small spot at one contrast and retained fraction.

    Returns
    -------
    tuple
        (fraction, contrast, peak_v, base_v, t, traces at the trace sites,
        summed peak excitatory and inhibitory conductance in nS)
    """
    fraction, contrast, center_xy, trace_sites = args
    model, fits = _MODEL
    g_exc, g_inh, t, all_v = run_spot(model, fits, SPOT_DIAMETER_SMALL, POLARITY, center_xy=center_xy,
                                      record_v=True, exc_scale=G_EXC, inh_scale=G_INH * fraction,
                                      contrast=contrast / 100)
    pre, post = (t >= START_MS - 10) & (t < START_MS), t >= START_MS
    step = max(1, int(round(1.0 / (t[1] - t[0]))))
    out = (fraction, contrast, all_v[:, post].max(axis=1).astype(np.float32),
           all_v[:, pre].mean(axis=1).astype(np.float32), t[::step].copy(),
           all_v[trace_sites][:, ::step].astype(np.float32), float(g_exc), float(g_inh))
    del all_v
    return out


# =============================================================================
# MAIN SCRIPT
# =============================================================================

def main():
    print("=" * 60)
    print(f"VG3 Amacrine Cell Model - Spread from a {SPOT_DIAMETER_SMALL:g} um Spot")
    print("=" * 60)
    started = time.time()
    if os.path.exists(os.path.join(OUTPUT_DIR, 'small_spot_summary.json')):
        print(f"\nResults already exist in {OUTPUT_DIR}")
        print("  Delete the folder to run the simulations again")
        sys.exit(0)
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("\n[1/4] Building the cell model")
    calibration = load_json(os.path.join(FIG4_DIR, 'calibration.json'))
    v_half, slope = calibration['calcium']['v_half_mV'], calibration['calcium']['slope_mV']
    source = int(calibration['focal']['source_segment'])
    model, fits = build_model()
    geometry = morphology_geometry(model)
    center_xy = model['node_coords'][source, :2].copy()
    distance = path_distance(model, source)
    sites, site_distances = traces_05.pick_distance_sites(distance,
                                                         traces_05.attenuation_curve_targets(N_TRACE_SITES))
    print(f"      Spot: {SPOT_DIAMETER_SMALL:g} um OFF spot centered on segment {source} "
          f"at {np.round(center_xy, 2)} um")
    print(f"      Calcium readout: midpoint {v_half:.2f} mV, slope {slope:.1f} mV")
    print("      Trace sites (um from the spot center): " + ", ".join(f"{d:.0f}" for d in site_distances))

    print(f"\n[2/4] Simulating {len(CONTRASTS) * len(FRACTIONS)} spots ({N_WORKERS} workers)")
    jobs = [(f, c, center_xy, sites) for f in FRACTIONS for c in CONTRASTS]
    with ProcessPoolExecutor(max_workers=N_WORKERS, initializer=_init_worker) as pool:
        results = list(pool.map(_spot_job, jobs))

    print("\n[3/4] Spatial spread")
    length = geometry['cable_lengths']
    maps, traces, rows = {**geometry, 'path_distance_from_source_um': distance,
                          'source_segment': np.asarray(source), 'spot_center_xy': center_xy}, {}, []
    for fraction, contrast, peak_v, base_v, t, site_traces, g_exc, g_inh in results:
        key = f"{'ctrl' if fraction == 1.0 else 'reduced'}_{contrast}"
        dv = peak_v.astype(float) - base_v.astype(float)
        ca = calcium_peak(peak_v.astype(float), base_v.astype(float), v_half, slope)
        voltage = fig4_run.spatial_metrics(dv / dv[source], distance, length) if dv[source] > 0 else None
        calcium = fig4_run.spatial_metrics(ca / ca[source], distance, length) if ca[source] > 1e-6 else None
        maps.update({f'{key}_dv': dv, f'{key}_ca': ca, f'{key}_v_peak': peak_v, f'{key}_v_base': base_v})
        traces[key] = site_traces
        rows.append(dict(condition=key, retained_fraction=fraction, contrast_pct=contrast,
                         total_g_exc_nS=g_exc, total_g_inh_nS=g_inh,
                         source_dv_mV=float(dv[source]), max_dv_mV=float(dv.max()),
                         min_trace_mV=float(site_traces.min()),
                         source_calcium_au=float(ca[source]), max_calcium_au=float(ca.max()),
                         voltage=voltage, calcium=calcium))
        v1e = voltage['path_1e_crossing_um'] if voltage else None
        c1e = calcium['path_1e_crossing_um'] if calcium else None
        print(f"      {key:>12}: gE {g_exc:5.2f} nS, gI {g_inh:5.2f} nS, source dV {dv[source]:6.2f} mV, "
              f"source Ca {ca[source]:.3f}, V 1/e {v1e if v1e is None else round(v1e, 1)} um, "
              f"Ca 1/e {c1e if c1e is None else round(c1e, 1)} um")

    print("\n[4/4] Saving results")
    np.savez_compressed(os.path.join(OUTPUT_DIR, 'small_spot_maps.npz'), **maps)
    np.savez_compressed(os.path.join(OUTPUT_DIR, 'small_spot_traces.npz'), t_ms=results[0][4],
                        sites=np.asarray(sites), site_distance_um=np.asarray(site_distances), **traces)
    write_json(os.path.join(OUTPUT_DIR, 'small_spot_summary.json'),
               dict(spot_diameter_um=SPOT_DIAMETER_SMALL, spot_center_xy=center_xy, source_segment=source,
                    calcium_readout=dict(v_half_mV=v_half, slope_mV=slope), contrasts_pct=CONTRASTS,
                    retained_fractions=FRACTIONS, runs=rows))

    print("\n" + "=" * 60)
    print("SIMULATIONS COMPLETE")
    print("=" * 60)
    print(f"\n  Elapsed: {time.time() - started:.0f} s")
    print(f"  Output dir: {OUTPUT_DIR}")
    print(f"\nNext step: Run 07_plot_small_spot_spread.py {SPOT_DIAMETER_SMALL:g}")


if __name__ == '__main__':
    main()
