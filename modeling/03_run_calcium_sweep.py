"""
03_run_calcium_sweep.py - Sweep the voltage-to-calcium transformation (Figures S6 and S7)

This script asks how much the Figure 4 result depends on the calcium activation
curve. For each of 2,601 pairs of activation midpoint and slope it:
- Converts the saved peak voltages into peak calcium for control and for 20%, 25% and 30% retained inhibition
- Measures how much reduced inhibition raises the response at low versus high contrast
- Matches each reduced response to an equivalent control contrast
- Measures how far calcium spreads from the focal input

No NEURON simulation is run. The voltages, geometry, excitation and inhibition are
the saved results of 01_run_fig4_simulations.py. Only the calcium transformation
changes between grid points. The activation is applied at each OFF contact before
averaging over the 75 OFF contacts. The readout is peak calcium (the Boltzmann of
the peak voltage, baseline-subtracted), not an AUC.

Prerequisites:
    - output/fig4/ from 01_run_fig4_simulations.py

Usage:
    python 03_run_calcium_sweep.py

Output:
    - output/calcium_sweep/calcium_sweep.npz : all metrics on the parameter grid
    - output/calcium_sweep/calcium_sweep_metrics.tsv : the same metrics as a table
    - output/calcium_sweep/summary.json : summary counts and ranges

Author: Pratyush Ramakrishna
Date: 2026
"""

import os
import sys
import csv
import time
import numpy as np
from scipy.special import expit
from scipy.interpolate import PchipInterpolator

from base_functions import load_json, write_json, SCRIPT_DIR

# =============================================================================
# CONFIGURATION
# =============================================================================

# Contrasts and retained inhibition fractions in the saved peaks
CONTRAST = np.array([0, 5, 10, 15, 20, 25, 30, 40, 50, 60, 70])  # %
FRACTIONS = np.array([1., .20, .25, .30])  # 1.0 is the control

# Ratios and normalizations are left undefined below this peak calcium
TECHNICAL_FLOOR = 1e-6

# Grid of activation parameters
HALF_GRID = np.arange(-45, -19.999, .5)                 # mV, midpoint
SLOPE_GRID = np.round(np.arange(1, 6.001, .1), 10)      # mV, slope

# Directory paths
FIG4_DIR = os.path.join(SCRIPT_DIR, 'output', 'fig4')
OUTPUT_DIR = os.path.join(SCRIPT_DIR, 'output', 'calcium_sweep')

# =============================================================================
# FUNCTIONS
# =============================================================================

def evaluate_peak(peaks, bases, pairs):
    """
    Mean peak calcium for every parameter pair and every condition.

    Parameters
    ----------
    peaks, bases : ndarray
        Absolute peak and baseline voltages, shape (n_fractions, n_contrasts, n_contacts)
    pairs : ndarray
        Activation (midpoint, slope) pairs, shape (n_pairs, 2)

    Returns
    -------
    ndarray
        Peak calcium with shape (n_pairs, n_fractions, n_contrasts)
    """
    n_f, n_c, n_k = peaks.shape
    vp, vb = peaks.reshape(-1), bases.reshape(-1)
    output = np.empty((len(pairs), n_f * n_c))
    for start in range(0, len(pairs), 128):
        p = pairs[start:start + 128]
        ca = (expit((vp[None, :] - p[:, 0, None]) / p[:, 1, None])
              - expit((vb[None, :] - p[:, 0, None]) / p[:, 1, None]))
        output[start:start + len(p)] = ca.reshape(len(p), n_f * n_c, n_k).mean(axis=2)
    return output.reshape(len(pairs), n_f, n_c)


def shape_metrics(peak, floor=TECHNICAL_FLOOR):
    """
    Contrast-response shape and the gain from reduced inhibition.

    Parameters
    ----------
    peak : ndarray
        Peak calcium with shape (n_pairs, n_fractions, n_contrasts)
    floor : float
        Peak calcium below which a ratio or normalization is left undefined

    Returns
    -------
    span : ndarray
        Peak calcium at 70% minus peak calcium at 15% contrast
    shape : ndarray
        Peak calcium normalized between the 15% and 70% responses of each condition
    ratio : ndarray
        Reduced / control peak calcium at 15%, 30% and 70% contrast (inf where control is
        silent and reduced inhibition responds)
    delta : ndarray
        Reduced - control peak calcium at 15%, 30% and 70% contrast
    pref : ndarray
        Gain at 15% and at 30% contrast divided by the gain at 70% contrast
    """
    il, im, ih = [int(np.flatnonzero(CONTRAST == c)[0]) for c in (15, 30, 70)]
    span = peak[:, :, ih] - peak[:, :, il]
    shape = np.divide(peak - peak[:, :, il, None], span[:, :, None],
                      out=np.full_like(peak, np.nan), where=span[:, :, None] > floor)
    key = peak[:, :, [il, im, ih]]
    ratio = np.divide(key[:, 1:, :], key[:, 0, None, :],
                      out=np.full_like(key[:, 1:, :], np.nan), where=key[:, 0, None, :] > floor)
    # Control silent but reduced inhibition responding: the gain is unbounded
    ratio[(key[:, 0, None, :] <= floor) & (key[:, 1:, :] > floor)] = np.inf
    delta = key[:, 1:, :] - key[:, 0, None, :]
    pref = np.divide(ratio[:, :, :2], ratio[:, :, 2, None],
                     out=np.full_like(ratio[:, :, :2], np.nan), where=np.isfinite(ratio[:, :, 2, None]))
    return span, shape, ratio, delta, pref


def shifts_for_peak(peak):
    """
    Match each reduced response to an equivalent control contrast.

    The control contrast-response curve is inverted with a monotone cubic
    interpolation (PCHIP), without extrapolation. Each equivalent contrast is
    divided by the contrast actually presented.

    Parameters
    ----------
    peak : ndarray
        Peak calcium with shape (n_pairs, n_fractions, n_contrasts)

    Returns
    -------
    shifts : ndarray
        Horizontal multiplier for each reduced fraction and contrast (NaN if not matched)
    counts : ndarray
        Number of matched contrasts
    mean : ndarray
        Mean multiplier over the matched contrasts
    """
    shifts = np.full((len(peak), 3, len(CONTRAST)), np.nan)
    for i, a in enumerate(peak):
        # Invert only the rising part of the control curve (it is flat at zero below threshold)
        start = int(np.flatnonzero(a[0] <= a[0][0] + TECHNICAL_FLOOR)[-1])
        ctrl, contrast = a[0][start:], CONTRAST[start:]
        if len(ctrl) < 2 or np.any(np.diff(ctrl) <= 0) or ctrl[-1] - ctrl[0] <= TECHNICAL_FLOOR:
            continue
        inverse = PchipInterpolator(ctrl, contrast, extrapolate=False)
        for j in range(3):
            reduced = a[j + 1]
            valid = (CONTRAST > 0) & (reduced > TECHNICAL_FLOOR) & (reduced >= ctrl[0]) & (reduced <= ctrl[-1])
            shifts[i, j, valid] = inverse(reduced[valid]) / CONTRAST[valid]

    counts = np.sum(np.isfinite(shifts), axis=2)
    mean = np.divide(np.nansum(shifts, axis=2), counts,
                     out=np.full(counts.shape, np.nan), where=counts > 0)
    return shifts, counts, mean


def spatial_metrics(maps, pairs):
    """
    Calcium spread from the focal input for every parameter pair.

    The peak voltage of the focal synapse (control and 25% retained) is
    converted to calcium at every segment (peak minus baseline activation),
    normalized to the source, and averaged in 5 um bins of path distance
    (weighted by cable length).

    Parameters
    ----------
    maps : dict
        Focal maps from 01_run_fig4_simulations.py
    pairs : ndarray
        Activation (midpoint, slope) pairs, shape (n_pairs, 2)

    Returns
    -------
    x : ndarray
        Mean path distance of each bin (um)
    profiles : ndarray
        Normalized calcium profiles, shape (n_pairs, 2, n_bins)
    cross : ndarray
        First distance where calcium falls to 1/e of the source (um), shape (n_pairs, 2)
    cable : ndarray
        Normalized calcium summed over the cable length
    source_ca : ndarray
        Calcium at the source
    """
    distance = maps['path_distance_from_source_um']
    length = maps['cable_lengths']
    source = int(maps['source_segment'])
    bins = [(distance >= low) & (distance < low + 5) for low in np.arange(0, distance.max() + 5, 5)]
    bins = [b for b in bins if np.any(b)]
    x = np.array([0.] + [float(np.average(distance[b], weights=length[b])) for b in bins])

    profiles = np.full((len(pairs), 2, len(x)), np.nan)
    cross = np.full((len(pairs), 2), np.nan)
    cable = np.full_like(cross, np.nan)
    source_ca = np.empty_like(cross)

    for j, condition in enumerate(('ctrl', 'reduced')):
        vf, vb = maps[f'local_{condition}_v_focal'], maps[f'local_{condition}_v_base']
        for start in range(0, len(pairs), 128):
            p = pairs[start:start + 128]
            ca = (expit((vf[None, :] - p[:, 0, None]) / p[:, 1, None])
                  - expit((vb[None, :] - p[:, 0, None]) / p[:, 1, None]))
            source_ca[start:start + len(p), j] = ca[:, source]
            norm = np.divide(ca, ca[:, source, None], out=np.full_like(ca, np.nan),
                             where=ca[:, source, None] > TECHNICAL_FLOOR)
            profile = np.column_stack([np.ones(len(p))] +
                                      [np.average(norm[:, b], weights=length[b], axis=1) for b in bins])
            bad = ~np.isfinite(norm[:, 0])
            profile[bad] = np.nan
            profiles[start:start + len(p), j] = profile
            cable[start:start + len(p), j] = norm @ length

            for k, y in enumerate(profile):
                at = np.flatnonzero((y[:-1] > 1 / np.e) & (y[1:] <= 1 / np.e))
                if len(at):
                    n = at[0]
                    cross[start + k, j] = x[n] + (1 / np.e - y[n]) * (x[n + 1] - x[n]) / (y[n + 1] - y[n])

    return x, profiles, cross, cable, source_ca


def extent(a):
    """Count, minimum, median and maximum of the finite values in an array."""
    v = np.asarray(a)
    v = v[np.isfinite(v)]
    return dict(n=int(v.size),
                min=float(v.min()) if len(v) else None,
                median=float(np.median(v)) if len(v) else None,
                max=float(v.max()) if len(v) else None)


def summarize(pairs, n_grid, peak, shape, ratio, delta, pref, count, shift_mean, cross):
    """
    Summarize the grid: counts and ranges for each retained fraction, and the two reference points.

    The reference points are the selected calcium parameters (the last pair) and
    the earlier default of -40 mV midpoint and 4 mV slope.
    """
    cases = []
    for j, r in enumerate(FRACTIONS[1:]):
        p = pref[:n_grid, j]
        valid = ~np.isnan(p).any(axis=1)
        interior = shape[:n_grid, j + 1, 6] - shape[:n_grid, 0, 6]
        cases.append(dict(
            retained_fraction=r,
            valid_preference_points=int(valid.sum()),
            both_low_intermediate_fractional_preference_count=int(np.sum(valid & np.all(p > 1, axis=1))),
            low_preference_count=int(np.sum(p[:, 0] > 1)),
            unbounded_low_preference_count=int(np.sum(np.isposinf(p[:, 0]))),
            unbounded_intermediate_preference_count=int(np.sum(np.isposinf(p[:, 1]))),
            intermediate_preference_count=int(np.sum(p[:, 1] > 1)),
            preference15_over70=extent(p[:, 0]),
            preference30_over70=extent(p[:, 1]),
            shape30_separation=extent(interior),
            positive_shape30_separation_count=int(np.sum(interior > 0)),
            gains={str(c): extent(ratio[:n_grid, j, k]) for k, c in enumerate((15, 30, 70))},
            absolute_increment={str(c): extent(delta[:n_grid, j, k]) for k, c in enumerate((15, 30, 70))},
            shift_mean=extent(shift_mean[:n_grid, j]),
            valid_shift_count_range=extent(count[:n_grid, j]),
            absolute_enhancement_all_nonzero_contrasts_count=int(
                np.sum(np.all(peak[:n_grid, j + 1, 1:] > peak[:n_grid, 0, 1:], axis=1)))))

    default = int(np.flatnonzero((pairs[:n_grid, 0] == -40) & np.isclose(pairs[:n_grid, 1], 4))[0])
    reference_points = []
    for i, label in [(n_grid, 'selected_exact'), (default, 'default_minus40_4')]:
        reference_points.append(dict(
            label=label, v_half_mV=pairs[i, 0], slope_mV=pairs[i, 1], mean_peak_ca=peak[i],
            shape30=shape[i, :, 6], gains_15_30_70=ratio[i], increments_15_30_70=delta[i],
            preference_15_30_over70=pref[i], calcium_only_shift_mean=shift_mean[i],
            shift_counts=count[i], spatial_1e_um=cross[i]))

    return dict(grid_points=n_grid, cases=cases, reference_points=reference_points,
                nominal_spatial=dict(control_1e_um=extent(cross[:n_grid, 0]),
                                     reduced_1e_um=extent(cross[:n_grid, 1]),
                                     ratio=extent(cross[:n_grid, 1] / cross[:n_grid, 0]),
                                     increase_count=int(np.sum(cross[:n_grid, 1] > cross[:n_grid, 0]))))


# =============================================================================
# MAIN SCRIPT
# =============================================================================

def main():
    print("=" * 60)
    print("VG3 Amacrine Cell Model - Calcium Transformation Sweep")
    print("=" * 60)
    started = time.time()

    if os.path.exists(os.path.join(OUTPUT_DIR, 'calcium_sweep.npz')):
        print(f"\nResults already exist in {OUTPUT_DIR}")
        print("  Delete the folder to run the sweep again")
        sys.exit(0)
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Load the saved peak voltages
    print("\n[1/5] Loading saved peak voltages")
    saved = np.load(os.path.join(FIG4_DIR, 'off_contact_peaks.npz'))
    peaks, bases = saved['peak_v'], saved['base_v']
    assert np.allclose(saved['retained_fractions'], FRACTIONS) and np.array_equal(saved['contrast_pct'], CONTRAST)
    maps = np.load(os.path.join(FIG4_DIR, 'nominal25', 'finite_local_maps.npz'))
    prior = load_json(os.path.join(FIG4_DIR, 'sensitivity_summary.json'))
    print(f"      {peaks.shape[0]} conditions x {peaks.shape[1]} contrasts at {peaks.shape[2]} OFF contacts")

    # Parameter grid, plus the selected calcium parameters as one extra point
    hh, kk = np.meshgrid(HALF_GRID, SLOPE_GRID, indexing='ij')
    grid = np.column_stack([hh.ravel(), kk.ravel()])
    n_grid = len(grid)
    calcium = prior['calcium']
    pairs = np.vstack([grid, [calcium['v_half_mV'], calcium['slope_mV']]])

    print(f"\n[2/5] Peak calcium at {n_grid} grid points and the selected point")
    peak = evaluate_peak(peaks, bases, pairs)

    print("\n[3/5] Computing contrast-response metrics")
    span, shape, ratio, delta, pref = shape_metrics(peak)
    shifts, count, shift_mean = shifts_for_peak(peak)
    targets = np.array(prior['imaging_OFF_normalized30_targets'])
    shape_error = (shape[:, 0, 6, None] - targets[0])**2 + (shape[:, 1:, 6] - targets[1])**2

    print("\n[4/5] Computing focal calcium spread")
    x, profiles, cross, cable, source_ca = spatial_metrics(maps, pairs)

    print("\n[5/5] Saving results")
    np.savez_compressed(os.path.join(OUTPUT_DIR, 'calcium_sweep.npz'),
                        v_half_grid_mV=HALF_GRID, slope_grid_mV=SLOPE_GRID, parameter_pairs=pairs,
                        retained_fractions=FRACTIONS, contrast_pct=CONTRAST, mean_peak_ca=peak,
                        peak_span_15_70=span, normalized_peak_ca=shape, gain_15_30_70=ratio,
                        increment_15_30_70=delta, preference_15_30_over70=pref, shape30_sse=shape_error,
                        calcium_shift_by_contrast=shifts, calcium_shift_count=count,
                        calcium_shift_mean=shift_mean, spatial_distance_um=x,
                        spatial_profiles_ctrl25=profiles, spatial_1e_ctrl25_um=cross,
                        spatial_weighted_cable_ctrl25_um=cable, spatial_source_ca_ctrl25=source_ca)

    fields = (['point_type', 'v_half_mV', 'slope_mV', 'retained_fraction']
              + [f'{p}{c}' for p in ['control_peak_', 'reduced_peak_', 'delta_peak_', 'gain_'] for c in [15, 30, 70]]
              + ['preference15_over70', 'preference30_over70', 'control_shape30', 'reduced_shape30',
                 'shape30_difference', 'shape30_sse', 'calcium_shift_mean', 'valid_shift_count',
                 'valid_shift_contrasts', 'control_spatial_1e_um', 'reduced25_spatial_1e_um'])
    with open(os.path.join(OUTPUT_DIR, 'calcium_sweep_metrics.tsv'), 'w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        for i, p in enumerate(pairs):
            for j, f in enumerate(FRACTIONS[1:]):
                row = dict(point_type='grid' if i < n_grid else 'selected_exact',
                           v_half_mV=p[0], slope_mV=p[1], retained_fraction=f,
                           preference15_over70=pref[i, j, 0], preference30_over70=pref[i, j, 1],
                           control_shape30=shape[i, 0, 6], reduced_shape30=shape[i, j + 1, 6],
                           shape30_difference=shape[i, j + 1, 6] - shape[i, 0, 6],
                           shape30_sse=shape_error[i, j], calcium_shift_mean=shift_mean[i, j],
                           valid_shift_count=count[i, j],
                           valid_shift_contrasts=','.join(map(str, CONTRAST[np.isfinite(shifts[i, j])])),
                           control_spatial_1e_um=cross[i, 0], reduced25_spatial_1e_um=cross[i, 1])
                for k, c in enumerate([15, 30, 70]):
                    ci = int(np.flatnonzero(CONTRAST == c)[0])
                    row.update({f'control_peak_{c}': peak[i, 0, ci], f'reduced_peak_{c}': peak[i, j + 1, ci],
                                f'delta_peak_{c}': delta[i, j, k], f'gain_{c}': ratio[i, j, k]})
                # Undefined values are left blank
                writer.writerow({k: '' if isinstance(v, (float, np.floating)) and not np.isfinite(v) else v
                                 for k, v in row.items()})

    summary = summarize(pairs, n_grid, peak, shape, ratio, delta, pref, count, shift_mean, cross)

    # How many grid points keep a valid preference if the floor is changed
    summary['technical_floor_sensitivity'] = []
    for floor in [1e-8, 1e-6, 1e-4]:
        _, _, _, _, pp = shape_metrics(peak, floor)
        summary['technical_floor_sensitivity'].append(dict(
            peak_floor=floor,
            valid_preference_points=np.sum(~np.isnan(pp[:n_grid]).any(axis=2), axis=0),
            both_preference_count=np.sum(~np.isnan(pp[:n_grid]).any(axis=2) & np.all(pp[:n_grid] > 1, axis=2), axis=0)))
    write_json(os.path.join(OUTPUT_DIR, 'summary.json'), summary)

    print("\n" + "=" * 60)
    print("SWEEP COMPLETE")
    print("=" * 60)
    print(f"\n  Elapsed: {time.time() - started:.0f} s")
    print(f"  Output dir: {OUTPUT_DIR}")
    print("\nNext step: Run 04_plot_calcium_sweep.py to plot the Figure S6 and S7 panels")


if __name__ == '__main__':
    main()
