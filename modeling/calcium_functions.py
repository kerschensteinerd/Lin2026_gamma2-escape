"""
calcium_functions.py - Voltage-to-calcium readout for the physiological VG3 model

This module contains functions for:
- Converting peak membrane voltage into a calcium proxy (Boltzmann activation)
- Averaging the readout over the OFF contacts (applied per contact, then averaged)
- Normalizing contrast responses between 15% and 70% contrast
- Fitting the activation midpoint and slope to the OFF imaging

The readout is a plain Boltzmann on PEAK voltage (no clearance, no AUC). Because
the Boltzmann is monotonic, peak calcium equals the Boltzmann of the peak voltage.
Values are activation fractions (0 to 1, baseline-subtracted), not calibrated
calcium concentration or dF/F.

Author: Pratyush Ramakrishna
Date: 2026
"""

import numpy as np
from scipy.special import expit

CTRL_TO_REDUCED_70 = 0.6     # prior on absolute control / reduced calcium at 70%


def activation(v_mV, v_half, slope):
    """
    Boltzmann activation of the calcium influx.

    Parameters
    ----------
    v_mV : ndarray
        Absolute membrane voltage (mV)
    v_half : float
        Voltage at half activation (mV)
    slope : float
        Slope of the activation curve (mV)

    Returns
    -------
    ndarray
        Activation between 0 and 1
    """
    return expit((np.asarray(v_mV) - v_half) / slope)


def calcium_peak(v_peak, v_base, v_half, slope):
    """
    Peak calcium proxy: activation at the peak voltage minus activation at baseline.

    Parameters
    ----------
    v_peak, v_base : ndarray
        Absolute peak and baseline voltage of each segment (mV)
    v_half, slope : float
        Activation midpoint and slope (mV)

    Returns
    -------
    ndarray
        Baseline-subtracted calcium proxy of each segment
    """
    return activation(v_peak, v_half, slope) - activation(v_base, v_half, slope)


def normalize_15_70(values, contrasts):
    """
    Normalize a contrast-response curve so its 15% value is 0 and its 70% value is 1.

    Parameters
    ----------
    values : ndarray
        Response at each contrast
    contrasts : sequence
        Contrasts in % (must contain 15 and 70)

    Returns
    -------
    ndarray
        Normalized responses (NaN if the 15-70% span is zero)
    """
    values = np.asarray(values, dtype=float)
    contrasts = list(contrasts)
    low, high = values[contrasts.index(15)], values[contrasts.index(70)]
    if high - low <= 0:
        return np.full_like(values, np.nan)
    return (values - low) / (high - low)


def fit_readout(peaks, bases, contrasts, targets30, targets15=(0.0, 0.02), weight15=0.25,
                ratio70=CTRL_TO_REDUCED_70, weight_ratio=0.5,
                half_grid=np.arange(-45.0, -20.0, 0.1), slope_grid=np.arange(2.5, 5.01, 0.1)):
    """
    Fit the activation midpoint and slope to the OFF imaging.

    For each condition, calcium is computed per OFF contact and averaged. The
    error is the squared distance of the 15-70% normalized response at 30% from
    the imaging target, plus a weak term on the raw 15% / 70% ratio (the imaging
    traces show almost no response at 15%), plus a term pulling the absolute
    control / reduced-inhibition ratio at 70% towards ratio70 (the imaging
    normalizes each genotype separately, so this ratio is a prior).

    Parameters
    ----------
    peaks, bases : dict
        condition -> array (n_contrasts, n_contacts) of absolute peak / baseline voltage,
        ordered control first, then reduced inhibition
    contrasts : sequence
        Contrasts in % (must contain 15, 30 and 70)
    targets30 : sequence
        Normalized 30% targets, ordered as the conditions in peaks
    targets15 : sequence
        Raw 15% / 70% ratio targets, ordered as the conditions in peaks
    weight15, weight_ratio : float
        Weights of the 15% term and of the 70% control / reduced term
    ratio70 : float
        Target absolute control / reduced-inhibition calcium at 70%
    half_grid, slope_grid : ndarray
        Candidate midpoints and slopes (mV). Slopes start at 2.5 mV: steeper
        readouts overfit the normalized imaging and are less like the voltage
        dependence of calcium channels.

    Returns
    -------
    v_half, slope : float
        Best midpoint and slope (mV)
    error : float
        Fit error at the best pair
    """
    contrasts = list(contrasts)
    i15, i30, i70 = contrasts.index(15), contrasts.index(30), contrasts.index(70)
    best = (np.inf, None, None)
    for slope in slope_grid:
        for v_half in half_grid:
            err, ca70 = 0.0, []
            for (cond, p), t30, t15 in zip(peaks.items(), targets30, targets15):
                ca = calcium_peak(p, bases[cond], v_half, slope).mean(axis=1)
                if ca[i70] - ca[i15] <= 0:
                    err = np.inf
                    break
                err += ((ca[i30] - ca[i15]) / (ca[i70] - ca[i15]) - t30) ** 2
                err += weight15 * (ca[i15] / ca[i70] - t15) ** 2
                ca70.append(ca[i70])
            if np.isfinite(err):
                err += weight_ratio * (ca70[0] / ca70[1] - ratio70) ** 2
            if err < best[0]:
                best = (err, float(v_half), float(slope))
    return best[1], best[2], best[0]
