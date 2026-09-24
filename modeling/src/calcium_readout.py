"""Explicit, post-hoc calcium readouts of absolute membrane voltage.

These functions do not add calcium channels to NEURON, change membrane voltage,
model channel kinetics, or predict calibrated calcium concentrations. They
combine an *instantaneous* Boltzmann activation with an optional calcium driving
factor, then optionally apply linear calcium clearance. All conditions share the
same driving-factor normalization at ``rest_v``; no response is normalized to its
own maximum. There is no indicator saturation or vesicle-release transformation.

``fixed`` isolates activation, ``ohmic`` uses a fixed calcium reversal potential,
and ``ghk`` uses the constant-field flux for a divalent ion at fixed intracellular
and extracellular concentrations. The GHK reversal therefore follows those
concentrations and temperature, rather than the separate ohmic ``eca_mV``.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.signal import lfilter
from scipy.special import expit


_R = 8.31446261815324  # J / (mol K)
_F = 96485.33212  # C / mol
_Z_CA = 2.0


def _voltage_array(v_mV: ArrayLike) -> NDArray[np.float64]:
    voltage = np.asarray(v_mV, dtype=float)
    if not np.all(np.isfinite(voltage)):
        raise ValueError("Voltages must be finite absolute membrane voltages in mV.")
    return voltage


def activation(
    v_mV: ArrayLike, v_half: float = -40.0, slope: float = 4.0
) -> NDArray[np.float64]:
    """Return instantaneous activation, 1 / (1 + exp(-(V-v_half)/slope)).

    Both voltage parameters are in mV. A single activation factor is used; its
    midpoint and slope are illustrative parameters, not a VG3 channel fit.
    """
    if not np.isfinite(v_half) or not np.isfinite(slope) or slope <= 0:
        raise ValueError("v_half must be finite and slope must be finite and positive.")
    return np.asarray(expit((_voltage_array(v_mV) - v_half) / slope))


def _bernoulli(x: NDArray[np.float64]) -> NDArray[np.float64]:
    """Numerically stable x/expm1(x), including its removable zero singularity."""
    result = np.empty_like(x)
    small = np.abs(x) < 1e-4
    xs = x[small]
    result[small] = 1.0 - xs / 2.0 + xs**2 / 12.0 - xs**4 / 720.0
    positive = (x >= 1e-4)
    xp = x[positive]
    # This form avoids overflow at strongly positive voltages.
    result[positive] = xp * np.exp(-xp) / (-np.expm1(-xp))
    negative = x <= -1e-4
    xn = x[negative]
    result[negative] = xn / np.expm1(xn)
    return result


def _ghk_inward(
    v_mV: NDArray[np.float64], cai_mM: float, cao_mM: float, temp_C: float
) -> NDArray[np.float64]:
    """Inward GHK flux with voltage-independent positive constants omitted.

    With x=zFV/(RT), this is Cao*x/(exp(x)-1) - Cai*x/(1-exp(-x)).
    Its V=0 limit is Cao-Cai, and it changes sign at the calcium Nernst
    potential. Concentrations stay fixed; flux does not feed back on them.
    """
    x = (_Z_CA * _F / (_R * (temp_C + 273.15))) * (v_mV / 1000.0)
    return cao_mM * _bernoulli(x) - cai_mM * _bernoulli(-x)


def driving_factor(
    v_mV: ArrayLike,
    mode: str = "fixed",
    rest_v: float = -45.0,
    eca_mV: float = 120.0,
    cai_mM: float = 0.0001,
    cao_mM: float = 2.0,
    temp_C: float = 33.0,
) -> NDArray[np.float64]:
    """Return a calcium driving factor normalized to 1 at the shared rest_v.

    ``fixed`` returns ones. ``ohmic`` returns (eca_mV-V)/(eca_mV-rest_v).
    ``ghk`` returns inward constant-field flux / inward flux at rest_v,
    using divalent calcium, concentrations in mM, and temperature in Celsius.
    ``eca_mV`` is used only in ``ohmic`` mode. The reference voltage must lie
    below the relevant reversal potential, so its inward reference is positive.

    No rectification is applied: above calcium reversal the signed factor is
    negative (outward calcium flux). Such voltages are outside the intended
    dendritic response regime, but preserving the reversal aids validation.
    """
    voltage = _voltage_array(v_mV)
    if not np.isfinite(rest_v):
        raise ValueError("rest_v must be finite.")
    if mode == "fixed":
        return np.ones_like(voltage)
    if mode == "ohmic":
        if not np.isfinite(eca_mV) or eca_mV <= rest_v:
            raise ValueError("Ohmic eca_mV must be finite and greater than rest_v.")
        return (eca_mV - voltage) / (eca_mV - rest_v)
    if mode != "ghk":
        raise ValueError("mode must be 'fixed', 'ohmic', or 'ghk'.")
    if (
        not np.isfinite(cai_mM)
        or not np.isfinite(cao_mM)
        or cai_mM <= 0
        or cao_mM <= 0
    ):
        raise ValueError("Calcium concentrations must be finite and positive.")
    if not np.isfinite(temp_C) or temp_C <= -273.15:
        raise ValueError("temp_C must be finite and above absolute zero.")
    reference = _ghk_inward(np.asarray(rest_v), cai_mM, cao_mM, temp_C)
    if reference <= 0:
        raise ValueError("rest_v must be below the GHK calcium reversal potential.")
    return _ghk_inward(voltage, cai_mM, cao_mM, temp_C) / reference


def influx_proxy(
    v_mV: ArrayLike,
    *,
    mode: str = "fixed",
    v_half: float = -40.0,
    slope: float = 4.0,
    rest_v: float = -45.0,
    eca_mV: float = 120.0,
    cai_mM: float = 0.0001,
    cao_mM: float = 2.0,
    temp_C: float = 33.0,
) -> NDArray[np.float64]:
    """Return activation(V) * driving_factor(V), in arbitrary influx units.

    Input is absolute voltage, not change from baseline. Resting influx is
    included. A common rest reference is used for every condition, and the
    result is never normalized by a condition's peak or baseline activation.
    """
    return activation(v_mV, v_half, slope) * driving_factor(
        v_mV, mode, rest_v, eca_mV, cai_mM, cao_mM, temp_C
    )


def calcium_trace(
    v_mV: ArrayLike,
    dt_ms: float,
    *,
    tau_ms: float = 50.0,
    mode: str = "fixed",
    v_half: float = -40.0,
    slope: float = 4.0,
    rest_v: float = -45.0,
    eca_mV: float = 120.0,
    cai_mM: float = 0.0001,
    cao_mM: float = 2.0,
    temp_C: float = 33.0,
) -> NDArray[np.float64]:
    """Filter full influx with linear clearance and subtract a common baseline.

    Time is the last array axis, regularly sampled at dt_ms. The underlying
    state obeys dC/dt = J(V) - C/tau_ms, with C[0] = tau_ms*J(rest_v), the
    steady resting state. Return C - tau_ms*J(rest_v) in arbitrary influx
    units * ms. This is neither a calibrated calcium concentration nor dF/F.

    For interval [t[i], t[i+1]), influx is held at J(V[i]). The exact update
    is C[i+1] = exp(-dt/tau)*C[i] + tau*(1-exp(-dt/tau))*J[i]. It is evaluated
    in baseline-subtracted coordinates for numerical stability, after full
    influx is computed. Therefore the first returned sample is zero even if
    the supplied voltage begins away from rest; any preceding equilibration
    should be included in the input trajectory. Negative changes are retained.

    Activation is instantaneous and clearance linear. There is no membrane
    feedback, calcium-dependent channel regulation, indicator saturation, or
    nonlinear calcium-to-release transformation.
    """
    voltage = _voltage_array(v_mV)
    if voltage.ndim == 0 or voltage.shape[-1] == 0:
        raise ValueError("v_mV must have a nonempty last axis representing time.")
    if not np.isfinite(dt_ms) or dt_ms <= 0:
        raise ValueError("dt_ms must be finite and positive.")
    if not np.isfinite(tau_ms) or tau_ms <= 0:
        raise ValueError("tau_ms must be finite and positive.")
    options = dict(
        mode=mode,
        v_half=v_half,
        slope=slope,
        rest_v=rest_v,
        eca_mV=eca_mV,
        cai_mM=cai_mM,
        cao_mM=cao_mM,
        temp_C=temp_C,
    )
    influx = influx_proxy(voltage, **options)
    resting_influx = influx_proxy(rest_v, **options)
    result = np.zeros_like(influx)
    if voltage.shape[-1] > 1:
        decay = np.exp(-dt_ms / tau_ms)
        gain = tau_ms * (-np.expm1(-dt_ms / tau_ms))
        result[..., 1:] = lfilter(
            [gain], [1.0, -decay], influx[..., :-1] - resting_influx, axis=-1
        )
    return result
