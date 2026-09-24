#!/usr/bin/env python3
"""Standalone numerical checks for the explicit Fig. 3 calcium readout.

Run from any directory with the repository virtualenv Python. This checks
mathematical implementation, not the biological suitability of its parameters.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
from numpy.testing import assert_allclose

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from calcium_readout import activation, calcium_trace, driving_factor, influx_proxy


def main() -> None:
    checks = []
    assert_allclose(activation([-1e4, -40.0, 1e4]), [0.0, 0.5, 1.0])
    checks.append("stable Boltzmann activation, including midpoint and extreme voltages")

    voltages = np.array([-45.0, -30.0, 0.0, 120.0])
    assert_allclose(driving_factor(voltages, "fixed"), 1.0)
    assert_allclose(driving_factor(voltages, "ohmic"), [1.0, 150 / 165, 120 / 165, 0.0])
    checks.append("shared rest normalization and analytic ohmic ratios/reversal")

    cai, cao, temp = 0.0001, 2.0, 33.0
    # Independent direct constant-field calculation at the nonzero reference.
    rt_over_zf_mV = 8.31446261815324 * (temp + 273.15) / (2 * 96485.33212) * 1000
    x_rest = -45.0 / rt_over_zf_mV
    raw_rest = x_rest * (cao * np.exp(-x_rest) - cai) / (-np.expm1(-x_rest))
    zero_limit = (cao - cai) / raw_rest
    around_zero = driving_factor(np.array([-1e-8, 0, 1e-8]), "ghk")
    assert_allclose(around_zero, zero_limit, rtol=5e-10)
    assert np.all(np.diff(around_zero) < 0)
    checks.append("GHK zero-voltage limit and continuity through zero")

    eca_ghk_mV = rt_over_zf_mV * np.log(cao / cai)
    around_reversal = driving_factor(np.array([eca_ghk_mV - 1, eca_ghk_mV, eca_ghk_mV + 1]), "ghk")
    assert_allclose(around_reversal[1], 0.0, atol=1e-18)
    assert around_reversal[0] > 0 > around_reversal[2]
    assert np.all(np.isfinite(driving_factor(np.array([-1e4, 1e4]), "ghk")))
    checks.append("GHK concentration/temperature-derived reversal, flux sign, and extreme-voltage stability")

    for mode in ("fixed", "ohmic", "ghk"):
        assert_allclose(calcium_trace(np.full((2, 150), -45.0), 0.1, mode=mode), 0.0, atol=0)
        assert_allclose(driving_factor(-45.0, mode), 1.0)
    checks.append("zero baseline response for all modes with steady-state initialization")

    dt, tau, n = 0.37, 43.0, 2000
    levels = np.array([-50.0, -30.0])
    voltage_steps = np.repeat(levels[:, None], n, axis=-1)
    times = np.arange(n) * dt
    for mode in ("fixed", "ohmic", "ghk"):
        delta_j = influx_proxy(levels, mode=mode) - influx_proxy(-45.0, mode=mode)
        expected = tau * delta_j[:, None] * (-np.expm1(-times / tau))
        actual = calcium_trace(voltage_steps, dt, tau_ms=tau, mode=mode)
        assert_allclose(actual, expected, atol=2e-12, rtol=2e-13)
        assert np.all(actual[:, 0] == 0)
    checks.append("last-axis filtering matches analytic positive/negative constant-input clearance")

    # Verify the implementation against the full-state recurrence on a varying
    # waveform; this catches incorrect baseline subtraction or sample alignment.
    trajectory = -45 + 18 * np.sin(np.linspace(0, 4 * np.pi, 300))
    full_j = influx_proxy(trajectory, mode="ghk")
    baseline_c = tau * influx_proxy(-45.0, mode="ghk")
    full_c = np.empty_like(trajectory)
    full_c[0] = baseline_c
    decay = np.exp(-dt / tau)
    for i in range(1, len(full_c)):
        full_c[i] = decay * full_c[i - 1] + tau * (1 - decay) * full_j[i - 1]
    assert_allclose(calcium_trace(trajectory, dt, tau_ms=tau, mode="ghk"), full_c - baseline_c, atol=2e-12)
    checks.append("time-varying trace equals independent full-influx recurrence")

    print(json.dumps({"status": "passed", "checks": checks, "ghk_reversal_mV": float(eca_ghk_mV), "ghk_zero_voltage_factor": float(zero_limit)}, indent=2))


if __name__ == "__main__":
    main()
