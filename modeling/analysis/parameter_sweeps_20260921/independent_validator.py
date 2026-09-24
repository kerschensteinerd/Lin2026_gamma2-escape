#!/usr/bin/env python3
"""Independent audit of the archived-voltage calcium parameter sweep.

The reference deliberately applies calcium_trace to complete local time series,
then trapezoidal integration and spatial averaging. It does not import the sweep
runner, its quadrature/adjoint implementation, or its metric helpers.
"""
from __future__ import annotations

import hashlib
import json
import sys
import csv
from pathlib import Path

import numpy as np
from scipy.special import expit

HERE = Path(__file__).resolve().parent
WORKSPACE = HERE.parents[1]
MODEL = WORKSPACE
ARCHIVE = MODEL / "outputs/fig4_v9_retained_inhibition_20260916"
sys.path.insert(0, str(MODEL / "src"))
from calcium_readout import calcium_trace

CONTRASTS = np.array([0, 5, 10, 15, 20, 25, 30, 40, 50, 60, 70])
REST = -45.0
TAU = 50.0
CHECKS = []


def check(name, passed, **details):
    CHECKS.append(dict(name=name, passed=bool(passed), **details))


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def archive_inputs():
    nominal = np.load(ARCHIVE / "nominal25/contrast_voltage_traces.npz")
    sensitivity = np.load(ARCHIVE / "sensitivity_voltage_traces.npz")
    t = nominal["t_ms"].astype(float)
    check("time_axis_exact_0_to_3052_by_2", np.array_equal(t, np.arange(0, 3053, 2)))
    check("time_axes_identical", np.array_equal(t, sensitivity["t_ms"]))
    traces = {"ctrl": [nominal[f"OFF_ctrl_{c}"].astype(float) for c in CONTRASTS]}
    for retained in [20, 25, 30]:
        traces[f"residual{retained}"] = [sensitivity[f"OFF_residual{retained}_{c}"].astype(float) for c in CONTRASTS]
    check("all_local_trace_shapes", all(v.shape == (75, len(t)) for seq in traces.values() for v in seq))
    check("all_voltage_values_finite", all(np.isfinite(v).all() for seq in traces.values() for v in seq))
    check("zero_contrast_exact_rest", all(np.array_equal(seq[0], np.full_like(seq[0], REST)) for seq in traces.values()))
    check("all_prestimulus_values_exact_rest", all(np.all(v[:, t <= 50] == REST) for seq in traces.values() for v in seq))
    check("nominal_reduced_archive_identical", all(np.array_equal(nominal[f"OFF_reduced_{c}"], sensitivity[f"OFF_residual25_{c}"]) for c in CONTRASTS))
    nominal.close()
    sensitivity.close()
    return t, traces


def full_trace_auc(t, voltage, midpoint, slope):
    ca = calcium_trace(voltage, float(t[1] - t[0]), tau_ms=TAU, mode="fixed", v_half=float(midpoint), slope=float(slope), rest_v=REST)
    window = (t >= 50 - 1e-6) & (t <= 3050 + 1e-6)
    return float(np.mean(np.trapezoid(ca[:, window], x=t[window], axis=1)))


def reference_auc_point(t, traces, midpoint, slope):
    return {condition: np.array([full_trace_auc(t, v, midpoint, slope) for v in seq]) for condition, seq in traces.items()}


def inverse_pchip_shift(control, reduced, floor=0.15):
    """Invert only strictly monotonic response nodes; never extrapolate."""
    if np.any(np.diff(control) <= 0) or control[-1] - control[0] <= floor:
        return dict(mean=None, n=0, equivalent_contrasts=[], contrasts=[], multipliers=[])
    take = (CONTRASTS > 0) & (reduced > floor) & (reduced >= control[0]) & (reduced <= control[-1])
    eq = np.asarray([scalar_pchip(control, CONTRASTS.astype(float), x) for x in reduced[take]], dtype=float)
    ratios = eq / CONTRASTS[take]
    return dict(mean=float(np.mean(ratios)) if len(ratios) else None, n=int(len(ratios)), equivalent_contrasts=eq.tolist(), contrasts=CONTRASTS[take].tolist(), multipliers=ratios.tolist())


def scalar_pchip(x, y, query):
    """Scalar cubic Hermite form, independent of scipy.interpolate."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if not x[0] <= query <= x[-1]:
        return float("nan")
    h = np.diff(x)
    secant = np.diff(y) / h
    assert np.all(h > 0)
    derivative = np.zeros(len(x))
    for i in range(1, len(x) - 1):
        if secant[i - 1] * secant[i] > 0:
            w1, w2 = 2 * h[i] + h[i - 1], h[i] + 2 * h[i - 1]
            derivative[i] = (w1 + w2) / (w1 / secant[i - 1] + w2 / secant[i])

    def endpoint(h0, h1, d0, d1):
        value = ((2 * h0 + h1) * d0 - h0 * d1) / (h0 + h1)
        if value * d0 <= 0:
            return 0.0
        if d0 * d1 < 0 and abs(value) > abs(3 * d0):
            return 3 * d0
        return value

    derivative[0] = endpoint(h[0], h[1], secant[0], secant[1])
    derivative[-1] = endpoint(h[-1], h[-2], secant[-1], secant[-2])
    i = min(len(x) - 2, int(np.searchsorted(x, query, side="right")) - 1)
    t = (query - x[i]) / h[i]
    return float((2 * t**3 - 3 * t**2 + 1) * y[i] + (t**3 - 2 * t**2 + t) * h[i] * derivative[i] + (-2 * t**3 + 3 * t**2) * y[i + 1] + (t**3 - t**2) * h[i] * derivative[i + 1])


def first_crossing(distance, length, normalized):
    x = [0.0]
    y = [1.0]
    for low in np.arange(0, distance.max() + 5, 5):
        use = (distance >= low) & (distance < low + 5)
        if use.any():
            x.append(float(np.average(distance[use], weights=length[use])))
            y.append(float(np.average(normalized[use], weights=length[use])))
    for i in range(len(y) - 1):
        if y[i] > 1 / np.e and y[i + 1] <= 1 / np.e:
            return x[i] + (1 / np.e - y[i]) * (x[i + 1] - x[i]) / (y[i + 1] - y[i])
    return float("nan")


def spatial_reference(midpoint, slope):
    with np.load(ARCHIVE / "nominal25/finite_local_maps.npz") as a:
        source = int(a["source_segment"])
        distance = a["path_distance_from_source_um"]
        lengths = a["cable_lengths"]
        result = {}
        for condition in ["ctrl", "reduced"]:
            calcium = TAU * (expit((REST + a[f"local_{condition}_finite_dv"] - midpoint) / slope) - expit((REST - midpoint) / slope))
            normalized = calcium / calcium[source]
            result[condition] = dict(first_1e_um=float(first_crossing(distance, lengths, normalized)), weighted_cable_um=float(np.sum(lengths * normalized)), source_calcium=float(calcium[source]))
        return result


def compare(name, observed, expected, tolerance=1e-10):
    observed, expected = np.asarray(observed), np.asarray(expected)
    mask_equal = np.array_equal(np.isnan(observed), np.isnan(expected))
    use = np.isfinite(observed) & np.isfinite(expected)
    max_error = float(np.max(np.abs(observed[use] - expected[use]) / np.maximum(1.0, np.abs(expected[use])))) if use.any() else 0.0
    check(name, mask_equal and max_error <= tolerance, max_error_scaled_by_max_abs_expected_or_one=max_error, tolerance=tolerance)


def extent(value):
    value = np.asarray(value)
    finite = value[np.isfinite(value)]
    return dict(n=int(len(finite)), min=float(finite.min()) if len(finite) else None, median=float(np.median(finite)) if len(finite) else None, max=float(finite.max()) if len(finite) else None)


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=HERE / "results")
    parser.add_argument("--report-dir", type=Path, default=HERE / "generated/validation")
    args = parser.parse_args()
    t, traces = archive_inputs()
    out = args.results.resolve()
    report_dir = args.report_dir.resolve()
    report_dir.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((out / "run_manifest.json").read_text())
    summary = json.loads((out / "summary.json").read_text())
    original_summary = json.loads((ARCHIVE / "sensitivity_summary.json").read_text())
    check("production_complete", manifest["status"] == "complete")
    check("declared_readout_and_integration", manifest["rest_mV"] == REST and manifest["tau_ms"] == TAU and manifest["readout"] == "fixed" and manifest["auc_window_ms"] == [50, 3050])
    check("technical_floor_disclosed", manifest["technical_auc_floor"] == 0.15)
    for name, digest in manifest["output_files"].items():
        check("output_hash_" + name, sha256(out / name) == digest)
    external_provenance = []
    checked_input_locations = []
    for name, digest in json.loads((out / "input_hashes.json").read_text()).items():
        if name.startswith("manuscript/MS/"):
            external_provenance.append(dict(path=name, sha256=digest,
                status="not_checked_external_provenance",
                reason="Historical manuscript file is not part of the modeling export or a numerical input."))
            continue
        prefix = "vg3_cid3_rf_simulation-main/"
        path = MODEL / name[len(prefix):] if name.startswith(prefix) else WORKSPACE / name
        historical_source = HERE / "historical_sources/run_calcium_sweep.py"
        if (name == "analysis/parameter_sweeps_20260921/run_calcium_sweep.py"
                and historical_source.is_file() and sha256(historical_source) == digest):
            # The original source is retained verbatim; the runnable copy only
            # changes path resolution, output defaults, and provenance scope.
            path = historical_source
        check("preserved_input_hash_" + name, path.is_file() and sha256(path) == digest)
        checked_input_locations.append(dict(historical_path=name, checked_path=str(path), sha256=digest))
    data = np.load(out / "calcium_sweep.npz")
    pairs = data["parameter_pairs"]
    check("complete_prescribed_grid", np.array_equal(pairs[:2601], np.array([(h, round(k / 10, 1)) for h in np.arange(-45, -19.999, 0.5) for k in range(10, 61)])))
    check("exact_reference_point_additional", pairs.shape == (2602, 2) and np.array_equal(pairs[-1], [original_summary["calcium"]["v_half_mV"], original_summary["calcium"]["slope_mV"]]))
    check("all_condition_axes", np.array_equal(data["retained_fractions"], [1, 0.2, 0.25, 0.3]) and np.array_equal(data["contrast_pct"], CONTRASTS))
    check("grid_axes_match", np.array_equal(data["v_half_grid_mV"], np.unique(pairs[:2601, 0])) and np.array_equal(data["slope_grid_mV"], np.unique(pairs[:2601, 1])))
    auc = data["mean_auc"]
    check("auc_shape_finite_nonnegative", auc.shape == (2602, 4, 11) and np.isfinite(auc).all() and np.min(auc) >= 0)
    check("every_zero_contrast_exact_zero", np.all(auc[:, :, 0] == 0))
    # Includes every corner, the original point, the historical point, the
    # nearest grid point, and several other positions not tested by the runner.
    point_indices = [0, 17, 50, 322, 540, 769, 998, 1300, 1733, 2257, 2550, 2600, 2601]
    direct_errors = []
    point_details = []
    for index in point_indices:
        midpoint, slope = pairs[index]
        reference = reference_auc_point(t, traces, midpoint, slope)
        reference_array = np.array([reference[k] for k in ["ctrl", "residual20", "residual25", "residual30"]])
        err = float(np.max(np.abs(auc[index] - reference_array) / np.maximum(1, np.abs(reference_array))))
        direct_errors.append(err)
        compare(f"full_time_series_auc_point_{index}", auc[index], reference_array)
        spatial = spatial_reference(midpoint, slope)
        compare(f"spatial_crossing_point_{index}", data["spatial_1e_ctrl25_um"][index], [spatial[k]["first_1e_um"] for k in ["ctrl", "reduced"]])
        compare(f"spatial_weighted_cable_point_{index}", data["spatial_weighted_cable_ctrl25_um"][index], [spatial[k]["weighted_cable_um"] for k in ["ctrl", "reduced"]])
        compare(f"spatial_source_calcium_point_{index}", data["spatial_source_ca_ctrl25"][index], [spatial[k]["source_calcium"] for k in ["ctrl", "reduced"]])
        point_details.append(dict(index=index, v_half_mV=float(midpoint), slope_mV=float(slope), max_full_filter_scaled_error=err, spatial=spatial))

    il, im, ih = [int(np.flatnonzero(CONTRASTS == c)[0]) for c in [15, 30, 70]]
    span = auc[:, :, ih] - auc[:, :, il]
    normalized = np.full_like(auc, np.nan)
    good_span = span > 0.15
    for c in range(len(CONTRASTS)):
        normalized[:, :, c][good_span] = ((auc[:, :, c] - auc[:, :, il]) / np.where(good_span, span, 1))[good_span]
    gain = np.full((2602, 3, 3), np.nan)
    delta = np.empty_like(gain)
    for j in range(3):
        for k, c in enumerate([il, im, ih]):
            good = auc[:, 0, c] > 0.15
            gain[good, j, k] = auc[good, j + 1, c] / auc[good, 0, c]
            delta[:, j, k] = auc[:, j + 1, c] - auc[:, 0, c]
    preference = gain[:, :, :2] / gain[:, :, 2, None]
    targets = original_summary["prior_OFF_shape30_targets"]
    sse = (normalized[:, 0, im, None] - targets[0])**2 + (normalized[:, 1:, im] - targets[1])**2
    for name, reference in [("auc_span_15_70", span), ("normalized_auc", normalized), ("gain_15_30_70", gain), ("increment_15_30_70", delta), ("preference_15_30_over70", preference), ("shape30_sse", sse)]:
        compare("every_grid_derived_" + name, data[name], reference)
    check("normalization_validity_mask_exact", np.array_equal(np.isfinite(data["normalized_auc"]), np.broadcast_to(good_span[:, :, None], auc.shape)))
    # Independent scalar inverse PCHIP at every parameter point, fraction, and
    # eligible contrast verifies the interpolation and non-extrapolation masks.
    shift = np.full((2602, 3, 11), np.nan)
    for i in range(len(pairs)):
        for j in range(3):
            reference = inverse_pchip_shift(auc[i, 0], auc[i, j + 1])
            for c, multiplier in zip(reference["contrasts"], reference["multipliers"]):
                shift[i, j, int(np.flatnonzero(CONTRASTS == c)[0])] = multiplier
    count = np.isfinite(shift).sum(axis=2)
    shift_mean = np.full((2602, 3), np.nan)
    nonzero = count > 0
    shift_mean[nonzero] = np.nansum(shift, axis=2)[nonzero] / count[nonzero]
    compare("all_inverse_pchip_matches_and_overlap_masks", data["calcium_shift_by_contrast"], shift)
    compare("all_shift_overlap_counts", data["calcium_shift_count"], count)
    compare("all_shift_means", data["calcium_shift_mean"], shift_mean)
    check("no_zero_contrast_equivalence", np.isnan(shift[:, :, 0]).all())
    distances = data["spatial_1e_ctrl25_um"][:2601]
    paired_crossings = np.isfinite(distances).all(axis=1)
    check("all_grid_resolved_spatial_crossings_broaden", np.all(distances[paired_crossings, 1] > distances[paired_crossings, 0]), valid_pairs=int(paired_crossings.sum()), control_no_crossing=int(np.isnan(distances[:, 0]).sum()), reduced_no_crossing=int(np.isnan(distances[:, 1]).sum()))
    compare("exact_reference_spatial_values", data["spatial_1e_ctrl25_um"][-1], [18.77733921143409, 34.337163202424904])

    with (out / "calcium_sweep_metrics.tsv").open() as stream:
        table = list(csv.DictReader(stream, delimiter="\t"))
    check("table_row_count", len(table) == len(pairs) * 3)
    table_ok = True
    for row_index, row in enumerate(table):
        i, j = divmod(row_index, 3)
        expected = dict(v_half_mV=pairs[i, 0], slope_mV=pairs[i, 1], retained_fraction=[0.2, 0.25, 0.3][j], preference15_over70=preference[i, j, 0], preference30_over70=preference[i, j, 1], control_shape30=normalized[i, 0, im], reduced_shape30=normalized[i, j + 1, im], shape30_difference=normalized[i, j + 1, im] - normalized[i, 0, im], shape30_sse=sse[i, j], calcium_shift_mean=shift_mean[i, j], valid_shift_count=count[i, j], control_spatial_1e_um=data["spatial_1e_ctrl25_um"][i, 0], reduced25_spatial_1e_um=data["spatial_1e_ctrl25_um"][i, 1])
        for k, c in enumerate([15, 30, 70]):
            ci = int(np.flatnonzero(CONTRASTS == c)[0])
            expected.update({f"control_auc_{c}": auc[i, 0, ci], f"reduced_auc_{c}": auc[i, j + 1, ci], f"delta_auc_{c}": delta[i, j, k], f"gain_{c}": gain[i, j, k]})
        for key, value in expected.items():
            table_ok &= (row[key] == "") if not np.isfinite(value) else np.isclose(float(row[key]), value, rtol=1e-12, atol=1e-12)
        table_ok &= row["point_type"] == ("grid" if i < 2601 else "original_exact")
        table_ok &= row["valid_shift_contrasts"] == ",".join(map(str, CONTRASTS[np.isfinite(shift[i, j])]))
    check("all_table_cells_reproduce_arrays_and_blank_masks", table_ok)

    for j, case in enumerate(summary["cases"]):
        p = preference[:2601, j]
        valid = np.isfinite(p).all(axis=1)
        interior = normalized[:2601, j + 1, im] - normalized[:2601, 0, im]
        expected = dict(valid_preference_points=int(valid.sum()), both_low_intermediate_fractional_preference_count=int(np.sum(valid & (p > 1).all(axis=1))), low_preference_count=int(np.sum(p[:, 0] > 1)), intermediate_preference_count=int(np.sum(p[:, 1] > 1)), positive_shape30_separation_count=int(np.sum(interior > 0)), absolute_enhancement_all_nonzero_contrasts_count=int(np.sum(np.all(auc[:2601, j + 1, 1:] > auc[:2601, 0, 1:], axis=1))))
        check(f"summary_counts_fraction_{j}", all(case[k] == v for k, v in expected.items()), **expected)
        ranges = dict(preference15_over70=extent(p[:, 0]), preference30_over70=extent(p[:, 1]), shape30_separation=extent(interior), shift_mean=extent(shift_mean[:2601, j]), valid_shift_count_range=extent(count[:2601, j]))
        check(f"summary_ranges_fraction_{j}", all(all((case[k][field] is None) if value[field] is None else np.isclose(case[k][field], value[field], rtol=1e-12, atol=1e-12) for field in value) for k, value in ranges.items()))
    # The full-precision historical summaries legitimately differ slightly from
    # recomputation from the saved float32 voltage trajectories.
    prior_data = json.loads((ARCHIVE / "nominal25/figure_data.json").read_text())
    compare("nominal_control_auc_matches_historical_float32_tolerance", auc[-1, 0], prior_data["contrast"]["ctrl"]["mean_ca_auc_au_ms2"], 1e-5)
    compare("nominal_reduced_auc_matches_historical_float32_tolerance", auc[-1, 2], prior_data["contrast"]["reduced"]["mean_ca_auc_au_ms2"], 1e-5)
    report = dict(status="pass" if all(c["passed"] for c in CHECKS) else "fail", checks_passed=sum(c["passed"] for c in CHECKS), checks_total=len(CHECKS), method="No runner import. Full existing calcium_trace at 13 parameter points, all 44 complete local time-series arrays; every saved derived metric/mask and TSV field independently recomputed; scalar Hermite inverse PCHIP for all parameter points; independent focal spatial measurements at the same 13 points.", max_full_filter_scaled_error=max(direct_errors), reference_points=point_details, checks=CHECKS, validation_source_sha256=sha256(__file__), production_manifest_sha256=sha256(out / "run_manifest.json"))
    report.update(external_provenance=external_provenance,
                  checked_input_locations=checked_input_locations,
                  scope="All required numerical input and output hashes and numerical checks remain mandatory. Historical manuscript hashes are recorded but are outside this export.")
    (report_dir / "independent_validation.json").write_text(json.dumps(report, indent=2) + "\n")
    lines = ["# Independent calcium-sweep validation", "", f"**{report['status'].upper()}: {report['checks_passed']}/{report['checks_total']} checks.** No NEURON simulation, parameter fit, manuscript edit, or production-output mutation was performed.", "", report["method"], "", f"The largest direct full-filter AUC error, scaled by max(1, reference AUC), was {max(direct_errors):.3g}. The 51 × 51 grid and additional exact original point, 20/25/30% retained-inhibition axes, 2-ms sampling, 50–3050-ms AUC integration, shared resting voltage and clearance, technical denominator masks, and non-extrapolated contrast matching agree with the declared implementation. All nominal reduced voltages are identical across their two archived copies. Input/manuscript hashes and recorded output hashes match.", "", "The 0.15 AUC floor is a numerical reporting threshold, not a measured biological detection limit. Grid coverage is not posterior probability or confidence. Spatial calculations cover only the archived original focal source and nominal 25% retained inhibition; the repeated spatial fields in each TSV fraction row are explicitly labeled reduced25. Calcium-only equivalence is not the manuscript's combined voltage/calcium marker. This audit verifies the calculation and saved data, not biological validity or robustness to alternative electrical parameters, morphologies, sources, or inhibitory distributions.", "", "Historical AUC values were computed before float32 voltage storage; comparison permits the documented 1e-5 relative tolerance. The sweep and full-trace validation both use the actual archived float32 values promoted to float64.", "", "Detailed checks and checked parameter coordinates are in `independent_validation.json`; reproduce with `vg3_cid3_rf_simulation-main/.venv/bin/python analysis/parameter_sweeps_20260921/independent_validator.py`.", ""]
    lines = [line.replace("Input/manuscript hashes and recorded output hashes match.",
                         "Required numerical input and output hashes are checked; historical manuscript hashes are external provenance and are not checked by this export.")
                 .replace("vg3_cid3_rf_simulation-main/.venv/bin/python analysis/parameter_sweeps_20260921/independent_validator.py",
                          "python analysis/parameter_sweeps_20260921/independent_validator.py") for line in lines]
    lines.extend(["", f"External provenance files not checked: {len(external_provenance)}. Details and original hashes are listed in the JSON report.", ""])
    (report_dir / "validation_independent.md").write_text("\n".join(lines))
    print(json.dumps({k: report[k] for k in ["status", "checks_passed", "checks_total", "max_full_filter_scaled_error"]}, indent=2))
    if report["status"] != "pass":
        print(json.dumps([c for c in CHECKS if not c["passed"]], indent=2))
        raise SystemExit(1)


if __name__ == "__main__":
    main()
