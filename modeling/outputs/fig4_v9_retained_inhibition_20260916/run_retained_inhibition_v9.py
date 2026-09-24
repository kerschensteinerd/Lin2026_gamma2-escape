#!/usr/bin/env python3
"""Figure 4 v9: fixed-readout sensitivity to retained effective inhibition.

Run with this repository's .venv/bin/python. Historical outputs are read-only.
The nominal residual is 25%; 20% and 30% are sensitivity cases, not fitted
transmitter fractions. All cases inherit the exact previous activation curve,
geometry, electrical parameters, excitation and focal drive without refitting.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
import shutil
import sys
import time

import numpy as np
import scipy
from scipy.special import expit
import neuron
from neuron import h

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts/fig3"))
from run_fig3_pilot import (build, load_named_config, load_fits, exc_eff_weight,
                           inh_eff_weight, gaussian_spot_weights, path_distance,
                           write_rows)
from run_anatomy_maps import morphology_geometry, array_hash
from run_anatomy_v3 import readout
from protocol_v3 import simulate_v3

PRIOR = ROOT / "outputs/fig3_v3"
PRIOR_SPATIAL = ROOT / "outputs/fig3_v3_spatial"
DEFAULT_OUT = ROOT / "outputs/fig4_v9_retained_inhibition_20260916"
FRACTIONS = (0.20, 0.25, 0.30)
KEY_CONTRASTS = (15, 30, 70)
DENSE = (0, 5, 10, 15, 20, 25, 30, 40, 50, 60, 70)


def file_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def spatial_metrics(y, distance, length):
    """Preserve the existing 5-µm, cable-weighted first-crossing definition."""
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


def simulate_effective(model, exc_nodes, ew, iw, record_nodes, prior_parameters,
                       fraction, contrast, control=False, dt=0.1, sample_dt=2.0,
                       tstop=3052.0):
    """Adapt the unchanged historical integrator to a total residual load.

    The legacy glycine_fraction argument is only a numerical adapter here:
    with gaba_remaining=0 its output is one total shunt of iw*fraction. No
    transmitter identity is assigned. Control uses its exact historical args.
    """
    legacy = copy.deepcopy(prior_parameters)
    if not control:
        legacy["glycine_fraction"] = fraction
    return simulate_v3(model, exc_nodes, ew, iw, record_nodes, legacy,
                       1.0 if control else 0.0, contrast / 100,
                       dt=dt, sample_dt=sample_dt, tstop=tstop)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    out = args.outdir.resolve()
    if out.exists() and any(out.iterdir()):
        raise SystemExit(f"Refusing to overwrite a nonempty output directory: {out}")
    out.mkdir(parents=True, exist_ok=True)
    nominal = out / "nominal25"
    nominal.mkdir()
    started = time.time()
    cal = json.loads((PRIOR / "calibration.json").read_text())
    prior_parameters = copy.deepcopy(cal["parameters"])
    calcium = copy.deepcopy(cal["calcium"])
    previous = json.loads((PRIOR / "figure_data.json").read_text())
    previous_metadata = json.loads((PRIOR / "metadata.json").read_text())
    focal_prior = json.loads((PRIOR_SPATIAL / "finite_local_summary.json").read_text())
    source = int(focal_prior["source_segment"])
    source_g = float(focal_prior["finite_sourceweight_uS"])
    checks = []

    def check(name, passed, **details):
        checks.append(dict(check=name, passed=bool(passed), **details))
        if not passed:
            write_json(out / "validation.json", checks)
            raise AssertionError(name)

    protected = [*sorted((ROOT / "data").glob("*")),
                 *sorted((ROOT / "data/fitted").glob("*.json")),
                 *sorted((ROOT / "src").glob("*.py")),
                 ROOT / "scripts/fig3/protocol_v3.py",
                 ROOT / "scripts/fig3/run_anatomy_v3.py",
                 ROOT / "scripts/fig3/run_anatomy_maps.py",
                 ROOT / "scripts/fig3/run_fig3_pilot.py",
                 ROOT / "scripts/fig3/calibrate_spatial_v3.py",
                 ROOT / "scripts/fig3/finite_local_spread_v3.py",
                 ROOT / "scripts/fig3/prepare_figure_v3.py",
                 *sorted(PRIOR.glob("*")),
                 PRIOR_SPATIAL / "finite_local_maps.npz",
                 PRIOR_SPATIAL / "finite_local_summary.json",
                 PRIOR_SPATIAL / "selected_spatial_maps.npz",
                 PRIOR_SPATIAL / "selected_parameters.json",
                 ROOT.parent / "Data_Figures 9.xlsx"]
    protected = list(dict.fromkeys(p for p in protected if p.is_file()))
    input_hashes = {str(p.relative_to(ROOT.parent)): file_hash(p) for p in protected}
    write_json(out / "input_hashes.json", input_hashes)
    shutil.copyfile(Path(__file__), out / Path(__file__).name)
    check("exact_frozen_calcium_parameters", calcium == focal_prior["calcium"] or
          all(calcium[k] == focal_prior["calcium"][k] for k in ("v_half_mV", "slope_mV", "tau_ms")))
    check("fixed_activation_mode", calcium["mode"] == "fixed")
    check("focal_drive_not_reselected", source == 4705 and source_g == 0.012785038085113628)
    cfg, _ = load_named_config("sep_g22")
    fits = load_fits()
    model = build(prior_parameters["g_pas_S_cm2"])
    geometry = morphology_geometry(model)
    prior_maps = np.load(PRIOR / "anatomy_maps.npz")
    for key, value in geometry.items():
        check("unchanged_geometry_" + key, key in prior_maps and np.array_equal(value, prior_maps[key]))
    rest = model["rest_v"]
    nodes, inh = model["off_nodes"], model["inh_nodes"]
    coords = model["node_coords"]
    ew = exc_eff_weight(fits, "OFF", 100) * cfg["G_exc"] * prior_parameters["excitation_scale"] * gaussian_spot_weights(coords[nodes, :2], model["soma_xy"], 100, 25)
    iw = inh_eff_weight(fits, "OFF", 100) * cfg["G_inh"] * prior_parameters["inhibition_scale"] * gaussian_spot_weights(coords[inh, :2], model["soma_xy"], 100, 25)
    check("control_background_fixed", np.isclose(iw.sum(), focal_prior["electrical_parameters"]["total_control_inhibition_uS"], rtol=1e-13, atol=0))
    immutable = dict(excitation_unit_contrast_sha256=array_hash(ew),
                     control_inhibition_sha256=array_hash(iw),
                     calcium_sha256=hashlib.sha256(json.dumps(calcium, sort_keys=True).encode()).hexdigest())
    parameters = {k: v for k, v in prior_parameters.items() if k != "glycine_fraction"}
    parameters.update(nominal_retained_effective_inhibition_fraction=0.25,
                      sensitivity_retained_effective_inhibition_fractions=list(FRACTIONS),
                      full_control_inhibition_fraction=1.0)
    provenance = dict(status="running", nominal_fraction=0.25,
                      change="Total effective inhibitory load only; no refit or new transmitter assignment.",
                      parameters=parameters, calcium=calcium,
                      prior_calibration="outputs/fig3_v3/calibration.json",
                      focal_source_segment=source, focal_conductance_uS=source_g,
                      total_control_inhibition_nS=float(iw.sum() * 1000),
                      total_unit_contrast_excitation_nS=float(ew.sum() * 1000),
                      n_model_reconstructions=1, n_off_sampling_sites=len(nodes),
                      n_segments=len(coords), immutable=immutable,
                      runtime=dict(python=sys.executable, numpy=np.__version__, scipy=scipy.__version__, neuron=neuron.__version__),
                      runner_sha256=file_hash(__file__),
                      historical_control="Historical numerical control arrays retained exactly; 70% control rerun verifies unchanged protocol.")
    write_json(out / "provenance.json", provenance)

    # Recompute spatial transfer at each load; never reuse reduced-condition
    # transfer arrays at a different residual fraction.
    distance = path_distance(model, source)
    length = geometry["cable_lengths"]
    areas = np.array([h.area(s.x, sec=s.sec) for s in model["segments"]])
    pernode = np.bincount(inh, weights=iw, minlength=len(coords))
    spatial_rows, spatial_arrays = [], {}
    previous_finite = np.load(PRIOR_SPATIAL / "finite_local_maps.npz")
    for fraction in (1.0, *FRACTIONS):
        for i, seg in enumerate(model["segments"]):
            seg.g_pas = prior_parameters["g_pas_S_cm2"] + 100 * pernode[i] * fraction / areas[i]
            seg.e_pas = rest
        h.finitialize(rest)
        imp = h.Impedance()
        ss = model["segments"][source]
        imp.loc(ss.x, sec=ss.sec)
        imp.compute(0.0)
        rin = float(imp.input(ss.x, sec=ss.sec))
        transfer = np.array([imp.transfer(s.x, sec=s.sec) for s in model["segments"]]) / rin
        source_dv = (model["p"]["exc_rev"] - rest) * source_g * rin / (1 + source_g * rin)
        dv = transfer * source_dv
        ca = (expit((rest + dv - calcium["v_half_mV"]) / calcium["slope_mV"]) - expit((rest - calcium["v_half_mV"]) / calcium["slope_mV"])) * calcium["tau_ms"]
        normalized_ca = ca / ca[source]
        row = dict(retained_effective_inhibition_fraction=fraction,
                   source_conductance_uS=source_g, source_Rin_MOhm=rin,
                   source_dv_mV=source_dv, source_calcium_au_ms=float(ca[source]),
                   voltage=spatial_metrics(transfer, distance, length),
                   calcium=spatial_metrics(normalized_ca, distance, length))
        spatial_rows.append(row)
        spatial_arrays[fraction] = dict(dv=dv, ca=ca, v_normalized=transfer, ca_normalized=normalized_ca)
        check(f"spatial_{fraction}_finite", np.isfinite(dv).all() and np.isfinite(ca).all() and ca.min() >= 0)
        if fraction in (1.0, 0.30):
            historical_key = "ctrl" if fraction == 1 else "gaba"
            for field, old_key in [("dv", "finite_dv"), ("ca", "ca"), ("v_normalized", "finite_v_normalized"), ("ca_normalized", "finite_ca_normalized")]:
                error = float(np.max(np.abs(spatial_arrays[fraction][field] - previous_finite[f"local_{historical_key}_{old_key}"])))
                check(f"spatial_{fraction}_{field}_historical_reproduction", error < 1e-8, maximum_absolute_error=error)
        print(f"Spatial retained={fraction:.2f}: V 1/e={row['voltage']['path_1e_crossing_um']:.6f}, Ca 1/e={row['calcium']['path_1e_crossing_um']:.6f} µm", flush=True)
    # Return passive membrane to its fixed baseline before point-conductance runs.
    for seg in model["segments"]:
        seg.g_pas = prior_parameters["g_pas_S_cm2"]
        seg.e_pas = rest
    finite_maps = {**geometry, "path_distance_from_source_um": distance,
                   "local_source_segment": np.asarray(source), "source_segment": np.asarray(source),
                   "local_source_xyz": coords[source]}
    for condition, fraction in [("ctrl", 1.0), ("reduced", 0.25)]:
        a = spatial_arrays[fraction]
        finite_maps.update({f"local_{condition}_ca": a["ca"],
                            f"local_{condition}_finite_dv": a["dv"],
                            f"local_{condition}_finite_v_normalized": a["v_normalized"],
                            f"local_{condition}_finite_ca_normalized": a["ca_normalized"]})
    # Preserve exact old control arrays in the deliverable, including profiles.
    for key in tuple(finite_maps):
        if key.startswith("local_ctrl_"):
            finite_maps[key] = previous_finite[key].copy()
    np.savez_compressed(nominal / "finite_local_maps.npz", **finite_maps)
    write_json(out / "spatial_sensitivity.json", spatial_rows)

    maps = {**geometry}
    for c in KEY_CONTRASTS:
        for suffix in ("dv", "ca", "ca_auc", "vpeak_mV"):
            key = f"spot_ctrl_{c}_{suffix}"
            maps[key] = prior_maps[key].copy()
    historical_traces = np.load(PRIOR / "contrast_voltage_traces.npz")
    t_ref = historical_traces["t_ms"]
    nominal_archive = {f"OFF_ctrl_{c}": historical_traces[f"OFF_ctrl_{c}"].copy() for c in DENSE}
    t, v = simulate_effective(model, nodes, ew, iw, np.arange(len(coords)), prior_parameters, 1.0, 70, control=True)
    dv, peak, auc = readout(v, t, calcium, rest)
    check("control_time_axis_unchanged", np.array_equal(t, t_ref))
    for field, observed, suffix in [("voltage", dv, "dv"), ("calcium", peak, "ca"), ("auc", auc, "ca_auc")]:
        error = float(np.max(np.abs(observed - prior_maps[f"spot_ctrl_70_{suffix}"])))
        check("control_70_" + field + "_reproduction", error < 1e-6, maximum_absolute_error=error)
    check("control_70_trace_reproduction", np.max(np.abs(v[nodes] - historical_traces["OFF_ctrl_70"])) < 5e-6)
    del v

    rows, trace_archive, input_checks = [], {}, []
    # Nominal full maps and dense curves are generated first for figure assembly.
    for fraction in (0.25, 0.20, 0.30):
        contrasts = (*KEY_CONTRASTS, *(c for c in DENSE if c not in KEY_CONTRASTS))
        for c in contrasts:
            full = fraction == 0.25 and c in KEY_CONTRASTS
            t, v = simulate_effective(model, nodes, ew, iw, np.arange(len(coords)) if full else nodes,
                                     prior_parameters, fraction, c)
            dv, peak, auc = readout(v, t, calcium, rest)
            if full:
                maps.update({f"spot_reduced_{c}_dv": dv, f"spot_reduced_{c}_ca": peak,
                             f"spot_reduced_{c}_ca_auc": auc, f"spot_reduced_{c}_vpeak_mV": dv + rest})
                dv, peak, auc, local_v = dv[nodes], peak[nodes], auc[nodes], v[nodes]
            else:
                local_v = v
            row = dict(retained_effective_inhibition_fraction=fraction, condition="reduced", polarity="OFF",
                       contrast_pct=c, n_model_sites=len(nodes), mean_peak_dv_mV=float(dv.mean()),
                       mean_peak_ca_au_ms=float(peak.mean()), mean_ca_auc_au_ms2=float(auc.mean()),
                       mean_peak_absolute_v_mV=float(dv.mean() + rest))
            rows.append(row)
            trace_archive[f"OFF_residual{round(fraction * 100)}_{c}"] = local_v.astype(np.float32)
            if fraction == 0.25:
                nominal_archive[f"OFF_reduced_{c}"] = local_v.astype(np.float32)
            baseline_error = float(np.max(np.abs(local_v[:, t < 49] - rest)))
            check(f"spot_{fraction}_{c}_finite_and_baseline", np.isfinite(v).all() and np.isfinite(peak).all() and np.isfinite(auc).all() and baseline_error < 1e-8,
                  baseline_max_error_mV=baseline_error)
            check(f"spot_{fraction}_{c}_same_time_axis", np.array_equal(t, t_ref))
            input_checks.append(dict(fraction=fraction, contrast_pct=c,
                                     excitation_sha256=array_hash(ew * c / 100),
                                     total_inhibition_sha256=array_hash(iw * fraction)))
            if fraction == 0.30:
                old = next(r for r in previous_metadata["contrast_rows"] if r["polarity"] == "OFF" and r["condition"] == "gaba" and r["contrast_pct"] == c)
                for field in ("mean_peak_dv_mV", "mean_peak_ca_au_ms", "mean_ca_auc_au_ms2"):
                    error = abs(row[field] - old[field])
                    check(f"residual30_{c}_{field}_historical_reproduction", error < 1e-6, absolute_error=error)
            print(f"Spot retained={fraction:.2f}, contrast={c}: ΔV={row['mean_peak_dv_mV']:.6f}, Ca={row['mean_peak_ca_au_ms']:.6f}, AUC={row['mean_ca_auc_au_ms2']:.6f}", flush=True)
            del v
        subset = [r for r in rows if r["retained_effective_inhibition_fraction"] == fraction]
        low = next(r["mean_ca_auc_au_ms2"] for r in subset if r["contrast_pct"] == 15)
        high = next(r["mean_ca_auc_au_ms2"] for r in subset if r["contrast_pct"] == 70)
        check(f"residual_{fraction}_positive_normalization_span", high > low)
        for row in subset:
            row["model_auc_shape"] = (row["mean_ca_auc_au_ms2"] - low) / (high - low)
        if fraction == 0.25:
            figure = dict(contrast={"ctrl": copy.deepcopy(previous["contrast"]["ctrl"])},
                          spatial={"ctrl": copy.deepcopy(previous["spatial"]["ctrl"])},
                          nominal_retained_effective_inhibition_fraction=0.25,
                          footer="25% effective inhibition retained; prior calcium activation fixed without refitting. Shared absolute scales; I endpoints imposed.")
            rr = sorted(subset, key=lambda r: r["contrast_pct"])
            figure["contrast"]["reduced"] = {
                "contrast_pct": [r["contrast_pct"] for r in rr],
                "mean_peak_dv_mV": [r["mean_peak_dv_mV"] for r in rr],
                "mean_peak_ca_au": [r["mean_peak_ca_au_ms"] for r in rr],
                "mean_ca_auc_au_ms2": [r["mean_ca_auc_au_ms2"] for r in rr],
                "auc_shape": [r["model_auc_shape"] for r in rr]}
            sp = next(r for r in spatial_rows if r["retained_effective_inhibition_fraction"] == 0.25)
            figure["spatial"]["reduced"] = dict(distance_um=sp["voltage"]["distance_um"],
                                                  transfer_v=sp["voltage"]["normalized_response"],
                                                  transfer_ca=sp["calcium"]["normalized_response"])
            np.savez_compressed(nominal / "anatomy_maps.npz", **maps)
            np.savez_compressed(nominal / "contrast_voltage_traces.npz", t_ms=t, **nominal_archive)
            write_json(nominal / "figure_data.json", figure)
            write_json(nominal / "metadata.json", dict(parameters=parameters, calcium=calcium, contrast_rows=subset,
                                                        spatial=sp, status="nominal_complete_sensitivity_running"))
            print(f"NOMINAL_READY {nominal}", flush=True)

    # One stronger-than-nominal response verifies time-step / sampling precision.
    tf, vf = simulate_effective(model, nodes, ew, iw, nodes, prior_parameters, 0.20, 70,
                               dt=0.05, sample_dt=1.0, tstop=3051.0)
    fine = readout(vf, tf, calcium, rest)
    ref = next(r for r in rows if r["retained_effective_inhibition_fraction"] == 0.20 and r["contrast_pct"] == 70)
    for field, value in zip(("mean_peak_dv_mV", "mean_peak_ca_au_ms", "mean_ca_auc_au_ms2"), fine):
        relative_error = abs(ref[field] - float(value.mean())) / abs(float(value.mean()))
        check("strongest_half_step_" + field, relative_error < 0.005, relative_error=relative_error)
    for c in KEY_CONTRASTS:
        r = [next(x for x in rows if x["retained_effective_inhibition_fraction"] == f and x["contrast_pct"] == c) for f in FRACTIONS]
        for field in ("mean_peak_dv_mV", "mean_peak_ca_au_ms", "mean_ca_auc_au_ms2"):
            check(f"monotonic_load_dependence_{c}_{field}", all(r[i][field] >= r[i + 1][field] for i in (0, 1)))
        hashes = [x["excitation_sha256"] for x in input_checks if x["contrast_pct"] == c]
        check(f"same_excitation_all_residuals_{c}", len(set(hashes)) == 1)
    for fraction in FRACTIONS:
        r = sorted([x for x in rows if x["retained_effective_inhibition_fraction"] == fraction], key=lambda x: x["contrast_pct"])
        for field in ("mean_peak_dv_mV", "mean_peak_ca_au_ms", "mean_ca_auc_au_ms2"):
            check(f"monotonic_contrast_{fraction}_{field}", all(a[field] <= b[field] for a, b in zip(r, r[1:])))
    check("frozen_calcium_unchanged_after_runs", calcium == cal["calcium"])
    check("excitation_weights_unchanged", array_hash(ew) == immutable["excitation_unit_contrast_sha256"])
    check("control_weights_unchanged", array_hash(iw) == immutable["control_inhibition_sha256"])
    check("all_figure_map_values_finite", all(np.isfinite(v).all() for v in maps.values()))
    check("published_control_contrast_exact", figure["contrast"]["ctrl"] == previous["contrast"]["ctrl"])
    check("published_control_spatial_exact", figure["spatial"]["ctrl"] == previous["spatial"]["ctrl"])
    final_input_hashes = {str(p.relative_to(ROOT.parent)): file_hash(p) for p in protected}
    check("historical_data_and_sources_unchanged", final_input_hashes == input_hashes)
    write_json(out / "stimulus_hashes.json", input_checks)
    write_rows(out / "contrast_sensitivity.tsv", rows)
    np.savez_compressed(out / "sensitivity_voltage_traces.npz", t_ms=t_ref, **trace_archive)

    targets = cal["targets_OFF_shape30"]
    ctrl30 = previous["contrast"]["ctrl"]["auc_shape"][previous["contrast"]["ctrl"]["contrast_pct"].index(30)]
    summary = []
    for fraction in FRACTIONS:
        rr = [r for r in rows if r["retained_effective_inhibition_fraction"] == fraction and r["contrast_pct"] in KEY_CONTRASTS]
        sp = next(r for r in spatial_rows if r["retained_effective_inhibition_fraction"] == fraction)
        shape30 = next(r["model_auc_shape"] for r in rr if r["contrast_pct"] == 30)
        row = dict(retained_effective_inhibition_fraction=fraction,
                   total_inhibition_nS=float(iw.sum() * fraction * 1000),
                   voltage_1e_um=sp["voltage"]["path_1e_crossing_um"],
                   calcium_1e_um=sp["calcium"]["path_1e_crossing_um"],
                   source_dv_mV=sp["source_dv_mV"], source_calcium_au_ms=sp["source_calcium_au_ms"],
                   normalized_auc30=shape30,
                   normalized_auc30_minus_experimental_target=shape30 - targets[1],
                   fixed_readout_two_group_shape_sse=(ctrl30 - targets[0]) ** 2 + (shape30 - targets[1]) ** 2,
                   contrasts=rr)
        summary.append(row)
    write_json(out / "sensitivity_summary.json", dict(status="complete", nominal_fraction=0.25,
               calcium=calcium, control_normalized_auc30=ctrl30, prior_OFF_shape30_targets=targets,
               control_spatial=spatial_rows[0], sensitivity=summary,
               interpretation="Residual fractions are prescribed, not fitted. Calcium activation was frozen. Report dependence and changed agreement rather than selecting the best-looking fraction."))
    provenance.update(status="complete", elapsed_seconds=time.time() - started,
                      checks_passed=sum(c["passed"] for c in checks), checks_total=len(checks))
    write_json(out / "provenance.json", provenance)
    nominal_metadata = json.loads((nominal / "metadata.json").read_text())
    nominal_metadata["status"] = "complete"
    write_json(nominal / "metadata.json", nominal_metadata)
    write_json(out / "validation.json", checks)
    (out / "README.md").write_text(
        "# Figure 4 v9 retained effective inhibition\n\n"
        "Nominal retained total effective inhibition is 25%; sensitivity cases are 20% and 30%. "
        "The control is unchanged. No activation parameter, focal conductance, or other parameter was refitted. "
        "These are deterministic samples from one reconstructed model, not biological replicates.\n\n"
        "`nominal25/figure_data.json`, `anatomy_maps.npz`, and `finite_local_maps.npz` supply A–E/I. "
        "Condition keys are `ctrl` and `reduced`. A–E have shared absolute scales; only I uses separate "
        "15%/70% endpoint normalization. `sensitivity_summary.json` and `contrast_sensitivity.tsv` report all "
        "prescribed cases; `provenance.json`, `input_hashes.json`, `stimulus_hashes.json`, the copied runner, "
        "and `validation.json` record reproducibility and targeted checks.\n\n"
        "The unchanged legacy protocol function uses an argument called `glycine_fraction` as a numerical "
        "adapter for the retained total load. This version does not assign a universal glycinergic identity "
        "to that load. Newly recomputed voltages at each fraction feed the exact prior calcium transformation.\n\n"
        "Reproduce with `.venv/bin/python scripts/fig4/run_retained_inhibition_v9.py --outdir <new-empty-directory>`. "
        "The runner refuses to overwrite nonempty outputs.\n")
    print(f"COMPLETE: {sum(c['passed'] for c in checks)}/{len(checks)} checks in {time.time()-started:.1f} s; {out}", flush=True)


if __name__ == "__main__":
    main()
