#!/usr/bin/env python3
"""Generate anatomical Fig. 3 maps using inherited VG3 voltage parameters.

No parameter is selected to amplify a desired visual effect. A 100-um OFF spot
uses 15%, 30%, and 70% contrast. A separate weak local-input experiment isolates
spatial transfer under a fixed pure-shunting inhibitory background. GABAergic
inhibition is removed while the assumed glycinergic component stays unchanged.
Calcium maps use instantaneous activation followed by linear clearance, without
a voltage-dependent driving term, indicator saturation, or release nonlinearity.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from neuron import h

from run_fig3_pilot import (
    ROOT, build, calcium_trace, exc_eff_weight, gaussian_spot_weights,
    inh_eff_weight, load_fits, load_named_config, path_distance, simulate,
)


CONTRAST_PCT = (15, 30, 70)
CONDITIONS = (("ctrl", 1.0), ("gaba", 0.0))
SAMPLE_DT_MS = 1.0
INTEGRATION_DT_MS = 0.1
GABA_FRACTION = 0.7
CHUNK_SIZE = 256


def morphology_geometry(model):
    """Preserve section polylines and split them at numerical segment borders.

    Every returned edge lies on an actual consecutive pt3d line piece. Splitting
    at segment boundaries gives exact color-to-compartment assignment without
    changing branch geometry or drawing chords between segment centers.
    """
    segments = model["segments"]
    sections = list(model["cell"].all)
    section_index = {sec.name(): i for i, sec in enumerate(sections)}
    segment_lookup = {
        (seg.sec.name(), min(seg.sec.nseg - 1, int(seg.x * seg.sec.nseg))): i
        for i, seg in enumerate(segments)
    }
    edges, edge_nodes, edge_sections = [], [], []
    for section_i, sec in enumerate(sections):
        n3d = int(h.n3d(sec=sec))
        xyz = np.array([
            [h.x3d(i, sec=sec), h.y3d(i, sec=sec), h.z3d(i, sec=sec)]
            for i in range(n3d)
        ])
        if len(xyz) < 2 or sec.L <= 0:
            continue
        arc = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(xyz, axis=0), axis=1))]
        boundaries = np.linspace(0.0, arc[-1], sec.nseg + 1)
        for j in range(n3d - 1):
            start, end = arc[j], arc[j + 1]
            if end <= start:
                continue
            inside = boundaries[(boundaries > start) & (boundaries < end)]
            cuts = np.r_[start, inside, end]
            for left, right in zip(cuts[:-1], cuts[1:]):
                if right - left <= 1e-12:
                    continue
                endpoints = xyz[j] + ((np.array([left, right]) - start) / (end - start))[:, None] * (xyz[j + 1] - xyz[j])
                local_segment = min(sec.nseg - 1, int(((left + right) / 2) / arc[-1] * sec.nseg))
                edges.append(endpoints)
                edge_nodes.append(segment_lookup[(sec.name(), local_segment)])
                edge_sections.append(section_i)
    return {
        "xyz_segments": np.asarray(model["node_coords"], dtype=float),
        "geometry_edges": np.asarray(edges, dtype=float),
        "edge_segment_index": np.asarray(edge_nodes, dtype=np.int32),
        "edge_section_index": np.asarray(edge_sections, dtype=np.int32),
        "cable_lengths": np.array([seg.sec.L / seg.sec.nseg for seg in segments]),
        "segment_section_index": np.array([section_index[seg.sec.name()] for seg in segments], dtype=np.int32),
        "segment_x": np.array([seg.x for seg in segments]),
        "section_nseg": np.array([sec.nseg for sec in sections], dtype=np.int32),
        "soma_xyz": np.asarray(model["cell"].soma_xyz, dtype=float),
        "soma_segment_index": np.asarray(model["soma_idx"], dtype=np.int32),
        "stimulus_center_xy": np.asarray(model["soma_xy"], dtype=float),
        "off_nodes": np.asarray(model["off_nodes"], dtype=np.int32),
        "inh_nodes": np.asarray(model["inh_nodes"], dtype=np.int32),
    }


def peak_calcium_chunked(voltage, rest_v, sample_dt_ms=SAMPLE_DT_MS):
    peaks = np.empty(voltage.shape[0], dtype=float)
    for start in range(0, len(peaks), CHUNK_SIZE):
        stop = min(len(peaks), start + CHUNK_SIZE)
        ca = calcium_trace(voltage[start:stop], sample_dt_ms, mode="fixed", rest_v=rest_v,
                           v_half=-40.0, slope=4.0, tau_ms=50.0)
        peaks[start:stop] = ca.max(axis=-1)
    return peaks


def peak_summary(delta_v, calcium, rest_v, off, lengths):
    absolute_peak_v = delta_v[off] + rest_v
    return {
        "mean_off_peak_dv_mV": float(np.mean(delta_v[off])),
        "length_weighted_off_peak_dv_mV": float(np.average(delta_v[off], weights=lengths[off])),
        "length_weighted_all_peak_dv_mV": float(np.average(delta_v, weights=lengths)),
        "mean_off_peak_ca_au_ms": float(np.mean(calcium[off])),
        "length_weighted_off_peak_ca_au_ms": float(np.average(calcium[off], weights=lengths[off])),
        "length_weighted_all_peak_ca_au_ms": float(np.average(calcium, weights=lengths)),
        "off_peak_absolute_v_mean_mV": float(absolute_peak_v.mean()),
        "off_peak_absolute_v_min_mV": float(absolute_peak_v.min()),
        "off_peak_absolute_v_max_mV": float(absolute_peak_v.max()),
        "all_peak_dv_min_mV": float(delta_v.min()),
        "all_peak_dv_max_mV": float(delta_v.max()),
        "all_peak_ca_min_au_ms": float(calcium.min()),
        "all_peak_ca_max_au_ms": float(calcium.max()),
    }


def array_hash(array):
    return hashlib.sha256(np.asarray(array, dtype=np.float64).tobytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=ROOT / "outputs/fig3_anatomy")
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    started = time.time()
    cfg, _ = load_named_config("sep_g22")
    model = build(cfg["g_pas"])
    fits = load_fits()
    rest = model["rest_v"]
    off, inh = model["off_nodes"], model["inh_nodes"]
    coords = model["node_coords"]
    all_nodes = np.arange(len(coords))
    eweights = exc_eff_weight(fits, "OFF", 100) * cfg["G_exc"] * gaussian_spot_weights(coords[off, :2], model["soma_xy"], 100, 25)
    iweights = inh_eff_weight(fits, "OFF", 100) * cfg["G_inh"] * gaussian_spot_weights(coords[inh, :2], model["soma_xy"], 100, 25)
    maps = morphology_geometry(model)
    lengths = maps["cable_lengths"]
    geometry_before = array_hash(maps["geometry_edges"])
    spot_summary, local_summary, checks = [], [], []
    input_hashes = {}

    for contrast_pct in CONTRAST_PCT:
        contrast = contrast_pct / 100.0
        for condition, gaba_remaining in CONDITIONS:
            t, v = simulate(model, off, eweights * contrast, iweights * contrast, all_nodes,
                            gaba_fraction=GABA_FRACTION, gaba_remaining=gaba_remaining,
                            dt=INTEGRATION_DT_MS, sample_dt=SAMPLE_DT_MS)
            key = f"spot_{condition}_{contrast_pct}"
            dv = v.max(axis=1) - rest
            ca = peak_calcium_chunked(v, rest)
            maps[f"{key}_dv"] = dv
            maps[f"{key}_ca"] = ca
            maps[f"{key}_vpeak_mV"] = dv + rest
            summary = dict(condition=condition, contrast_pct=contrast_pct)
            summary.update(peak_summary(dv, ca, rest, off, lengths))
            spot_summary.append(summary)
            baseline_error = float(np.abs(v[:, t < 49.0] - rest).max())
            checks.append(dict(check=f"{key}_pre_stimulus_baseline", passed=baseline_error < 1e-9,
                               max_abs_error_mV=baseline_error))
            checks.append(dict(check=f"{key}_all_trace_samples_finite", passed=bool(np.isfinite(v).all())))
            input_hashes[key] = {
                "excitation": array_hash(eweights * contrast),
                "retained_glycine": array_hash(iweights * contrast * (1.0 - GABA_FRACTION)),
            }
            print(f"{key}: OFF mean peak dV={summary['mean_off_peak_dv_mV']:.5f} mV, Ca={summary['mean_off_peak_ca_au_ms']:.5f}", flush=True)
            del v

    pilot_metadata_path = ROOT / "outputs/fig3_pilot/metadata.json"
    pilot_sources = json.loads(pilot_metadata_path.read_text()).get("sources", []) if pilot_metadata_path.exists() else []
    source = int(pilot_sources[1]) if len(pilot_sources) >= 2 else 4705
    if source not in off:
        raise ValueError("Inherited local source must be a fixed OFF input segment.")
    maps["local_source_segment"] = np.asarray(source, dtype=np.int32)
    maps["local_source_xyz"] = coords[source]
    maps["path_distance_from_source_um"] = path_distance(model, source)
    local_weight = exc_eff_weight(fits, "OFF", 100) * cfg["G_exc"] * 0.1
    for condition, gaba_remaining in CONDITIONS:
        t, v = simulate(model, [source], [local_weight], iweights, all_nodes,
                        tonic=True, e_inh=rest, gaba_fraction=GABA_FRACTION,
                        gaba_remaining=gaba_remaining, dt=INTEGRATION_DT_MS, sample_dt=SAMPLE_DT_MS)
        dv = v.max(axis=1) - rest
        ca = peak_calcium_chunked(v, rest)
        source_peak_index = int(v[source].argmax())
        source_peak_time = float(t[source_peak_index])
        same_time_dv = v[:, source_peak_index] - rest
        maps[f"local_{condition}_dv"] = dv
        maps[f"local_{condition}_ca"] = ca
        maps[f"local_{condition}_normalized"] = dv / dv[source]
        maps[f"local_{condition}_at_source_peak_dv"] = same_time_dv
        maps[f"local_{condition}_at_source_peak_normalized"] = same_time_dv / dv[source]
        summary = dict(condition=condition, source_segment=source,
                       source_peak_dv_mV=float(dv[source]), source_peak_time_ms=source_peak_time,
                       peak_transfer_weighted_cable_um=float(np.sum(lengths * dv / dv[source])),
                       same_time_transfer_weighted_cable_um=float(np.sum(lengths * same_time_dv / dv[source])))
        summary.update(peak_summary(dv, ca, rest, off, lengths))
        local_summary.append(summary)
        baseline_error = float(np.abs(v[:, t < 49.0] - rest).max())
        checks.append(dict(check=f"local_{condition}_pre_stimulus_baseline", passed=baseline_error < 1e-9,
                           max_abs_error_mV=baseline_error))
        checks.append(dict(check=f"local_{condition}_all_trace_samples_finite", passed=bool(np.isfinite(v).all())))
        input_hashes[f"local_{condition}"] = {
            "excitation": array_hash([local_weight]),
            "retained_glycine": array_hash(iweights * (1.0 - GABA_FRACTION)),
        }
        print(f"local_{condition}: source peak={dv[source]:.5f} mV, source-normalized cable={summary['peak_transfer_weighted_cable_um']:.3f} um", flush=True)
        del v

    # Test map sampling against the pilot's finer 0.5-ms readout grid at fixed
    # input sites. These are peak-metric checks; integration dt remains 0.1 ms.
    for condition, gaba_remaining in CONDITIONS:
        t, v = simulate(model, off, eweights * .7, iweights * .7, off,
                        gaba_fraction=GABA_FRACTION, gaba_remaining=gaba_remaining,
                        dt=INTEGRATION_DT_MS, sample_dt=.5)
        reference = {"dv": v.max(axis=1) - rest, "ca": peak_calcium_chunked(v, rest, .5)}
        for quantity in ("dv", "ca"):
            map_values = maps[f"spot_{condition}_70_{quantity}"][off]
            scale = max(float(np.mean(np.abs(reference[quantity]))), 1e-12)
            error = float(np.max(np.abs(map_values - reference[quantity])) / scale)
            checks.append(dict(check=f"spot_{condition}_70_{quantity}_sampling_convergence",
                               passed=error < .001, max_abs_difference_over_mean_reference=error,
                               tolerance=.001))
        del v

    edge_lengths = np.linalg.norm(np.diff(maps["geometry_edges"], axis=1)[:, 0, :], axis=1)
    lengths_from_edges = np.bincount(maps["edge_segment_index"], weights=edge_lengths, minlength=len(coords))
    length_error = float(np.max(np.abs(lengths_from_edges - lengths)))
    checks.append(dict(check="every_segment_has_geometry", passed=bool(np.all(lengths_from_edges > 0))))
    checks.append(dict(check="geometry_preserves_cable_length_per_segment", passed=length_error < 1e-6,
                       max_abs_error_um=length_error))
    checks.append(dict(check="identical_geometry_all_conditions", passed=geometry_before == array_hash(morphology_geometry(model)["geometry_edges"])))
    for contrast_pct in CONTRAST_PCT:
        checks.append(dict(check=f"spot_{contrast_pct}_identical_excitation_and_retained_glycine",
                           passed=input_hashes[f"spot_ctrl_{contrast_pct}"] == input_hashes[f"spot_gaba_{contrast_pct}"]))
    checks.append(dict(check="local_identical_excitation_and_retained_glycine",
                       passed=input_hashes["local_ctrl"] == input_hashes["local_gaba"]))
    checks.append(dict(check="all_saved_arrays_finite", passed=all(np.isfinite(a).all() for a in maps.values())))
    checks.append(dict(check="stored_maps_are_unclipped_raw_peaks", passed=True,
                       explanation="Raw maxima are stored without clipping, rescaling, or condition-specific normalization; explicitly named local normalized arrays divide by each local source peak."))

    metadata = {
        "status": "exploratory_inherited_voltage_parameters_not_fit_to_fig3",
        "n_segments": len(coords), "n_geometry_edges": len(maps["geometry_edges"]),
        "n_off_input_sites": len(off), "n_inhibitory_contacts": len(inh),
        "geometry": "Actual NEURON pt3d section polylines, split at numerical segment boundaries; XYZ coordinates unchanged; inherited model uniform cable diameter 0.5 um.",
        "coordinate_units": "um", "map_dv_units": "mV above shared rest", "map_ca_units": "arbitrary influx units * ms",
        "rest_v_mV": rest, "dt_ms": INTEGRATION_DT_MS, "sample_dt_ms": SAMPLE_DT_MS,
        "tstop_ms": 650.0, "stimulus_onset_ms": 50.0,
        "spot": {"diameter_um": 100.0, "gaussian_blur_sigma_um": 25.0, "contrasts_pct": list(CONTRAST_PCT),
                 "polarity": "OFF", "inhibitory_reversal_mV": -60.0,
                 "protocol": "One inherited sustained conductance event; excitation and evoked inhibition scale linearly with contrast."},
        "inhibition": {"assumed_gaba_fraction": GABA_FRACTION, "retained_glycine_fraction": 1.0 - GABA_FRACTION,
                       "ctrl_gaba_remaining": 1.0, "gaba_condition_gaba_remaining": 0.0,
                       "identities": "Transmitter identities are unmeasured; effective conductance at every inherited inhibitory contact is split into co-located GABA/glycine components with inherited identical kinetics."},
        "calcium": {"mode": "fixed", "v_half_mV": -40.0, "slope_mV": 4.0, "tau_ms": 50.0,
                    "description": "Instantaneous activation saturation followed by linear clearance; common resting influx baseline, no condition-specific peak normalization; no driving-force, indicator, release, or membrane-feedback transformation."},
        "local_activation": {"source_segment": source, "source_xyz_um": coords[source].tolist(),
                             "exc_weight_uS": float(local_weight), "fraction_of_unblurred_100pct_input": 0.1,
                             "inhibitory_reversal_mV": rest,
                             "protocol": "A single weak local excitatory conductance with fixed tonic pure-shunting inhibitory load; separate isolation protocol from the evoked contrast experiment.",
                             "normalized_maps": "Per-segment peak response / source peak response; additional at_source_peak maps are simultaneous values at the source peak time."},
        "summary_conventions": "Mean of local peak responses at the same 75 OFF input sites, plus length-weighted means at those sites and over the whole cable; calcium transformed locally before peak and spatial averaging.",
        "spot_summary": spot_summary, "local_summary": local_summary,
        "input_hashes": input_hashes, "config": cfg,
        "locked_parameters": model["p"], "elapsed_seconds": time.time() - started,
        "hashes": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                   for path in [Path(__file__), ROOT / "scripts/fig3/run_fig3_pilot.py", ROOT / "src/calcium_readout.py",
                                ROOT / "data/cid3.swc", ROOT / "data/sm2nrn_cid3.mat"]},
    }
    ctrl_local, gaba_local = local_summary
    metadata["local_effects"] = {
        "source_peak_gain": gaba_local["source_peak_dv_mV"] / ctrl_local["source_peak_dv_mV"],
        "peak_transfer_weighted_cable_gain": gaba_local["peak_transfer_weighted_cable_um"] / ctrl_local["peak_transfer_weighted_cable_um"],
        "same_time_transfer_weighted_cable_gain": gaba_local["same_time_transfer_weighted_cable_um"] / ctrl_local["same_time_transfer_weighted_cable_um"],
    }
    metadata["spot_effects"] = []
    for contrast_pct in CONTRAST_PCT:
        ctrl = next(r for r in spot_summary if r["condition"] == "ctrl" and r["contrast_pct"] == contrast_pct)
        gaba = next(r for r in spot_summary if r["condition"] == "gaba" and r["contrast_pct"] == contrast_pct)
        metadata["spot_effects"].append({
            "contrast_pct": contrast_pct,
            "off_voltage_peak_gain": gaba["mean_off_peak_dv_mV"] / ctrl["mean_off_peak_dv_mV"],
            "off_calcium_peak_gain": gaba["mean_off_peak_ca_au_ms"] / ctrl["mean_off_peak_ca_au_ms"],
            "off_voltage_peak_increment_mV": gaba["mean_off_peak_dv_mV"] - ctrl["mean_off_peak_dv_mV"],
            "off_calcium_peak_increment_au_ms": gaba["mean_off_peak_ca_au_ms"] - ctrl["mean_off_peak_ca_au_ms"],
        })
    np.savez_compressed(args.outdir / "anatomy_maps.npz", **maps)
    (args.outdir / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    validation = {"passed": bool(all(check["passed"] for check in checks)), "checks": checks}
    (args.outdir / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    print(f"Saved {args.outdir}; all checks passed={validation['passed']}; elapsed={metadata['elapsed_seconds']:.1f}s", flush=True)
    if not validation["passed"]:
        raise SystemExit("One or more numerical or protocol checks failed; inspect validation.json.")


if __name__ == "__main__":
    main()
