# Modeling parameter sweeps

This directory publishes the calcium, electrical, and spatial sensitivity analyses used with the Figure 4 model. The source was `g2_loom_cg/analysis/parameter_sweeps_20260921`, alongside the original `vg3_cid3_rf_simulation-main` model. Plotting scripts produce individual panels from the archived numerical results.

Use the Python environment described in [the modeling README](../README.md). Run the following commands from `modeling/`. All paths resolve relative to the scripts, so commands also work from another working directory when given the appropriate script path. The numerical requirements are NumPy, SciPy, Matplotlib, and NEURON; the validated environment used Python 3.11 and NumPy 2.

## Calcium analysis and panel plots

Validate the archived calcium grid without rerunning NEURON:

```bash
python analysis/parameter_sweeps_20260921/independent_validator.py
```

This independently evaluates the full calcium filter at 13 parameter coordinates for all 44 condition/contrast trajectories, recomputes derived metrics and table fields, and checks required numerical input/output hashes. The new report is written under `parameter_sweeps_20260921/generated/validation/`; the historical reports remain unchanged. Eight historical manuscript-file hashes are explicitly reported as external provenance outside this export, rather than counted as successful checks. Missing or changed required numerical files still fail validation. The old production-runner hash is checked against its exact copy in `historical_sources/`; the runnable copy contains portability edits.

Plot individual calcium-sweep panels or review visualizations:

```bash
python analysis/parameter_sweeps_20260921/calcium_supplement/plot_calcium_panels.py
python analysis/parameter_sweeps_20260921/build_sweep_outputs.py
```

These write separate panel plots, tables, and manifests under `generated/`; they leave the numerical reference archives unchanged. Related heatmaps retain shared color limits and include their own color bars and keys.

To recompute the 2,601-point calcium grid and exact reference point from the preserved Figure 4 voltage archives:

```bash
python analysis/parameter_sweeps_20260921/run_calcium_sweep.py --outdir analysis/parameter_sweeps_20260921/reproduction
python analysis/parameter_sweeps_20260921/independent_validator.py --results analysis/parameter_sweeps_20260921/reproduction --report-dir analysis/parameter_sweeps_20260921/generated/reproduction_validation
```

The runner refuses a nonempty output directory. Builders accept `--results` and `--outdir` to use newly generated results. This sweep changes the calcium readout while keeping the archived electrical responses fixed.

## Electrical and spatial analyses

The numerical reference NPZ/TSV files, protocol, summary reports, and historical validation reports are included. The large intermediate voltage-prefix/state caches are regenerable and excluded. Do not attempt to resume the historical run with those caches absent. Start a fresh run using the commands below; the default output directory is `reproduction/` beneath each analysis and the sources remain unchanged.

```bash
python analysis/parameter_sweeps_20260921/electrical/validate_tail.py
python analysis/parameter_sweeps_20260921/electrical/run_electrical.py --workers 6
python analysis/parameter_sweeps_20260921/electrical/validate_focal.py
python analysis/parameter_sweeps_20260921/electrical/build_electrical_outputs.py --results analysis/parameter_sweeps_20260921/electrical/reproduction

python analysis/parameter_sweeps_20260921/spatial/run_spatial.py --stage focal
python analysis/parameter_sweeps_20260921/spatial/run_spatial.py --stage spot
python analysis/parameter_sweeps_20260921/spatial/build_outputs.py --results analysis/parameter_sweeps_20260921/spatial/reproduction
```

The electrical tail validation must run first: the grid's final manifest records its validation file. The electrical and spatial runners/validators accept `--outdir`; use the same directory for their dependent stages. Electrical runs resume only when their input hashes match. The spatial spot stage requires the focal stage's preserved-input manifest. The spatial runner imports the included `electrical/passive_tail.py` helper.

To rebuild review plots directly from the included numerical references, run either electrical or spatial builder without `--results`. Their default output is a separate directory under `generated/`; the required numerical inputs are copied there before derived exports are produced.

These full electrical/spatial grids were not rerun during the repository export. The existing `execution_audit/` reports document the original complete audit. Its unchanged source is at `historical_sources/execution_audit/independent_execution_validation.py`; it requires original intermediate caches, manuscript/artwork inputs, and the original workspace layout, and is therefore an archival record rather than a portable entry point.

## Provenance and scope

`parameter_sweeps_20260921/protocol.md`, `calcium_report_20260921.md`, `completed_sweeps_report.md`, the electrical/spatial READMEs, and saved JSON manifests are historical records. Their original paths and claims about the state of that earlier workspace are preserved. Use this README for current commands and packaging scope. Numerical arrays, result tables, protocol, and historical hash/validation records are copied unchanged. `historical_sources/` holds original numerical scripts; the `run_spatial_*_snapshot.py` files are also historical execution snapshots, not current entry points.

Portable analysis entry points adapt model-path resolution, output-directory arguments/defaults, and provenance/report handling. Plotting entry points export each visualization separately. Parameter grids, numerical formulas, simulation protocols, and archived results are unchanged. The original model-folder prefix in historical hashes is mapped to `modeling/` by the calcium validator.

Excluded items are virtual environments, bytecode, logs, manuscript drafts, generated plot exports, QA screenshots, and intermediate `states/`, `transient_traces/`, `focal_cache/`, and `spot_voltage_cache/` directories. Plotting scripts regenerate derived output. Historical full-audit coverage of excluded files is not asserted for this export.
