# VGluT3/cid3 compartmental modeling

This folder contains the local modeling snapshot used for the Lin2026 gamma2
escape study: a passive multicompartment NEURON model of one reconstructed
VGluT3 amacrine cell, its morphology and synapse inputs, fitted parameters,
the Figure 4 numerical pipeline, and accompanying parameter sweeps.

The original modeling code is by Pratyush Ramakrishna and retains its
[MIT license](LICENSE). The original repository README is preserved as
[README_original.md](README_original.md); its example settings and figure
references predate the current paper analysis. No Git history was present in
the local source snapshot.

## Install

From the top level of this repository:

```bash
cd modeling
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-tested.txt
```

The recorded environment was Python 3.11.15, NEURON 9.0.2, NumPy 2.4.6,
SciPy 1.17.1, and Matplotlib 3.11.1. The paper scripts use `numpy.trapezoid`
and therefore require NumPy 2 or later. `requirements.txt` also retains the
original unpinned dependency list. On a headless machine, set `MPLBACKEND=Agg`
before running plot-producing scripts.

## Figure 4 model

The paper uses full effective inhibition for control and 25% retained
effective inhibition for the reduced condition, with 20% and 30% sensitivity
cases. These fractions specify total inhibitory load; they do not assign
transmitter identities. The activation curve, electrical parameters,
geometry, and focal conductance remain fixed across these conditions.

Run from `modeling/` with the environment activated:

```bash
python scripts/fig3/validate_calcium_readout.py
python scripts/fig4/run_retained_inhibition_v9.py \
  --outdir outputs/reproduction_fig4
```

The second command runs NEURON simulations and can take substantial time.
It refuses to overwrite a nonempty output directory. Use a fresh output
directory for each run. It also checks reproduction of the historical control,
geometry, calcium readout, input weights, and time-step convergence.

Required historical calibration inputs are included in `outputs/fig3_v3/`
and `outputs/fig3_v3_spatial/`. The names retain their original lineage;
the current paper analysis is Figure 4. Numerical reference results are in
[`outputs/fig4_v9_retained_inhibition_20260916/`](outputs/fig4_v9_retained_inhibition_20260916/README.md):

- `nominal25/figure_data.json`: contrast curves and focal spatial profiles.
- `nominal25/anatomy_maps.npz` and `finite_local_maps.npz`: compartment maps.
- `nominal25/contrast_voltage_traces.npz` and
  `sensitivity_voltage_traces.npz`: local voltage traces for calcium analyses.
- `provenance.json`, `validation.json`, and the sensitivity tables: recorded
  settings, checks, and reference results.

Saved provenance records retain historical filenames and workspace paths.
Those records describe the original runs. `IMPORT_MANIFEST.json` records the
files in this portable import and their source hashes.

## Parameter sweeps

See [analysis/parameter_sweeps_20260921/README.md](analysis/parameter_sweeps_20260921/README.md)
for calcium activation, electrical-parameter, focal-site, and inhibitory-layout
analyses. The calcium analyses supply the numerical basis for Figures S6/S7.
The package includes archived numerical results so they can be inspected
without rerunning the full simulation grids.

## Original examples and historical figure scripts

The original example and inhibition-removal sweep remain available:

```bash
python scripts/example/run_example.py --stim spot --outdir outputs/example
python scripts/sweep/run_inhibition_sweep.py --outdir outputs/sweep
```

These examples use their original protocols; use the Figure 4 command above
for the current paper settings. `src/` contains the shared model and calcium
implementation, and `data/` contains the SWC morphology, synapse mapping,
synapse classification lists, and fitted parameters.

Historical `scripts/fig3/` and `scripts/fig4/` figure assembly, experimental
panel, and Illustrator helpers are retained as source history. They reference
separate manuscript workbooks, figure files, or Adobe Illustrator and do not
form a self-contained artwork pipeline here. The numerical Figure 4 runner
and the documented sweep analyses use the bundled model inputs.

Local virtual environments, caches, temporary logs, obsolete candidate runs,
and manuscript artwork are excluded from this import. Model parameters and
archived numerical arrays are preserved.

## Import validation

On 24 September 2026, the copied model passed the seven existing calcium
readout checks and 52 additional input, geometry, and reproduction checks
while running from outside the repository directory. Its 12,398 compartments
and 30 ON, 75 OFF, and 242 inhibitory contacts matched the source inputs.
Archived geometry and focal voltage/calcium maps were reproduced; sampled
70% contrast traces at all OFF sites through 148 ms agreed with the archived
float32 traces within 1.90 × 10⁻⁶ mV.

These checks reused the recorded dependency environment. A fresh dependency
installation, the full NEURON sweep grids, and historical Illustrator assembly
were not rerun for this import.
