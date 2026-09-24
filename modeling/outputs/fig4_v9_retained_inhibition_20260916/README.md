# Figure 4 v9 retained effective inhibition

Nominal retained total effective inhibition is 25%; sensitivity cases are 20% and 30%. The control is unchanged. No activation parameter, focal conductance, or other parameter was refitted. These are deterministic samples from one reconstructed model, not biological replicates.

`nominal25/figure_data.json`, `anatomy_maps.npz`, and `finite_local_maps.npz` supply A–E/I. Condition keys are `ctrl` and `reduced`. A–E have shared absolute scales; only I uses separate 15%/70% endpoint normalization. `sensitivity_summary.json` and `contrast_sensitivity.tsv` report all prescribed cases; `provenance.json`, `input_hashes.json`, `stimulus_hashes.json`, the copied runner, and `validation.json` record reproducibility and targeted checks.

The unchanged legacy protocol function uses an argument called `glycine_fraction` as a numerical adapter for the retained total load. This version does not assign a universal glycinergic identity to that load. Newly recomputed voltages at each fraction feed the exact prior calcium transformation.

Reproduce with `.venv/bin/python scripts/fig4/run_retained_inhibition_v9.py --outdir <new-empty-directory>`. The runner refuses to overwrite nonempty outputs.
