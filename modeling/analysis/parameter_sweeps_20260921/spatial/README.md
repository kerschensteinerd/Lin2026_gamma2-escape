# Spatial-control results — 21 September 2026

Completed the prespecified site/strength and conductance-distribution controls on the unchanged Figure 4 v9 model. All calcium parameters and excitation inputs remain shared within each control–reduced pair. No calibration or fitting was performed.

Across 675 reduced-versus-control site/strength comparisons, all 669 observed calcium 1/e crossings were broader. Six reduced conditions had no observed 1/e crossing within the arbor and are explicitly censored. The normalized calcium cable integral increased in all 675/675 comparisons, including those censored cases.

At nominal 25% retained inhibition and the original focal strength, calcium 1/e distance increased at 75/75 sites: median 16.75 µm, range 2.13–34.43 µm. These are sitewise descriptors from one arbor, not biological replicates.

## Matched-total inhibition layouts

For each 20%, 25%, and 30% retained fraction, tested uniform weakening, 20 random retention orders (PCG64 seeds 0–19), proximal-first retention, and distal-first retention. Nonuniform patterns retain full original contact weights until the conductance budget is reached, partially retain only the boundary contact, and set the remainder to zero. No contact is strengthened; the total is matched exactly. Contact counts need not match the retained conductance fraction.

| Retained | Layouts with R15/R70 > 1 and R30/R70 > 1 | R15/R70 range | R30/R70 range | Mean calcium shift range | Valid contrasts per layout |
|---|---:|---:|---:|---:|---:|
| 20% | 23/23 | 18.22–27.17 | 7.63–7.75 | 3.951–5.984 | 2–3 |
| 25% | 23/23 | 11.35–22.09 | 7.07–7.40 | 3.267–4.840 | 3–4 |
| 30% | 23/23 | 7.37–17.63 | 6.34–6.79 | 2.784–4.221 | 3–4 |

Here R15, R30, and R70 are reduced/control calcium AUC ratios at 15%, 30%, and 70% contrast. All 69 reduced layouts enhanced absolute calcium AUC at every nonzero contrast. This is consistent with the monotone calcium readout and pointwise reduced inhibitory conductances; the preferential fractional-enhancement metrics and spatial profiles provide the more informative sensitivity checks.

The horizontal multiplier uses inverse PCHIP of the full-control calcium AUC curve without extrapolation. Each mean summarizes only its valid sampled overlap; full matches, exclusions, and contrast support are saved. It is neither an amplitude ratio nor the earlier pooled voltage/calcium Figure 7 marker.

For nominal 25% retained inhibition, uniform weakening yields a calcium-only mean horizontal multiplier of 3.268536 over 5%, 10%, 15%, and 20% contrast. The smallest layout mean is 3.267108 (random seed 10); distal-first retention gives 4.840160 over 5%, 10%, and 15% contrast. These ranges therefore compare descriptive means with different overlap support, not fits of a common shift over a fixed contrast domain.

## Numerical validation

Focal calculations use the passive DC impedance solution with a finite excitatory conductance. Six independently specified full 3050-ms NEURON simulations matched the complete voltage maps within 2.33e-10 mV. Historical canonical-source descriptors were reproduced.

The 770 spot cases preserve the original 0.1-ms NEURON onset transient, 2-ms sampling, local nonlinear calcium transformation, and 50–3050-ms trapezoidal AUC. A late voltage tail is extended only after a conservative all-compartment passive contraction bound is satisfied; the calcium clearance recurrence is evaluated through the complete stimulus. Maximum declared relative AUC bound was 1.31e-10. Eight independent full 3052-ms original-protocol comparisons had maximum mean AUC relative difference 0. This is a documented numerical acceleration of the protocol, not a replacement of transient AUC by steady activation.

## Files and reproducibility

- `run_spatial.py`: rerunnable focal/spot stages; exact executed snapshots are saved beside it.
- `layout_summary.tsv`, `layout_contact_weights.tsv`, `source_sites.tsv`: complete input identities and weights.
- `focal_all_sites.tsv`, `focal_layout_controls.tsv`, `focal_binned_profiles.tsv`: absolute and normalized focal descriptors.
- `focal_site_paired_effects.tsv`, `focal_site_summary.tsv`, `focal_layout_paired_effects.tsv`: paired changes and censoring.
- `spot_responses.tsv`, `spot_local_responses.tsv`, `spot_layout_metrics.tsv`, `spot_contrast_effects.tsv`, `horizontal_matches.tsv`: all model responses and comparisons.
- `focal_cache/`, `spot_voltage_cache/`: reproducible transfer maps/input resistances and voltage prefixes.
- `focal_validation.json`, `spot_validation.json`, `spot_tail_diagnostics.tsv`, `source_hashes_before.json`: validation and provenance.
- `spatial_sites.pdf/.svg/.png`, `spatial_layouts.pdf/.svg/.png`: editable and preview figures.

Rerun with the repository `.venv/bin/python`, then `run_spatial.py --stage focal`, `run_spatial.py --stage spot`, and `build_outputs.py`. The spot stage imports the preserved `../electrical/passive_tail.py` helper and verifies its hash. Original sources, calibration, figures, and manuscript remain unchanged.

## Interpretation

The specified electrical and calcium settings support the spatial effect across all sampled input locations, and the qualitative contrast-gain redistribution survives the tested residual-contact distributions. The magnitudes depend on site, input strength, and inhibitory layout. Source-normalized 1/e distances are first-crossing descriptors, not fitted exponential length constants. Calcium diffusion and active membrane currents are not modeled. Spatial/gain covariance is not a causal mediation analysis. The 75 sites, 20 random seeds, and parameter conditions do not provide biological confidence intervals or evidence across independent reconstructions.
