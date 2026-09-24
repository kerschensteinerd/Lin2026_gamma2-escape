# Calcium sweep results and recommended interpretation

The completed sweep supports broad but conditional dependence of contrast gain on the voltage-to-calcium transformation. The qualitative preference survives the older, shallower sigmoid, so it does not require the exact selected curve. The preference weakens and reverses in identifiable parts of the explored domain. The approximately 3.3-fold horizontal contrast shift is much less sensitive to this transformation than fractional contrast preference is.

## Work completed

- Saved the full [sweep protocol](protocol.md) before executing the calcium grid.
- Evaluated 2,601 regular combinations: midpoint −45 to −20 mV in 0.5-mV steps, slope 1 to 6 mV in 0.1-mV steps. Evaluated the exact original pair separately.
- Reprocessed the archived local voltages for control and 20%, 25%, and 30% retained inhibition at all 11 contrasts. This produced 114,444 regular-grid condition–contrast mean AUCs, plus 44 at the exact reference point.
- Recomputed focal calcium spread for control and nominal 25% retained inhibition using the available complete local voltage maps.
- Generated four review figures, vector SVGs, machine-readable results, a full-contrast table with matching/exclusion reasons, and an eight-panel supplemental figure plan.
- Preserved the manuscript, historical source code, electrical simulations, calibration, and existing figures. No electrical simulation or calcium refit was performed.

## Main results

Define gain at contrast C as the reduced/control ratio of absolute model calcium AUC. Define P15 = gain(15%)/gain(70%) and P30 = gain(30%)/gain(70%). Values above one indicate preferential fractional enhancement at the lower contrast. These definitions do not imply a larger absolute response increment at low contrast.

| Retained inhibition | Both P15 and P30 > 1 | P30 > 1 | Positive normalized 30% separation | Calcium-only mean shift range |
|---|---:|---:|---:|---:|
| 20% | 1,692/2,525 | 2,013/2,555 | 2,147/2,594 | 3.773–4.329 |
| 25% | 1,561/2,525 | 1,835/2,555 | 1,966/2,594 | 3.115–3.486 |
| 30% | 1,457/2,525 | 1,708/2,555 | 1,829/2,594 | 2.694–2.965 |

Denominators count technically usable grid points for that metric. All rows refer to the same specified rectangular parameter grid. These are deterministic grid-coverage summaries, not biological frequencies, probabilities, confidence intervals, or estimates of physiological parameter distributions. Strict counts are unchanged when a 1e−10 directional margin is required; the nearest fractional-preference boundary lies more than 1e−4 from equality.

At nominal 25% retained inhibition, 1,561/2,525 usable combinations satisfy both preference criteria (61.8% of usable points; 60.0% of the entire grid, with 76 technically unavailable). Intermediate preference alone persists at 1,835/2,555 usable points. All 2,601 combinations increase absolute calcium AUC at every nonzero contrast; this is expected from stronger depolarization followed by a monotone readout and provides limited independent evidence for contrast selectivity.

### Dependence on the selected activation curve

| Calcium setting | Gain at 15% | Gain at 30% | Gain at 70% | P15 | P30 |
|---|---:|---:|---:|---:|---:|
| Exact original (−37.479 mV, 1.398 mV) | 17.995 | 11.427 | 1.584 | 11.362 | 7.215 |
| Historical reference (−40 mV, 4 mV) | 3.412 | 2.722 | 1.581 | 2.158 | 1.723 |
| Domain corner (−20 mV, 6 mV) | 4.309 | 5.194 | 5.997 | 0.719 | 0.866 |

The historical curve retains the qualitative mechanism but produces appreciably smaller low/intermediate fractional gain than the selected steep curve. The late, broad activation corner reverses both preference measures while retaining positive absolute enhancement. The separate corner (−20 mV, 1 mV) has technically weak control responses; it is not a usable ratio counterexample. The example sheet displays all four domain corners as a complete diagnostic set, alongside the exact and historical references, rather than selecting only supporting examples. This corner-set addition was made after the full maps revealed the reversal region; it did not alter the protocol, grid, or calculations.

### Spatial spread

At nominal retained inhibition, both first 1/e crossings are observed at 2,574 grid points, and reduced inhibition broadens calcium spread at all of them. The reduced/control crossing ratio ranges from 1.261 to 9.673 across those points. At 27 additional points the reduced profile never crosses 1/e within the sampled arbor; these are recorded as unobserved crossings, not zero, failed broadening, or a maximum-distance value.

The exact reference reproduces 18.777 → 34.337 µm. The historical curve gives 24.980 → 41.068 µm. These calcium descriptors change while the underlying electrical fields remain unchanged. Source normalization and calcium saturation make the distance amplitude-dependent; it is not a newly fitted cable length constant. Spatial results apply only to the original focal site, stimulus strength, and 25% retained inhibition. They establish no multi-site or multi-morphology robustness.

### Contrast equivalence

The nominal calcium-only mean horizontal shift spans 3.115–3.486 over the 2,590 points with an available match. The exact reference gives 3.269 and the historical curve gives 3.314. The mean may use different sampled contrasts as overlap changes, so the complete match support is saved. Restricting the summary to a common set of 5%, 10%, and 15% contrasts yields 3.179–3.486 over 2,528 grid points, with median 3.283.

These calcium-only summaries differ by definition from the manuscript's approximately 3.3 combined voltage/calcium marker. Their stability is consistent with an electrical horizontal rescaling transmitted through a shared monotone transformation. It is not an independent validation of the calcium mechanism or the behavioral comparison.

## Numerical and reporting safeguards

The calcium activation is applied locally before averaging. Clearance remains 50 ms, rest −45 mV, and the AUC window is 50–3050 ms. Backward linear-filter weights and aggregation of exactly repeated stored voltage values accelerate the calculation without a steady-state approximation or voltage binning. All source voltages are promoted from archived float32 to float64 for postprocessing.

All 12 production checks and 114 independent checks passed. Independent full time-series filtering at 13 parameter coordinates across all 44 condition/contrast traces agrees with the accelerated AUC to a maximum scaled error of 1.435×10⁻¹⁴. Original higher-precision output summaries are reproduced within 2.54×10⁻⁷ scaled error, consistent with float32 voltage storage. Independent checks cover all derived masks and scalar Hermite/PCHIP matches, the 7,806 parameter/fraction rows, sampled focal metrics, and source/manuscript hashes. See [independent validation](validation_independent.md).

Ratios and endpoint normalizations use an AUC denominator/span threshold of 0.15 arbitrary influx units × ms² (10⁻⁶ × clearance × stimulus duration). This is a technical reporting threshold, not a measured detection limit. Raw AUCs are retained. Nominal counts for both preferences are 1,561/2,583 at threshold 0.0015, 1,561/2,525 at threshold 0.15, and 1,551/2,355 at threshold 15. No choice of threshold is used to optimize the conclusion.

Spatial source normalization uses 0.00005 arbitrary influx units × ms. All source responses pass this threshold. Unobserved 1/e crossings remain missing with an explicit reason. Full-contrast inverse PCHIP never extrapolates, never forces monotonicity with jitter, and records exclusion reasons. No prior calibration objective or high-contrast penalty was applied during this sweep; discrepancy from old targets is saved only as descriptive context.

The execution report adds full-contrast tables, activation at resting/maximal archived voltage, absolute cable-weighted calcium, and common-support shift summaries to the initially saved core output. Those are postprocessing of the fixed completed grid. The exact original point is excluded from regular-grid coverage counts.

## Recommended supplemental presentation

Use [the overview](review/calcium_sweep_overview.png) as the current review sheet and [the detailed figure plan](supplemental_figure_plan.md) to assemble the eventual supplement. The plan reserves its calcium section for completed results and clearly marks electrical, focal-site, and inhibitory-layout sections as future work.

The result worth emphasizing is that preferential fractional enhancement occurs across a continuous activation regime, persists with the historical shallower curve, and has explicit boundaries. The figures should show that regime and its counterexamples. Keep the absolute-response sheet available because large gain ratios can reflect very small control signals.

A possible Results sentence, subject to the final figure layout, is:

> With electrical responses held fixed, preferential fractional enhancement at low and intermediate contrasts persisted across a range of shared calcium activation functions, including a shallower reference function, but weakened or reversed as the activation midpoint and slope changed. Thus, the contrast dependence did not require the exact calibrated transformation, although its magnitude and direction depended on the activation regime.

This text does not imply robustness to electrical parameters, additional morphologies, or input distributions. The planned electrical and spatial studies remain necessary for those claims.

## Files and reproducibility

- [Protocol](protocol.md): completed calcium protocol and exact future electrical/spatial designs.
- [Supplemental figure plan](supplemental_figure_plan.md): panel architecture, actual findings, and draft caption.
- [Overview SVG](review/calcium_sweep_overview.svg) and [PNG](review/calcium_sweep_overview.png).
- [Residual-inhibition sensitivity SVG](review/calcium_sweep_residual_sensitivity.svg) and [PNG](review/calcium_sweep_residual_sensitivity.png).
- [Absolute responses SVG](review/calcium_sweep_absolute_responses.svg) and [PNG](review/calcium_sweep_absolute_responses.png).
- [Original settings and all domain corners SVG](review/calcium_sweep_examples.svg) and [PNG](review/calcium_sweep_examples.png).
- [Core NPZ arrays](results/calcium_sweep.npz), [parameter metrics TSV](results/calcium_sweep_metrics.tsv), [summary](results/summary.json), [run manifest](results/run_manifest.json), and [input hashes](results/input_hashes.json).
- [Full contrast responses and matching reasons](review/contrast_responses_and_matches.tsv), [activation and spatial metadata](review/activation_spatial_metadata.tsv), and [review manifest](review/review_manifest.json).
- [Production runner](run_calcium_sweep.py), [plot/export builder](build_sweep_outputs.py), and [independent validator](independent_validator.py).

From the workspace root, reproduce core data into a new directory (the runner refuses a nonempty output directory):

```bash
vg3_cid3_rf_simulation-main/.venv/bin/python analysis/parameter_sweeps_20260921/run_calcium_sweep.py --outdir analysis/parameter_sweeps_20260921/reproduction_01
```

The plot/export builder and independent validator read the canonical `results` directory. Existing validated production data do not need to be recomputed to rerender the plots.
