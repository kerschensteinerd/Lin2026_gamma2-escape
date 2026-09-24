# Independent calcium-sweep validation

**PASS: 114/114 checks.** No NEURON simulation, parameter fit, manuscript edit, or production-output mutation was performed.

No runner import. Full existing calcium_trace at 13 parameter points, all 44 complete local time-series arrays; every saved derived metric/mask and TSV field independently recomputed; scalar Hermite inverse PCHIP for all parameter points; independent focal spatial measurements at the same 13 points.

The largest direct full-filter AUC error, scaled by max(1, reference AUC), was 1.43e-14. The 51 × 51 grid and additional exact original point, 20/25/30% retained-inhibition axes, 2-ms sampling, 50–3050-ms AUC integration, shared resting voltage and clearance, technical denominator masks, and non-extrapolated contrast matching agree with the declared implementation. All nominal reduced voltages are identical across their two archived copies. Input/manuscript hashes and recorded output hashes match.

The 0.15 AUC floor is a numerical reporting threshold, not a measured biological detection limit. Grid coverage is not posterior probability or confidence. Spatial calculations cover only the archived original focal source and nominal 25% retained inhibition; the repeated spatial fields in each TSV fraction row are explicitly labeled reduced25. Calcium-only equivalence is not the manuscript's combined voltage/calcium marker. This audit verifies the calculation and saved data, not biological validity or robustness to alternative electrical parameters, morphologies, sources, or inhibitory distributions.

Historical AUC values were computed before float32 voltage storage; comparison permits the documented 1e-5 relative tolerance. The sweep and full-trace validation both use the actual archived float32 values promoted to float64.

Detailed checks and checked parameter coordinates are in `independent_validation.json`; reproduce with `vg3_cid3_rf_simulation-main/.venv/bin/python analysis/parameter_sweeps_20260921/independent_validator.py`.

## Review-export addendum

The separate review builder and regenerated exports also pass a scoped independent check. All **85,866 long-form rows** (2,602 parameter points × 3 retained fractions × 11 contrasts) have correct parameter/condition coordinates, AUCs, multiplier blanks, and match statuses: 29,345 valid, 48,142 outside the control response range, 7,806 zero-contrast exclusions, 363 weak reduced responses, and 210 weak control spans. All 2,602 activation/spatial metadata rows have correct activation values, crossing statuses, and shift extrema. The global maximum across both control and sensitivity voltage archives is **−20.679054260253906 mV**.

Direct whole-arbor sums at seven parameter points reproduce the exported absolute calcium cable integrals to 2.55 × 10⁻¹⁶ scaled error. These quantities weight the steady calcium proxy by cable length; they are not volume-weighted calcium content. The common-support descriptor correctly averages only valid 5%, 10%, and 15% contrast matches, retaining 2,536/2,528/2,522 grid points for 20%/25%/30% residual inhibition. Thus it uses a common contrast set but different eligible parameter subsets. The original exact point is excluded from these grid summaries.

Code inspection confirms correct grid transposition and missing-value masking, a shared color range within each metric across retained-inhibition fractions, and all four domain corners alongside the exact original and historical examples. Gray regions denote the stated technical masks or unresolved spatial crossing, not zero response. All 12 recorded review-export hashes, builder hash, and source-archive hash match. This addendum does not repeat the core numerical validation or alter its JSON record; image layout inspection is handled separately by the main task.
