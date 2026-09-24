# Electrical robustness sweep

Completed all **10,692 condition–contrast responses (972 contrast-response curves) and 972 focal cases**, using 324 electrical states and 19 distinct excitation products per state. Duplicate products were reused exactly; all planned design rows remain in the outputs.

The maximum voltage was -4.472577 mV; every result remained within passive reversal bounds. All late-voltage approximation bounds passed. Transient prefixes ended at 68–118 ms; the largest AUC error bound divided by max(1, computed AUC) was 6.54e-10.

The historical 0.1-ms NEURON integration and 2-ms calcium sampling are retained. Exact resting baseline is skipped (validated sample-for-sample). A conservative all-compartment passive contraction bound permits replacement of the late voltage tail by its last recorded value; the original local sigmoid, 50-ms clearance recurrence, and trapezoidal AUC are then evaluated through 3050 ms. This is a declared numerical approximation, not substitution of a steady-state response for stimulus AUC. Five full-duration nominal/corner checks agree to at most 8.28×10⁻¹⁴ in scaled AUC; focal impedance calculations have separate point-process checks.

| Shared excitation | Retained inhibition | Both P15 and P30 > 1 | Resolved calcium spread pairs | P30 range |
|---|---|---|---|---|
| 0.5× | 20% | 45/81 | 75/81 | 0.3979–8.66 |
| 0.5× | 25% | 42/81 | 78/81 | 0.4241–8.09 |
| 0.5× | 30% | 36/81 | 81/81 | 0.4675–7.11 |
| 1× | 20% | 63/81 | 67/81 | 0.7725–7.95 |
| 1× | 25% | 63/81 | 72/81 | 0.6086–7.25 |
| 1× | 30% | 63/81 | 76/81 | 0.5541–6.56 |
| 2× | 20% | 81/81 | 60/81 | 1.0013–7.76 |
| 2× | 25% | 81/81 | 67/81 | 1.0013–6.71 |
| 2× | 30% | 81/81 | 70/81 | 1.0013–5.57 |

All nonzero-contrast reduced-inhibition AUCs exceed their paired controls. This expected monotone-readout effect is distinct from preferential fractional enhancement, which depends on electrical settings. Missing 1/e crossings mean the profile did not reach that level within this arbor; those cases remain in the source data and are never assigned zero spread. Complete per-contrast horizontal-equivalence matches include their support and exclusion reasons.

These exploratory fourfold ranges are not physiological confidence intervals. Strongly depolarized corners approach the 0-mV excitation reversal; they remain visible without claiming physiological plausibility. The sweep keeps morphology, reversal potentials, contact placement, calcium activation, and all within-pair excitation settings fixed. It does not test active calcium feedback, another morphology, or measured biological distributions. Absolute cable-calcium integrals are length-weighted proxy responses, not calcium content.

Reproduce: run `run_electrical.py --workers 6`, then `validate_focal.py`, then `build_electrical_outputs.py`, with the model `.venv/bin/python`. Completed states are resumable only with identical hashed inputs. Raw exact transient prefixes are retained in `transient_traces/`; full-arbor normalized focal voltage maps are retained in `states/`.

Black contours mark zero effect in signed maps: P = 1 in preference panels and zero normalized-shape difference in the shape panel. No contour is drawn if the boundary is absent from the sampled range.
