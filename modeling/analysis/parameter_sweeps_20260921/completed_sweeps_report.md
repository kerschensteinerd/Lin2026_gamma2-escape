# Completed electrical, spatial, and calcium-supplement work

21 September 2026. The electrical and spatial controls in sections B and C of the saved protocol are now completed. Their results strengthen the spatial-broadening conclusion and show that preferential fractional enhancement depends on the electrical operating regime. The calcium supplement is ready for review, with provisional figure numbering.

## Electrical sweep

The complete design contains 10,692 spot condition–contrast responses and 972 focal conditions. The 9 × 9 grid varies full-control inhibition and passive leak independently from 0.25× to 4× their reference values. It crosses shared excitation strengths of 0.5×, 1×, and 2× with control and 20%, 25%, and 30% retained inhibition. Geometry, calcium parameters, and all other settings are fixed. No parameters were refitted.

Define P15 and P30 as the reduced/control calcium-AUC ratio at 15% or 30% contrast divided by that ratio at 70%. The table gives electrical grid points where **both** preferences exceed 1. All 81 points in each row/column combination have usable denominators.

| Shared excitation | 20% retained inhibition | 25% retained inhibition | 30% retained inhibition |
|---|---:|---:|---:|
| 0.5× | 45/81 | 42/81 | 36/81 |
| 1× | 63/81 | 63/81 | 63/81 |
| 2× | 81/81 | 81/81 | 81/81 |

At the reference excitation and 25% retained inhibition, the effect persists at 63/81 electrical points. P15 ranges from 0.344 to 18.383 and P30 from 0.609 to 7.254: the broad parameter domain includes reversals. At 2× excitation, the smallest P30 is only about 1.0013, so a consistent sign does not imply a large effect everywhere.

All 646 electrical-grid comparisons with paired observed calcium 1/e crossings broadened. The other 83 reduced-condition crossings were not observed within the available arbor; all control crossings were observed. At reference excitation and 25% retained inhibition, this is 72/72 observed comparisons plus nine censored reduced crossings. Cable-integral descriptors retain information for censored cases. Absolute calcium AUC increased at every nonzero contrast, as expected for this passive disinhibition model with a monotone calcium transform.

The nominal point reproduces the original result: calcium first-crossing distance increases from 18.777 to 34.337 µm; P15 = 11.362 and P30 = 7.215. Calcium-only horizontal equivalence is 3.2685 at that point. Across the nominal-excitation, 25%-retained grid, the mean equivalence ranges from 2.600 to 3.780, with three to five valid sampled contrast matches. Those means have unequal support and are not fitted common gains.

The full exploratory grid reaches a maximum local voltage of −4.473 mV. This is inside the model's reversal bounds but emphasizes that the sweep spans strongly depolarized regimes; the grid is not a measured physiological distribution.

Detailed results, figures, data, and reproduction instructions: [electrical/README.md](electrical/README.md).

## Spatial controls

The site analysis contains 900 focal conditions: all 75 unique OFF input sites, three shared conductance strengths, and four inhibitory conditions. Among the resulting 675 reduced/control comparisons, **669/669 observed calcium 1/e crossings broadened**. Six reduced crossings were censored (five at 20% and one at 25% retained inhibition, all with 2× focal drive). The normalized calcium cable integral increased in **675/675** comparisons, including these six cases.

At the original drive and 25% retained inhibition, all 75 sites broadened; the median first-crossing increase was 16.75 µm, with a range of 2.13–34.43 µm. These are descriptive differences across locations on one reconstruction, not independent cell measurements.

The layout analysis contains 770 spot responses and 210 canonical-source focal conditions. At each retained fraction, it compares uniform weakening, 20 seeded random retention orders, proximal-first retention, and distal-first retention. Each layout preserves exactly the same total inhibitory conductance. Original contact weights are never increased; nonuniform layouts have at most one partially retained boundary contact.

All **69/69 reduced layouts** retain both lower-contrast preferences at the reference electrical settings. At 25% retained inhibition, P15 ranges from 11.35 to 22.09 and P30 from 7.07 to 7.40. Thus uniform weakening is not necessary for the qualitative contrast preference, although the magnitude depends on layout. The canonical focal calcium crossing also broadened in all 69 layouts at the original focal drive. The multi-site and multi-layout controls are separate designs; they do not constitute a full site × layout factorial study.

Detailed results and figures: [spatial/README.md](spatial/README.md).

## Calcium supplement

The package contains an eight-panel main figure and a companion page comparing 20%, 25%, and 30% retained inhibition. It presents absolute responses, preferential fractional enhancement, normalized contrast-response shape, focal calcium spread, calcium-only horizontal equivalence, and a usable parameter regime where contrast preference reverses. The main figure uses the existing 2,601-point calcium grid plus separately evaluated exact selected parameters; it does not mix the new electrical sweep into the calcium maps.

- [Main figure, PDF](calcium_supplement/Figure_Sx_calcium_sweep.pdf)
- [Retained-inhibition companion, PDF](calcium_supplement/Figure_Sx_calcium_retained_inhibition.pdf)
- [Full figure legend](calcium_supplement/legend.md)
- [Proposed Results, Methods, and Discussion text](calcium_supplement/proposed_manuscript_text.md)
- [Source data, exports, and reproduction guide](calcium_supplement/README.md)

SVG and PNG exports accompany both PDFs. The selected Figure 4 colors and Arial typography are retained. Final supplemental numbering remains to be assigned. The existing manuscript and existing figures were not edited.

## Validation and numerical execution

The numerical acceleration is documented in [execution_audit/execution_scope.md](execution_audit/execution_scope.md). Each spot run retains the original NEURON onset transient at 0.1-ms integration and 2-ms recording. Only a converged late voltage tail is continued at its last value, after satisfying a conservative bound measured across every compartment. The local nonlinear calcium transform and clearance recurrence continue through the original 50–3050-ms AUC window.

Five electrical spot reference checks and eight layout spot reference checks used the full original simulation duration. The largest observed scaled electrical AUC discrepancy was 8.28 × 10⁻¹⁴; layout checks had zero mean-AUC discrepancy at numerical precision. Five direct electrical focal checks and six full-duration spatial focal checks verified the impedance-derived voltage maps; all maximum voltage errors were below 3.2 × 10⁻¹⁰ mV. The exact model trajectory and numerical bounds are retained in the output archives.

All **244 independent audit checks passed**, and all **82 protected original files** retained their starting hashes. The independent root audit reconstructs every inhibitory layout and every exported spatial focal metric, checks all matrix dimensions and original-source hashes, confirms cross-analysis agreement at shared conditions, and independently recomputes calcium clearance from saved NEURON prefixes. Its executable record and current result are [independent_execution_validation.py](execution_audit/independent_execution_validation.py) and [independent_execution_validation.json](execution_audit/independent_execution_validation.json). Both calcium PDFs were rendered and visually inspected; their 22 plotting checks passed.

## Suggested manuscript interpretation

“Across a sweep of inhibitory load and membrane leak, reduced inhibition broadened modeled calcium spread wherever first crossings were observed, while preferential fractional enhancement of weaker inputs depended on the shared excitation strength and electrical operating regime. At the reference excitation strength, the lower-contrast preference persisted at 63 of 81 electrical parameter combinations. The spatial effect also persisted across all 75 tested input locations, and the contrast preference survived all 69 tested residual-inhibition layouts at matched total conductance. These controls support the proposed mechanism across multiple model settings while identifying the conditions under which its contrast preference weakens or reverses.”

Keep the calcium figure as Figure Sx and develop the electrical maps plus compact site/layout panels as a separate related supplement, using the completed review plots as source material. Grid coverage, input sites, and random layouts describe sensitivity within one deterministic reconstruction; they are not biological probabilities, confidence intervals, or additional biological replicates.
