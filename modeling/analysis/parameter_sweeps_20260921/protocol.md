# Figure 4 compartmental-model robustness protocol

Protocol specified 21 September 2026, before the production calcium sweep. This is a computational sensitivity analysis of one deterministic reconstruction, not a biological parameter-estimation study. The calcium sweep is the analysis to execute now; the electrical and spatial studies below are planned follow-up work. Any implementation departures must be identified in the run report rather than silently replacing this protocol.

## Questions and scope

1. Does a shared voltage-to-calcium transformation yield preferential fractional enhancement at low/intermediate contrast across a broad range of activation midpoints and slopes?
2. Does reduced inhibition broaden electrical and calcium spread and increase contrast sensitivity across different inhibitory loads, membrane leaks, and shared excitation strengths?
3. How much do the results depend on the illustrated focal input site, its strength, and uniform weakening of all inhibitory contacts?

The key precedent from the Ramakrishna study is to display parameter dependence and structural controls together. Its practice of adjusting synaptic strengths during a length-constant sweep addresses robustness under an amplitude constraint. The primary analysis here instead keeps excitation fixed within each comparison and does not recalibrate calcium parameters or rescale input to rescue responses. A later amplitude-matched analysis, if performed, must be labeled separately.

## Frozen reference model and pairing

Use `vg3_cid3_rf_simulation-main/outputs/fig4_v9_retained_inhibition_20260916/` and its provenance; historical `fig3` filenames describe source lineage, not the current figure number.

| Quantity | Reference setting |
|---|---|
| Morphology | One cid3 EM reconstruction; 12,398 compartments; 1,546.07 µm cable; uniform 0.5-µm diameter |
| Sites | Same 75 OFF excitatory contacts and 242 effective inhibitory contacts |
| Passive membrane | `g_pas = 1e-4 S/cm²`, axial resistivity 100 Ω·cm, capacitance 1 µF/cm² |
| Reversals | Rest/leak and inhibition −45 mV; excitation 0 mV |
| Full-control inhibition | 193.66351679974264 nS, with the archived contact weights |
| Reduced inhibition | 20%, 25%, or 30% of full-control total conductance; 25% is nominal |
| Excitation at unit contrast | Total spatially weighted OFF conductance 49.164017292788145 nS |
| Spot stimulus | 100-µm diameter, σ=25-µm blur, anatomical soma center, sustained from 50 to 3050 ms |
| Contrast grid (%) | 0, 5, 10, 15, 20, 25, 30, 40, 50, 60, 70 |
| Original electrical integration | 0.1-ms time step; archived voltage sampled every 2 ms, through 3052 ms |
| Selected calcium parameters | Midpoint −37.47892660844856 mV; slope 1.3984142481316566 mV; clearance 50 ms; fixed driving factor |
| Focal source | Compartment 4705; conductance 0.012785038085113628 µS |

The retained-inhibition fractions are fractions of effective conductance, not measured fractions of GABAergic or glycinergic synapses. Each control–reduced pair shares geometry, excitatory contacts and weights, leak, axial resistance, capacitance, reversals, timing, and calcium parameters. Only the specified inhibitory manipulation differs within a pair. A single shared control is reused across the three retained fractions when its other parameters are identical. Changing an electrical parameter creates a new shared control; it does not justify reusing the original control response.

Hash all source arrays and scripts before and after execution. Write new outputs only to this dated analysis directory. Do not modify manuscript, historical model, calibration, or published figure files.

## A. Calcium transformation sweep — execute now

### Exact design

Cross 51 activation midpoints, `−45 + 0.5*j mV` for `j=0,…,50`, with 51 slopes, `1 + 0.1*k mV` for `k=0,…,50`. This is 2,601 parameter combinations over −45 to −20 mV and 1 to 6 mV. Evaluate each combination at every archived contrast for control and each of 20%, 25%, and 30% retained inhibition. There are 114,444 condition–contrast aggregate responses on the regular grid.

Evaluate the exact selected pair separately because it is not a grid point. Also identify the historical illustrative pair (−40 mV, 4 mV), which is on the regular grid. Do not replace the exact selected point with its nearest grid point.

These ranges are broad exploratory choices that cover the selected and historical transformations and sample both early saturation and weak activation. They are not measured priors, confidence bounds, or a claim that every combination is physiologically plausible. Report activation at rest and at the observed voltage extrema to show what the chosen domain implies. The sweep changes the post hoc calcium transformation; it leaves the saved electrical responses unchanged.

### Transformation and computation

For each local absolute voltage trace, use

`a(V) = 1/[1 + exp(−(V−Vhalf)/slope)]`

`dq/dt = a(V) − a(−45 mV) − q/50 ms`, with `q(0)=0`.

Use the same activation function, baseline reference, clearance, and fixed driving factor in all conditions. Preserve signed baseline-subtracted values. Transform each local voltage before averaging across the 75 OFF contacts. Do not transform the spatially averaged voltage.

The primary response is local calcium AUC over the inclusive archived 50–3050-ms window, integrated by the same trapezoid rule as the original readout, then averaged over the 75 contacts. Units are arbitrary influx units × ms². Its absolute value is a model proxy, not a calcium concentration or fluorescence signal.

The production calculation may exploit the linearity of clearance and integration: derive backward quadrature weights for the exact historical recurrence, apply those weights to the baseline-subtracted local sigmoid values, and aggregate identical archived float32 voltages. This must reproduce the original `calcium_trace` recursion on representative traces within a recorded numerical tolerance. The nonlinear activation must still precede spatial averaging. No steady-state approximation may replace the transient AUC calculation.

### Prespecified primary outputs

Let `A_f(C)` be the mean calcium AUC for retained fraction `f`, with full control `f=1`.

1. **Absolute responses and enhancement.** Save `A_f(C)` at every contrast, and `D_f(C)=A_f(C)−A_1(C)`. Summarize the sign and magnitude at 15%, 30%, and 70%; retain the complete curves. Positive enhancement alone is not evidence for preferential low-contrast gain, because it follows readily from increased depolarization plus a monotone activation transform.
2. **Fractional response enhancement.** Compute `R_f(C)=A_f(C)/A_1(C)` for nonzero contrasts with technically usable denominators. Primary contrast-selectivity measures are `P15=R_f(15)/R_f(70)` and `P30=R_f(30)/R_f(70)`. Values greater than 1 indicate larger fractional enhancement at 15% or 30% than at 70%. Map the continuous values and the `P=1` boundary; report whether both exceed 1. Do not reinterpret this as a larger absolute increment at low contrast.
3. **Endpoint-normalized shape.** For each condition independently, `S_f(C)=[A_f(C)−A_f(15)]/[A_f(70)−A_f(15)]`. Save all normalized curves and the interior differences `S_f(C)−S_1(C)` at 20%, 25%, 30%, 40%, 50%, and 60%. Show the 30% difference prominently. The endpoints at 15% and 70% are imposed, not independent predictions. Values outside the endpoint range at other contrasts are retained, not clipped.
4. **Calibration context.** At 30%, save the unweighted sum of squared deviations from the previously used OFF targets: `[S_1(30)−0.0345741]² + [S_f(30)−0.7036578]²`, or their full-precision stored equivalents if available. Label this descriptive discrepancy from prior calibration targets, not a likelihood, independent validation, or a new fitting objective. Do not select a replacement parameter set by its minimum.
5. **Calcium-only horizontal equivalence.** Construct inverse PCHIP from the full control curve `A_1(C)` to contrast using the unchanged contrast grid, without extrapolation. At each nonzero reduced-condition sampled contrast, calculate `C_equiv/C` only when its response lies within the sampled control response range and the inverse is technically defined. Save each match, valid flag, exclusion reason, count, supported contrasts, mean, minimum, and maximum. Do not pool calcium with voltage to make this sweep's shift metric. The earlier Figure 7 pooled voltage/calcium marker can be reproduced separately as a historical check, but is not the primary sweep output.
6. **Spatial calcium spread at nominal retained inhibition.** Recompute steady local calcium from the saved full-arbor focal voltage maps for control and 25% retained inhibition, using `q_inf=50*[a(V)−a(−45)]`. Normalize to each condition's source response solely for spatial-profile descriptors; retain absolute source response and absolute cable-weighted response. Reuse the exact historical source location, cable lengths, path distances, and 5-µm bins. Record the first interpolated 1/e crossing and cable-weighted integral of the source-normalized profile, and their reduced-minus-control differences and ratios. Voltage profiles are invariant across this sweep. Do not derive full-arbor calcium maps for 20% or 30% from binned spatial summary curves; those saved voltage maps are not available in the current v9 archive.

The regular-grid maps and grid-area summaries use only the 2,601 regular points. Separately evaluated selected parameters must not change those denominators.

### Technical flags, not physiological exclusions

Preserve all finite raw AUCs. Use a primary ratio/normalization technical threshold of `1e−6 × 50 ms × 3000 ms = 0.15` arbitrary influx units × ms²: mark a ratio denominator or an endpoint span at or below this value as technically weak. Repeat the relevant counts with factors `1e−8` and `1e−4` (thresholds 0.0015 and 15). These are numerical/display safeguards, not measured detection limits. Report both the number technically evaluable and the number satisfying each directional criterion; show masked locations separately from effects that fail the criterion. Do not substitute masked values with zero or wins.

For focal calcium, use the corresponding peak/steady-state technical scale `1e−6 × 50 = 0.00005` arbitrary influx units × ms for source normalization, and record the threshold explicitly. Preserve unnormalized responses even when normalization is unstable. Missing first crossings are censored/not observed within the arbor, not zero distance; do not replace them with the maximum path length.

PCHIP requires finite strictly increasing control responses. Record nonmonotonicity, numerical ties, absent overlap, and technically weak spans as distinct reasons for an unavailable equivalence. Do not jitter or silently deduplicate to force a fit. Report all individual matches and matched contrast support, because a mean shift can change when the overlap domain changes. Treat comparisons of shift maps with unequal support cautiously.

Directional classifications use strict signs beyond a stated floating-point tolerance. They are descriptive criteria, not significance tests. Retain continuous effects so a tiny positive effect is visibly distinguishable from a substantial one. A regular-grid percentage describes this chosen coordinate box and sampling measure; it is not the probability that the biological mechanism is correct and will change under reparameterization or another range.

### Validation and reporting

- Confirm archived contrast grids, trace shapes, uniform sampling, baseline voltage, and finite values. Record source storage dtype and hash each input.
- Reproduce original selected-parameter AUCs and nominal spatial descriptors, allowing and quantifying archived float32 versus historical higher-precision differences.
- Independently compare the optimized AUC implementation with the full original calcium recurrence at the exact selected pair, historical pair, four sweep corners, and a prespecified interior point (−32.5 mV, 3.5 mV), across all four electrical conditions and available contrasts. Report maximal absolute/relative differences and the tolerance used; investigate discrepancies rather than relaxing it to fit results.
- Confirm zero-contrast responses, shared-control equality across retained-fraction comparisons, finite outputs, and source normalization. Do not regard universally positive enhancement as a surprising test result for a monotone readout.
- Save executable analysis, run configuration, dependency versions, timings, source hashes, machine-readable grid results, per-contrast results/matches, validation results, and vector/raster figures. Describe unresolved numerical issues before biological interpretation.

## B. Electrical sweep — planned, not executed here

### Primary two-dimensional map

Use common multiplier grid `M={2^(k/2): k=−4,…,4}` = approximately `{0.25, 0.353553, 0.5, 0.707107, 1, 1.414214, 2, 2.828427, 4}` on each axis:

- Full-control inhibitory load `G_I = m_I × 193.66351679974264 nS`, i.e. 48.4158792–774.6540672 nS, preserving relative contact weights.
- Passive leak `g_pas = m_L × 1e−4 S/cm²`, i.e. 2.5e−5–4e−4 S/cm², uniformly across compartments.

Start with the 81 combinations at original excitation strength and nominal retained fraction 25%, paired to a newly calculated full-control condition at each point. Run all 11 contrast values and the fixed focal conductance. This is 1,782 condition–contrast responses, plus focal/DC calculations. The nominal pair must reproduce the reference. Axial resistance, diameter, capacitance, morphology, contact locations, reversals, and calcium parameters are fixed at the exact selected settings.

These are exploratory fourfold ranges in both directions, not physiological bounds. A range that reaches a mechanism boundary is informative. Do not delete weak or strongly saturated regions for disagreement with the preferred result.

### Shared excitation and retained-inhibition extension

At the same 81 electrical combinations, cross shared spot-excitation multipliers `{0.5, 1, 2}` with the full control and retained fractions `{0.20, 0.25, 0.30}`. The total complete design comprises `81 × 3 × 4 × 11 = 10,692` condition–contrast responses; the primary map is a subset, not additional duplicate runs. Both members of a pair receive the same excitation multiplier. Cache each shared control once.

For focal spread, use corresponding shared focal-conductance multipliers `{0.5, 1, 2}` on the original 12.785038085113628-nS drive. Do not readjust the drive to recover a target depolarization as leak or inhibition changes. In this passive model the steady transfer/input-resistance approach can avoid repeated long transients, but must be checked against direct steady stimulation at the nominal point and selected corners.

At each electrical pair save voltage and calcium contrast curves, peak local and mean depolarization, original AUC metrics, source input resistance, focal source amplitude, 1/e and cable-integral spread, rest/baseline values, and technical masks. Plot calcium spread change versus contrast-sensitivity change across electrical combinations, with nominal point marked. This covariance is a consistency result within a model, not proof that spread mediates all gain.

Record absolute voltage ranges and activation saturation. Without an independently established biological voltage range, flag values transparently rather than calling them physiologically accepted/rejected. A passive model with reversals −45 and 0 mV should remain between those reversals; violations indicate an implementation problem.

No fit to imaging data or behavior occurs at each point. Do not force a 3.3-fold horizontal match. Report changes in the mechanism and horizontal-equivalence support across the full box. Electrical and calcium sensitivities are separate planned factorial studies; the current calcium sweep does not establish electrical robustness.

## C. Spatial robustness — planned, not executed here

### Focal sites and strength

Stimulate each of the same 75 OFF input contacts in turn, retaining compartment 4705 as the originally illustrated site. Use shared focal-conductance multipliers `{0.5, 1, 2}` and conditions `{control, 20%, 25%, 30%}` at the reference electrical/calcium parameters: `75 × 3 × 4 = 900` focal cases. Recompute path distance from each source. Use the same finite conductance at a given multiplier in both members of every pair; do not tune conductance site by site to a fixed source voltage.

Save paired source amplitude, input resistance, first 1/e crossing, cable integral, and absolute cable-weighted calcium response. Show every site's paired change, with site order or color indicating soma path distance and branch location. Report median, full range, and descriptive interquartile range, but no SEM, biological confidence interval, or paired biological significance test. If an OFF contact is repeated at one compartment, retain contact identity for trace comparability and additionally report the number of unique physical sites.

### Residual-inhibition distribution at matched total conductance

At each retained fraction `{0.20, 0.25, 0.30}`, compare 23 prespecified distributions over the existing 242 inhibitory contacts:

1. Uniform weakening: every contact weight multiplied by the fraction (current model).
2. Twenty random retention patterns, seeded `0,…,19`, with `numpy.random.Generator(PCG64(seed))` generating a permutation of contact indices.
3. Proximal-first retention and distal-first retention, ordering contacts by path distance from the anatomically identified soma, with compartment index then contact index as tie-breakers.

For each nonuniform pattern, traverse the prespecified ordering and retain each original contact's full weight until the target budget `f × total_control_weight` is reached; use only the remaining budget at the boundary contact and set the rest to zero. This matches total retained conductance exactly, never increases a contact above its full-control weight, and permits at most one partially retained contact. Do not normalize a randomly selected fixed contact count upward; conductance fraction and contact fraction need not agree. Save the contact-level weights, retained count, partial contact, and seed/order.

Keep excitation and full-control inhibition identical for all comparisons. Run the full spot contrast grid at reference electrical parameters (69 distinct reduced layouts plus one shared control = 770 condition–contrast responses). At the canonical focal site run the three shared focal strengths for the same 70 layouts (210 focal cases). Compare individual patterns and their deterministic ensemble ranges with uniform weakening. Random seeds represent sensitivity to layout, not additional biological cells or a null distribution for inference.

If this analysis exposes a strong site-by-layout interaction, a subsequent full site × layout crossing can be specified and labeled an extension. It should not be silently substituted for the prespecified primary spatial studies. A claim of robustness across different reconstructed cells still requires additional morphologies; none is supplied by resampling sites on cid3.

## D. Secondary model-assumption controls — planned priorities

These are distinct controls rather than silently enlarging the main grids:

1. **Readout physics:** at the exact selected pair and historical pair, compare fixed driving with the existing ohmic option (`E_Ca=+120 mV`) and GHK option (`[Ca]i=0.0001 mM`, `[Ca]o=2 mM`, 33 °C), each normalized to its common −45-mV reference. Cross clearance `{20, 50, 100, 200} ms` and all retained fractions on the archived voltage traces. Clearance strongly affects absolute integrated units; under sustained stimulation normalized shapes may change little. Do not confuse this with active calcium currents feeding back on voltage.
2. **Inhibitory reversal:** rerun reference electrical comparisons at `E_I={−55, −50, −45} mV`. For hyperpolarizing cases, equilibrate each condition to its actual prestimulus baseline and subtract that baseline in readouts, while preserving the same absolute-voltage activation curve across conditions. Record the condition-dependent baseline instead of artificially re-centering all activation to −45 mV. This requires fresh simulations; it cannot be obtained by shifting the existing voltage traces.
3. **Temporal protocol:** at reference electrical parameters test sustained excitation durations `{50, 200, 1000, 3000} ms`, with a prestimulus baseline and at least five clearance time constants after offset. Report AUC during stimulation and over the same specified post-offset duration, separately. Keep inhibition tonic initially. Contrast-dependent or stimulus-locked inhibition would require an independently justified temporal profile and should be specified before testing.
4. **Cable assumptions:** if the main electrical effects are robust, examine axial-resistivity multipliers `{0.5, 1, 2}` and uniform-diameter multipliers `{0.5, 1, 2}` one at a time at reference leak/inhibition. These controls change transfer and input resistance jointly; they do not replace independent morphologies.

## Interpretation and figure commitments

Display where effects persist, weaken, reverse, or become technically uninterpretable. Mark the selected and historical calcium pairs. Show absolute responses beside fractional or endpoint-normalized summaries. Report all three retained fractions, or show nominal maps with aligned sensitivity maps in an additional page.

The calcium sweep answers whether the qualitative output depends tightly on the originally selected transformation. It cannot validate voltage-gated calcium biophysics, compensate for the model's missing active currents, show physiological parameter identification, or establish a quantitative link to looming behavior. A broad parameter region supporting a sign is stronger than one fitted example, but uniform positive calcium enhancement under a monotone transformation has limited independent evidential value. Spatial broadening plus preferential fractional enhancement is the more informative combined diagnostic.

One supplementary figure should ultimately combine calcium maps, electrical maps, and spatial controls. The current deliverable can contain completed calcium evidence and a clearly labeled layout plan for the electrical/site/layout panels. Do not fill unperformed panels with simulated-looking placeholders or describe the proposed studies as completed.

## Local source chain

- `analysis/fig4_results_v9_20260916/model_update_summary.md`
- `analysis/fig4_results_v9_20260916/check_model_shift_v9.py`
- `analysis/fig4_results_v8_20260916/model_sources.md`
- `vg3_cid3_rf_simulation-main/outputs/fig4_v9_retained_inhibition_20260916/provenance.json`
- `vg3_cid3_rf_simulation-main/scripts/fig4/run_retained_inhibition_v9.py`
- `vg3_cid3_rf_simulation-main/scripts/fig3/protocol_v3.py`
- `vg3_cid3_rf_simulation-main/scripts/fig3/run_anatomy_v3.py`
- `vg3_cid3_rf_simulation-main/src/calcium_readout.py`
- `vg3_cid3_rf_simulation-main/scripts/fig3/finite_local_spread_v3.py`
- `vg3_cid3_rf_simulation-main/scripts/fig3/calibrate_spatial_v3.py`

Paths in this source list are relative to the workspace root `/Users/kerschensteinerd/Documents/g2_loom_cg`.
