# Calcium analysis panels

`plot_calcium_panels.py` reads the archived calcium transformation sweep and produces 14 standalone visualizations in PNG and SVG, with their numerical source data. It does not run simulations or fit parameters.

From `modeling/`, run:

```sh
python analysis/parameter_sweeps_20260921/calcium_supplement/plot_calcium_panels.py
```

The default input is `../results/calcium_sweep.npz` relative to this directory; outputs go to `../generated/calcium_panels/`. Use `--results PATH` to select a results directory and `--outdir PATH` to select an output directory.

## Outputs

- `calcium_response_selected` and `calcium_response_domain_corner`: calcium AUC across contrast for the exact selected transformation and the domain corner at midpoint −20 mV and slope 6 mV. Control is gray and 25% retained inhibition is orange; both curves are solid. Markers identify the 15%, 30%, and 70% contrasts. The two plots retain their different, labeled y-axis ranges.
- `retained_{20,25,30}_log2_preference15` and `retained_{20,25,30}_log2_preference30`: fractional preference maps, where gain(C) = AUC_reduced(C) / AUC_control(C), P15 = gain(15%) / gain(70%), and P30 = gain(30%) / gain(70%). Positive log2 preference indicates preferential fractional enhancement at the lower contrast. All maps share the −7 to 7 color range; black contours mark zero.
- `retained_{20,25,30}_shape30_difference`: reduced-minus-control difference at 30% contrast after each response curve is separately normalized between its 15% and 70% endpoints. These maps share the −1 to 1 color range; black contours mark zero.
- `spatial_spread_difference_um`: change in the first calcium 1/e crossing after inhibition reduction at the original focal input site, using the 0–80 µm color range. Missing reduced-condition crossings are censored, shown in gray, and never assigned zero or a maximum distance.
- `calcium_spatial_profile_selected`: source-normalized calcium profiles for the selected transformation, with the 1/e reference and both conditions identified.
- `calcium_contrast_equivalence`: mean equivalent-control/reduced-contrast multiplier, using the 3.1–3.5 color range. Dotted contours mark valid-match count boundaries.

The three 25% retained-inhibition preference/shape maps are emitted once each. Every map includes its own colorbar and marker/missing-value key. Stars identify the exact selected parameter pair; open squares identify the domain corner. The full prescribed domain is retained, including directional reversals and unavailable metrics.

`source_data/` contains exact plotted arrays in `plotted_arrays.npz`, parameter-map and example-curve/profile TSV files, and `plot_mapping.json`. Raw responses and individual contrast-equivalence matches are included so plotted summaries can be reconstructed. Array names describe metrics and retained-inhibition conditions. `qa/plot_validation.json` records exact artist-array comparisons, text-bound checks, and verification that the source archive is unchanged. The generated `manifest.json` records output hashes and dependency versions.

## Interpretation

The grid contains 2,601 midpoint/slope pairs: −45 to −20 mV in 0.5-mV steps and 1 to 6 mV in 0.1-mV steps. The exact selected pair is evaluated separately. Calcium conversion uses a fixed driving factor and 50-ms clearance and precedes averaging across 75 OFF input contacts on one reconstruction. These sites and parameter combinations are not independent biological replicates. Electrical responses, geometry, excitation, inhibitory conductance, and input locations remain fixed within each archived condition; no parameters are recalibrated.

AUC is integrated over the 3-s stimulus, expressed in arbitrary influx units × ms², and divided by 1,000 in the response plots. It is a model proxy for calcium, not a calibrated concentration or fluorescence signal. Separately endpoint-normalizing each condition does not imply equal absolute responses. Unavailable preference/shape metrics reflect a control denominator or normalization span at or below the 0.15 arbitrary-unit × ms² safeguard.

Spatial profiles use cable-length-weighted 5-µm path-distance bins. First crossings are interpolated between adjacent bins, not fitted exponential length constants. The analysis uses the original focal site (compartment 4705) with identical focal conductance in both conditions and tests the calcium transformation, not robustness across input locations. At 25% retained inhibition, 27 grid pairs have no reduced-condition crossing within the sampled arbor.

Calcium-only horizontal equivalence inverts each full control AUC curve using PCHIP without extrapolation. Only technically usable reduced-condition responses within the simulated control range contribute. The valid contrast points can differ across parameter pairs; means and parameter-grid ranges are descriptive, not uncertainty intervals or a fitted common biological gain.
