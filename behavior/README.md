# Looming behavior analysis

`looming_behavior_analysis.m` analyzes precomputed speed traces for WT and KO mice at 22%, 64% and 91% looming contrast. It produces per-mouse speed heatmaps, reaction-time and maximum-speed summaries, and escape probabilities. Supply a merged MAT input separately; experimental speed data are not distributed in this source folder.

## Run and input format

Use MATLAB R2018b or later in a desktop session. The script uses base MATLAB plotting, interpolation and summary functions, with a local colormap helper; no additional toolbox is required. Add this directory to the path and run:

```matlab
looming_behavior_analysis
```

Choose the merged MAT file in the file dialog. It must contain `HeatmapData`, a struct array with one entry for each genotype/contrast combination:

| Field | Required content |
| --- | --- |
| `Genotype` | `'WT'` or `'KO'` |
| `Contrast` | Numeric value 22, 64 or 91 |
| `heatmap` | Numeric `[mice frames]` matrix of speed, with one row per mouse and a common time axis |

Each of the six genotype/contrast combinations must occur exactly once. Frame 1 is time zero. Defaults assume 30 frames/s and a trace extending through 10 seconds (301 samples at 30 Hz). Speeds must use the same calibrated units as the threshold of 35. Supply precomputed speed traces whose rows represent the intended independent mouse-level observations.

## Analysis and settings

Edit the settings inside the script to match the supplied traces. All four sections use the same loaded dataset and run in this order:

1. **Speed heatmaps.** Plot genotype × contrast panels, interpolate the display time axis fivefold with shape-preserving cubic interpolation (`pchip`), and show looming markers at 2, 2.8, 3.6, 4.4 and 5.2 seconds. Default displayed times are 1–8 seconds. The color scale saturates at 90 and uses a background shift of 10. These display operations do not alter the traces used for the numerical analyses.
2. **Reaction time.** Find the first sample with speed strictly greater than 35 within the inclusive 2–10-second interval and subtract the 2-second trigger time. Nonresponding rows receive `NaN`. Plot mean ± SEM across responding mice only.
3. **Maximum speed.** Find each mouse's maximum speed over the inclusive 2–10-second interval and plot mean ± SEM across mice, omitting `NaN` results.
4. **Escape probability.** Count a mouse as escaping if any speed sample is strictly greater than 35 during the same 2–10-second interval; divide by all rows in that genotype/contrast group. Print numerator, denominator and probability to the Command Window, and plot probabilities.

`speedThreshold` and `analysisWindow` are repeated in their relevant sections; keep them consistent when changing the response definition. `fps`, `contrastList` and `genotypeList` are shared settings. A mouse with a missing/NaN-only response has no threshold crossing and remains in the escape-probability denominator. Decide missing-data handling upstream before running the script.

The script leaves arrays such as `KO_RT_all`, `WT_RT_all`, maximum-speed arrays, group means/SEMs and escape counts/probabilities in the workspace. It opens figures and prints escape counts. Export figures or save selected workspace variables explicitly as needed, outside the source checkout.

Escape probabilities in these plots use the speed-threshold definition described above.
