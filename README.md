# Lin2026 gamma2 escape study code

Code for compartmental modeling, two-photon image processing, and looming
behavior analysis accompanying the Lin 2026 gamma2 escape study.

| Directory | Contents | Start here |
| --- | --- | --- |
| `modeling/` | Reconstructed VGluT3/cid3 NEURON model, morphology and synapse map, fitted parameters, Figure 4 simulations, calcium-readout sweeps, and small-spot simulations | [Model installation and workflow](modeling/README.md) |
| `two_photon_imaging/` | MATLAB channel extraction, translation registration, ROI preparation, interactive tdTomato/activity masking, stimulus timing, and spot-size response analysis | [Imaging inputs and execution order](two_photon_imaging/README.md) |
| `behavior/` | MATLAB speed heatmaps, reaction times, maximum speeds, and speed-threshold escape probabilities from precomputed mouse speed traces | [Behavior input schema and analysis](behavior/README.md) |

## Running the code

The model is self-contained: its required reconstruction, synapse map, and
parameter files are included. Use Python 3.11 and the dependencies in
`modeling/requirements-tested.txt`, compile the NEURON mechanism with
`nrnivmodl mod` from `modeling/`, then run the numbered scripts in the order
shown in the model README. Arrays and figures are generated locally.

Imaging requires MATLAB R2022b or later; behavior requires R2018b or later.
Imaging additionally requires the Image Processing and Signal Processing Toolboxes, a desktop
session for interactive threshold selection, the original interleaved TIFF
recordings, and matching stimulus-order metadata. Run `deinterleave_channels`,
`register_channels`, `prepare_roi_projections`, and
`segment_tdtomato_dendrites` in sequence; run `extract_stimulus_timing` after
channel extraction, then `analyze_spot_size_responses`. Edit the acquisition,
crop, timing, and filter settings to match the recordings before execution.
The included spot-order MAT file contains stimulus sizes only.

For behavior, run `looming_behavior_analysis` and select a MAT file containing
`HeatmapData`, with one struct entry for each WT/KO and 22/64/91 contrast
combination. Each entry contains `Genotype`, `Contrast`, and a `[mice frames]`
speed matrix. Defaults use 30 Hz, a 2-second stimulus trigger, and an inclusive
2–10-second analysis window with speed threshold 35 in the input's calibrated
units. The script displays results and leaves numerical summaries in the
MATLAB workspace.

## Reproduction scope

Experimental recordings, tracking outputs, and merged behavior speed data
must be supplied separately. The imaging workflow includes preprocessing and a tdTomato-guided
spot-size example. Behavior analysis starts from calibrated speed traces.
Each folder README gives exact inputs, outputs, settings, and reprocessing
behavior.

The model workflow was executed using NEURON 9.0.2 on macOS arm64. MATLAB
scripts were checked and exercised with synthetic inputs in MATLAB R2026a.
These checks validate software execution and file flow; experimental figure
reproduction also requires the original recordings and analysis settings.

## License and attribution

The modeling code retains its original [MIT license](modeling/LICENSE) and
scientific attribution in the model README. That license applies to the
modeling directory; no additional license is asserted for the MATLAB analysis
scripts. Preserve the supplied notices when redistributing the model.
