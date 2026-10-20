# Two-photon imaging analysis

MATLAB scripts for deinterleaving, motion correction, tdTomato-guided selection of active dendritic pixels, stimulus-marker timing, and spot-size response analysis. Supply the imaging recordings and their acquisition settings separately. The included `example_spot_order.mat` contains only the five stimulus sizes `[400; 50; 200; 100; 25]` in presentation order; use it only when that order matches the recording.

## Requirements and configuration

Use MATLAB with the Image Processing Toolbox (`imregcorr`, `imwarp`, `medfilt2`, `imshow`) and Signal Processing Toolbox (`sgolayfilt`). Use MATLAB R2022b or later: the registration script reads the `Translation` property of the transform returned by `imregcorr`. Validation was performed with MATLAB R2026a Update 3; use that release to match the tested registration behavior. The script uses the release default correlation method, which changed in R2024b. See the [MathWorks `imregcorr` documentation](https://www.mathworks.com/help/images/ref/imregcorr.html) for version behavior. A desktop session is needed for folder selection and interactive thresholding. Add this directory to the MATLAB path and run each script from the base workspace, for example:

```matlab
addpath('/path/to/Lin2026_gamma2-escape/two_photon_imaging');
deinterleave_channels
```

Edit each script's settings block before running. Defaults describe one acquisition layout and are not inferred from TIFF metadata. In particular, confirm channel order, image dimensions, registration/crop windows, marker thresholds, frame rate, stimulus order, and repeat timing. Full recordings and intermediate movies are held in memory; channel movies are saved as `-v7.3` MAT files.

## Execution order and file flow

Use a common recording stem (`<base>`) across channels and stimulus metadata. Select the folders listed below when prompted.

| Order | Script | Select as input | Files written |
| --- | --- | --- | --- |
| 1 | `deinterleave_channels.m` | Folder containing page-interleaved `.tif`/`.tiff` recordings | `extract_4channel/<base>_Gcamp.mat`, `_Tdtomato.mat`, `_TransmittedLight.mat`, `_Digit.mat` |
| 2 | `register_channels.m` | The extracted-channel folder | `registered_crop_centerReg/<base>_<channel>_registered.mat` |
| 3 | `prepare_roi_projections.m` | The registered-channel folder | `ROI_detection/<base>_ROIprep.mat` |
| 4 | `segment_tdtomato_dendrites.m` | The `ROI_detection` folder | `Segmentation/<base>_Segmentation.mat` and `Segmentation/Figures/<base>/` mask PNGs |
| 5 | `extract_stimulus_timing.m` | The extracted-channel folder, containing `_Digit.mat` | `StimTiming/<base>_StimTiming.mat` and a timing QC PNG |
| 6 | `analyze_spot_size_responses.m` | Select the `Segmentation` folder, `StimTiming` folder, then a spot-order metadata folder | `Segmentation/SpotAnalysis/<base>_SpotAnalysis.mat` and a mean-response PNG when pixels qualify |

Stimulus timing can run any time after deinterleaving. Registration and image preparation preserve the number and order of frames, so marker indices and extracted fluorescence traces share a time axis. Spatial `[row col]` positions in segmentation outputs refer to the cropped ROI-preparation images, not the raw TIFF coordinates.

### Channel extraction and motion correction

`deinterleave_channels` defaults to four interleaved pages per time point in the order `Gcamp`, `Tdtomato`, `TransmittedLight`, `Digit`. Variables inside extracted MAT files have those same names and size `[height width frames]`. The script preserves source pixel type. Any incomplete trailing channel cycle is dropped; confirm complete, correctly ordered input acquisitions.

`register_channels` keeps raw-image coordinates `8:122` in both dimensions, producing 115 × 115 images for the default geometry. It estimates translation on cropped coordinates `20:105` against the mean of the first 20 transmitted-light frames (or all frames if fewer). It applies the same translation to GCaMP and the optional tdTomato channel, using zero fill outside the image. Saved files include `shifts` (`[dy dx]` per frame) and crop/registration windows. The reference output variable is `registered_Ref`; the other outputs are `registered_Gcamp` and `registered_Tdtomato` by default.

Confirm equal channel dimensions/frame counts and spatial alignment. Inspect shifts and images to ensure the later crop removes zero-filled borders. Reference and second-channel names are configurable, but `Gcamp` is the required activity-channel name. A missing optional second channel produces a warning; subsequent tdTomato segmentation requires it.

### Preparation and manual pixel selection

`prepare_roi_projections` keeps coordinates `20:110` within each registered image (91 × 91 pixels by default). With `APPLY_HIGHPASS_FILTER = true`, it subtracts a 50-frame moving mean from each GCaMP pixel, then applies a 3 × 3 spatial median filter to each frame with symmetric padding. Setting that flag to `false` omits the moving-mean subtraction. The second channel receives only the spatial median filter. This flag affects both the activity-selection image and exported traces and must match the intended analysis.

Outputs include `Fcrop` (cropped, unfiltered GCaMP), `Ffilt` (processed GCaMP), `stdImg = std(Ffilt,0,3)`, `meanImg = mean(Fcrop,3)`, and the second channel's filtered mean image `meanImg_Red`. Crop windows, moving-mean window, median-filter size, high-pass flag and second-channel name are saved. QC projections are displayed; their limits are guarded for flat/sparse images without changing saved data.

`segment_tdtomato_dendrites` calls the included `thresholdGUI.m` twice per recording: first on the second-channel mean image, then on the GCaMP standard-deviation image. Click **Save Threshold** in each window. The accepted thresholds and separate/intersection masks are saved. Closing a threshold window without saving skips the recording. The helper writes `threshold` into the base workspace, so run the script from that workspace.

Use finite images and inspect the final saved masks. The supplied segmentation rule sets values below threshold to zero and then all remaining nonzero values to one; at a zero threshold it excludes zero-valued pixels, whereas the GUI preview uses `image >= threshold`. Positive thresholds on finite images give matching preview and saved masks.

`RedGreen_pixel_resp` contains the selected red-and-green pixel traces and `Green_pixel_resp` all green-mask traces, each `[frames pixels]`. Their columns match the rows of `zip_RedGreen` and `zip_Green`. Both contain the selected pixels' `Ffilt` values. The tdTomato-guided segmentation script uses `meanImg_Red` and skips recordings without it.

### Marker timing and spot responses

`extract_stimulus_timing` computes the spatial standard deviation of each Digit-channel frame. Defaults use frame 200 as the baseline, replace marker spikes greater than baseline + 60 with the following value, and threshold at baseline + 5. It expects five complete marker square waves and warns on other counts. `Stim_ON` is the first ON frame and `Stim_OFF` the last ON frame, both 1-based. Incomplete events at recording boundaries are discarded. Saved `stim_durations` are `Stim_OFF - Stim_ON`; these are not inclusive frame counts. Inspect each timing QC PNG and change thresholds/baseline frame for the actual marker signal.

Each spot-order MAT file must contain a vector `spot_order`, with one unique size per detected marker square wave. The analysis first looks for `<base>_spot_order.mat`; otherwise, a single `*spot_order*.mat` in the selected folder is used for all recordings. Multiple possible fallback files cause the recording to be skipped. Spot-size units must match the experiment's metadata; the script does not perform a spatial calibration.

`analyze_spot_size_responses` uses these defaults:

| Setting | Value and meaning |
| --- | --- |
| `baselineWindow` | `[-34 0]`, inclusive offsets from each `Stim_ON`; median of 35 samples, including the ON frame |
| `nRepeats` | 3 repeat windows per marker square wave |
| `repeatOffset` | First repeat starts 33 frames after `Stim_ON` |
| `repeatPeriod`, `repeatLength` | Repeat starts are 67 frames apart; each window is 67 frames long |
| `sgolayOrder`, `sgolayFrameLen` | Fifth-order Savitzky–Golay smoothing over 13 frames |
| `rrThreshold` | Strictly `RR > 0.4` for at least one spot size |
| `fps`, `offWindow` | 16.67 Hz for plotting; shade seconds `[2 4]` after the repeat-window start as stimulus OFF |

The baseline window must fit within the recording and start after the preceding square wave ends. The default includes the ON frame: confirm that this reflects the marker timing before using it as a prestimulus estimate. The three repeat windows are defined by fixed offsets/lengths, not by detecting individual repeats; verify that they fit each marker event as well as the recording.

For each pixel and size, repeat reliability is the mean of pairwise Pearson correlations among the three smoothed repeat traces. A pixel qualifies when any size exceeds the threshold; constant repeats can yield undefined correlations and do not qualify on those values. Responses for each size are averaged across qualified pixels and repeats. The output includes `repeatResp` (`[time repeats sizes pixels]`), `RR` (`[pixels sizes]`), `isQualified`, pixel positions, per-wave baseline levels, `meanResp` (`[time sizes]`), sorted sizes, presentation order, marker timing, and analysis/plot settings. These responses are baseline-subtracted fluorescence in arbitrary units.

## Reprocessing and scope

Extraction, registration and preparation overwrite same-named outputs. Segmentation and marker timing skip existing result MAT files; remove stale outputs before deliberately rerunning those stages after changing inputs or settings. Spot analysis overwrites its result files. Keep all raw and generated data outside the source checkout.

These scripts provide the preprocessing and spot-size analysis pipeline described above. Synthetic execution checks establish software data flow. Apply the pipeline using the recordings, acquisition calibration, actual settings and manually chosen masks for the experiment.
