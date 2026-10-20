# VGluT3 amacrine cell compartmental model

This directory reproduces the model panels for Figure 4 and the calcium-readout
analyses for Figures S6/S7 in the Lin 2026 gamma2 escape study. It contains a
multicompartment NEURON model of the reconstructed VGluT3 amacrine cell cid3,
its morphology and synapse map, fitted size-tuning parameters, simulation
scripts, and plotting scripts. The modeling code is by Pratyush Ramakrishna
and retains the [MIT license](LICENSE).

## Installation

Use Python 3.11 with NEURON, NumPy, SciPy, and Matplotlib. From the repository
root, create an environment and compile the potassium-channel mechanism for
your operating system:

```bash
cd modeling
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-tested.txt
nrnivmodl mod
```

On Windows, activate the virtual environment with `.venv\Scripts\activate`
and compile `mod/khh_vg3.mod` using NEURON's `mknrndll` utility from the
`modeling` directory. A C/C++ compiler compatible with the installed NEURON
version is required. Compiled mechanism libraries are platform specific and
are generated locally. The scripts load the compiled mechanism from this
directory, including when invoked from another working directory.

The recorded validation environment uses Python 3.11.15, NEURON 9.0.2,
NumPy 2.4.6, SciPy 1.17.1, and Matplotlib 3.11.1 on macOS arm64.
`requirements.txt` provides the unpinned dependency list. Six worker processes
are used for the larger simulations; adjust `N_WORKERS` in the simulation
scripts for machines with fewer available cores or limited memory.

## Required inputs and configuration

All inputs needed for the simulation are included; raw imaging or behavioral
recordings are not required.

| Input | Purpose |
| --- | --- |
| `data/cid3.swc` | Reconstructed cid3 morphology |
| `data/sm2nrn_cid3.mat` | EM synapse positions, types, and presynaptic identities |
| `data/vg3_onbip_list.txt`, `data/vg3_offbip_list.txt` | ON/OFF bipolar classification lists |
| `data/fitted/fitted_off_epsc_params.json`, `fitted_on_epsc_params.json`, `fitted_inh_params.json` | Fitted excitatory and inhibitory size-tuning gains |
| `data/locked_params.json` | Passive, synaptic, stimulus, and contrast-gain parameters |
| `data/published/focal_source.json` | Focal source coordinates, initial conductance estimate, target source depolarization, and OFF imaging targets |
| `mod/khh_vg3.mod` | Sustained potassium-current mechanism source |

`src/vg3_base_functions.py` supplies defaults, morphology construction, and
synapse mapping. `src/run_functions.py` combines the locked parameters with
the fitted size-tuning gains. `base_functions.py` defines the excitation and
inhibition scale factors, the Figure 4 stimulus, and the focal-input protocol.
`calcium_functions.py` defines the calcium proxy, normalization, and readout
fit. Simulation grids and worker counts are near the top of each entry script.

## Execution order

Run the following commands from `modeling/` with the environment activated:

```bash
python 00_plot_contrast_gains.py
python 01_run_fig4_simulations.py
python 02_plot_fig4_panels.py
python 03_run_calcium_sweep.py
python 04_plot_calcium_sweep.py
python 05_plot_example_voltage_traces.py
python 06_run_small_spot_spread.py 50
python 07_plot_small_spot_spread.py 50
python calibrate.py
```

| Script | Prerequisite | Main outputs |
| --- | --- | --- |
| `00_plot_contrast_gains.py` | Compiled model and input data | Conductance-versus-contrast panel in `figures/fig4_panels/` |
| `01_run_fig4_simulations.py` | Compiled model and input data | `output/fig4/`: per-contact voltage peaks, calibration, contrast tables, and focal spatial profiles; maps and voltage traces in `nominal25/` |
| `02_plot_fig4_panels.py` | Script 01 | Individual Figure 4 panels as SVG, PDF, and PNG in `figures/fig4_panels/` |
| `03_run_calcium_sweep.py` | Script 01 | `output/calcium_sweep/`: parameter-grid arrays, TSV metrics, and JSON summary |
| `04_plot_calcium_sweep.py` | Script 03 | Sweep panels in `figures/calcium_sweep_panels/`, with plotted arrays and TSV source data in `source_data/` |
| `05_plot_example_voltage_traces.py` | Script 01 | Focal voltage traces at five path distances in `figures/example_traces/` |
| `06_run_small_spot_spread.py [diameter_um]` | Script 01 | Voltage/calcium maps, traces, and summary in `output/small_spot_<diameter>um/` |
| `07_plot_small_spot_spread.py [diameter_um]` | Scripts 01 and 06 with matching diameter | Small-spot maps, attenuation profiles, and traces in `figures/small_spot_<diameter>um/` |
| `calibrate.py [--solve]` | Compiled model and input data | Printed somatic size tuning, input resistance, clamp conductances, and potassium-current checks; `--solve` also estimates potassium density without changing input files |

Scripts 01, 03, and 06 skip execution when their completion output already
exists. To recompute after changing inputs or parameters, move or remove the
corresponding generated output directory and then rerun dependent steps.
Plotting scripts overwrite their generated panels. Generated arrays, figures,
logs, caches, and compiled libraries are excluded from the release.

## Model and analysis definitions

Excitatory and inhibitory inputs are single-event `Exp2Syn` conductances at
mapped EM contacts. The Figure 4 protocol uses a 100 µm OFF spot with 25 µm
Gaussian blur, 100 ms onset, a 0.1 ms fixed integration step, and 200 ms total
duration. The cell has a passive leak and a sustained Kv current. Control
retains all inhibitory weights; the nominal reduced condition retains 25%
of every inhibitory contact's weight, with 20% and 30% sensitivity cases.
These fractions describe effective inhibitory conductance, without assigning
transmitter identities. The code does not remove synapses in these runs.

Excitation and inhibition use separate saturating contrast gains normalized
to one at 100% contrast. The excitation half-maximal contrast and exponent
are 35% and 2.0; the inhibitory values are 15% and 1.4. These contrast gains
and the Kv activation kinetics are modeling assumptions. Excitatory rise/decay
time constants are 10/300 ms; inhibitory values are 20/300 ms. The passive
leak is 1 mS/cm², axial resistivity 100 Ω·cm, capacitance 1 µF/cm², and cable
diameter 0.5 µm. Rest, excitatory reversal, and inhibitory reversal are −45,
0, and −60 mV. The soma is identified from native morphology before imposing
the uniform cable diameter.

The calcium proxy is baseline-subtracted Boltzmann activation of each site's
peak absolute voltage. It has no calcium-clearance or indicator kinetics and
is expressed as activation fraction, rather than calibrated concentration or
ΔF/F. Contrast responses apply the readout separately at the 75 OFF contacts
and then average. Each condition's normalized curve sets its 15% response to
zero and its 70% response to one. The readout fit searches midpoint and slope,
with slope at least 2.5 mV, using the supplied OFF imaging targets plus weak
constraints on low-contrast calcium and an assumed control/reduced ratio of
0.6 at 70%. This ratio is a prior, since the imaging targets are normalized
separately by condition.

The focal protocol applies one excitatory synapse at the segment nearest the
specified coordinates, solves its conductance for 15 mV source depolarization
in control, and combines it with phasic inhibition from the 100 µm spot at
30% contrast. Source-normalized voltage and calcium profiles use cable-length
weighted means in 5 µm path-distance bins. The first interpolated 1/e crossing
quantifies spread. Script 05 extends the focal traces to 600 ms; the Figure 4
peak maps use the 100–200 ms interval. Script 06 instead stimulates both local
excitation and local inhibition with a small spot and records through 600 ms.

The sweep evaluates 2,601 midpoint/slope pairs: −45 to −20 mV in 0.5 mV
steps and 1 to 6 mV in 0.1 mV steps, plus the fitted pair. It reuses the
simulated voltages, calculates low/high-contrast gain preferences and
normalized response shape, finds equivalent control contrasts, and measures
focal calcium spread. Hatched map cells mark ratios with near-zero control
responses; undefined or censored values are represented separately in the
saved arrays and plot legends.

## Calibration and scope

The fitted size-tuning inputs and `calibrate.py` refer to Kim et al. (2015),
Figure 1C/F. The sustained Kv current magnitude is checked against Grimes,
Seal, Oesch, Edwards, and Diamond (2011), *Visual Neuroscience* 28:381–392,
Figure 3E. Passive and cable parameters follow the Friedrichsen (2024) VG3
model; the excitatory contrast half-maximal value is motivated by
Poleg-Polsky and Diamond (2016). Kv activation midpoint, slope, and time
constant are explicit model choices rather than measured VG3 values.

The model is phenomenological. Its somatic input resistance and clamp
conductances differ from the cited measurements, and its calcium fit is
approximate. The scripts expose these quantities through calibration outputs;
they do not assert that every measured property is reproduced.

## References

- Kim T, Soto F, Kerschensteiner D (2015). An excitatory amacrine cell detects
  object motion and provides feature-selective input to ganglion cells in the
  mouse retina. *eLife* 4:e08025. [doi:10.7554/eLife.08025](https://doi.org/10.7554/eLife.08025).
- Grimes WN, Seal RP, Oesch N, Edwards RH, Diamond JS (2011). Genetic targeting
  and physiological features of VGLUT3+ amacrine cells. *Visual Neuroscience*
  28:381–392. [doi:10.1017/S0952523811000290](https://doi.org/10.1017/S0952523811000290).
- Friedrichsen et al. (2024). Subcellular pathways through VGluT3-expressing
  mouse amacrine cells provide locally tuned object-motion-selective signals
  in the retina. *Nature Communications*. [doi:10.1038/s41467-024-46996-0](https://doi.org/10.1038/s41467-024-46996-0).
- Poleg-Polsky A, Diamond JS (2016). Retinal Circuitry Balances Contrast Tuning
  of Excitation and Inhibition to Enable Reliable Computation of Direction
  Selectivity. *Journal of Neuroscience* 36:5861–5876.
  [doi:10.1523/JNEUROSCI.4013-15.2016](https://doi.org/10.1523/JNEUROSCI.4013-15.2016).
