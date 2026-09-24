# VG3 (cid3) Retinal Amacrine Cell — Receptive Field Simulation

A multi-compartmental NEURON model of a VGluT3 (cid3) retinal amacrine cell,
fit to Kim et al. 2015 electrophysiology. The model reproduces ON/OFF spatial
integration and shows how inhibition shapes response gain and confines signal
spread across the dendritic arbor.

* **Morphology:** cid3 EM reconstruction (SWC), synapse positions mapped from
  the same reconstruction.
* **Biophysics:** passive dendrites, fitted OFF/ON bipolar EPSC and GABAergic
  IPSC size-tuning.
* **Synaptic input:** ON- and OFF-bipolar synapses activated by a Gaussian
  spot or a rectangular bar, plus inhibitory synapses that can be partially
  removed to test their role in shaping the response.
* **Output:** overlaid voltage traces, 3D peak-depolarization spread along the
  skeleton, and contrast/gain tuning curves, comparing full inhibition ("Ctrl")
  against 75% of inhibitory synapses removed ("−75% inh.").

<p align="center">
  <img src="figures/skeleton_3d_example.png" width="600" alt="Peak-voltage spread example">
  <br>
  <em>Figure: example output from run_example.py (generate by running it once).</em>
</p>

## How to run

```
./setup.sh
```
Installs NEURON and the other dependencies. See Requirements below if it can't
install NEURON automatically on your platform.

### 1. Example
* **Script:** `scripts/example/run_example.py`
* **Function:** runs one stimulus (`--stim spot` or `--stim bar`) at full vs
  75%-removed inhibition; produces the overlaid trace, 3D spread, and contrast
  tuning curve.
* **Output:** `outputs/example/` — `trace_overlay.png`, `skeleton_3d.png`,
  `contrast_tuning.png`, `contrast_tuning_data.txt`.

### 2. Gain sweep
* **Script:** `scripts/sweep/run_gain_sweep.py`
* **Function:** sweeps a single gain multiplier over the fitted excitatory and
  inhibitory weights, at full vs 75%-removed inhibition.
* **Output:** `outputs/sweep/` — `gain_sweep_tuning.png`, `gain_sweep_data.txt`.

## Model parameters

| Parameter | Value | Description |
| :--- | :--- | :--- |
| `Ra` | 100 Ω·cm | Axial resistivity |
| `Cm` | 1 µF/cm² | Membrane capacitance |
| `g_pas` | 0.001 S/cm² | Passive leak conductance |
| `E_leak` (rest) | −45 mV | Resting / leak reversal |
| `E_exc` | 0 mV | Excitatory (AMPA-like) reversal |
| `E_inh` | −60 mV | Inhibitory (GABA/Cl⁻) reversal |
| `exc_weight` (OFF) | 8.59e-6 µS | OFF-bipolar EPSC weight, fit to Kim 2015 |
| `d50_exc` / `n_exc` (OFF) | 137.7 µm / 1.86 | OFF center-excitation Hill |
| `d50_supp` / `n_supp` (OFF) | 157.8 µm / 2.32 | OFF surround-suppression Hill |
| `exc_weight` (ON) | 3.15e-6 µS | ON-bipolar EPSC weight, fit to Kim 2015 |
| `on_scale` | 2.0 | ON count-asymmetry compensation |
| `d50_exc` / `n_exc` (ON) | 400 µm / 0.84 | ON center-excitation Hill |
| `d50_supp` / `n_supp` (ON) | 263.9 µm / 3.69 | ON surround-suppression Hill |
| `inh_weight` | 5.61e-6 µS | Inhibitory synapse weight, fit to Kim 2015 |
| `A_off` / `A_on` | 1.0 / 0.57 | Per-polarity inhibition amplitude |
| `G_exc` / `G_inh` | 400 / 600 | Voltage-scale gains, fit to Kim 2015 Fig 1C |

## Requirements

* Python 3.8+
* NEURON with Python support — `pip install neuron` on macOS/Linux/WSL; native
  Windows needs the official binary installer
  (https://github.com/neuronsimulator/nrn/releases). `setup.sh` detects your
  platform and prints the right instructions.
* NumPy, SciPy, Matplotlib

## License

MIT — see [LICENSE](LICENSE).
