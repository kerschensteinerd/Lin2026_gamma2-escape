# Electrical and spatial sweep execution

The user requested execution of protocol sections B and C and preparation of the calcium supplement on 21 September 2026. The original protocol and calcium results remain immutable. New files reside in `electrical/`, `spatial/`, and `calcium_supplement/`.

## Authorized matrices

- Electrical: 81 inhibition/leak pairs, three excitation multipliers, four inhibitory conditions, and 11 contrasts = 10,692 design rows. The primary nominal-excitation, 25%-retained map is a subset. Focal spread uses the same 81 pairs, three focal strengths, and four inhibitory conditions = 972 cases.
- Spatial: 75 OFF sources, three focal strengths, and four inhibitory conditions = 900 cases. The 69 reduced-inhibition layouts and one control yield 770 spot design rows and 210 focal layout/strength cases.
- Calcium supplement: figure and companion sensitivity plots from the validated existing calcium data, with legend and proposed manuscript text. No current manuscript or existing figure is replaced.

## Numerical execution clarification

Full-length reference NEURON runs retain the historical fixed 0.1-ms integration and 2-ms recording, sustained step, local calcium transformation, and 50–3050-ms calcium-AUC definition. Production may shorten only the electrically settled tail after validating a quantitative error bound against full-length NEURON. The initial transient is retained, and the calcium response continues through the full original response window. This is an explicit numerical acceleration of the specified protocol, not a steady-state replacement for the whole AUC.

The passive membrane permits a bound on unobserved convergence. With uniform capacitance and positive leak/shunts, the backward-Euler voltage propagation contracts in the maximum norm at least as fast as `q=(1+dt/tau_max)^(-Delta/dt)`, where `tau_max=cm/(1000*g_leak)` ms. The maximum voltage change over a checked block across all compartments bounds the remaining voltage error at the current integration time by `q/(1-q)` times that change. Production conservatively uses `1/(1-q)` instead, covering the lag between the last 2-ms recorded sample and the current integration time. The shared sigmoid is Lipschitz with constant `1/(4*k)`; multiplying this voltage bound by `tau * remaining_duration/(4*k)` bounds calcium-AUC error. The 10-ms block check stops only once the bound is at most `1e-7 * max(1, accumulated_prefix_mean_AUC)`. The accumulated prefix is no larger than the final nonnegative AUC, making the criterion conservative. Actual bounds and full-length validation results are recorded for each run. No biological parameter is retuned.

Repeated excitation-multiplier × contrast products can be cached because they prescribe identical point conductances. Every requested design row must still be emitted, with a documented link to its unique simulation. Uniform zero-contrast, pure-shunt baselines are exactly −45 mV and yield zero calcium response. Skipping an exactly quiescent prestimulus period is permissible only after demonstrating unchanged onset sampling against the original implementation.

Original source and manuscript preservation is verified against `input_hashes.json`. Root review independently checked output completeness, paired conditions, nominal reproduction, spatial-layout conductance budgets, numerical validations, and final figures. Its executable audit and results are saved alongside this note. Grid points, compartments, sites, and random layouts remain computational sensitivity samples rather than biological replicates.
