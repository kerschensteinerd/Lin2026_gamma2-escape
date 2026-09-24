#!/usr/bin/env python3
"""Create a small, auditable figure-data adapter from validated v3 simulations."""
from pathlib import Path
import json
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / 'outputs/fig3_v3'
SPATIAL = ROOT / 'outputs/fig3_v3_spatial'


def main():
    meta = json.loads((DATA/'metadata.json').read_text())
    selected = json.loads((SPATIAL/'selected_parameters.json').read_text())
    for field in ('g_pas_S_cm2', 'inhibition_scale', 'glycine_fraction'):
        assert np.isclose(meta['parameters'][field], selected[field]), f'Spatial/spot parameter mismatch: {field}'
    archive = np.load(SPATIAL/'selected_spatial_maps.npz')
    finite = np.load(SPATIAL/'finite_local_maps.npz')
    distance, length = archive['path_distance_from_source_um'], archive['cable_lengths']
    output = {'contrast': {}, 'spatial': {}}
    for key in ('ctrl', 'gaba'):
        rr = sorted([r for r in meta['contrast_rows'] if r['polarity'] == 'OFF' and r['condition'] == key],
                    key=lambda r: r['contrast_pct'])
        output['contrast'][key] = {
            'contrast_pct': [r['contrast_pct'] for r in rr],
            'mean_peak_dv_mV': [r['mean_peak_dv_mV'] for r in rr],
            'mean_peak_ca_au': [r['mean_peak_ca_au_ms'] for r in rr],
            'mean_ca_auc_au_ms2': [r['mean_ca_auc_au_ms2'] for r in rr],
            'auc_shape': [r['model_auc_shape'] for r in rr]}
        x, yca, yv = [0.], [1.], [1.]
        response_ca = finite[f'local_{key}_finite_ca_normalized']
        response_v = finite[f'local_{key}_finite_v_normalized']
        for low in np.arange(0, distance.max()+5, 5):
            use = (distance >= low) & (distance < low+5)
            if use.any():
                x.append(float(np.average(distance[use], weights=length[use])))
                yca.append(float(np.average(response_ca[use], weights=length[use])))
                yv.append(float(np.average(response_v[use], weights=length[use])))
        output['spatial'][key] = {'distance_um': x, 'transfer_ca': yca, 'transfer_v': yv}
    gly = 100*meta['parameters']['glycine_fraction']
    ca = meta['calcium']
    output['footer'] = (f'Illustrative calibration; fixed inhibition, {gly:g}% glycine retained. '
                        'Shared scales in D-H; I endpoints imposed.')
    output['companion_footer'] = (f'Same model as Figure 3 v3; three-second excitation with fixed shunting inhibition.\n'
        f'{gly:g}% glycinergic conductance retained; common activation midpoint {ca["v_half_mV"]:.2f} mV, slope {ca["slope_mV"]:.2f} mV.')
    (DATA/'figure_data.json').write_text(json.dumps(output, indent=2)+'\n')
    print(DATA/'figure_data.json')


if __name__ == '__main__':
    main()
