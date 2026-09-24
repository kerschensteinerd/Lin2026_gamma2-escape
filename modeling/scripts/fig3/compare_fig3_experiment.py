#!/usr/bin/env python3
"""Unfitted, condition-normalized shape comparison with the Fig. 3 workbook.

Run using the repository .venv. Twelve independent single-event simulations
retain the inherited voltage model. The experimental 3 s ON / 3 s OFF protocol
is represented only by separate 3 s response integration windows; the inherited
10/300 ms synaptic kernel is not an exact stimulus or calcium-trace match.
"""

from pathlib import Path
import argparse
import csv
import hashlib
import json
import time

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

from run_fig3_pilot import (
    ROOT, COLORS, MODES, build, load_fits, load_named_config,
    exc_eff_weight, inh_eff_weight, gaussian_spot_weights, simulate, write_rows,
)
from calcium_readout import calcium_trace


CONTRASTS = np.array([.15, .30, .70])
READOUT = dict(v_half=-40., slope=4., tau_ms=50., eca_mV=120.,
               cai_mM=.0001, cao_mM=2., temp_C=33.)


SPATIAL_SAMPLING = dict(
    model='ON/OFF bipolar input-contact samples; imaging sublaminae are not registered',
    experiment='Imaging sublaminae at approximately 65% ON and 35% OFF IPL depth',
    model_ON_contact_depth_percent=dict(median=51.0, minimum=37.9, maximum=61.4,
                                       below_50_percent=12, within_65_plusminus_5_percent=1,
                                       n=30),
    model_OFF_contact_depth_percent=dict(median=40.3, all_below_50_percent=True, n=75),
    depth_mapping='Inherited 40 um IPL mapping; audit summary, not image-to-morphology registration',
    inference='Shape differences apply to a provisional contact-based readout; they do not isolate a discrepancy in the measured sublaminae.',
)
READOUT_DESCRIPTION = ('Mean across fixed ON/OFF bipolar input-contact samples of local '
                       'baseline-subtracted calcium AUC; signed values retained; '
                       'imaging sublaminae are not registered')


def plot_comparison(rows, out):
    plt.rcParams.update({'font.size': 9, 'axes.spines.top': False,
                         'axes.spines.right': False, 'svg.fonttype': 'none'})
    fig, axes = plt.subplots(2, 3, figsize=(11.6, 7.4))
    titles = {'fixed': 'Activation only', 'ohmic': 'Activation × ohmic driving term',
              'ghk': 'Activation × GHK permeation'}
    for row_idx, polarity in enumerate(['ON', 'OFF']):
        for col_idx, mode in enumerate(MODES):
            ax = axes[row_idx, col_idx]
            for condition, color in COLORS.items():
                subset = [r for r in rows if r['polarity'] == polarity and
                          r['condition'] == condition and r['mode'] == mode]
                ax.errorbar(CONTRASTS * 100, [r['experiment_shape'] for r in subset],
                            yerr=[r['experiment_shape_sem'] for r in subset],
                            color=color, marker='o', ms=4, lw=1.7, capsize=2, zorder=3)
                ax.plot(CONTRASTS * 100, [r['model_shape'] for r in subset], '--s',
                        color=color, ms=4, mfc='white', lw=1.6, zorder=4)
            ax.set(xlim=(10, 75), ylim=(-.12, 1.26), xticks=[15, 30, 70],
                   xlabel='Contrast (%)', ylabel=f'{polarity}: provisional response shape')
            if row_idx == 0:
                ax.set_title(titles[mode], fontsize=10)
            # Endpoints are imposed separately per condition, so only 30% is a
            # free shape comparison; those endpoints do not validate the model.
            ax.axvline(30, color='#dddddd', lw=.6, zorder=0)
            ax.grid(axis='y', color='#eeeeee', lw=.6)
    handles = [
        Line2D([], [], color=COLORS['control'], lw=2, label='Control'),
        Line2D([], [], color=COLORS['gaba_removed'], lw=2, label='cKO / modeled GABA removal'),
        Line2D([], [], color='black', marker='o', lw=1.5, label='Experiment (mean ± workbook SEM)'),
        Line2D([], [], color='black', marker='s', mfc='white', ls='--', lw=1.5, label='Unfitted model'),
    ]
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(.5, .91),
               ncol=2, frameon=False, fontsize=9)
    fig.suptitle('Fig. 3 shape comparison • 100 µm spot • 3 s response window', y=.985, fontsize=13)
    fig.text(.5, .94, 'Model: ON/OFF input-contact samples; imaging sublaminae are not registered',
             ha='center', fontsize=9)
    fig.text(.5, .035,
             'Each condition scaled as [AUC(C) − AUC(15%)] / [AUC(70%) − AUC(15%)]; '
             'endpoints 0 and 1 are imposed.\n'
             'Model: single event at 50 ms, inherited 10/300 ms kernel; '
             '70% GABA / 30% glycine, with glycine retained.\n'
             'Illustrative voltage-to-calcium readouts; no fit, indicator model, '
             'or exact experimental trace match. Spatial sampling differs; discrepancies are provisional.',
             ha='center', va='bottom', fontsize=8.5, linespacing=1.5)
    fig.subplots_adjust(left=.075, right=.985, top=.80, bottom=.18, hspace=.37, wspace=.29)
    for extension in ['png', 'svg']:
        fig.savefig(out / f'experimental_shape_comparison.{extension}', dpi=180, facecolor='white')
    plt.close(fig)



def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--experiment', type=Path, default=
                        ROOT.parent / 'analysis/fig3_model/experimental_fig3.json')
    parser.add_argument('--outdir', type=Path, default=ROOT / 'outputs/fig3_pilot')
    parser.add_argument('--plot-only', action='store_true',
                        help='Replot cached comparison TSV and update spatial-readout labels; no NEURON runs.')
    parser.add_argument('--dt', type=float, default=.05)
    parser.add_argument('--sample-dt', type=float, default=1.)
    args = parser.parse_args()
    if args.plot_only:
        out = args.outdir
        with (out / 'experiment_comparison.tsv').open() as stream:
            rows = list(csv.DictReader(stream, delimiter='\t'))
        text_fields = {'polarity', 'condition', 'mode', 'experiment_source_range'}
        rows = [{key: value if key in text_fields else float(value)
                 for key, value in row.items()} for row in rows]
        plot_comparison(rows, out)
        metadata_path = out / 'experiment_comparison_metadata.json'
        metadata = json.loads(metadata_path.read_text())
        metadata['spatial_sampling'] = SPATIAL_SAMPLING
        metadata['readout'] = READOUT_DESCRIPTION
        metadata['interpretation'] = ('Only 30 percent compares normalized shape; endpoints are imposed. '
                                      'Discrepancies are provisional because model input-contact samples '
                                      'and experimental imaging sublaminae are not registered.')
        metadata['cached_replot_only'] = True
        metadata['hashes'][str(Path(__file__).resolve())] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        metadata_path.write_text(json.dumps(metadata, indent=2) + '\n')
        print(f'Replotted cached comparison with explicit spatial-sampling limitations: {out}', flush=True)
        return
    if args.dt <= 0 or args.sample_dt <= 0:
        parser.error('Time steps must be positive.')
    experiment = json.loads(args.experiment.read_text())
    if not np.allclose(experiment['contrasts'], CONTRASTS):
        raise ValueError('Experimental contrasts do not match the planned comparison.')
    out = args.outdir
    out.mkdir(parents=True, exist_ok=True)
    started = time.time()
    cfg, _ = load_named_config('sep_g22')
    model = build(cfg['g_pas'])
    fits = load_fits()
    rest = model['rest_v']
    coords = model['node_coords']
    inh = model['inh_nodes']
    rows = []
    archive = {}
    recorded_nodes = {}

    for polarity in ['ON', 'OFF']:
        nodes = model[f'{polarity.lower()}_nodes']
        recorded_nodes[polarity] = nodes.tolist()
        eweights = (exc_eff_weight(fits, polarity, 100.) * cfg['G_exc'] *
                    gaussian_spot_weights(coords[nodes, :2], model['soma_xy'], 100., 25.))
        iweights = (inh_eff_weight(fits, polarity, 100.) * cfg['G_inh'] *
                    gaussian_spot_weights(coords[inh, :2], model['soma_xy'], 100., 25.))
        for condition, remaining, experimental_label in [
                ('control', 1., 'Ctrl'), ('gaba_removed', 0., 'cKO')]:
            group = experiment['groups'][f'{experimental_label}_{polarity}']
            for contrast in CONTRASTS:
                t, voltage = simulate(
                    model, nodes, eweights * contrast, iweights * contrast, nodes,
                    gaba_fraction=.7, gaba_remaining=remaining,
                    e_inh=cfg['inh_rev'], dt=args.dt, sample_dt=args.sample_dt,
                    tstop=3051.,
                )
                archive[f'{polarity}_{condition}_{contrast:.2f}'] = voltage.astype(np.float32)
                window = (t >= 50. - 1e-7) & (t <= 3050. + 1e-7)
                if not np.isclose(t[window][0], 50.) or not np.isclose(t[window][-1], 3050.):
                    raise ValueError('Recording grid must contain both integration endpoints.')
                recorded = group['contrasts'][str(float(contrast))]
                peak_dv = float((voltage - rest).max(axis=1).mean())
                for mode in MODES:
                    calcium = calcium_trace(voltage, args.sample_dt, mode=mode,
                                            rest_v=rest, **READOUT)
                    # The readout is already baseline-subtracted. Negative values
                    # remain signed; transform each site before spatial averaging.
                    local_auc = np.trapezoid(calcium[:, window], x=t[window], axis=1)
                    rows.append(dict(
                        polarity=polarity, condition=condition, mode=mode,
                        contrast=float(contrast), n_model_sites=len(nodes),
                        model_mean_local_peak_dV_mV=peak_dv,
                        model_mean_local_calcium_auc_au_ms2=float(local_auc.mean()),
                        model_min_local_calcium_auc_au_ms2=float(local_auc.min()),
                        model_max_local_calcium_auc_au_ms2=float(local_auc.max()),
                        model_shape=float('nan'),
                        experiment_mean_auc=float(recorded['mean']),
                        experiment_sem_auc=float(recorded['sem']),
                        experiment_n=int(recorded['n']),
                        experiment_shape=float(recorded['normalized_mean']),
                        experiment_shape_sem=float(recorded['normalized_sem']),
                        model_minus_experiment_shape=float('nan'),
                        experiment_source_range=recorded['range'],
                    ))
                print(f'{polarity} {condition} contrast={contrast:.2f}: '
                      f'mean local peak dV={peak_dv:.4f} mV', flush=True)

    for polarity in ['ON', 'OFF']:
        for condition in COLORS:
            for mode in MODES:
                subset = [r for r in rows if r['polarity'] == polarity and
                          r['condition'] == condition and r['mode'] == mode]
                low = subset[0]['model_mean_local_calcium_auc_au_ms2']
                span = subset[-1]['model_mean_local_calcium_auc_au_ms2'] - low
                if not np.isfinite(span) or abs(span) < 1e-10:
                    raise ValueError('Cannot normalize a zero/nonfinite model response span.')
                for row in subset:
                    row['model_shape'] = (row['model_mean_local_calcium_auc_au_ms2'] - low) / span
                    row['model_minus_experiment_shape'] = row['model_shape'] - row['experiment_shape']
    write_rows(out / 'experiment_comparison.tsv', rows)
    np.savez_compressed(out / 'experiment_comparison_voltage_traces.npz', t_ms=t,
                        ON_record_nodes=model['on_nodes'], OFF_record_nodes=model['off_nodes'],
                        **archive)

    plot_comparison(rows, out)

    mismatches = [dict(polarity=r['polarity'], condition=r['condition'], mode=r['mode'],
                       model_at_30_percent=r['model_shape'],
                       experiment_at_30_percent=r['experiment_shape'],
                       signed_difference=r['model_minus_experiment_shape'])
                  for r in rows if np.isclose(r['contrast'], .3)]
    metadata = dict(
        status='unfitted_exploratory_shape_comparison', experimental_source=str(args.experiment.resolve()),
        experimental_protocol='100 um spot; 3 s ON then 3 s OFF; contrasts 15, 30, 70 percent',
        model_protocol='Separate polarity simulations; one event at 50 ms; inherited 10/300 ms synaptic kernels',
        n_simulations=12, config=cfg, dt_ms=args.dt, sample_dt_ms=args.sample_dt,
        tstop_ms=3051., integration_window_ms=[50., 3050.], spot_diameter_um=100., spot_blur_sigma_um=25.,
        readout=READOUT_DESCRIPTION, spatial_sampling=SPATIAL_SAMPLING,
        recorded_segment_indices=recorded_nodes,
        calcium=dict(**READOUT, rest_v=rest, activation='instantaneous Boltzmann',
                     feedback_to_voltage=False, indicator_model=False, release_model=False),
        inhibition=dict(gaba_fraction=.7, glycine_fraction=.3, gaba_remaining_control=1.,
                        gaba_remaining_removal=0., glycine_remaining_both=1.,
                        transmitter_assignment='effective co-located split; identities unmeasured'),
        shape_normalization='(AUC(C)-AUC(.15))/(AUC(.70)-AUC(.15)), independently for each condition and readout',
        interpretation='Only 30 percent compares normalized shape; imposed endpoints do not show equal absolute maxima or fit quality. Discrepancies are provisional because model input-contact samples and experimental imaging sublaminae are not registered.',
        experiment_error_bars='Workbook normalized SEM: raw SEM divided by condition span; normalization uncertainty not propagated.',
        mismatch_at_30_percent=mismatches, elapsed_seconds=time.time()-started,
        numerical_validation='Parent performs separate integration-step and sampling checks.',
        hashes={str(path.resolve()): hashlib.sha256(path.read_bytes()).hexdigest()
                for path in [args.experiment, Path(__file__),
                             Path(__file__).with_name('run_fig3_pilot.py'), ROOT / 'src/calcium_readout.py']},
    )
    (out / 'experiment_comparison_metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
    print(json.dumps(mismatches, indent=2), flush=True)
    print(f'Wrote comparison to {out}; {time.time()-started:.1f} s', flush=True)


if __name__ == '__main__':
    main()
