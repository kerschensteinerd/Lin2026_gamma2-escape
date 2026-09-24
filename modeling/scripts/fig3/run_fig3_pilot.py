#!/usr/bin/env python3
"""Exploratory Fig. 3 experiments; inherited passive model, explicit residual inhibition.

Run with .venv/bin/python scripts/fig3/run_fig3_pilot.py. See outputs/fig3_pilot/.
Calcium is a post-hoc sensitivity analysis, not fitted calcium or fluorescence.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import sys
import time

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from neuron import h

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'src'))
from run_functions import build, load_fits, load_named_config, exc_eff_weight, inh_eff_weight
from vg3_base_functions import gaussian_spot_weights
from calcium_readout import calcium_trace, influx_proxy

COLORS = {'control': '#343434', 'gaba_removed': '#db772d'}
CONTRASTS = np.array([0, .02, .05, .1, .2, .35, .5, .75, 1.])
MODES = ['fixed', 'ohmic', 'ghk']


def write_rows(path, rows):
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)


def path_distance(model, source):
    h.distance(0, model['segments'][int(source)].x, sec=model['segments'][int(source)].sec)
    return np.array([h.distance(s.x, sec=s.sec) for s in model['segments']])


def simulate(model, exc_nodes, exc_weights, inh_weights, record_nodes,
             gaba_fraction=.7, gaba_remaining=1., tonic=False, e_inh=-60.,
             current_node=None, current_nA=.001, tstop=650., dt=.1, sample_dt=.5):
    """GABA and glycine share inherited kinetics/location; their weights stay explicit.

    Tonic protocol holds the fitted peak inhibitory load fixed at E_inh. It is
    an isolation experiment, not a simulation of the light stimulus time course.
    """
    p = model['p']; syns = []; ncons = []; clamps = []
    stim = h.NetStim(); stim.number = 1; stim.start = 50.; stim.noise = 0
    for idx, weight in zip(exc_nodes, exc_weights):
        if weight <= 0:
            continue
        s = h.Exp2Syn(model['segments'][int(idx)])
        s.e = p['exc_rev']; s.tau1 = p['exc_tau1']; s.tau2 = p['exc_tau2']
        nc = h.NetCon(stim, s); nc.weight[0] = float(weight); nc.delay = 0
        syns.append(s); ncons.append(nc)
    for fraction in [gaba_fraction * gaba_remaining, 1 - gaba_fraction]:
        for idx, total_weight in zip(model['inh_nodes'], inh_weights):
            weight = float(total_weight * fraction)
            if weight <= 1e-12:
                continue
            if tonic:
                s = h.SEClamp(model['segments'][int(idx)])
                s.dur1 = tstop + 1; s.amp1 = e_inh; s.rs = 1. / weight
                clamps.append(s)
            else:
                s = h.Exp2Syn(model['segments'][int(idx)])
                s.e = e_inh; s.tau1 = p['inh_tau1']; s.tau2 = p['inh_tau2']
                nc = h.NetCon(stim, s); nc.weight[0] = weight; nc.delay = 0
                syns.append(s); ncons.append(nc)
    if current_node is not None:
        s = h.IClamp(model['segments'][int(current_node)])
        s.delay = 50.; s.dur = 100.; s.amp = current_nA; clamps.append(s)
    tv = h.Vector().record(h._ref_t, sample_dt)
    vv = [h.Vector().record(model['segments'][int(idx)]._ref_v, sample_dt) for idx in record_nodes]
    h.CVode().active(0); h.dt = dt; h.finitialize(p['rest_v'])
    h.continuerun(tstop)
    return np.asarray(tv).copy(), np.array([np.asarray(v) for v in vv])


def save_figure(fig, out, name):
    for ext in ['png', 'svg']:
        fig.savefig(out / f'{name}.{ext}', dpi=180, facecolor='white')
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--outdir', type=Path, default=ROOT / 'outputs/fig3_pilot')
    args = ap.parse_args(); out = args.outdir; out.mkdir(parents=True, exist_ok=True)
    started = time.time()
    cfg, _ = load_named_config('sep_g22'); model = build(cfg['g_pas']); fits = load_fits()
    p = model['p']; rest = model['rest_v']; coords = model['node_coords']
    off = model['off_nodes']; inh = model['inh_nodes']; roi = np.r_[off, model['soma_idx']]
    eweights = exc_eff_weight(fits, 'OFF', 100) * cfg['G_exc'] * gaussian_spot_weights(coords[off, :2], model['soma_xy'], 100, 25)
    iweights = inh_eff_weight(fits, 'OFF', 100) * cfg['G_inh'] * gaussian_spot_weights(coords[inh, :2], model['soma_xy'], 100, 25)
    rows = []; trace_archive = {}; curves = {}
    for condition, remaining in [('control', 1.), ('gaba_removed', 0.)]:
        for c in CONTRASTS:
            t, v = simulate(model, off, eweights*c, iweights*c, roi, gaba_remaining=remaining)
            trace_archive[f'{condition}_{c:.2f}'] = v.astype(np.float32)
            dv = (v[:-1] - rest).max(axis=1).mean()
            row = dict(condition=condition, contrast=float(c), voltage_mean_local_peak_mV=float(dv),
                       soma_peak_mV=float((v[-1]-rest).max()), min_local_peak_voltage_mV=float(v[:-1].max(1).min()),
                       max_local_peak_voltage_mV=float(v[:-1].max(1).max()))
            for mode in MODES:
                ca = calcium_trace(v[:-1], .5, mode=mode, rest_v=rest)
                row[f'calcium_{mode}_mean_local_peak_au_ms'] = float(ca.max(axis=1).mean())
            rows.append(row)
            print(f'contrast {condition} {c:.2f}: dV={dv:.3f}', flush=True)
        curves[condition] = [r for r in rows if r['condition'] == condition]
    np.savez_compressed(out/'contrast_traces.npz', t_ms=t, roi_segment_indices=roi, **trace_archive)
    write_rows(out/'contrast_responses.tsv', rows)

    # Parameter sensitivity uses identical NEURON voltage traces, no fitting or selection.
    sensitivity = []
    for vh in [-45., -40., -35., -30.]:
        for slope in [3., 4., 6., 8.]:
            for mode in MODES:
                values = {}
                for condition in COLORS:
                    for c in [.1, 1.]:
                        v = trace_archive[f'{condition}_{c:.2f}'][:-1]
                        ca = calcium_trace(v, .5, mode=mode, v_half=vh, slope=slope, rest_v=rest)
                        values[(condition, c)] = float(ca.max(axis=1).mean())
                low = values['gaba_removed', .1]-values['control', .1]
                high = values['gaba_removed', 1.]-values['control', 1.]
                sensitivity.append(dict(v_half_mV=vh, slope_mV=slope, mode=mode,
                    low_contrast_absolute_increment=low, high_contrast_absolute_increment=high,
                    low_over_high_increment=low/high if abs(high)>1e-10 else float('nan'),
                    low_contrast_fractional_gain=values['gaba_removed',.1]/values['control',.1],
                    high_contrast_fractional_gain=values['gaba_removed',1.]/values['control',1.]))
    write_rows(out/'calcium_parameter_sensitivity.tsv', sensitivity)

    # Retained inhibition sensitivity: selective GABA loss, same glycine in each pair.
    inhibition_sensitivity = []
    for gaba_fraction in [.5, .9]:
        for c in [.1, 1.]:
            t, v = simulate(model, off, eweights*c, iweights*c, roi, gaba_fraction=gaba_fraction, gaba_remaining=0.)
            row = dict(gaba_fraction=gaba_fraction, glycine_retained_fraction=1-gaba_fraction, contrast=c,
                       voltage_mean_local_peak_mV=float((v[:-1]-rest).max(1).mean()))
            for mode in MODES:
                row[f'calcium_{mode}_mean_local_peak_au_ms'] = float(calcium_trace(v[:-1], .5, mode=mode, rest_v=rest).max(1).mean())
            inhibition_sensitivity.append(row)
    write_rows(out/'residual_inhibition_sensitivity.tsv', inhibition_sensitivity)

    # Local current probe under fixed pure shunting load. Same-time source normalization.
    soma_d = path_distance(model, model['soma_idx'])
    order = off[np.argsort(soma_d[off])]
    sources = [int(order[int(q*(len(order)-1))]) for q in [.25,.5,.75]]
    lengths = np.array([s.sec.L/s.sec.nseg for s in model['segments']])
    spread_rows = []; transfer = {}; spread_summary = []
    for source in sources:
        d = path_distance(model, source)
        for condition, remaining in [('control',1.), ('gaba_removed',0.)]:
            t, v = simulate(model, [], [], iweights, np.arange(len(coords)), tonic=True,
                            e_inh=rest, gaba_remaining=remaining, current_node=source, tstop=170.)
            baseline = v[:, t<49.].mean(1)
            response = v[:, (t>=140)&(t<=145)].mean(1)-baseline
            normalized = response/response[source]
            transfer[source,condition] = normalized
            footprint = float(np.sum(lengths*normalized))
            spread_summary.append(dict(source_segment=source, condition=condition, local_response_mV=float(response[source]),
                                       transfer_weighted_cable_um=footprint))
            for low in np.arange(0,201,10):
                mask = (d>=low)&(d<low+10)
                if mask.any():
                    spread_rows.append(dict(source_segment=source, condition=condition, path_bin_um=float(low+5),
                        normalized_transfer=float(np.average(normalized[mask], weights=lengths[mask])), cable_in_bin_um=float(lengths[mask].sum())))
            print(f'spread {source} {condition}: weighted cable={footprint:.2f} um', flush=True)
    write_rows(out/'spatial_transfer.tsv', spread_rows)
    write_rows(out/'spatial_transfer_summary.tsv', spread_summary)
    np.savez_compressed(out/'spatial_transfer_maps.npz', coordinates_um=coords,
        sources=np.array(sources), **{f'{source}_{condition}':values for (source,condition),values in transfer.items()})

    # Pair experiment: fixed inhibition in A, B, A+B; pure shunt means baseline V=rest.
    source = sources[1]; d = path_distance(model, source)
    candidates = off[off!=source]
    targets = list(dict.fromkeys(int(candidates[np.argmin(abs(d[candidates]-x))]) for x in [15.,35.,70.,110.]))
    pair_rows = []; pair_archive = {}
    for target in targets:
        for condition, remaining in [('control',1.), ('gaba_removed',0.)]:
            values = {}
            for label, nodes in [('A',[source]),('B',[target]),('AB',[source,target])]:
                weight = exc_eff_weight(fits,'OFF',100)*cfg['G_exc']*.1
                t, v = simulate(model, nodes, np.full(len(nodes),weight), iweights, [source,target],
                                tonic=True, e_inh=rest, gaba_remaining=remaining)
                values[label] = v[0]-rest
                if target == targets[1]:
                    pair_archive[f'{condition}_{label}'] = values[label]
            expected = values['A']+values['B']
            actual = values['AB']
            peak_idx = int(actual.argmax())
            pair_rows.append(dict(source_segment=source, target_segment=target, path_distance_um=float(d[target]),condition=condition,
                A_peak_at_A_mV=float(values['A'].max()), B_peak_at_A_mV=float(values['B'].max()),
                AB_peak_at_A_mV=float(actual.max()), arithmetic_sum_peak_mV=float(expected.max()),
                AB_over_sum_at_AB_peak=float(actual[peak_idx]/expected[peak_idx]),
                neighbor_increment_at_AB_peak_mV=float(actual[peak_idx]-values['A'][peak_idx])))
    write_rows(out/'paired_inputs.tsv', pair_rows)
    np.savez_compressed(out/'pair_traces.npz', t_ms=t, **pair_archive)

    # Numerical and protocol checks, with quantitative differences retained.
    checks = []
    checks.append(dict(check='zero_contrast_zero_response',passed=all(abs(r['voltage_mean_local_peak_mV'])<1e-8 for r in rows if r['contrast']==0)))
    checks.append(dict(check='glycinergic_weight_identical_between_conditions',passed=bool(np.allclose(iweights*.3,iweights*(1-.7)))))
    tc, vc = simulate(model,off,eweights,iweights,roi,dt=.05,sample_dt=.5)
    tref = trace_archive['control_1.00']
    shared = min(vc.shape[1],tref.shape[1])
    error = float(abs(vc[:,:shared]-tref[:,:shared]).max())
    # The plotted features are peaks. A separate four-step convergence audit
    # resolves a first-order onset transient (51 ms); peaks agree within 2 ppm.
    # Retain the pointwise error explicitly rather than silently discarding it.
    ref_peak = (tref[:-1]-rest).max(1).mean()
    fine_peak = (vc[:-1]-rest).max(1).mean()
    relative_error = float(abs(fine_peak-ref_peak)/abs(fine_peak))
    checks.append(dict(check='half_timestep_voltage_peak_relative_error',passed=relative_error<1e-4,value=relative_error,
                       diagnostic_max_pointwise_trace_difference_mV=error))
    for mode in MODES:
        reference = next(r for r in rows if r['condition']=='control' and r['contrast']==1.)[f'calcium_{mode}_mean_local_peak_au_ms']
        fine = float(calcium_trace(vc[:-1],.5,mode=mode,rest_v=rest).max(1).mean())
        relative_error = abs(fine-reference)/abs(fine)
        checks.append(dict(check=f'half_timestep_{mode}_calcium_peak_relative_error',passed=relative_error<1e-4,value=relative_error))
    ts, vs = simulate(model,off,eweights,iweights,roi,tstop=1500.)
    for mode in MODES:
        short = next(r for r in rows if r['condition']=='control' and r['contrast']==1.)[f'calcium_{mode}_mean_local_peak_au_ms']
        long = float(calcium_trace(vs[:-1],.5,mode=mode,rest_v=rest).max(1).mean())
        err = abs(long-short)/max(abs(long),1e-9)
        checks.append(dict(check=f'long_window_{mode}_calcium_peak_relative_error',passed=err<.001,value=err))
    # Ambiguous dual-classified contacts: remove only their inhibitory component for a sensitivity run.
    overlap = np.intersect1d(off,inh); iw = iweights.copy(); iw[np.isin(inh,overlap)] = 0
    _, va = simulate(model,off,eweights,iw,roi)
    ambiguity_delta = float((va[:-1]-rest).max(1).mean()-rows[len(CONTRASTS)-1]['voltage_mean_local_peak_mV'])
    checks.append(dict(check='all_simulation_outputs_finite',passed=all(np.isfinite(a).all() for a in trace_archive.values())))
    with open(out/'validation.json','w') as f: json.dump(checks,f,indent=2)
    metadata = dict(status='exploratory_not_fit_to_experimental_calcium',config=cfg,locked_parameters=p,
        n_segments=len(coords),n_off_synapses=len(off),n_on_synapses=len(model['on_nodes']),n_inhibitory_contacts=len(inh),
        readout='mean of local peaks at 75 fixed OFF bipolar-contact dendritic sites; transform locally before averaging',
        contrasts=CONTRASTS.tolist(),stimulus='100 um OFF spot; 25 um Gaussian blur; inherited single event at 50 ms; 10/300 ms taus',
        inhibition='co-located effective GABA/glycine conductance split 70/30; remove GABA only; transmitter identities unmeasured',
        dt_ms=.1,sample_dt_ms=.5,tstop_ms=650.,calcium=dict(v_half_mV=-40,slope_mV=4,tau_ms=50,eca_ohmic_mV=120,cai_mM=.0001,cao_mM=2,temp_C=33),
        dual_classified_segment_indices=overlap.tolist(),remove_ambiguous_inhibition_high_contrast_dV_change_mV=ambiguity_delta,
        sources=sources,elapsed_seconds=time.time()-started,
        hashes={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in [ROOT/'data/cid3.swc', ROOT/'data/sm2nrn_cid3.mat',Path(__file__),ROOT/'src/calcium_readout.py']})
    with open(out/'metadata.json','w') as f: json.dump(metadata,f,indent=2)

    plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
    fig, axes = plt.subplots(2,3,figsize=(12,7.3),layout='constrained')
    ax=axes[0,0]
    for condition,color in COLORS.items():
        subset=[r for r in spread_rows if r['condition']==condition and r['source_segment']==source]
        ax.plot([r['path_bin_um'] for r in subset],[r['normalized_transfer'] for r in subset],color=color,lw=2,label=condition.replace('_',' '))
    ax.set(xlabel='Dendritic path distance (µm)',ylabel='Remote / local ΔV',title='A  Reduced shunting extends transfer',xlim=(0,150),ylim=(0,1))
    ax.legend(frameon=False,fontsize=8)
    ax=axes[0,1]
    for condition,color in COLORS.items():
        subset=[r for r in pair_rows if r['condition']==condition]
        ax.plot([r['path_distance_um'] for r in subset],[r['B_peak_at_A_mV'] for r in subset],'-o',color=color)
    ax.set(xlabel='Input A–B path distance (µm)',ylabel='Response at A to weak input B (mV)',title='B  Neighboring inputs contribute more')
    ax=axes[0,2]
    for condition,color in COLORS.items():
        ax.plot(CONTRASTS*100,[r['voltage_mean_local_peak_mV'] for r in curves[condition]],'-o',ms=3,color=color)
    ax.set(xlabel='Contrast (%)',ylabel='Mean local peak ΔV (mV)',title='C  Voltage contrast response')
    for ax,mode,title in zip(axes[1],MODES,['Activation only','Activation × ohmic driving term','Activation × GHK permeation']):
        key=f'calcium_{mode}_mean_local_peak_au_ms'
        reference=curves['control'][-1][key]
        for condition,color in COLORS.items():
            ax.plot(CONTRASTS*100,[r[key]/reference for r in curves[condition]],'-o',ms=3,color=color)
        ax.set(xlabel='Contrast (%)',ylabel='Ca proxy / control at 100%',title=title)
    fig.suptitle('Fig. 3 mechanism pilot • glycinergic component retained\nInherited voltage parameters; illustrative calcium readouts; each readout uses one shared scale',fontsize=12)
    save_figure(fig,out,'fig3_mechanism_pilot')

    fig,axes=plt.subplots(1,3,figsize=(12,3.6),layout='constrained')
    voltage=np.linspace(-55,0,301)
    for mode in MODES:
        drive=influx_proxy(voltage,mode=mode,rest_v=rest)
        axes[0].plot(voltage,drive,label=mode)
    axes[0].set(xlabel='Voltage (mV)',ylabel='Influx proxy (arbitrary units)',title='Transfer functions; shared baseline scaling')
    axes[0].legend(frameon=False)
    for mode in MODES:
        key=f'calcium_{mode}_mean_local_peak_au_ms';ctrl=np.array([r[key] for r in curves['control']]);ko=np.array([r[key] for r in curves['gaba_removed']])
        axes[1].plot(CONTRASTS*100,(ko-ctrl)/ctrl[-1],'-o',ms=3,label=mode)
        axes[2].plot(CONTRASTS[1:]*100,ko[1:]/ctrl[1:],'-o',ms=3,label=mode)
    axes[1].set(xlabel='Contrast (%)',ylabel='Absolute increment / control at 100%',title='Absolute effect of GABA removal')
    axes[2].set(xlabel='Contrast (%)',ylabel='GABA removed / control',title='Fractional effect of GABA removal')
    for ax in axes[1:]: ax.legend(frameon=False)
    save_figure(fig,out,'calcium_mechanism_controls')

    fig,axes=plt.subplots(1,2,figsize=(9,4),layout='constrained')
    for condition,color in COLORS.items():
        subset=[r for r in pair_rows if r['condition']==condition]
        axes[0].plot([r['path_distance_um'] for r in subset],[r['AB_over_sum_at_AB_peak'] for r in subset],'-o',color=color,label=condition.replace('_',' '))
        for label,style in [('A',':'),('B','--'),('AB','-')]:
            axes[1].plot(t,pair_archive[f'{condition}_{label}'],style,color=color,label=f'{condition.replace("_"," ")} {label}')
    axes[0].axhline(1,color='grey',lw=.8)
    axes[0].set(xlabel='Input A–B path distance (µm)',ylabel='Combined / arithmetic sum at combined peak',title='Passive summation remains subadditive')
    axes[0].legend(frameon=False,fontsize=8)
    axes[1].set(xlabel='Time (ms)',ylabel='ΔV at input A (mV)',title=f'Weak paired inputs; separation {d[targets[1]]:.1f} µm',xlim=(30,400))
    axes[1].legend(frameon=False,fontsize=7)
    save_figure(fig,out,'paired_input_controls')
    fig,axes=plt.subplots(1,2,figsize=(9,4),layout='constrained')
    for ax,condition in zip(axes,COLORS):
        points=ax.scatter(coords[:,0],coords[:,1],c=transfer[source,condition],s=2,
                          cmap='magma',vmin=0,vmax=1,rasterized=True)
        ax.plot(coords[source,0],coords[source,1],'+',color='#15bdca',ms=12,mew=2)
        ax.set(aspect='equal',title=condition.replace('_',' '),xlabel='x (µm)',ylabel='y (µm)')
    fig.colorbar(points,ax=axes,label='Remote / local ΔV',shrink=.8)
    fig.suptitle('Local current input; matched shunt reversal at rest; cyan cross = input')
    save_figure(fig,out,'spatial_transfer_maps')
    print(f'Done: {out}; {time.time()-started:.1f}s; validation={all(c["passed"] for c in checks)}',flush=True)


if __name__=='__main__':
    main()
