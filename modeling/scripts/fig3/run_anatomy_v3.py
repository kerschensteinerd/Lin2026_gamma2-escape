#!/usr/bin/env python3
"""Generate shared-scale v3 anatomy maps and dense contrast curves.

Uses independently selected spatial parameters and a shared, phenomenologically
calibrated activation curve. OFF and ON contacts are fixed across conditions.
"""
from pathlib import Path
import argparse, json, time, hashlib
import numpy as np
from run_fig3_pilot import (ROOT, build, load_named_config, load_fits, simulate,
    exc_eff_weight, inh_eff_weight, gaussian_spot_weights, write_rows)
from run_anatomy_maps import morphology_geometry, array_hash
from calcium_readout import calcium_trace
from protocol_v3 import simulate_v3

CONDITIONS=(("ctrl",1.,"Ctrl"),("gaba",0.,"cKO"))
DENSE=(0,5,10,15,20,25,30,40,50,60,70)
MAP_CONTRASTS=(15,30,70)


def readout(voltage,t,calcium,rest):
    peak=np.empty(len(voltage));auc=np.empty(len(voltage));win=(t>=50-1e-6)&(t<=3050+1e-6)
    for start in range(0,len(voltage),256):
        end=min(start+256,len(voltage))
        ca=calcium_trace(voltage[start:end],t[1]-t[0],rest_v=rest,mode='fixed',
                         v_half=calcium['v_half_mV'],slope=calcium['slope_mV'],tau_ms=calcium['tau_ms'])
        peak[start:end]=ca.max(axis=1)
        auc[start:end]=np.trapezoid(ca[:,win],x=t[win],axis=1)
    return voltage.max(axis=1)-rest,peak,auc


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--outdir',type=Path,default=ROOT/'outputs/fig3_v3')
    ap.add_argument('--skip-on',action='store_true')
    args=ap.parse_args();out=args.outdir;out.mkdir(parents=True,exist_ok=True)
    report=json.loads((out/'calibration.json').read_text());params=report['parameters'];ca_params=report['calcium']
    cfg,_=load_named_config('sep_g22');fits=load_fits();model=build(params['g_pas_S_cm2'])
    rest=model['rest_v'];coords=model['node_coords'];allnodes=np.arange(len(coords));maps=morphology_geometry(model)
    archive={}; rows=[]; checks=[]; hashes={};started=time.time()
    experiment=json.loads((ROOT.parent/'analysis/fig3_model/experimental_fig3.json').read_text())
    # Cached calibration traces at measured contrasts save redundant OFF runs.
    cache=np.load(out/'calibration_voltage_traces.npz')
    for pol in ['OFF']+([] if args.skip_on else ['ON']):
        nodes=model[f'{pol.lower()}_nodes'];inh=model['inh_nodes']
        ew=exc_eff_weight(fits,pol,100)*cfg['G_exc']*params['excitation_scale']*gaussian_spot_weights(coords[nodes,:2],model['soma_xy'],100,25)
        iw=inh_eff_weight(fits,pol,100)*cfg['G_inh']*params['inhibition_scale']*gaussian_spot_weights(coords[inh,:2],model['soma_xy'],100,25)
        for cond,remaining,explabel in CONDITIONS:
            for c in tuple(MAP_CONTRASTS)+tuple(c for c in DENSE if c not in MAP_CONTRASTS) if pol=='OFF' else MAP_CONTRASTS:
                key=f'{pol}_{cond}_{c}'
                full=pol=='OFF' and c in MAP_CONTRASTS
                t,v=simulate_v3(model,nodes,ew,iw,allnodes if full else nodes,params,remaining,c/100,dt=.1,sample_dt=2.,tstop=3052.)
                dv,peak,auc=readout(v,t,ca_params,rest)
                if full:
                    maps[f'spot_{cond}_{c}_dv']=dv;maps[f'spot_{cond}_{c}_ca']=peak
                    maps[f'spot_{cond}_{c}_ca_auc']=auc;maps[f'spot_{cond}_{c}_vpeak_mV']=dv+rest
                    dv,peak,auc=dv[nodes],peak[nodes],auc[nodes]
                    vv=v[nodes]
                    err=float(np.max(np.abs(vv-cache[key])))
                    checks.append(dict(check=f'{key}_full_map_matches_calibration_trace',passed=err<5e-6,max_abs_error_mV=err))
                else: vv=v
                archive[key]=vv.astype(np.float32)
                if full: np.savez_compressed(out/'anatomy_maps.npz',**maps)
                baseline_err=float(np.max(np.abs(vv[:,t<49]-rest)))
                checks.append(dict(check=f'{key}_baseline_and_finite',passed=baseline_err<1e-8 and bool(np.isfinite(v).all()),baseline_max_abs_error_mV=baseline_err))
                row=dict(polarity=pol,condition=cond,contrast_pct=c,n_model_sites=len(nodes),
                    mean_peak_dv_mV=float(dv.mean()),mean_peak_ca_au_ms=float(peak.mean()),
                    mean_ca_auc_au_ms2=float(auc.mean()),mean_peak_absolute_v_mV=float(dv.mean()+rest),
                    minimum_local_peak_absolute_v_mV=float(dv.min()+rest),maximum_local_peak_absolute_v_mV=float(dv.max()+rest),
                    model_auc_shape=float('nan'),experiment_auc_shape=float('nan'),experiment_auc_shape_sem=float('nan'))
                if c in MAP_CONTRASTS:
                    group=experiment['groups'][f'{explabel}_{pol}']['contrasts'][str(c/100)]
                    row['experiment_auc_shape']=group['normalized_mean'];row['experiment_auc_shape_sem']=group['normalized_sem']
                rows.append(row)
                hashes[key]=dict(excitation=array_hash(ew*c/100),glycine=array_hash(iw*(1 if params.get('inhibition_protocol')=='tonic_pure_shunt' else c/100)*params['glycine_fraction']))
                print(f'{key}: ΔV={dv.mean():.4f}; Ca peak={peak.mean():.4f}; Ca AUC={auc.mean():.4f}',flush=True)
                del v
    for pol in sorted({r['polarity'] for r in rows}):
        for cond,_,_ in CONDITIONS:
            rr=[r for r in rows if r['polarity']==pol and r['condition']==cond]
            lo=next(r['mean_ca_auc_au_ms2'] for r in rr if r['contrast_pct']==15)
            hi=next(r['mean_ca_auc_au_ms2'] for r in rr if r['contrast_pct']==70)
            for row in rr: row['model_auc_shape']=(row['mean_ca_auc_au_ms2']-lo)/(hi-lo)
        for c in DENSE if pol=='OFF' else MAP_CONTRASTS:
            checks.append(dict(check=f'{pol}_{c}_identical_excitation_and_retained_glycine',passed=hashes[f'{pol}_ctrl_{c}']==hashes[f'{pol}_gaba_{c}']))
    # One fine-step/fine-sample run checks the strongest response quantitatively.
    nodes=model['off_nodes'];inh=model['inh_nodes']
    ew=exc_eff_weight(fits,'OFF',100)*cfg['G_exc']*params['excitation_scale']*gaussian_spot_weights(coords[nodes,:2],model['soma_xy'],100,25)
    iw=inh_eff_weight(fits,'OFF',100)*cfg['G_inh']*params['inhibition_scale']*gaussian_spot_weights(coords[inh,:2],model['soma_xy'],100,25)
    tf,vf=simulate_v3(model,nodes,ew,iw,nodes,params,0.,.7,dt=.05,sample_dt=1.,tstop=3051.)
    measures=readout(vf,tf,ca_params,rest)
    ref=next(r for r in rows if r['polarity']=='OFF' and r['condition']=='gaba' and r['contrast_pct']==70)
    for key,fine in zip(['mean_peak_dv_mV','mean_peak_ca_au_ms','mean_ca_auc_au_ms2'],measures):
        error=abs(ref[key]-float(fine.mean()))/max(abs(float(fine.mean())),1e-12)
        checks.append(dict(check=f'half_step_and_sampling_{key}',passed=error<.005,relative_error=error,tolerance=.005))
    lengths=maps['cable_lengths'];edges=maps['geometry_edges'];edgeindex=maps['edge_segment_index']
    lengths_edges=np.bincount(edgeindex,weights=np.linalg.norm(edges[:,1]-edges[:,0],axis=1),minlength=len(lengths))
    geometry_error=float(np.max(np.abs(lengths_edges-lengths)))
    checks.append(dict(check='geometry_cable_length_preserved',passed=geometry_error<1e-6,max_abs_error_um=geometry_error))
    checks.append(dict(check='all_map_arrays_finite',passed=all(bool(np.isfinite(v).all()) for v in maps.values())))
    checks.append(dict(check='glycine_nonzero',passed=params['glycine_fraction']>0))
    np.savez_compressed(out/'anatomy_maps.npz',**maps)
    np.savez_compressed(out/'contrast_voltage_traces.npz',t_ms=t,**archive)
    write_rows(out/'contrast_responses.tsv',rows)
    gains=[]
    for c in MAP_CONTRASTS:
        r0=next(r for r in rows if r['polarity']=='OFF' and r['condition']=='ctrl' and r['contrast_pct']==c)
        r1=next(r for r in rows if r['polarity']=='OFF' and r['condition']=='gaba' and r['contrast_pct']==c)
        gains.append(dict(contrast_pct=c,**{f'{k}_{op}':r1[k]/r0[k] if op=='gain' else r1[k]-r0[k] for k in ['mean_peak_dv_mV','mean_peak_ca_au_ms','mean_ca_auc_au_ms2'] for op in ['gain','increment']}))
    metadata=dict(status='illustrative_calibrated_OFF_shape_not_identified_biophysics',parameters=params,calcium=ca_params,
        geometry='Unmodified real section pt3d polylines split at segment boundaries.',n_segments=len(coords),n_off_input_sites=len(model['off_nodes']),n_on_input_sites=len(model['on_nodes']),
        dt_ms=.1,sample_dt_ms=2.,tstop_ms=3052.,stimulus_onset_ms=50.,auc_window_ms=[50,3050],
        spot=dict(diameter_um=100.,gaussian_blur_sigma_um=25.,polarity='OFF',inhibitory_reversal_mV=rest if params.get('inhibition_protocol')=='tonic_pure_shunt' else -60.,protocol='Sustained 3-s excitation scaled with contrast against fixed pure-shunt inhibition.' if params.get('excitation_protocol')=='sustained_3s' else 'One inherited 10/300-ms conductance event.'),
        spatial_sampling='Mean of local responses at fixed 75 OFF or 30 ON bipolar-input contacts; not registered imaging sublaminae.',
        interpretation='OFF intermediate-contrast shape was calibrated; ON is held out. All primary response magnitudes use a shared scale. Normalization to condition-specific 15/70 endpoints is only for the labeled experimental shape comparison.',
        map_ca_units='arbitrary influx units * ms, local peak',map_ca_auc_units='arbitrary influx units * ms^2, baseline-subtracted signed AUC',
        spot_summary=[dict(r,mean_off_peak_dv_mV=r['mean_peak_dv_mV'],mean_off_peak_ca_au_ms=r['mean_peak_ca_au_ms']) for r in rows if r['polarity']=='OFF' and r['contrast_pct'] in MAP_CONTRASTS],
        contrast_rows=rows,gains=gains,elapsed_seconds=time.time()-started,
        checks_passed=sum(c['passed'] for c in checks),checks_total=len(checks),
        hashes={str(Path(__file__).resolve()):hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),str(out/'calibration.json'):hashlib.sha256((out/'calibration.json').read_bytes()).hexdigest()})
    def json_safe(value):
        if isinstance(value,dict): return {k:json_safe(v) for k,v in value.items()}
        if isinstance(value,list): return [json_safe(v) for v in value]
        if isinstance(value,float) and not np.isfinite(value): return None
        return value
    (out/'metadata.json').write_text(json.dumps(json_safe(metadata),indent=2,allow_nan=False)+'\n')
    (out/'validation.json').write_text(json.dumps(checks,indent=2)+'\n')
    print(f"Complete: {metadata['checks_passed']}/{len(checks)} checks; {time.time()-started:.1f}s",flush=True)

if __name__=='__main__':main()
