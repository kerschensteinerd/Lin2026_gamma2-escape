#!/usr/bin/env python3
"""Reproducible fixed-voltage calcium sensitivity; no NEURON simulation or fitting.

AUC is a linear functional of baseline-subtracted influx. Its backward weights
exactly reproduce the archived left-held filter plus trapezoidal integration.
Repeated *identical* float32 voltage samples are grouped without quantization.
The nonlinear local activation is applied before any spatial averaging.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
import time
import numpy as np
import scipy
from scipy.special import expit
from scipy.interpolate import PchipInterpolator

ROOT = Path(__file__).resolve().parents[2]
MODEL = ROOT
SOURCE = MODEL / 'outputs/fig4_v9_retained_inhibition_20260916'
sys.path.insert(0, str(MODEL / 'src'))
from calcium_readout import calcium_trace
CONTRAST = np.array([0,5,10,15,20,25,30,40,50,60,70])
FRACTIONS = np.array([1.,.20,.25,.30])
REST, TAU = -45., 50.
AUC_SCALE = TAU * 3000.
TECHNICAL_FLOOR = 1e-6 * AUC_SCALE


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def json_write(path, obj):
    def clean(x):
        if isinstance(x, dict): return {str(k):clean(v) for k,v in x.items()}
        if isinstance(x, (list,tuple)): return [clean(v) for v in x]
        if isinstance(x, np.ndarray): return clean(x.tolist())
        if isinstance(x, (np.integer,)): return int(x)
        if isinstance(x, (np.bool_,)): return bool(x)
        if isinstance(x, (float,np.floating)): return float(x) if np.isfinite(x) else None
        return x
    path.write_text(json.dumps(clean(obj),indent=2,allow_nan=False)+'\n')


def weights_for_auc(t):
    use=np.flatnonzero((t>=50-1e-6)&(t<=3050+1e-6))
    trap=np.zeros(len(t))
    delta=np.diff(t[use])
    trap[use[:-1]]+=delta/2
    trap[use[1:]]+=delta/2
    decay=np.exp(-(t[1]-t[0])/TAU)
    gain=TAU*(-np.expm1(-(t[1]-t[0])/TAU))
    adj=np.zeros(len(t))
    adj[-1]=trap[-1]
    for i in range(len(t)-2,-1,-1): adj[i]=trap[i]+decay*adj[i+1]
    weights=np.zeros(len(t)); weights[:-1]=gain*adj[1:]
    return weights


def voltage_groups(v, weights):
    vv=v.astype(float).reshape(-1)
    unique,inverse=np.unique(vv,return_inverse=True)
    totals=np.bincount(inverse,weights=np.broadcast_to(weights,v.shape).reshape(-1)/v.shape[0])
    keep=totals>0
    return unique[keep],totals[keep]


def evaluate_auc(groups, pairs):
    output=np.empty((len(pairs),len(groups)))
    for j,(v,w) in enumerate(groups):
        for start in range(0,len(pairs),128):
            p=pairs[start:start+128]
            influx=expit((v[None,:]-p[:,0,None])/p[:,1,None])-expit((REST-p[:,0])/p[:,1])[:,None]
            output[start:start+len(p),j]=influx@w
    return output.reshape(len(pairs),len(FRACTIONS),len(CONTRAST))


def shape_metrics(auc, floor=TECHNICAL_FLOOR):
    il,im,ih=[int(np.flatnonzero(CONTRAST==c)[0]) for c in (15,30,70)]
    span=auc[:,:,ih]-auc[:,:,il]
    shape=np.divide(auc-auc[:,:,il,None],span[:,:,None],out=np.full_like(auc,np.nan),where=span[:,:,None]>floor)
    key=auc[:,:,[il,im,ih]]
    ratio=np.divide(key[:,1:,:],key[:,0,None,:],out=np.full_like(key[:,1:,:],np.nan),where=key[:,0,None,:]>floor)
    delta=key[:,1:,:]-key[:,0,None,:]
    pref=np.divide(ratio[:,:,:2],ratio[:,:,2,None],out=np.full_like(ratio[:,:,:2],np.nan),where=np.isfinite(ratio[:,:,2,None]))
    return span,shape,ratio,delta,pref


def shifts_for_auc(auc):
    shifts=np.full((len(auc),3,len(CONTRAST)),np.nan)
    for i,a in enumerate(auc):
        ctrl=a[0]
        if np.any(np.diff(ctrl)<=0) or ctrl[-1]-ctrl[0]<=TECHNICAL_FLOOR: continue
        inverse=PchipInterpolator(ctrl,CONTRAST,extrapolate=False)
        for j in range(3):
            reduced=a[j+1]
            valid=(CONTRAST>0)&(reduced>TECHNICAL_FLOOR)&(reduced>=ctrl[0])&(reduced<=ctrl[-1])
            shifts[i,j,valid]=inverse(reduced[valid])/CONTRAST[valid]
    counts=np.sum(np.isfinite(shifts),axis=2)
    mean=np.divide(np.nansum(shifts,axis=2),counts,out=np.full(counts.shape,np.nan),where=counts>0)
    return shifts,counts,mean


def spatial_metrics(maps,pairs):
    distance=maps['path_distance_from_source_um'];length=maps['cable_lengths']
    source=int(maps['source_segment'])
    bins=[(distance>=low)&(distance<low+5) for low in np.arange(0,distance.max()+5,5)]
    bins=[b for b in bins if np.any(b)]
    x=np.array([0.]+[float(np.average(distance[b],weights=length[b])) for b in bins])
    profiles=np.full((len(pairs),2,len(x)),np.nan)
    cross=np.full((len(pairs),2),np.nan)
    cable=np.full_like(cross,np.nan);source_ca=np.empty_like(cross)
    for j,condition in enumerate(('ctrl','reduced')):
        v=REST+maps[f'local_{condition}_finite_dv']
        for start in range(0,len(pairs),128):
            p=pairs[start:start+128]
            ca=TAU*(expit((v[None,:]-p[:,0,None])/p[:,1,None])-expit((REST-p[:,0])/p[:,1])[:,None])
            source_ca[start:start+len(p),j]=ca[:,source]
            norm=np.divide(ca,ca[:,source,None],out=np.full_like(ca,np.nan),where=ca[:,source,None]>1e-6*TAU)
            profile=np.column_stack([np.ones(len(p))]+[np.average(norm[:,b],weights=length[b],axis=1) for b in bins])
            bad=~np.isfinite(norm[:,0]);profile[bad]=np.nan
            profiles[start:start+len(p),j]=profile
            cable[start:start+len(p),j]=norm@length
            for k,y in enumerate(profile):
                at=np.flatnonzero((y[:-1]>1/np.e)&(y[1:]<=1/np.e))
                if len(at):
                    n=at[0];cross[start+k,j]=x[n]+(1/np.e-y[n])*(x[n+1]-x[n])/(y[n+1]-y[n])
    return x,profiles,cross,cable,source_ca


def summaries(auc,shape,ratio,delta,pref,shifts,count,shift_mean,cross,n_grid):
    def extent(a):
        v=np.asarray(a);v=v[np.isfinite(v)]
        return dict(n=int(v.size),min=float(v.min()) if len(v) else None,median=float(np.median(v)) if len(v) else None,max=float(v.max()) if len(v) else None)
    cases=[]
    for j,r in enumerate(FRACTIONS[1:]):
        p=pref[:n_grid,j]; valid=np.isfinite(p).all(axis=1)
        interior=shape[:n_grid,j+1,6]-shape[:n_grid,0,6]
        cases.append(dict(retained_fraction=r,valid_preference_points=int(valid.sum()),
             both_low_intermediate_fractional_preference_count=int(np.sum(valid&np.all(p>1,axis=1))),
             low_preference_count=int(np.sum(p[:,0]>1)),intermediate_preference_count=int(np.sum(p[:,1]>1)),
             preference15_over70=extent(p[:,0]),preference30_over70=extent(p[:,1]),
             shape30_separation=extent(interior),positive_shape30_separation_count=int(np.sum(interior>0)),
             gains={str(c):extent(ratio[:n_grid,j,k]) for k,c in enumerate((15,30,70))},
             absolute_increment={str(c):extent(delta[:n_grid,j,k]) for k,c in enumerate((15,30,70))},
             shift_mean=extent(shift_mean[:n_grid,j]),valid_shift_count_range=extent(count[:n_grid,j]),
             absolute_enhancement_all_nonzero_contrasts_count=int(np.sum(np.all(auc[:n_grid,j+1,1:]>auc[:n_grid,0,1:],axis=1)))))
    nominal=[]
    for i,label in [(n_grid,'original_exact'),(int(np.flatnonzero((PARAM_PAIRS[:n_grid,0]==-40)&np.isclose(PARAM_PAIRS[:n_grid,1],4))[0]),'historical_minus40_4')]:
        nominal.append(dict(label=label,v_half_mV=PARAM_PAIRS[i,0],slope_mV=PARAM_PAIRS[i,1],mean_auc=auc[i],shape30=shape[i,:,6],gains_15_30_70=ratio[i],increments_15_30_70=delta[i],preference_15_30_over70=pref[i],calcium_only_shift_mean=shift_mean[i],shift_counts=count[i],spatial_1e_um=cross[i]))
    return dict(grid_points=n_grid,cases=cases,reference_points=nominal,
                nominal_spatial=dict(control_1e_um=extent(cross[:n_grid,0]),reduced_1e_um=extent(cross[:n_grid,1]),
                ratio=extent(cross[:n_grid,1]/cross[:n_grid,0]),increase_count=int(np.sum(cross[:n_grid,1]>cross[:n_grid,0]))))


def main():
    global PARAM_PAIRS
    ap=argparse.ArgumentParser();ap.add_argument('--outdir',type=Path,default=Path(__file__).parent/'reproduction');args=ap.parse_args()
    out=args.outdir.resolve()
    if out.exists() and any(out.iterdir()): raise SystemExit('Refusing nonempty output directory')
    protocol=Path(__file__).with_name('protocol.md')
    if not protocol.is_file(): raise SystemExit('Save protocol before production run')
    out.mkdir(parents=True,exist_ok=True);start=time.time()
    files=[SOURCE/'nominal25/contrast_voltage_traces.npz',SOURCE/'sensitivity_voltage_traces.npz',SOURCE/'nominal25/finite_local_maps.npz',SOURCE/'nominal25/figure_data.json',SOURCE/'contrast_sensitivity.tsv',SOURCE/'sensitivity_summary.json',SOURCE/'provenance.json',MODEL/'src/calcium_readout.py',protocol,Path(__file__)]
    # Publication export: manuscript hashes remain historical provenance only.
    protected=files
    hashes={str(p.relative_to(ROOT)):sha(p) for p in protected}
    json_write(out/'input_hashes.json',hashes)
    nominal=np.load(files[0]); sensitivity=np.load(files[1]); maps=np.load(files[2])
    t=nominal['t_ms'];assert np.array_equal(t,sensitivity['t_ms'])
    weights=weights_for_auc(t)
    traces=[];labels=[];groups=[]
    for f in FRACTIONS:
        for c in CONTRAST:
            key=f'OFF_ctrl_{c}' if f==1 else f'OFF_residual{round(f*100)}_{c}'
            v=nominal[key] if f==1 else sensitivity[key]
            assert v.shape==(75,len(t)) and np.isfinite(v).all()
            assert np.max(np.abs(v[:,t<49]-REST))<1e-8
            traces.append(v);labels.append(key);groups.append(voltage_groups(v,weights))
    half=np.arange(-45,-19.999,.5);slope=np.round(np.arange(1,6.001,.1),10)
    hh,kk=np.meshgrid(half,slope,indexing='ij');grid=np.column_stack([hh.ravel(),kk.ravel()]);n_grid=len(grid)
    prior=json.loads(files[5].read_text());ca=prior['calcium']
    PARAM_PAIRS=np.vstack([grid,[ca['v_half_mV'],ca['slope_mV']]])
    json_write(out/'run_manifest.json',dict(status='running',created_utc=__import__('datetime').datetime.now(__import__('datetime').timezone.utc).isoformat(),n_grid=n_grid,half_grid=half,slope_grid=slope,fractions=FRACTIONS,contrasts=CONTRAST,rest_mV=REST,tau_ms=TAU,readout='fixed',auc_window_ms=[50,3050],technical_auc_floor=TECHNICAL_FLOOR,protocol_sha256=sha(protocol),source_voltage_dtypes=sorted(set(str(v.dtype) for v in traces)),method='Exact backward AUC weights, identical-value aggregation; activation before spatial averaging',unique_voltage_counts=dict(zip(labels,[len(g[0]) for g in groups])),runtime=dict(python=sys.executable,numpy=np.__version__,scipy=scipy.__version__)))
    print(f'Evaluating {n_grid} grid points plus exact reference; {len(traces)} condition/contrast traces',flush=True)
    auc=evaluate_auc(groups,PARAM_PAIRS)
    checks=[]
    def check(name,passed,**details):
        checks.append(dict(check=name,passed=bool(passed),**details))
        if not passed:
            json_write(out/'validation.json',checks);raise AssertionError(name)
    check('finite_nonnegative_auc',np.isfinite(auc).all() and auc.min()>=-1e-10)
    check('exact_zero_contrast',np.max(np.abs(auc[:,:,0]))<1e-10)
    win=(t>=50-1e-6)&(t<=3050+1e-6)
    validate_points=[0,50,1300,2550,2600,n_grid]
    for pi in validate_points:
        errors=[]
        for ti,v in enumerate(traces):
            q=calcium_trace(v,float(t[1]-t[0]),tau_ms=TAU,mode='fixed',v_half=PARAM_PAIRS[pi,0],slope=PARAM_PAIRS[pi,1],rest_v=REST)
            direct=float(np.trapezoid(q[:,win],x=t[win],axis=1).mean())
            calculated=auc[pi].reshape(-1)[ti]
            errors.append(abs(direct-calculated)/max(abs(direct),1.))
        check(f'adjoint_vs_full_filter_point{pi}',max(errors)<1e-10,max_scaled_error=max(errors))
    reference=json.loads(files[3].read_text())
    with files[4].open() as s: oldrows=list(csv.DictReader(s,delimiter='\t'))
    oldauc=np.empty((4,11));oldauc[0]=reference['contrast']['ctrl']['mean_ca_auc_au_ms2']
    for j,f in enumerate(FRACTIONS[1:],start=1):
        oldauc[j]=[float(next(r for r in oldrows if float(r['retained_effective_inhibition_fraction'])==f and int(r['contrast_pct'])==c)['mean_ca_auc_au_ms2']) for c in CONTRAST]
    err=np.abs(auc[-1]-oldauc)/np.maximum(abs(oldauc),1.)
    check('archived_exact_parameters_reproduce_original_auc_float32_tolerance',err.max()<1e-5,max_scaled_error=float(err.max()),explanation='Archived voltages float32; original summary computed before float32 storage')
    span,shape,ratio,delta,pref=shape_metrics(auc)
    shifts,count,shift_mean=shifts_for_auc(auc)
    x,profiles,cross,cable,source_ca=spatial_metrics(maps,PARAM_PAIRS)
    check('spatial_original_reference',np.max(np.abs(cross[-1]-[18.77733921143409,34.337163202424904]))<1e-9,observed_um=cross[-1])
    check('finite_all_spatial_profiles',np.isfinite(profiles).all())
    targets=np.array(prior['prior_OFF_shape30_targets'])
    shape_error=(shape[:,0,6,None]-targets[0])**2+(shape[:,1:,6]-targets[1])**2
    np.savez_compressed(out/'calcium_sweep.npz',v_half_grid_mV=half,slope_grid_mV=slope,parameter_pairs=PARAM_PAIRS,retained_fractions=FRACTIONS,contrast_pct=CONTRAST,mean_auc=auc,auc_span_15_70=span,normalized_auc=shape,gain_15_30_70=ratio,increment_15_30_70=delta,preference_15_30_over70=pref,shape30_sse=shape_error,calcium_shift_by_contrast=shifts,calcium_shift_count=count,calcium_shift_mean=shift_mean,spatial_distance_um=x,spatial_profiles_ctrl25=profiles,spatial_1e_ctrl25_um=cross,spatial_weighted_cable_ctrl25_um=cable,spatial_source_ca_ctrl25=source_ca)
    with (out/'calcium_sweep_metrics.tsv').open('w') as stream:
        fields=['point_type','v_half_mV','slope_mV','retained_fraction']+[f'{p}{c}' for p in ['control_auc_','reduced_auc_','delta_auc_','gain_'] for c in [15,30,70]]+['preference15_over70','preference30_over70','control_shape30','reduced_shape30','shape30_difference','shape30_sse','calcium_shift_mean','valid_shift_count','valid_shift_contrasts','control_spatial_1e_um','reduced25_spatial_1e_um']
        writer=csv.DictWriter(stream,fieldnames=fields,delimiter='\t');writer.writeheader()
        for i,p in enumerate(PARAM_PAIRS):
            for j,f in enumerate(FRACTIONS[1:]):
                row=dict(point_type='grid' if i<n_grid else 'original_exact',v_half_mV=p[0],slope_mV=p[1],retained_fraction=f,preference15_over70=pref[i,j,0],preference30_over70=pref[i,j,1],control_shape30=shape[i,0,6],reduced_shape30=shape[i,j+1,6],shape30_difference=shape[i,j+1,6]-shape[i,0,6],shape30_sse=shape_error[i,j],calcium_shift_mean=shift_mean[i,j],valid_shift_count=count[i,j],valid_shift_contrasts=','.join(map(str,CONTRAST[np.isfinite(shifts[i,j])])) ,control_spatial_1e_um=cross[i,0],reduced25_spatial_1e_um=cross[i,1])
                for k,c in enumerate([15,30,70]):
                    ci=int(np.flatnonzero(CONTRAST==c)[0]);row.update({f'control_auc_{c}':auc[i,0,ci],f'reduced_auc_{c}':auc[i,j+1,ci],f'delta_auc_{c}':delta[i,j,k],f'gain_{c}':ratio[i,j,k]})
                writer.writerow({k:'' if isinstance(v,(float,np.floating)) and not np.isfinite(v) else v for k,v in row.items()})
    summary=summaries(auc,shape,ratio,delta,pref,shifts,count,shift_mean,cross,n_grid)
    summary['technical_floor_sensitivity']=[]
    for relative in [1e-8,1e-6,1e-4]:
        _,ss,rr,dd,pp=shape_metrics(auc,relative*AUC_SCALE)
        summary['technical_floor_sensitivity'].append(dict(relative_to_tau_duration=relative,auc_floor=relative*AUC_SCALE,valid_preference_points=np.sum(np.isfinite(pp[:n_grid]).all(axis=2),axis=0),both_preference_count=np.sum(np.isfinite(pp[:n_grid]).all(axis=2)&np.all(pp[:n_grid]>1,axis=2),axis=0)))
    summary['interpretation']='Fractions of this declared uniform grid are descriptive coverage, not probabilities, posterior mass, or biological confidence. Calcium-only horizontal shift differs from prior combined voltage+calcium marker. No fitting or screening for agreement was performed.'
    json_write(out/'summary.json',summary)
    finalhash={str(p.relative_to(ROOT)):sha(p) for p in protected}
    check('numerical_inputs_and_sources_unchanged',hashes==finalhash)
    json_write(out/'validation.json',checks)
    manifest=json.loads((out/'run_manifest.json').read_text());manifest.update(status='complete',elapsed_seconds=time.time()-start,checks_passed=sum(c['passed'] for c in checks),output_files={p.name:sha(p) for p in out.iterdir() if p.is_file() and p.name!='run_manifest.json'})
    json_write(out/'run_manifest.json',manifest)
    print(f'COMPLETE {len(checks)} checks, {time.time()-start:.1f}s. {out}',flush=True)

if __name__=='__main__': main()
