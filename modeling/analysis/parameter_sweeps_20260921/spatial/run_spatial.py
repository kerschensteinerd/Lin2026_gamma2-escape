#!/usr/bin/env python3
"""Prespecified Figure 4 spatial controls, without fitting or source mutation."""
from __future__ import annotations
import argparse, csv, hashlib, json, sys, time
from pathlib import Path
import numpy as np
import scipy
from scipy.special import expit
from scipy.interpolate import PchipInterpolator
import neuron
from neuron import h

OUT = Path(__file__).resolve().parent
SCRIPT_DIR = OUT
WORK = SCRIPT_DIR.parents[2]
REPO = WORK
sys.path.insert(0, str(REPO / 'scripts/fig3'))
from run_fig3_pilot import build, load_named_config, load_fits, exc_eff_weight, inh_eff_weight, gaussian_spot_weights, path_distance
from run_anatomy_maps import morphology_geometry
from run_anatomy_v3 import readout
from protocol_v3 import simulate_v3

PRIOR = REPO/'outputs/fig4_v9_retained_inhibition_20260916'
FRACTIONS=(.20,.25,.30)
STRENGTHS=(.5,1.,2.)
CONTRASTS=(0,5,10,15,20,25,30,40,50,60,70)
CANONICAL=4705
SOURCE_G=.012785038085113628
THRESHOLD=.15

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def arrsha(a): return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()
def json_write(p,x): Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def tsv(p,rows):
    with Path(p).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def loadtsv(p):
    with Path(p).open() as f:return list(csv.DictReader(f,delimiter='\t'))

def make_model():
    provenance=json.loads((PRIOR/'provenance.json').read_text())
    model=build(provenance['parameters']['g_pas_S_cm2'])
    cfg,_=load_named_config('sep_g22');fits=load_fits()
    xy=model['node_coords'];off=model['off_nodes'];inh=model['inh_nodes'];params=provenance['parameters']
    ew=exc_eff_weight(fits,'OFF',100)*cfg['G_exc']*params['excitation_scale']*gaussian_spot_weights(xy[off,:2],model['soma_xy'],100,25)
    iw=inh_eff_weight(fits,'OFF',100)*cfg['G_inh']*params['inhibition_scale']*gaussian_spot_weights(xy[inh,:2],model['soma_xy'],100,25)
    assert abs(iw.sum()*1000-provenance['total_control_inhibition_nS'])<1e-10
    assert abs(ew.sum()*1000-provenance['total_unit_contrast_excitation_nS'])<1e-10
    return model,provenance,ew,iw

def protected():
    files=[*sorted((REPO/'data').rglob('*')),*sorted((REPO/'src').glob('*.py')),
        *[REPO/'scripts/fig3'/n for n in ('run_fig3_pilot.py','run_anatomy_maps.py','run_anatomy_v3.py','protocol_v3.py')],
        REPO/'scripts/fig4/run_retained_inhibition_v9.py',SCRIPT_DIR.parent/'protocol.md',
        *sorted(PRIOR.rglob('*'))]
    return {str(p.relative_to(WORK)):sha(p) for p in files if p.is_file()}

def configure_shunts(model,iw,leak):
    areas=np.array([h.area(s.x,sec=s.sec) for s in model['segments']])
    pernode=np.bincount(model['inh_nodes'],weights=iw,minlength=len(model['segments']))
    for i,s in enumerate(model['segments']):s.g_pas=leak+100*pernode[i]/areas[i];s.e_pas=model['rest_v']
    h.finitialize(model['rest_v'])

def clear_shunts(model,leak):
    for s in model['segments']:s.g_pas=leak;s.e_pas=model['rest_v']

def matched_weights(iw,order,fraction):
    result=np.zeros_like(iw);remaining=float(iw.sum()*fraction)
    for j in order:
        take=min(float(iw[j]),remaining)
        if take<=0:break
        result[j]=take;remaining-=take
    assert abs(result.sum()-iw.sum()*fraction)<1e-13
    assert np.all(result>=0) and np.all(result<=iw)
    return result

def make_layouts(model,iw):
    soma_dist=path_distance(model,model['soma_idx'])
    nodes=model['inh_nodes'];idx=np.arange(len(nodes))
    layouts=[dict(layout_id='control',family='control',fraction=1.,seed=None,weights=iw.copy(),order=None)]
    for fraction in FRACTIONS:
        prefix=f'f{round(fraction*100)}'
        layouts.append(dict(layout_id=f'{prefix}_uniform',family='uniform',fraction=fraction,seed=None,weights=iw*fraction,order=None))
        for seed in range(20):
            order=np.random.Generator(np.random.PCG64(seed)).permutation(len(iw))
            layouts.append(dict(layout_id=f'{prefix}_random{seed:02d}',family='random',fraction=fraction,seed=seed,weights=matched_weights(iw,order,fraction),order=order))
        for family,sign in [('proximal',1),('distal',-1)]:
            order=np.lexsort((idx,nodes,sign*soma_dist[nodes]))
            layouts.append(dict(layout_id=f'{prefix}_{family}',family=family,fraction=fraction,seed=None,weights=matched_weights(iw,order,fraction),order=order))
    assert len(layouts)==70
    return layouts,soma_dist

def export_layouts(model,iw,layouts,soma_dist):
    summary=[];contacts=[]
    for lay in layouts:
        weights=lay['weights'];partial=(weights>0)&(weights<iw*(1-1e-12))
        summary.append(dict(layout_id=lay['layout_id'],family=lay['family'],retained_fraction=lay['fraction'],seed=lay['seed'],
            total_inhibition_nS=float(weights.sum()*1000),target_inhibition_nS=float(iw.sum()*lay['fraction']*1000),
            budget_error_nS=float((weights.sum()-iw.sum()*lay['fraction'])*1000),n_contacts_nonzero=int((weights>0).sum()),
            n_contacts_partial=int(partial.sum()),partial_contact_indices=','.join(str(i) for i in np.flatnonzero(partial)),
            weights_sha256=arrsha(weights),weighted_mean_soma_path_um=float(np.average(soma_dist[model['inh_nodes']],weights=weights))))
        rank=np.empty(len(iw),int)
        if lay['order'] is not None:rank[lay['order']]=np.arange(len(iw))
        for j,(node,w,w0) in enumerate(zip(model['inh_nodes'],weights,iw)):
            contacts.append(dict(layout_id=lay['layout_id'],contact_index=j,segment_index=int(node),
                soma_path_um=float(soma_dist[node]),control_conductance_nS=float(w0*1000),retained_conductance_nS=float(w*1000),
                retention_ratio=float(w/w0) if w0 else None,selection_rank=int(rank[j]) if lay['order'] is not None else None))
    tsv(OUT/'layout_summary.tsv',summary);tsv(OUT/'layout_contact_weights.tsv',contacts)
    np.savez_compressed(OUT/'layout_weights.npz',layout_ids=np.array([x['layout_id'] for x in layouts]),weights_uS=np.stack([x['weights'] for x in layouts]),
        inh_nodes=model['inh_nodes'],control_weights_uS=iw)
    return summary

def profile_metrics(y,d,length):
    bin_id=np.floor(d/5).astype(int)
    w=np.bincount(bin_id,weights=length);nonzero=w>0
    distance=np.r_[0.,np.bincount(bin_id,weights=length*d)[nonzero]/w[nonzero]]
    mean=np.r_[1.,np.bincount(bin_id,weights=length*y)[nonzero]/w[nonzero]]
    cross=np.flatnonzero((mean[:-1]>1/np.e)&(mean[1:]<=1/np.e))
    value=None
    if len(cross):
        i=cross[0];value=float(distance[i]+(1/np.e-mean[i])*(distance[i+1]-distance[i])/(mean[i+1]-mean[i]))
    return dict(crossing_1e_um=value,crossing_observed=value is not None,cable_integral_um=float(np.sum(length*y)),
        max_path_um=float(d.max()),profile_distance_um=distance,profile_response=mean)

def finite_focal(model,provenance,lay,source,transfer,rin,distance,length,strength,source_index=None,soma_distance=None):
    rest=model['rest_v'];ca=provenance['calcium'];g=SOURCE_G*strength
    source_dv=(model['p']['exc_rev']-rest)*g*rin/(1+g*rin);dv=transfer*source_dv
    q=ca['tau_ms']*(expit((rest+dv-ca['v_half_mV'])/ca['slope_mV'])-expit((rest-ca['v_half_mV'])/ca['slope_mV']))
    vmet=profile_metrics(transfer,distance,length);qsource=float(q[source]);technical=qsource>1e-6*ca['tau_ms']
    qmet=profile_metrics(q/qsource,distance,length) if technical else None
    row=dict(layout_id=lay['layout_id'],family=lay['family'],retained_fraction=lay['fraction'],seed=lay['seed'],
        source_index=source_index,source_segment=int(source),source_soma_path_um=float(soma_distance[source]) if soma_distance is not None else None,
        focal_strength_multiplier=strength,focal_conductance_nS=g*1000,source_Rin_MOhm=rin,source_dv_mV=float(source_dv),source_calcium_au_ms=qsource,
        voltage_1e_um=vmet['crossing_1e_um'],voltage_crossing_observed=vmet['crossing_observed'],voltage_cable_integral_um=vmet['cable_integral_um'],
        calcium_1e_um=qmet['crossing_1e_um'] if technical else None,calcium_crossing_observed=qmet['crossing_observed'] if technical else False,
        calcium_cable_integral_um=qmet['cable_integral_um'] if technical else None,calcium_source_technical_ok=technical,
        calcium_absolute_cable_integral_au_ms_um=float(np.sum(length*q)),max_path_um=vmet['max_path_um'])
    assert np.isfinite(dv).all() and np.isfinite(q).all() and q.min()>=-1e-10
    assert np.max(dv)<=source_dv*(1+1e-10) and abs(transfer[source]-1)<1e-9
    return row,vmet,qmet,dv

def focal_stage():
    started=time.time();model,prov,ew,iw=make_model();layouts,soma_d=make_layouts(model,iw)
    if not (OUT/'source_hashes_before.json').exists():json_write(OUT/'source_hashes_before.json',protected())
    else:assert json.loads((OUT/'source_hashes_before.json').read_text())==protected()
    lsummary=export_layouts(model,iw,layouts,soma_d)
    geometry=morphology_geometry(model);length=geometry['cable_lengths'];off=model['off_nodes'];n=len(off)
    assert n==75 and CANONICAL in off
    distance=np.stack([path_distance(model,int(source)) for source in off]);canon_i=list(off).index(CANONICAL)
    np.savez_compressed(OUT/'source_geometry.npz',**geometry,source_path_distances_um=distance,soma_path_distances_um=soma_d)
    tsv(OUT/'source_sites.tsv',[dict(source_index=i,segment_index=int(s),soma_path_um=float(soma_d[s]),
        x_um=float(model['node_coords'][s,0]),y_um=float(model['node_coords'][s,1]),z_um=float(model['node_coords'][s,2]),
        is_original_source=bool(s==CANONICAL)) for i,s in enumerate(off)])
    allsites=[];canonical=[];profiles=[];cache={};checks=[];validation_cases=[]
    uniform=[lay for lay in layouts if lay['family'] in ('control','uniform')]
    for li,lay in enumerate(layouts):
        configure_shunts(model,lay['weights'],prov['parameters']['g_pas_S_cm2'])
        sources=range(n) if lay in uniform else [canon_i]
        laymaps=[];laysource=[];layrin=[]
        for si in sources:
            source=int(off[si]);imp=h.Impedance();s=model['segments'][source];imp.loc(s.x,sec=s.sec);imp.compute(0.)
            rin=float(imp.input(s.x,sec=s.sec));transfer=np.array([imp.transfer(seg.x,sec=seg.sec) for seg in model['segments']])/rin
            laymaps.append(transfer);laysource.append(source);layrin.append(rin)
            for strength in STRENGTHS:
                row,vm,cm,dv=finite_focal(model,prov,lay,source,transfer,rin,distance[si],length,strength,si,soma_d)
                if lay in uniform:allsites.append(row)
                if source==CANONICAL:canonical.append(row)
                for j,(x,y,z) in enumerate(zip(vm['profile_distance_um'],vm['profile_response'],cm['profile_response'])):
                    profiles.append(dict(layout_id=lay['layout_id'],source_segment=source,focal_strength_multiplier=strength,
                        path_bin_index=j,mean_path_um=float(x),source_normalized_voltage=float(y),source_normalized_calcium=float(z)))
                # Direct validations chosen by input coordinates, not observed outcomes.
                key=(lay['layout_id'],source,strength)
                requested=[('control',CANONICAL,1.),('f25_uniform',CANONICAL,1.),
                    ('f20_uniform',int(off[np.argmin(soma_d[off])]),.5),('f30_uniform',int(off[np.argmax(soma_d[off])]),2.),
                    ('f25_random00',CANONICAL,1.),('f25_distal',CANONICAL,2.)]
                if key in requested:validation_cases.append((lay,source,strength,dv.copy(),row))
        np.savez_compressed(OUT/'focal_cache'/f'{lay["layout_id"]}.npz',sources=np.array(laysource),source_Rin_MOhm=np.array(layrin),normalized_voltage_transfer=np.stack(laymaps))
        print(f'focal layout {li+1}/70 {lay["layout_id"]}: {len(sources)} sources; {time.time()-started:.1f}s',flush=True)
    assert len(allsites)==900 and len(canonical)==210
    tsv(OUT/'focal_all_sites.tsv',allsites);tsv(OUT/'focal_layout_controls.tsv',canonical);tsv(OUT/'focal_binned_profiles.tsv',profiles)
    prior=json.loads((PRIOR/'spatial_sensitivity.json').read_text())
    for old in prior:
        f=old['retained_effective_inhibition_fraction'];lid='control' if f==1 else f'f{round(f*100)}_uniform'
        row=next(r for r in canonical if r['layout_id']==lid and r['focal_strength_multiplier']==1.)
        for field,oldv in [('voltage_1e_um',old['voltage']['path_1e_crossing_um']),('calcium_1e_um',old['calcium']['path_1e_crossing_um']),('source_dv_mV',old['source_dv_mV'])]:
            error=abs(row[field]-oldv);checks.append(dict(check=f'archived_{lid}_{field}',passed=bool(error<1e-8),absolute_error=error,tolerance=1e-8))
    clear_shunts(model,prov['parameters']['g_pas_S_cm2'])
    p=dict(prov['parameters'],glycine_fraction=1.)
    for k,(lay,source,strength,expected,row) in enumerate(validation_cases):
        tick=time.time();t,v=simulate_v3(model,[source],[SOURCE_G*strength],lay['weights'],np.arange(len(length)),p,0.,1.,dt=.1,sample_dt=10.,tstop=3050.)
        observed=v[:,-1]-model['rest_v'];err=float(np.max(np.abs(observed-expected)))
        check=dict(check=f'direct_3050ms_focal_{lay["layout_id"]}_{source}_{strength}',passed=bool(err<1e-7),max_absolute_dv_error_mV=err,tolerance_mV=1e-7,
            maximum_baseline_error_mV=float(np.max(np.abs(v[:,t<49]-model['rest_v']))),elapsed_seconds=time.time()-tick)
        checks.append(check);assert check['maximum_baseline_error_mV']<1e-8
        print(f'direct focal {k+1}/{len(validation_cases)}: error {err:.3g}mV; {time.time()-tick:.1f}s',flush=True)
        del v
    checks.append(dict(check='all_source_geometry_count',passed=bool(len(allsites)==900 and len(canonical)==210),all_source_cases=len(allsites),layout_cases=len(canonical),n_unique_source_compartments=int(len(np.unique(off)))))
    checks.append(dict(check='source_files_unchanged',passed=protected()==json.loads((OUT/'source_hashes_before.json').read_text())))
    json_write(OUT/'focal_validation.json',checks)
    json_write(OUT/'focal_run.json',dict(status='complete' if all(c['passed'] for c in checks) else 'validation_failed',elapsed_seconds=time.time()-started,
        all_source_cases=len(allsites),layout_cases=len(canonical),source_count=n,unique_source_count=int(len(np.unique(off))),layout_count=len(layouts),
        direct_full_duration_validations=len(validation_cases),python=sys.executable,numpy=np.__version__,scipy=scipy.__version__,neuron=neuron.__version__,
        calcium=prov['calcium'],runner_sha256=sha(__file__),crossing_definition='First crossing of cable-weighted5um path bins, linearly interpolated; none means no crossing observed within arbor.',
        baseline='All inhibitory/leak reversals are -45mV; pure shunt preserves uniform resting voltage.',
        source_normalization_threshold=5e-5))
    assert all(c['passed'] for c in checks)
    print(f'FOCAL_COMPLETE {time.time()-started:.1f}s',flush=True)

def spot_stage():
    """Run all770layout/contrast cases with independently checked passive tail."""
    sys.path.insert(0,str(SCRIPT_DIR.parent/'electrical'))
    from passive_tail import simulate_spot
    helper=SCRIPT_DIR.parent/'electrical/passive_tail.py';helper_hash=sha(helper)
    started=time.time();model,prov,ew,iw=make_model();layouts,soma_d=make_layouts(model,iw)
    assert json.loads((OUT/'source_hashes_before.json').read_text())==protected()
    rows=[];localrows=[];diagnostics=[];checks=[];validation_cache={}
    requested=[('control',70),('f25_uniform',30),('f20_random00',5),('f25_random07',30),
               ('f30_random19',70),('f25_proximal',5),('f25_distal',70),('f20_proximal',70)]
    (OUT/'spot_voltage_cache').mkdir(exist_ok=True)
    for li,lay in enumerate(layouts):
        archive={};layout_rows=[]
        for c in CONTRASTS:
            tick=time.time();res=simulate_spot(model,ew,lay['weights'],contrast=c/100,tolerance=1e-7,return_trace=False)
            row=dict(layout_id=lay['layout_id'],family=lay['family'],retained_fraction=lay['fraction'],seed=lay['seed'],
                contrast_pct=c,n_model_sites=len(model['off_nodes']),mean_peak_dv_mV=res['mean_peak_dv'],
                mean_peak_ca_au_ms=res['mean_peak_ca'],mean_ca_auc_au_ms2=res['mean_ca_auc'],
                min_local_peak_dv_mV=res['min_peak_dv'],max_local_peak_dv_mV=res['max_peak_dv'],
                tail_start_ms=res['tail_start_ms'],tail_voltage_error_bound_mV=res['tail_voltage_error_bound_mV'],
                tail_auc_error_bound=res['tail_auc_error_bound'],model_auc_shape=None)
            rows.append(row);layout_rows.append(row)
            for j,(a,v,q) in enumerate(zip(res['local_auc'],res['local_peak_dv'],res['local_peak_ca'])):
                localrows.append(dict(layout_id=lay['layout_id'],contrast_pct=c,source_index=j,segment_index=int(model['off_nodes'][j]),
                    peak_dv_mV=float(v),peak_calcium_au_ms=float(q),calcium_auc_au_ms2=float(a)))
            diagnostics.append(dict(layout_id=lay['layout_id'],contrast_pct=c,elapsed_seconds=time.time()-tick,
                integration_end_ms=res['end_integration_ms'],tail_start_ms=res['tail_start_ms'],nblocks=res['nblocks'],
                tail_voltage_error_bound_mV=res['tail_voltage_error_bound_mV'],tail_auc_error_bound=res['tail_auc_error_bound'],
                bound_scale=res['bound_scale'],relative_tolerance=res['requested_relative_tolerance'],
                relative_auc_bound=res['tail_auc_error_bound']/max(1,res['mean_ca_auc']),
                tau_max_ms=res['tau_max_ms'],contraction_per_block=res['contraction_per_block'],
                baseline_mV=res['baseline_mV'],min_voltage_mV=res['min_voltage'],max_voltage_mV=res['max_voltage']))
            assert res['min_voltage']>=-45-1e-8 and res['max_voltage']<=1e-8
            assert res['tail_auc_error_bound']<=1e-7*res['bound_scale']*(1+1e-8)
            if c==0:assert abs(res['mean_ca_auc'])<1e-8 and abs(res['mean_peak_dv'])<1e-8
            archive[f't_ms_c{c}']=res['recorded_t_ms'];archive[f'voltage_mV_c{c}']=res['recorded_v_mV']
            archive[f'local_auc_c{c}']=res['local_auc'];archive[f'local_peak_ca_c{c}']=res['local_peak_ca'];archive[f'local_peak_dv_c{c}']=res['local_peak_dv']
            if (lay['layout_id'],c) in requested:validation_cache[lay['layout_id'],c]=(res,lay)
        span=next(r['mean_ca_auc_au_ms2'] for r in layout_rows if r['contrast_pct']==70)-next(r['mean_ca_auc_au_ms2'] for r in layout_rows if r['contrast_pct']==15)
        lo=next(r['mean_ca_auc_au_ms2'] for r in layout_rows if r['contrast_pct']==15)
        for r in layout_rows:r['model_auc_shape']=(r['mean_ca_auc_au_ms2']-lo)/span if span>THRESHOLD else None
        np.savez_compressed(OUT/'spot_voltage_cache'/f'{lay["layout_id"]}.npz',**archive)
        print(f'spot layout {li+1}/70 {lay["layout_id"]}: elapsed{time.time()-started:.1f}s',flush=True)
        tsv(OUT/'spot_responses.tsv',rows)
    assert len(rows)==770 and len(validation_cache)==8
    tsv(OUT/'spot_local_responses.tsv',localrows);tsv(OUT/'spot_tail_diagnostics.tsv',diagnostics)
    # Validate original3052-ms protocol at prespecified parameter coordinates.
    p=dict(prov['parameters'],glycine_fraction=1.)
    for lid,c in requested:
        res,lay=validation_cache[lid,c];tick=time.time()
        t,v=simulate_v3(model,model['off_nodes'],ew,lay['weights'],model['off_nodes'],p,0.,c/100,dt=.1,sample_dt=2.,tstop=3052.)
        dv,q,a=readout(v,t,prov['calcium'],model['rest_v'])
        keep=(t>=50-1e-6)&(t<=3050+1e-6)
        fullv=v[:,keep];fastv=np.empty_like(fullv);prefix=res['recorded_v_mV'];fastv[:,:prefix.shape[1]]=prefix;fastv[:,prefix.shape[1]:]=prefix[:,-1,None]
        trace_err=float(np.max(np.abs(fullv-fastv)));auc_err=float(np.max(np.abs(a-res['local_auc'])))
        mean_rel=float(abs(a.mean()-res['mean_ca_auc'])/max(1,abs(a.mean())))
        check=dict(check=f'full_original_spot_{lid}_{c}',passed=bool(auc_err<=res['tail_auc_error_bound']+1e-6 and trace_err<=res['tail_voltage_error_bound_mV']+1e-8 and mean_rel<1e-7),
            max_local_auc_error=auc_err,mean_auc_relative_error=mean_rel,max_voltage_trace_error_mV=trace_err,
            reported_auc_error_bound=res['tail_auc_error_bound'],reported_voltage_error_bound_mV=res['tail_voltage_error_bound_mV'],
            maximum_local_peak_ca_error=float(np.max(np.abs(q-res['local_peak_ca']))),maximum_local_peak_dv_error=float(np.max(np.abs(dv-res['local_peak_dv']))),
            full_duration_s=time.time()-tick,numerical_allowance_auc=1e-6,numerical_allowance_mV=1e-8)
        checks.append(check);print(f'full spot {lid} C{c}: relativeAUCerror{mean_rel:.3g}; elapsed{time.time()-tick:.1f}s',flush=True)
    # Historical v9 means derive from the original full-duration runs.
    old=json.loads((PRIOR/'nominal25/figure_data.json').read_text())['contrast']['ctrl']
    oldrows=loadtsv(PRIOR/'contrast_sensitivity.tsv')
    for row in rows:
        if row['family'] not in ('control','uniform'):continue
        c=row['contrast_pct']
        if row['family']=='control':
            i=old['contrast_pct'].index(c);ref={name:old[name][i] for name in ['mean_peak_dv_mV','mean_ca_auc_au_ms2']}
        else:ref=next(r for r in oldrows if float(r['retained_effective_inhibition_fraction'])==row['retained_fraction'] and float(r['contrast_pct'])==c)
        for name in ('mean_peak_dv_mV','mean_ca_auc_au_ms2'):
            err=abs(row[name]-float(ref[name]));rel=err/max(1,abs(float(ref[name])))
            checks.append(dict(check=f'historical_{row["layout_id"]}_{c}_{name}',passed=bool(rel<1e-7),absolute_error=err,scaled_relative_error=rel,tolerance=1e-7))
    checks.append(dict(check='passive_tail_helper_unchanged',passed=sha(helper)==helper_hash,sha256=helper_hash))
    checks.append(dict(check='source_files_unchanged',passed=protected()==json.loads((OUT/'source_hashes_before.json').read_text())))
    json_write(OUT/'spot_validation.json',checks)
    json_write(OUT/'spot_run.json',dict(status='complete' if all(c['passed'] for c in checks) else 'validation_failed',elapsed_seconds=time.time()-started,
        case_count=len(rows),layout_count=len(layouts),n_local_records=len(localrows),calcium=prov['calcium'],runner_sha256=sha(__file__),
        helper_path=str(helper),helper_sha256=helper_hash,full_original_protocol_checks=8,
        max_tail_auc_relative_bound=max(r['relative_auc_bound'] for r in diagnostics),max_tail_start_ms=max(r['tail_start_ms'] for r in diagnostics),
        max_direct_mean_auc_relative_error=max(c['mean_auc_relative_error'] for c in checks if 'mean_auc_relative_error' in c),
        archive_definition='Per-layout saved exactNEURONprefix at2ms from50ms; constant last-voltage extension through3050ms; baselinebefore50ms identically-45mV; same original localcalciumrecurrence/trapezoid; allboundsstored.',
        interpretation='Numerically bounded latevoltagetail acceleration; no steady-state substitution for stimulus onset or calcium AUC.'))
    assert all(c['passed'] for c in checks)
    print(f'SPOT_COMPLETE {time.time()-started:.1f}s',flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--outdir',type=Path,default=SCRIPT_DIR/'reproduction');ap.add_argument('--stage',choices=['focal','spot'],required=True);args=ap.parse_args()
    OUT=args.outdir.resolve();OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'focal_cache').mkdir(exist_ok=True)
    if args.stage=='focal':focal_stage()
    elif args.stage=='spot':spot_stage()
