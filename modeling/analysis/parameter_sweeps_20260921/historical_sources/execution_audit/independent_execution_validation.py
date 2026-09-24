"""Independent, read-only validation of completed electrical/spatial outputs.

Uses exported arrays, never imports the production sweep runners.
"""
from pathlib import Path
import os
for variable in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[variable] = '1'
import csv
import hashlib
import itertools
import json
import numpy as np
from scipy.special import expit
from scipy.signal import lfilter

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
WORK = BASE.parents[1]
CHECKS = []


def check(name, passed, **details):
    CHECKS.append(dict(check=name, passed=bool(passed), **details))


def rows(name):
    with Path(name).open() as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def near(a, b, rtol=1e-10, atol=1e-9):
    return np.allclose(a, b, rtol=rtol, atol=atol, equal_nan=True)


def optional_number(v):
    return float(v) if v not in ('', None) else np.nan


def crossing(distance, response, length):
    # Independently form cable-weighted bins (no production helper imported).
    bins = np.floor(distance/5).astype(int)
    xs, ys = [0.], [1.]
    for k in np.unique(bins):
        select = bins == k
        xs.append(float(np.average(distance[select], weights=length[select])))
        ys.append(float(np.average(response[select], weights=length[select])))
    for i in range(1, len(xs)):
        if ys[i-1] > 1/np.e >= ys[i]:
            return xs[i-1] + (xs[i]-xs[i-1])*(ys[i-1]-1/np.e)/(ys[i-1]-ys[i])
    return np.nan


def spatial():
    folder = BASE/'spatial'
    layouts = rows(folder/'layout_summary.tsv')
    arrays = np.load(folder/'layout_weights.npz')
    geom = np.load(folder/'source_geometry.npz')
    control = arrays['control_weights_uS']
    index = {str(x): i for i, x in enumerate(arrays['layout_ids'])}
    nodes = arrays['inh_nodes']
    d = geom['soma_path_distances_um'][nodes]
    ids = np.arange(len(nodes))
    check('spatial_design_70_layouts_242_contacts', len(layouts)==70 and arrays['weights_uS'].shape==(70,242))
    check('spatial_soma_distance_origin', geom['soma_path_distances_um'][int(geom['soma_segment_index'])] == 0)
    max_budget_error = 0.
    max_reconstruction_error = 0.
    for row in layouts:
        fraction = float(row['retained_fraction'])
        actual = arrays['weights_uS'][index[row['layout_id']]]
        target = fraction*sum(control)
        max_budget_error = max(max_budget_error, abs(sum(actual)-target))
        if row['family'] in ('control', 'uniform'):
            expected = control*fraction
        else:
            if row['family'] == 'random':
                order = np.random.Generator(np.random.PCG64(int(row['seed']))).permutation(len(control))
            else:
                sign = 1 if row['family'] == 'proximal' else -1
                order = np.array(sorted(ids, key=lambda j: (sign*d[j],nodes[j],j)))
            # The cumulative-budget formula is independent of the greedy loop.
            previous = np.r_[0, np.cumsum(control[order])[:-1]]
            expected = np.zeros(len(control))
            expected[order] = np.minimum(control[order], np.maximum(target-previous, 0.))
            partial = (actual > 1e-14) & (actual < control-1e-14)
            check('at_most_one_partial_'+row['layout_id'], partial.sum() <= 1)
        max_reconstruction_error = max(max_reconstruction_error, float(np.max(abs(expected-actual))))
        check('layout_weights_'+row['layout_id'], near(expected,actual,rtol=0,atol=1e-13) and np.all(actual>=0) and np.all(actual<=control))
    check('matched_total_inhibitory_budgets', max_budget_error<1e-13, maximum_error_uS=max_budget_error)
    check('independent_layout_reconstruction', max_reconstruction_error<1e-13, maximum_error_uS=max_reconstruction_error)
    allsites = rows(folder/'focal_all_sites.tsv')
    canonical = rows(folder/'focal_layout_controls.tsv')
    source_nodes = [int(x) for x in geom['off_nodes']]
    expected = set(itertools.product(['control','f20_uniform','f25_uniform','f30_uniform'], source_nodes, [.5,1.,2.]))
    keys = [(r['layout_id'],int(r['source_segment']),float(r['focal_strength_multiplier'])) for r in allsites]
    check('900_unique_focal_site_conditions', len(keys)==900 and len(set(keys))==900 and set(keys)==expected)
    expected = set(itertools.product(index,[4705],[.5,1.,2.]))
    keys2 = [(r['layout_id'],int(r['source_segment']),float(r['focal_strength_multiplier'])) for r in canonical]
    check('210_unique_focal_layout_conditions', len(keys2)==210 and len(set(keys2))==210 and set(keys2)==expected)
    caches = {}
    for lid in index:
        with np.load(folder/'focal_cache'/f'{lid}.npz') as cached:
            caches[lid] = {k: cached[k] for k in cached.files}
    err = {k:0. for k in ['source_dv_mV','source_calcium_au_ms','voltage_1e_um','calcium_1e_um','voltage_cable_integral_um','calcium_cable_integral_um','calcium_absolute_cable_integral_au_ms_um']}
    match_missing = True
    for row in allsites+canonical:
        source = int(row['source_segment'])
        cache = caches[row['layout_id']]
        j = list(cache['sources']).index(source)
        transfer = cache['normalized_voltage_transfer'][j]
        rin = float(cache['source_Rin_MOhm'][j])
        g = float(row['focal_conductance_nS'])/1000
        dv0 = 45*g*rin/(1+g*rin)
        q = 50*(expit((-45+dv0*transfer+37.47892660844856)/1.3984142481316566)-expit((-45+37.47892660844856)/1.3984142481316566))
        dist = geom['source_path_distances_um'][source_nodes.index(source)]
        length = geom['cable_lengths']
        qn = q/q[source]
        calc = dict(source_dv_mV=dv0, source_calcium_au_ms=q[source], voltage_1e_um=crossing(dist,transfer,length), calcium_1e_um=crossing(dist,qn,length), voltage_cable_integral_um=np.dot(length,transfer), calcium_cable_integral_um=np.dot(length,qn), calcium_absolute_cable_integral_au_ms_um=np.dot(length,q))
        for k,value in calc.items():
            observed = optional_number(row[k])
            if np.isnan(value) or np.isnan(observed):
                match_missing = match_missing and np.isnan(value) and np.isnan(observed)
            else:
                err[k] = max(err[k],abs(value-observed))
    for field,error in err.items():
        check('recomputed_focal_'+field, error<1e-8, maximum_absolute_error=error)
    check('independent_focal_censor_flags', match_missing)
    for name in ('focal_validation.json','spot_validation.json'):
        p = folder/name
        if p.exists():
            data=json.loads(p.read_text())
            if isinstance(data,list):
                check('production_'+name,all(x['passed'] for x in data),checks=len(data))
    spot_path=folder/'spot_responses.tsv'
    if spot_path.exists():
        spot=rows(spot_path)
        keys=[(r['layout_id'],int(r['contrast_pct'])) for r in spot]
        expected=set(itertools.product(index,[0,5,10,15,20,25,30,40,50,60,70]))
        check('spatial_770_unique_spot_conditions',len(keys)==770 and len(set(keys))==770 and set(keys)==expected)
        check('spatial_zero_contrast',all(float(r['mean_ca_auc_au_ms2'])==0 and float(r['mean_peak_dv_mV'])==0 for r in spot if int(r['contrast_pct'])==0))
        diag=rows(folder/'spot_tail_diagnostics.tsv')
        bound=max(float(r['tail_auc_error_bound'])/float(r['bound_scale']) for r in diag)
        check('spatial_all_tail_bounds',bound<=1e-7,maximum_relative_bound=bound)
        check('spatial_all_baselines_and_reversals',all(float(r['baseline_mV'])==-45 and float(r['min_voltage_mV'])>=-45-1e-9 and float(r['max_voltage_mV'])<=1e-9 for r in diag))


def preservation():
    saved=json.loads((HERE/'input_hashes.json').read_text())['files']
    changed=[]
    for name,digest in saved.items():
        p=WORK/name
        if not p.exists() or hashlib.sha256(p.read_bytes()).hexdigest()!=digest:
            changed.append(name)
    check('protected_original_inputs_unchanged',not changed,files=len(saved),changed=changed)


def electrical():
    folder=BASE/'electrical'
    path=folder/'electrical_sweep.npz'
    if not path.exists():
        return
    with np.load(path) as archive:
        a={k:archive[k] for k in archive.files}
    m=2.**(np.arange(-4,5)/2)
    check('electrical_prescribed_axes',near(a['inhibition_multipliers'],m) and near(a['leak_multipliers'],m) and near(a['excitation_multipliers'],[.5,1,2]) and near(a['retained_fractions'],[1,.2,.25,.3]) and near(a['contrast_pct'],[0,5,10,15,20,25,30,40,50,60,70]))
    design=rows(folder/'spot_responses.tsv')
    focal=rows(folder/'focal_responses.tsv')
    keyfields=['inhibition_multiplier','leak_multiplier','excitation_multiplier','retained_fraction','contrast_pct']
    keys=[tuple(float(r[k]) for k in keyfields) for r in design]
    check('electrical_10692_unique_spot_conditions',len(keys)==10692 and len(set(keys))==10692 and a['mean_auc'].shape==(9,9,3,4,11))
    keys=[tuple(float(r[k]) for k in keyfields[:-1]) for r in focal]
    check('electrical_972_unique_focal_conditions',len(keys)==972 and len(set(keys))==972 and a['focal_source_ca'].shape==(9,9,3,4))
    for k in ('mean_auc','local_auc','mean_peak_dv_mV','local_peak_dv_mV','mean_peak_ca','tail_auc_bound','tail_voltage_bound_mV','tail_bound_scale'):
        check('electrical_finite_nonnegative_'+k,np.isfinite(a[k]).all() and a[k].min()>=0)
    check('electrical_zero_contrast',np.max(abs(a['local_auc'][...,0,:]))<1e-10 and np.max(abs(a['local_peak_dv_mV'][...,0,:]))<1e-10)
    check('electrical_baseline_and_reversal_bounds',near(a['baseline_mV'],-45) and a['min_voltage_mV'].min()>=-45-1e-9 and a['max_voltage_mV'].max()<=1e-9)
    check('electrical_all_tail_bounds',np.all(a['tail_auc_bound']<=1e-7*a['tail_bound_scale']+1e-10),maximum_relative_bound=float((a['tail_auc_bound']/a['tail_bound_scale']).max()))
    for aggregate,local in [('mean_auc','local_auc'),('mean_peak_dv_mV','local_peak_dv_mV'),('mean_peak_ca','local_peak_ca')]:
        check('electrical_local_then_mean_'+aggregate,near(a[aggregate],a[local].mean(axis=-1)))
    check('electrical_contrast_monotonicity',np.min(np.diff(a['mean_auc'],axis=-1))>=-1e-8)
    check('electrical_reduced_inhibition_increases_absolute_auc',np.min(a['mean_auc'][:,:,:,1:,1:]-a['mean_auc'][:,:,:,:1,1:])>=-1e-8)
    original=np.load(BASE/'results/calcium_sweep.npz')['mean_auc'][-1]
    observed=a['mean_auc'][4,4,1]
    error=float(np.max(abs(observed-original)/np.maximum(1,abs(original))))
    check('electrical_nominal_matches_calcium_archive',error<5e-7,maximum_scaled_error=error,tolerance=5e-7)
    # Duplicate products must agree for every electrical state and fraction.
    products={}
    for ei,e in enumerate(a['excitation_multipliers']):
        for ci,c in enumerate(a['contrast_pct']):
            product=round(float(e*c/100),12)
            if product in products:
                ej,cj=products[product]
                check(f'electrical_identical_drive_E{ei}_C{ci}',np.array_equal(a['local_auc'][:,:,ei,:,ci],a['local_auc'][:,:,ej,:,cj]))
            else: products[product]=(ei,ci)
    # Independently apply a linear IIR filter to recorded local activation.
    # This avoids importing the production calcium implementation.
    for ii,ll,ff,ei,ci in [(0,0,1,0,1),(0,0,1,2,10),(4,4,0,1,10),(4,4,2,1,6),(8,8,0,0,1),(8,0,3,2,3)]:
        trace=np.load(folder/'transient_traces'/f'i{ii}_l{ll}_f{ff}.npz')
        product=round(float(a['excitation_multipliers'][ei]*a['contrast_pct'][ci]/100),12)
        k=str(product).replace('.','p')
        ts=trace['t_'+k];prefix=trace['v_'+k]
        check(f'electrical_recording_alignment_{ii}_{ll}_{ff}_{ei}_{ci}',near(ts,np.arange(50,50+2*len(ts),2)))
        v=np.pad(prefix,((0,0),(0,1501-prefix.shape[1])),mode='edge')
        influx=expit((v+37.47892660844856)/1.3984142481316566)-expit((-45+37.47892660844856)/1.3984142481316566)
        decay=np.exp(-2/50)
        ca=lfilter([0.,50*(1-decay)],[1.,-decay],influx,axis=1)
        auc=np.trapezoid(ca,dx=2,axis=1)
        ref=a['local_auc'][ii,ll,ei,ff,ci]
        rel=float(np.max(abs(auc-ref)/np.maximum(1,abs(ref))))
        check(f'electrical_independent_clearance_{ii}_{ll}_{ff}_{ei}_{ci}',rel<1e-12,maximum_scaled_error=rel)
    # Independently executed spatial-control nominal overlaps must agree.
    layoutrows=rows(BASE/'spatial/focal_layout_controls.tsv')
    fmap={'control':0,'f20_uniform':1,'f25_uniform':2,'f30_uniform':3}
    for r in layoutrows:
        if r['layout_id'] not in fmap: continue
        fi=fmap[r['layout_id']];ei=[.5,1.,2.].index(float(r['focal_strength_multiplier']))
        for field,key in [('calcium_1e_um','focal_calcium_1e_um'),('voltage_1e_um','focal_voltage_1e_um'),('source_dv_mV','focal_source_dv_mV')]:
            check(f'electrical_spatial_overlap_{fi}_{ei}_{field}',near(optional_number(r[field]),a[key][4,4,ei,fi]))
    p=BASE/'spatial/spot_responses.tsv'
    if p.exists():
        spot=rows(p)
        err=0.
        for r in spot:
            if r['layout_id'] not in fmap: continue
            fi=fmap[r['layout_id']];ci=list(a['contrast_pct']).index(int(r['contrast_pct']))
            err=max(err,abs(float(r['mean_ca_auc_au_ms2'])-a['mean_auc'][4,4,1,fi,ci]))
        check('electrical_spatial_spot_overlap',err<1e-8,maximum_absolute_error=err)
    for name in ('tail_validation.json','focal_validation.json'):
        p=folder/name
        if p.exists():
            data=json.loads(p.read_text())
            check('production_electrical_'+name,all(x['passed'] for x in data),checks=len(data))


def supplement():
    folder=BASE/'calcium_supplement'
    manifest=json.loads((folder/'manifest.json').read_text())
    mismatches=[]
    for name,digest in manifest['artifact_hashes'].items():
        p=folder/name
        if not p.exists() or hashlib.sha256(p.read_bytes()).hexdigest()!=digest:
            mismatches.append(name)
    check('calcium_supplement_artifact_hashes',not mismatches,files=len(manifest['artifact_hashes']),mismatches=mismatches)
    source=hashlib.sha256((BASE/'results/calcium_sweep.npz').read_bytes()).hexdigest()
    check('calcium_supplement_frozen_numerical_source',source==manifest['source_archive_sha256'])


if __name__ == '__main__':
    spatial()
    electrical()
    supplement()
    preservation()
    output=dict(status='passed' if all(c['passed'] for c in CHECKS) else 'failed', checks_passed=sum(c['passed'] for c in CHECKS), checks_total=len(CHECKS), checks=CHECKS)
    (HERE/'independent_execution_validation.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({k:v for k,v in output.items() if k!='checks'}))
    for item in CHECKS:
        if not item['passed']: print(item)
    assert output['status']=='passed'
