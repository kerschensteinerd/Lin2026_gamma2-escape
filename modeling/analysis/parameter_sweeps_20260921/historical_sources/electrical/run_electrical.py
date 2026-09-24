#!/usr/bin/env python3
"""Full protocol-B matrix with exact duplicate-drive caching and bounded tails."""
from __future__ import annotations
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[name]='1'
from pathlib import Path
import sys,json,hashlib,time,argparse,multiprocessing as mp,csv
import numpy as np
from scipy.special import expit
from neuron import h
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];MODEL=ROOT/'vg3_cid3_rf_simulation-main'
sys.path.insert(0,str(MODEL/'scripts/fig3'))
from run_fig3_pilot import build,load_fits,load_named_config,exc_eff_weight,inh_eff_weight,gaussian_spot_weights,path_distance
from passive_tail import simulate_spot,VHALF,SLOPE,TAU

M=2.**(np.arange(-4,5)/2);E=np.array([.5,1.,2.]);F=np.array([1.,.2,.25,.3]);C=np.array([0,5,10,15,20,25,30,40,50,60,70])
SOURCE=4705;FOCAL=.012785038085113628
SPOT={'mean_auc':'mean_ca_auc','mean_peak_dv_mV':'mean_peak_dv','mean_peak_ca':'mean_peak_ca','max_peak_dv_mV':'max_peak_dv','min_peak_dv_mV':'min_peak_dv','min_voltage_mV':'min_voltage','max_voltage_mV':'max_voltage','baseline_mV':'baseline_mV','tail_auc_bound':'tail_auc_error_bound','tail_voltage_bound_mV':'tail_voltage_error_bound_mV','tail_start_ms':'tail_start_ms','tail_bound_scale':'bound_scale','max_activation':'max_activation'}
LOCAL={'local_auc':'local_auc','local_peak_dv_mV':'local_peak_dv','local_peak_ca':'local_peak_ca'}
FOCAL_NAMES=['source_rin_mohm','source_dv_mV','source_ca','voltage_1e_um','calcium_1e_um','voltage_weighted_cable_um','calcium_weighted_cable_um','absolute_calcium_cable','max_voltage_mV']
STATE={}

def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def write_json(path,value):
    def clean(x):
        if isinstance(x,dict):return {str(k):clean(v) for k,v in x.items()}
        if isinstance(x,(list,tuple)):return [clean(v) for v in x]
        if isinstance(x,np.ndarray):return clean(x.tolist())
        if isinstance(x,np.integer):return int(x)
        if isinstance(x,np.bool_):return bool(x)
        if isinstance(x,(float,np.floating)):return float(x) if np.isfinite(x) else None
        return x
    path.write_text(json.dumps(clean(value),indent=2,allow_nan=False)+'\n')

def initialize():
    m=build(1e-4);cfg,_=load_named_config('sep_g22');fits=load_fits();xy=m['node_coords'];off=m['off_nodes'];inh=m['inh_nodes'];segs=m['segments']
    ew=exc_eff_weight(fits,'OFF',100)*cfg['G_exc']*gaussian_spot_weights(xy[off,:2],m['soma_xy'],100,25)
    iw=inh_eff_weight(fits,'OFF',100)*cfg['G_inh']*3*gaussian_spot_weights(xy[inh,:2],m['soma_xy'],100,25)
    distance=path_distance(m,SOURCE);length=np.array([s.sec.L/s.sec.nseg for s in segs]);areas=np.array([h.area(s.x,sec=s.sec) for s in segs])
    bins=[(distance>=lo)&(distance<lo+5) for lo in np.arange(0,distance.max()+5,5)];bins=[b for b in bins if b.any()]
    x=np.array([0.]+[np.average(distance[b],weights=length[b]) for b in bins])
    STATE.update(model=m,ew=ew,iw=iw,distance=distance,length=length,areas=areas,bins=bins,x=x)

def profile(values):return np.r_[1.,[np.average(values[b],weights=STATE['length'][b]) for b in STATE['bins']]]
def crossing(y):
    hits=np.flatnonzero((y[:-1]>1/np.e)&(y[1:]<=1/np.e))
    if not len(hits):return np.nan
    i=hits[0];x=STATE['x'];return x[i]+(1/np.e-y[i])*(x[i+1]-x[i])/(y[i+1]-y[i])

def focal(leak,iw):
    m=STATE['model'];segs=m['segments'];pernode=np.bincount(m['inh_nodes'],weights=iw,minlength=len(segs))
    for i,s in enumerate(segs):s.g_pas=leak+100*pernode[i]/STATE['areas'][i]
    h.finitialize(-45.);imp=h.Impedance();source=segs[SOURCE];imp.loc(source.x,sec=source.sec);imp.compute(0.)
    rin=float(imp.input(source.x,sec=source.sec));transfer=np.array([imp.transfer(s.x,sec=s.sec) for s in segs])/rin
    vp=profile(transfer);data={k:np.empty(3) for k in FOCAL_NAMES};cp=[]
    for j,e in enumerate(E):
        g=FOCAL*e;source_dv=45*g*rin/(1+g*rin);dv=transfer*source_dv
        ca=TAU*(expit((-45+dv-VHALF)/SLOPE)-expit((-45-VHALF)/SLOPE));cn=ca/ca[SOURCE];p=profile(cn);cp.append(p)
        values=[rin,source_dv,ca[SOURCE],crossing(vp),crossing(p),np.sum(transfer*STATE['length']),np.sum(cn*STATE['length']),np.sum(ca*STATE['length']),float((-45+dv).max())]
        for k,v in zip(FOCAL_NAMES,values):data[k][j]=v
    for s in segs:s.g_pas=leak
    data.update(voltage_profile=vp,calcium_profiles=np.array(cp),normalized_voltage_map=transfer)
    return data

def run_state(task):
    ii,ll,ff=task;path=HERE/'states'/f'i{ii}_l{ll}_f{ff}.npz'
    if path.is_file():return dict(task=task,cached=True)
    start=time.time();leak=M[ll]*1e-4;m=STATE['model'];m['g_pas']=leak;m['p']['g_pas']=leak
    for seg in m['segments']:seg.g_pas=leak;seg.e_pas=-45.
    iw=STATE['iw']*M[ii]*F[ff]
    products=sorted({round(float(e*c/100),12) for e in E for c in C})
    results={};traces={}
    for product in products:
        r=simulate_spot(m,STATE['ew'],iw,contrast=product)
        if not(r['min_voltage']>=-45-1e-9 and r['max_voltage']<=1e-9):raise AssertionError('Reversal bounds')
        if not r['tail_auc_error_bound']<=1e-7*r['bound_scale']+1e-10:raise AssertionError('Tail bound')
        results[product]=r
        key=str(product).replace('.','p');traces[f't_{key}']=r['recorded_t_ms'];traces[f'v_{key}']=r['recorded_v_mV']
    out={}
    for save,key in SPOT.items():out[save]=np.array([[results[round(float(e*c/100),12)][key] for c in C] for e in E])
    for save,key in LOCAL.items():out[save]=np.array([[results[round(float(e*c/100),12)][key] for c in C] for e in E])
    foc=focal(leak,iw)
    out.update({f'focal_{k}':v for k,v in foc.items()});out.update(mi=M[ii],ml=M[ll],retained_fraction=F[ff],indices=np.array(task),unique_products=np.array(products),elapsed_seconds=time.time()-start)
    np.savez_compressed(HERE/'transient_traces'/path.name,**traces)
    temporary=path.with_suffix('.writing.npz');np.savez_compressed(temporary,**out);temporary.replace(path)
    return dict(task=task,elapsed_s=time.time()-start,unique=len(products),max_bound=float(out['tail_auc_bound'].max()),max_tail_ms=float(out['tail_start_ms'].max()))

def aggregate():
    shape=(9,9,3,4,11);out={key:np.empty(shape) for key in SPOT};out.update({key:np.empty((*shape,75)) for key in LOCAL})
    out.update({f'focal_{key}':np.empty((9,9,3,4)) for key in FOCAL_NAMES})
    nbin=len(STATE['x']);out['focal_voltage_profile']=np.empty((9,9,4,nbin));out['focal_calcium_profiles']=np.empty((9,9,3,4,nbin));state_seconds=[]
    for ii in range(9):
        for ll in range(9):
            for ff in range(4):
                with np.load(HERE/'states'/f'i{ii}_l{ll}_f{ff}.npz') as a:
                    for key in [*SPOT,*LOCAL]:out[key][ii,ll,:,ff]=a[key]
                    for key in FOCAL_NAMES:out[f'focal_{key}'][ii,ll,:,ff]=a[f'focal_{key}']
                    out['focal_voltage_profile'][ii,ll,ff]=a['focal_voltage_profile'];out['focal_calcium_profiles'][ii,ll,:,ff]=a['focal_calcium_profiles'];state_seconds.append(float(a['elapsed_seconds']))
    out.update(inhibition_multipliers=M,leak_multipliers=M,excitation_multipliers=E,retained_fractions=F,contrast_pct=C,path_bin_um=STATE['x'],path_distance_um=STATE['distance'],cable_lengths_um=STATE['length'],off_nodes=STATE['model']['off_nodes'])
    np.savez_compressed(HERE/'electrical_sweep.npz',**out)
    with (HERE/'spot_responses.tsv').open('w') as file:
        fields=['inhibition_multiplier','leak_multiplier','excitation_multiplier','retained_fraction','contrast_pct',*SPOT];writer=csv.DictWriter(file,fieldnames=fields,delimiter='\t');writer.writeheader()
        for ii,ll,ee,ff,cc in np.ndindex(shape):writer.writerow(dict(inhibition_multiplier=M[ii],leak_multiplier=M[ll],excitation_multiplier=E[ee],retained_fraction=F[ff],contrast_pct=C[cc],**{key:out[key][ii,ll,ee,ff,cc] for key in SPOT}))
    with (HERE/'focal_responses.tsv').open('w') as file:
        fields=['inhibition_multiplier','leak_multiplier','excitation_multiplier','retained_fraction',*FOCAL_NAMES];writer=csv.DictWriter(file,fieldnames=fields,delimiter='\t');writer.writeheader()
        for ii,ll,ee,ff in np.ndindex((9,9,3,4)):
            row=dict(inhibition_multiplier=M[ii],leak_multiplier=M[ll],excitation_multiplier=E[ee],retained_fraction=F[ff],**{key:out['focal_'+key][ii,ll,ee,ff] for key in FOCAL_NAMES});writer.writerow({k:'' if isinstance(v,float) and not np.isfinite(v) else v for k,v in row.items()})
    return dict(design_spot_curves=10692,design_focal_cases=972,unique_state_files=324,unique_drive_products=len(sorted({round(float(e*c/100),12) for e in E for c in C})),total_state_compute_s=sum(state_seconds),max_tail_ms=float(out['tail_start_ms'].max()),max_tail_auc_bound=float(out['tail_auc_bound'].max()),maximum_voltage=float(out['max_voltage_mV'].max()))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--workers',type=int,default=6);ap.add_argument('--aggregate-only',action='store_true');args=ap.parse_args()
    for sub in ['states','transient_traces']: (HERE/sub).mkdir(exist_ok=True)
    sources=list((MODEL/'src').glob('*.py'))+list((MODEL/'data').rglob('*.json'))+[MODEL/'data/cid3.swc',MODEL/'data/sm2nrn_cid3.mat',MODEL/'data/vg3_onbip_list.txt',MODEL/'data/vg3_offbip_list.txt',MODEL/'scripts/fig3/run_fig3_pilot.py',MODEL/'scripts/fig3/protocol_v3.py',MODEL/'scripts/fig3/run_anatomy_v3.py',HERE.parent/'protocol.md',HERE/'passive_tail.py',Path(__file__)]
    hashes={str(p.relative_to(ROOT)):digest(p) for p in sources}
    if (HERE/'input_hashes.json').exists():
        prior=json.loads((HERE/'input_hashes.json').read_text())
        if prior!=hashes:raise RuntimeError('Inputs changed since resumable run began')
    else:write_json(HERE/'input_hashes.json',hashes)
    start=time.time();tasks=list(np.ndindex((9,9,4)));missing=[t for t in tasks if not (HERE/'states'/f'i{t[0]}_l{t[1]}_f{t[2]}.npz').exists()]
    write_json(HERE/'progress.json',dict(status='running',completed=324-len(missing),total_states=324,workers=args.workers))
    if missing and args.aggregate_only:raise RuntimeError('Missing states')
    if missing:
        with mp.get_context('spawn').Pool(args.workers,initializer=initialize) as pool:
            for number,result in enumerate(pool.imap_unordered(run_state,missing),start=324-len(missing)+1):
                write_json(HERE/'progress.json',dict(status='running',completed=number,total_states=324,last=result,elapsed_s=time.time()-start))
                print(json.dumps(dict(completed=number,**result)),flush=True)
    initialize();summary=aggregate()
    if hashes!={str(p.relative_to(ROOT)):digest(p) for p in sources}:raise RuntimeError('Protected inputs modified')
    write_json(HERE/'run_manifest.json',dict(status='complete',elapsed_wall_s=time.time()-start,workers=args.workers,summary=summary,method='Exact duplicate effective-excitation caching; original NEURON dt=0.1ms/sample=2ms transient; all-arbor passive contraction bound <=1e-7 times accumulated AUC scale, followed by historical calcium recurrence under bounded late-voltage continuation. Baseline time skip exactly validated.',axes_order=['inhibition_multiplier','leak_multiplier','excitation_multiplier','retained_fraction','contrast'],numpy_version=np.__version__,neuron_version=__import__('neuron').__version__,input_hashes=hashes,primary_output_hashes={name:digest(HERE/name) for name in ['electrical_sweep.npz','spot_responses.tsv','focal_responses.tsv','tail_validation.json']}))
    write_json(HERE/'progress.json',dict(status='complete',completed=324,total_states=324,summary=summary));print(json.dumps(summary),flush=True)

if __name__=='__main__':main()
