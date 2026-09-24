#!/usr/bin/env python3
"""Calibrate one shared activation curve to two OFF calcium shape constraints.

The two 30% points are the only independent shape targets: per-condition 15/70%
normalization imposes both endpoints. Absolute AUCs are retained, reported, and
never rescaled separately. The primary protocol uses a sustained 3-s excitatory
conductance against fixed pure shunts; transient alternatives are diagnostics.
No measured synaptic waveform or VG3 calcium-channel parameters are identified.
"""
from pathlib import Path
import argparse, json, hashlib, time
import numpy as np
from scipy.optimize import least_squares
from scipy.special import expit
from scipy.signal import lfilter
from run_fig3_pilot import (ROOT, build, load_named_config, load_fits,
    exc_eff_weight, inh_eff_weight, gaussian_spot_weights, simulate)
from calcium_readout import calcium_trace
from protocol_v3 import simulate_v3

CONTRASTS=(15,30,70)
CONDITIONS=(("ctrl",1.,"Ctrl"),("gaba",0.,"cKO"))
BOUNDS=np.array([[-45.,1.],[-20.,6.]])
TAU=50.
SAMPLE_DT=2.
TSTOP=3052.


def shared_readout(v, t, vhalf, slope, rest=-45., tau=TAU):
    """AUC of spatial mean equals mean of local AUCs (linear clearance)."""
    j=(expit((v-vhalf)/slope)-expit((rest-vhalf)/slope)).mean(axis=0)
    dt=t[1]-t[0]; decay=np.exp(-dt/tau); gain=tau*(-np.expm1(-dt/tau))
    ca=np.zeros_like(j); ca[1:]=lfilter([gain],[1.,-decay],j[:-1])
    window=(t>=50-1e-6)&(t<=3050+1e-6)
    return float(np.trapezoid(ca[window],x=t[window]))


def fit_activation(traces, t, targets, rest=-45., shared_scale_compromise=False):
    def values(p):
        return np.array([[shared_readout(traces[f'OFF_{cond}_{c}'],t,*p,rest=rest)
                          for c in CONTRASTS] for cond,_,_ in CONDITIONS])
    def shape(y): return (y[:,1]-y[:,0])/(y[:,2]-y[:,0])
    def residual(p):
        y=values(p);r=shape(y)-targets
        if shared_scale_compromise:
            r=np.r_[r,max(0.,y[1,2]/y[0,2]-1.3)*.3]
        return r
    options=[]
    for initial in [(-40,2),(-35,2),(-30,2),(-25,1.5),(-35,5)]:
        result=least_squares(residual,initial,bounds=BOUNDS,ftol=1e-10,xtol=1e-10,gtol=1e-10,max_nfev=200)
        y=values(result.x); sh=shape(y); delta=y[1]-y[0]
        options.append(dict(v_half_mV=float(result.x[0]),slope_mV=float(result.x[1]),
                            shape30=sh.tolist(),squared_shape_error=float(np.sum((sh-targets)**2)),objective_sse=float(np.sum(residual(result.x)**2)),
                            absolute_auc=y.tolist(),gaba_over_control=(y[1]/y[0]).tolist(),
                            absolute_auc_increment=delta.tolist(),
                            low_over_high_absolute_increment=float(delta[0]/delta[2]),
                            boundary_hits=[name for i,name in enumerate(['v_half','slope']) if min(abs(result.x[i]-BOUNDS[:,i]))<1e-4],
                            optimizer_success=bool(result.success)))
    options.sort(key=lambda r:r['objective_sse'])
    return options[0],options


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--parameters',type=Path,required=True)
    ap.add_argument('--outdir',type=Path,default=ROOT/'outputs/fig3_v3')
    ap.add_argument('--excitation-scale',type=float,default=1.)
    ap.add_argument('--fit-only',action='store_true')
    ap.add_argument('--shared-scale-compromise',action='store_true',help='Add a disclosed soft penalty on 70% calcium gain above 1.3; a qualitative model criterion, not a measured bound.')
    args=ap.parse_args(); out=args.outdir;out.mkdir(parents=True,exist_ok=True)
    parameters=json.loads(args.parameters.read_text())
    parameters['excitation_scale']=args.excitation_scale
    cfg,_=load_named_config('sep_g22'); fits=load_fits()
    cache=out/'calibration_voltage_traces.npz'
    if args.fit_only:
        a=np.load(cache); t=a['t_ms']; traces={k:a[k] for k in a.files if k!='t_ms'}
        cached=json.loads((out/'voltage_parameters.json').read_text())
        if cached!=parameters: raise ValueError('Cached voltage parameters differ.')
        rest=-45.
    else:
        model=build(parameters['g_pas_S_cm2']); coords=model['node_coords'];rest=model['rest_v']
        traces={}; started=time.time()
        for polarity in ['OFF']:
            nodes=model[f'{polarity.lower()}_nodes'];inh=model['inh_nodes']
            ew=exc_eff_weight(fits,polarity,100)*cfg['G_exc']*parameters['excitation_scale']*gaussian_spot_weights(coords[nodes,:2],model['soma_xy'],100,25)
            iw=inh_eff_weight(fits,polarity,100)*cfg['G_inh']*parameters['inhibition_scale']*gaussian_spot_weights(coords[inh,:2],model['soma_xy'],100,25)
            for cond,remaining,_ in CONDITIONS:
                for c in CONTRASTS:
                    t,v=simulate_v3(model,nodes,ew,iw,nodes,parameters,remaining,c/100,dt=.1,sample_dt=SAMPLE_DT,tstop=TSTOP)
                    traces[f'{polarity}_{cond}_{c}']=v.astype(np.float32)
                    print(f'{polarity} {cond} {c}: mean peak ΔV {np.mean(v.max(1)-rest):.5f}; max V {v.max():.3f}',flush=True)
        np.savez_compressed(cache,t_ms=t,**traces)
        (out/'voltage_parameters.json').write_text(json.dumps(parameters,indent=2)+'\n')
        print(f'Voltage runs: {time.time()-started:.1f}s',flush=True)
    exppath=ROOT.parent/'analysis/fig3_model/experimental_fig3.json'
    experiment=json.loads(exppath.read_text())
    targets=np.array([experiment['groups'][f'{exp}_OFF']['contrasts']['0.3']['normalized_mean'] for _,_,exp in CONDITIONS])
    best,options=fit_activation(traces,t,targets,rest,args.shared_scale_compromise)
    baseline=np.array([[shared_readout(traces[f'OFF_{cond}_{c}'],t,-40,4,rest=rest) for c in CONTRASTS] for cond,_,_ in CONDITIONS])
    baseline_shape=(baseline[:,1]-baseline[:,0])/(baseline[:,2]-baseline[:,0])
    v2=np.load(ROOT/'outputs/fig3_pilot/experiment_comparison_voltage_traces.npz')
    v2traces={f'OFF_{cond}_{c}':v2[f"OFF_{'control' if cond=='ctrl' else 'gaba_removed'}_{c/100:.2f}"] for cond,_,_ in CONDITIONS for c in CONTRASTS}
    v2values=np.array([[shared_readout(v2traces[f'OFF_{cond}_{c}'],v2['t_ms'],-40,4,rest=rest) for c in CONTRASTS] for cond,_,_ in CONDITIONS])
    v2shape=(v2values[:,1]-v2values[:,0])/(v2values[:,2]-v2values[:,0])
    report=dict(status='illustrative_OFF_shape_calibration_only_not_channel_identification',parameters=parameters,
        selection_objective='OFF shape SSE + 0.09*max(0, high-contrast Ca gain-1.3)^2; 1.3 is an illustrative soft compression target, not a measured bound.' if args.shared_scale_compromise else 'OFF shape SSE only',
        calcium=dict(v_half_mV=best['v_half_mV'],slope_mV=best['slope_mV'],tau_ms=TAU,mode='fixed'),
        bounds=dict(v_half_mV=BOUNDS[:,0].tolist(),slope_mV=BOUNDS[:,1].tolist(),tau_ms='fixed at 50'),
        targets_OFF_shape30=targets.tolist(),selected=best,multistart_results=options,
        same_voltage_v2_activation_shape30=baseline_shape.tolist(),
        original_v2_shape30=v2shape.tolist(),original_v2_squared_shape_error=float(np.sum((v2shape-targets)**2)),
        protocols=dict(dt_ms=.1,sample_dt_ms=SAMPLE_DT,tstop_ms=TSTOP,auc_window_ms=[50,3050],
                       readout='Mean of signed, local, baseline-subtracted calcium AUC over 75 fixed OFF contacts.',
                       limitations='Fixed pure-shunting inhibition and sustained 3-s excitation isolate the mechanism; no measured synaptic time course is fitted; model contact samples unregistered to imaging sublaminae.' if parameters.get('excitation_protocol')=='sustained_3s' else 'Single inherited 10/300-ms conductance event, not literal 3-s square stimulus; model contact samples unregistered to imaging sublaminae.'),
        interpretation='Only two OFF interior shape values calibrated. 15 and 70 percent endpoints are imposed separately per condition only for shape comparison. Absolute AUC remains shared scale. No indicator, release, driving-force, or membrane feedback transformation.',
        hashes={str(exppath):hashlib.sha256(exppath.read_bytes()).hexdigest(),str(cache):hashlib.sha256(cache.read_bytes()).hexdigest(),str(Path(__file__).resolve()):hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),str(Path(__file__).with_name('protocol_v3.py')):hashlib.sha256(Path(__file__).with_name('protocol_v3.py').read_bytes()).hexdigest()})
    (out/'calibration.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['selected','original_v2_shape30','original_v2_squared_shape_error']},indent=2),flush=True)

if __name__=='__main__':main()
