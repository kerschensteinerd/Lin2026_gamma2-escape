from pathlib import Path
import sys,time,json
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];MODEL=ROOT
import argparse
ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument('--outdir',type=Path,default=HERE/'reproduction')
args=ap.parse_args()
HERE=args.outdir.resolve();HERE.mkdir(parents=True,exist_ok=True)
sys.path.insert(0,str(MODEL/'scripts/fig3'))
from run_fig3_pilot import build,load_fits,load_named_config,exc_eff_weight,inh_eff_weight,gaussian_spot_weights
from protocol_v3 import simulate_v3
from run_anatomy_v3 import readout
from passive_tail import simulate_spot
cfg,_=load_named_config('sep_g22');fits=load_fits();checks=[]
cases=[(1,1,1,.25,.3),(1,1,1,1,.7),(.25,.25,.5,.2,.05),(.25,.25,2,.2,.7),(4,4,.5,1,.05)]
for mi,ml,me,f,c in cases:
 m=build(ml*1e-4);n=m['off_nodes'];ins=m['inh_nodes'];xy=m['node_coords'];ew=exc_eff_weight(fits,'OFF',100)*cfg['G_exc']*gaussian_spot_weights(xy[n,:2],m['soma_xy'],100,25);iw=inh_eff_weight(fits,'OFF',100)*cfg['G_inh']*3*mi*f*gaussian_spot_weights(xy[ins,:2],m['soma_xy'],100,25)
 start=time.time();fast=simulate_spot(m,ew,iw,c,me,return_trace=True);fasttime=time.time()-start
 p=dict(glycine_fraction=1.,excitation_protocol='sustained_3s',inhibition_protocol='tonic_pure_shunt')
 start=time.time();t,v=simulate_v3(m,n,ew*me,iw,n,p,0,c,tstop=3052);fulltime=time.time()-start
 ca=dict(v_half_mV=-37.47892660844856,slope_mV=1.3984142481316566,tau_ms=50,mode='fixed')
 dv,peak,auc=readout(v,t,ca,-45)
 use=(t>=50-1e-6)&(t<=3050+1e-6)
 maxtrace=float(np.max(np.abs(fast['v_mV']-v[:,use])))
 exactprefix=float(np.max(np.abs(fast['recorded_v_mV']-v[:,use][:,:len(fast['recorded_t_ms'])])))
 error=float(abs(fast['mean_ca_auc']-auc.mean()));relative=error/max(abs(float(auc.mean())),1)
 check=dict(mi=mi,ml=ml,me=me,f=f,contrast=c,fast_s=fasttime,full_s=fulltime,tail_start_ms=fast['tail_start_ms'],prefix_max_mV=exactprefix,max_trace_mV=maxtrace,voltage_bound_mV=fast['tail_voltage_error_bound_mV'],auc_error=error,auc_scaled_error=relative,auc_bound=fast['tail_auc_error_bound'],mean_peak_dv_error=float(abs(fast['mean_peak_dv']-dv.mean())),mean_peak_ca_error=float(abs(fast['mean_peak_ca']-peak.mean())),passed=relative<1e-7 and exactprefix<1e-8 and error<=fast['tail_auc_error_bound']+1e-7 and maxtrace<=fast['tail_voltage_error_bound_mV']+1e-8)
 checks.append(check);print(json.dumps(check),flush=True)
 (HERE/'tail_validation.json').write_text(json.dumps(checks,indent=2)+'\n')
 if not check['passed']:raise AssertionError('tail validation failed')
