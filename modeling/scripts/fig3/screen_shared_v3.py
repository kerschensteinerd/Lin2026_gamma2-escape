"""Steady-state screening to expose calcium shape/amplitude tradeoffs."""
import json
import numpy as np
from scipy.special import expit
from scipy.optimize import least_squares
from run_fig3_pilot import *
from protocol_v3 import simulate_v3

out=ROOT/'outputs/fig3_v3';cfg,_=load_named_config('sep_g22');m=build(.0001);f=load_fits();xy=m['node_coords'];o=m['off_nodes'];i=m['inh_nodes']
ew=exc_eff_weight(f,'OFF',100)*cfg['G_exc']*gaussian_spot_weights(xy[o,:2],m['soma_xy'],100,25)
iw=inh_eff_weight(f,'OFF',100)*cfg['G_inh']*3*gaussian_spot_weights(xy[i,:2],m['soma_xy'],100,25)
rows=[];archive={};target=np.array([.0345740518597,.7036578295329])
for gly in [.1,.2,.3]:
 for escale in [1.,3.,5.]:
  pp=dict(g_pas_S_cm2=.0001,inhibition_scale=3.,glycine_fraction=gly,excitation_scale=escale,inhibition_protocol='tonic_pure_shunt',excitation_protocol='sustained_3s')
  v=[]
  for rem in [1.,0.]:
   vc=[]
   for c in [.15,.3,.7]:
    t,vv=simulate_v3(m,o,ew*escale,iw,o,pp,rem,c,tstop=100.,sample_dt=1.)
    vc.append(vv[:,-1])
   v.append(vc)
  v=np.array(v);archive[f'gly{gly}_exc{escale}']=v
  def evaluate(p):
   act=(expit((v-p[0])/p[1])-expit((-45-p[0])/p[1])).mean(-1)
   shape=(act[:,1]-act[:,0])/(act[:,2]-act[:,0]);delta=act[1]-act[0]
   return dict(v_half_mV=float(p[0]),slope_mV=float(p[1]),shape30=shape.tolist(),shape_sse=float(np.sum((shape-target)**2)),ca=act.tolist(),gain=(act[1]/act[0]).tolist(),delta=delta.tolist(),low_over_high=float(delta[0]/delta[2]))
  def objective(p):
   r=evaluate(p)
   return np.r_[np.array(r['shape30'])-target, max(0.,r['gain'][2]-1.3)*.3]
  candidates=[]
  for start in [(-40,2),(-32,2),(-25,2)]:
   result=least_squares(objective,start,bounds=([-45,1],[-20,6]),ftol=1e-9,xtol=1e-9,gtol=1e-9,max_nfev=100)
   r=evaluate(result.x);r['objective_sse']=float(np.sum(objective(result.x)**2));candidates.append(r)
  r=min(candidates,key=lambda x:x['objective_sse']);r['parameters']=pp;r['mean_peak_dv']=np.mean(v+45,axis=-1).tolist()
  rows.append(r);print(json.dumps(r),flush=True)
np.savez_compressed(out/'shared_screen_voltage.npz',**archive)
(out/'shared_screen.json').write_text(json.dumps(rows,indent=2))
