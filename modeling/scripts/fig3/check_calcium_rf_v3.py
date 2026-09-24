#!/usr/bin/env python3
"""Steady local calcium RF check with the SAME v3 shunts and spot activation.

NEURON transfer resistances + exact synaptic conductance-load correction give the
DC membrane solution. Same sigmoid activation acts locally before Gaussian RF
fits. Finite stimulus kinetics / experimental pixel AUC are not represented.
"""
import argparse,json
import numpy as np
from scipy.special import expit
from calibrate_spatial_v3 import *
from run_fig3_pilot import exc_eff_weight

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--v-half',type=float,default=-35.);ap.add_argument('--slope',type=float,default=2.);ap.add_argument('--exc-scale',type=float,default=1.);ap.add_argument('--contrast',type=float,default=1.);a=ap.parse_args()
 m=np.load(OUT/'selected_spatial_maps.npz');summary=json.loads((OUT/'selected_spatial_summary.json').read_text());selection=summary['selection'];fits=load_fits();cfg,_=load_named_config('sep_g22')
 off=m['off_nodes'];positions=m['bar_positions_um'];bw=m['bar_weights'];excweight=exc_eff_weight(fits,'OFF',100)*cfg['G_exc']*a.exc_scale*a.contrast;g=bw*excweight;out={};arrays={};checks=[]
 for cond in ['ctrl','gaba']:
  z=m[f'rf_{cond}_transfer_mohm'];dv=np.stack([np.linalg.solve(np.eye(len(off))+z*gj[None,:],z@gj*45.) for gj in g.T],axis=1)
  ca=expit((-45+dv-a.v_half)/a.slope)-expit((-45-a.v_half)/a.slope)
  ctr,r2=rf_centers(ca,positions);vctr,vr2=rf_centers(dv,positions)
  out[cond]=dict(calcium_regression=regress_binned(m['off_x_um'],ctr),voltage_regression=regress_binned(m['off_x_um'],vctr),calcium_gaussian_r2_min=float(r2.min()),calcium_gaussian_r2_mean=float(r2.mean()),absolute_v_range_mV=[float(-45+dv.min()),float(-45+dv.max())])
  arrays[f'{cond}_voltage_dv']=dv;arrays[f'{cond}_calcium_activation']=ca;arrays[f'{cond}_rf_center_um']=ctr
  checks.append(dict(check=f'{cond}_finite_positive_calcium',passed=bool(np.isfinite(ca).all() and np.min(ca)>=-1e-10)))
 # Independent point-process conductances validate the matrix load correction.
 model=build(selection['g_pas_S_cm2']);segments=model['segments'];coords=model['node_coords'];inh=model['inh_nodes']
 iw=inh_eff_weight(fits,'OFF',100)*cfg['G_inh']*selection['inhibition_scale']*gaussian_spot_weights(coords[inh,:2],model['soma_xy'],100,25)
 middle=int(np.argmin(np.abs(positions)))
 for cond,retained in [('ctrl',1.),('gaba',selection['glycine_fraction'])]:
  clamps=[]
  for nodes,weights,reversal in [(inh,iw*retained,-45.),(off,g[:,middle],0.)]:
   for idx,weight in zip(nodes,weights):
    if weight<=0:continue
    c=h.SEClamp(segments[int(idx)]);c.dur1=151.;c.amp1=reversal;c.rs=1./float(weight);clamps.append(c)
  h.CVode().active(0);h.dt=.1;h.finitialize(-45.);h.continuerun(150.)
  observed=np.array([segments[int(i)].v+45. for i in off]);expected=arrays[f'{cond}_voltage_dv'][:,middle];error=float(np.max(np.abs(observed-expected)))
  checks.append(dict(check=f'{cond}_DC_matrix_matches_NEURON_point_process',passed=error<1e-6,max_error_mV=error))
  del c,clamps
 result=dict(status='steady-state calcium RF consistency diagnostic, not a fit to flashed imaging',parameters=vars(a),validation=checks,electrical_selection=selection,summary=out,limitations=['Bar geometry and Gaussian center fit match the spatial diagnostic, but the readout samples 75 OFF contacts on one EM reconstruction; it does not match imaging planes or pixel/retina replicates.','DC membrane potentials use exact excitation conductance loading. Calcium is steady activation above a common resting baseline; linear clearance would multiply all RF responses by the same factor.','No finite channel/inactivation kinetics, fluorescence indicator saturation, or finite flashed-stimulus AUC; these RFs are a diagnostic, not an experimentally matched calcium fit.','The retained glycinergic conductance is an illustrative fraction distributed at the same contacts, not a measured fraction.'])
 (OUT/'calcium_rf_check.json').write_text(json.dumps(result,indent=2));np.savez_compressed(OUT/'calcium_rf_check.npz',bar_positions_um=positions,off_x_um=m['off_x_um'],**arrays)
 print(json.dumps(result,indent=2),flush=True)
 if not all(c['passed'] for c in checks):raise SystemExit('calcium RF consistency validation failed')
if __name__=='__main__':main()
