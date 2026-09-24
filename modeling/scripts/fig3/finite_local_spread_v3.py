#!/usr/bin/env python3
"""Finite local excitation with shared v3 shunts and calcium activation.

A fixed local excitatory conductance at the same source produces 15 mV in
control; its exact strength stays identical after GABA removal. The full steady
voltage solution includes excitation conductance loading. Calcium is local
steady activation above baseline times the shared 50 ms linear clearance time.
This is a focal drive illustration, not a measured unitary synaptic event.
"""
import argparse,json,hashlib
import numpy as np
from scipy.special import expit
from calibrate_spatial_v3 import *

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--v-half',type=float,default=-37.47893);ap.add_argument('--slope',type=float,default=1.39843);ap.add_argument('--tau',type=float,default=50.);ap.add_argument('--calibration',type=Path);args=ap.parse_args()
 if args.calibration:
  cal=json.loads(args.calibration.read_text())['calcium'];args.v_half=cal['v_half_mV'];args.slope=cal['slope_mV'];args.tau=cal['tau_ms']
 archive=np.load(OUT/'selected_spatial_maps.npz');maps={k:archive[k] for k in archive.files};saved=json.loads((OUT/'selected_spatial_summary.json').read_text());selection=saved['selection'];source=int(maps['local_source_segment']);d=maps['path_distance_from_source_um'];lengths=maps['cable_lengths'];rest=-45.;erev=0.;baseline=expit((rest-args.v_half)/args.slope);rin_ctrl=saved['summary']['ctrl']['source_Rin_MOhm'];rows=[];checks=[]
 for target in [10.,12.5,15.]:
  g=target/(rin_ctrl*(erev-rest-target))
  for cond in ['ctrl','gaba']:
   transfer=maps[f'local_{cond}_normalized'];rin=saved['summary'][cond]['source_Rin_MOhm'];source_dv=(erev-rest)*g*rin/(1+g*rin);dv=transfer*source_dv;ca=(expit((rest+dv-args.v_half)/args.slope)-baseline)*args.tau;norm_ca=ca/ca[source]
   row=dict(control_source_target_dv_mV=target,condition=cond,source_conductance_uS=g,source_dv_mV=source_dv,source_calcium_au_ms=float(ca[source]),voltage=metrics(dv/source_dv,d,lengths),calcium=metrics(norm_ca,d,lengths),absolute_ca_weighted_cable_au_ms_um=float(np.sum(ca*lengths)))
   rows.append(row)
   checks.append(dict(check=f'{cond}_{target:g}mV_finite_and_bounded',passed=bool(np.isfinite(dv).all() and np.isfinite(ca).all() and ca.min()>=-1e-12 and ca.max()<=args.tau+1e-12)))
   if target==15.:
    maps[f'local_{cond}_ca']=ca;maps[f'local_{cond}_finite_dv']=dv;maps[f'local_{cond}_finite_ca_normalized']=norm_ca;maps[f'local_{cond}_finite_v_normalized']=dv/source_dv
 primary=[r for r in rows if r['control_source_target_dv_mV']==15.];g=primary[0]['source_conductance_uS'];model=build(selection['g_pas_S_cm2']);segs=model['segments'];coords=model['node_coords'];inh=model['inh_nodes'];fits=load_fits();cfg,_=load_named_config('sep_g22');iw=inh_eff_weight(fits,'OFF',100)*cfg['G_inh']*selection['inhibition_scale']*gaussian_spot_weights(coords[inh,:2],model['soma_xy'],100,25)
 for cond,frac in [('ctrl',1.),('gaba',selection['glycine_fraction'])]:
  clamps=[]
  for idx,w in zip(inh,iw*frac):
   if w<=0:continue
   clamp=h.SEClamp(segs[int(idx)]);clamp.dur1=151.;clamp.amp1=rest;clamp.rs=1./float(w);clamps.append(clamp)
  exc=h.SEClamp(segs[source]);exc.dur1=151.;exc.amp1=erev;exc.rs=1./g;clamps.append(exc)
  h.CVode().active(0);h.dt=.1;h.finitialize(rest);h.continuerun(150.)
  observed=np.array([s.v-rest for s in segs]);error=float(np.max(np.abs(observed-maps[f'local_{cond}_finite_dv'])));checks.append(dict(check=f'{cond}_finite_voltage_matches_independent_NEURON_point_conductance',passed=error<1e-6,max_error_mV=error))
  del clamp,exc,clamps
 checks.append(dict(check='same_local_conductance_between_conditions',passed=primary[0]['source_conductance_uS']==primary[1]['source_conductance_uS']))
 np.savez_compressed(OUT/'finite_local_maps.npz',**maps)
 meta=dict(status='finite focal-drive illustration with shared calcium activation; not fitted unitary excitation or functional length constant',electrical_parameters=selection,calcium=dict(v_half_mV=args.v_half,slope_mV=args.slope,tau_ms=args.tau,mode='fixed driving, activation plus linear clearance'),primary_control_source_target_dv_mV=15.,finite_sourceweight_uS=g,finite_sourceweight_nS=g*1000,source_segment=source,source_xyz_um=maps['local_source_xyz'].tolist(),source_definition='single effective focal conductance at fixed OFF-contact compartment4705; amplitude selected to place control source15mVabove rest; not measured unitary synaptic strength',primary=primary,sensitivity=rows,validation=checks,limitations=['The same conductance is used in both conditions and at the same source. It is a finite focal excitation chosen for illustration, not an independently fitted stimulus or a unitary synapse.','Calcium maps use a shared absolute scale of0to50arbitrary units times milliseconds; no per-condition normalization is applied to anatomical calcium maps. Profile normalization is explicit and source-specific.','Steady activation times linear clearance is the DC calcium concentration proxy, not fluorescence or a finite-stimulus calcium AUC. The spot simulations supply the corresponding finite-time integration.','Ca1/e distances are first crossings of cable-length-weighted5um path bins after source normalization. They are amplitude-dependent descriptive distances on one branched arbor, not the published functional calcium lambda.','No calcium-current feedback, calcium driving-force dependence, indicator saturation, or release model is included.'],hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),ROOT/'outputs/fig3_v3_spatial/selected_parameters.json']})
 (OUT/'finite_local_summary.json').write_text(json.dumps(meta,indent=2));print(json.dumps(dict(primary=primary,validation=checks),indent=2),flush=True)
 if not all(c['passed'] for c in checks):raise SystemExit('finite local validation failed')
if __name__=='__main__':main()
