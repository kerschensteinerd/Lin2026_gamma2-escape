#!/usr/bin/env python3
"""Select an illustrative paper-leak model retaining30%glycinergic shunt.

An optional exact two-slope voltage-surrogate calibration is diagnostic only.
Neither analysis is a fitted or independently validated calcium imaging assay.
"""
import hashlib,json,time
import numpy as np
from scipy.optimize import brentq
from neuron import h
from calibrate_spatial_v3 import *
from run_fig3_pilot import simulate
from run_anatomy_maps import morphology_geometry

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--exact-surrogate',action='store_true');parser.add_argument('--glycine-fraction',type=float,default=.3);args=parser.parse_args();prefix='exact_surrogate_' if args.exact_surrogate else '';start=time.time();cfg,_=load_named_config('sep_g22');leak=1e-4;model=build(leak);fits=load_fits();segs=model['segments'];coords=model['node_coords'];off=model['off_nodes'];inh=model['inh_nodes'];areas=np.array([h.area(s.x,sec=s.sec) for s in segs]);lengths=np.array([s.sec.L/s.sec.nseg for s in segs]);source=4705;d=path_distance(model,source)
 iw=inh_eff_weight(fits,'OFF',100)*cfg['G_inh']*gaussian_spot_weights(coords[inh,:2],model['soma_xy'],100,25)
 base=np.bincount(inh,weights=iw,minlength=len(segs));positions=np.arange(-150,151,10.);bw=bar_weights(coords[off,:2],model['soma_xy'],positions);x=coords[off,0]-model['soma_xy'][0]
 cache={}
 def evaluate(strength):
  if strength in cache:return cache[strength]
  apply_shunts(model,leak,base*strength,areas);z=spatial_transfer(model,off,off);response=z@bw;centers,r2=rf_centers(response,positions);reg=regress_binned(x,centers)
  cache[strength]=(reg,z,response,centers,r2);return cache[strength]
 strengths=({c:float(brentq(lambda q:evaluate(q)[0]['slope']-target,a,b,xtol=1e-7)) for c,target,a,b in [('ctrl',.877,3.,10.),('gaba',.162,.03,.3)]} if args.exact_surrogate else {'ctrl':3.,'gaba':3.*args.glycine_fraction})
 frac=strengths['gaba']/strengths['ctrl'];selection=dict(status=('exactly calibrated small-signal voltage surrogate; diagnostic only' if args.exact_surrogate else 'illustrative shared-parameter model;30%glycinergic conductance retained; not a calcium imaging fit'),g_pas_S_cm2=leak,inhibition_scale=strengths['ctrl'],glycine_fraction=frac,G_inh=cfg['G_inh']*strengths['ctrl'],G_exc_inherited=cfg['G_exc'],inh_rev_spot_mV=-45.,inh_rev_local_mV=model['rest_v'],effective_retained_inhibition_scale=strengths['gaba'],target_rf_center_slopes={'ctrl':.877,'gaba':.162},passive_unbranched_lambda_um=float(1e4*np.sqrt(.5e-4/(4*100*leak))),source_segment=source,total_control_inhibition_uS=float(iw.sum()*strengths['ctrl']),total_retained_glycine_uS=float(iw.sum()*strengths['gaba']),total_passive_leak_uS=float(leak*areas.sum()*.01),local_probe='infinitesimal NEURON DC current transfer under fixed pure shunt; normalized at source')
 (OUT/(prefix+'selected_parameters.json')).write_text(json.dumps(selection,indent=2));print(json.dumps(selection,indent=2),flush=True)
 maps=morphology_geometry(model);maps.update(path_distance_from_source_um=d,local_source_segment=np.asarray(source),local_source_xyz=coords[source],bar_positions_um=positions,bar_weights=bw,off_x_um=x)
 checks=[];summary={}
 for condition,strength in strengths.items():
  reg,z,response,centers,r2=evaluate(strength);apply_shunts(model,leak,base*strength,areas);lm,rin=local_map(model,source)
  maps[f'local_{condition}_normalized']=lm;maps[f'local_{condition}_dc_dv_mV']=lm*rin*.001;maps[f'rf_{condition}_centers_um']=centers;maps[f'rf_{condition}_response_mohm']=response;maps[f'rf_{condition}_transfer_mohm']=z
  summary[condition]=dict(**reg,**metrics(lm,d,lengths),source_Rin_MOhm=rin,gaussian_rf_r2_min=float(r2.min()),gaussian_rf_r2_mean=float(r2.mean()),inhibition_scale=strength,total_inhibition_uS=float(iw.sum()*strength),remote_at_100um_cable_weighted=float(np.average(lm[d>=90],weights=lengths[d>=90])))
  checks.append(dict(check=f'{condition}_reciprocal_transfer_matrix',passed=bool(np.allclose(z,z.T,rtol=1e-9,atol=1e-9)),max_abs_mohm=float(np.max(np.abs(z-z.T)))))
  checks.append(dict(check=f'{condition}_RF_center_slope_finite',passed=bool(np.isfinite(reg['slope'])),slope=reg['slope'],target=selection['target_rf_center_slopes'][condition],target_difference=reg['slope']-selection['target_rf_center_slopes'][condition]))
  checks.append(dict(check=f'{condition}_bounded_positive_normalized_transfer',passed=bool(np.all((lm>=0)&(lm<=1+1e-9)))))
  # Independent original point-process shunts and current step, no distributed shunts.
  apply_shunts(model,leak,np.zeros_like(base),areas)
  t,v=simulate(model,[],[],iw*strengths['ctrl'],np.arange(len(segs)),gaba_fraction=1-frac,gaba_remaining=1. if condition=='ctrl' else 0.,tonic=True,e_inh=model['rest_v'],current_node=source,current_nA=.001,tstop=150.,sample_dt=1.)
  step=v[:,(t>=140)&(t<=145)].mean(1)-model['rest_v'];normal=step/step[source];err=float(np.max(np.abs(normal-lm)));rinerr=float(abs(step[source]/.001-rin)/rin)
  maps[f'local_{condition}_currentstep_dv_mV']=step;maps[f'local_{condition}_currentstep_normalized']=normal
  checks.append(dict(check=f'{condition}_DC_agrees_original_point_shunt_step',passed=err<1e-4 and rinerr<1e-4,max_absolute_transfer_error=err,relative_Rin_error=rinerr))
  del v
 # Distribution uncertainty: same calibrated total inhibitory load fractions from1%to30%.
 residual_rows=[]
 for r in [.01,.03,.1,.3,frac]:
  reg,_,_,_,_=evaluate(strengths['ctrl']*r);residual_rows.append(dict(glycine_fraction=r,rf_center_slope=reg['slope']))
 checks.append(dict(check='glycinergic_conductance_positive_and_identical_between_conditions',passed=0<frac<1,total_glycine_uS=selection['total_retained_glycine_uS']))
 checks.append(dict(check='source_is_fixed_OFF_contact',passed=bool(source in off)))
 checks.append(dict(check='finite_maps',passed=all(np.isfinite(a).all() for a in maps.values())))
 np.savez_compressed(OUT/(prefix+'selected_spatial_maps.npz'),**maps)
 result=dict(selection=selection,summary=summary,retained_fraction_sensitivity=residual_rows,calibration_protocol=dict(stimulus='50 x80um rectangular bars at x offsets-150..150 in10umsteps, with25um Gaussian blur; y centered at soma',drive='infinitesimal equal current injection at75OFFcontacts weighted by bar overlap',background='fixed inherited100um OFFspot inhibitory conductance spatialweights, pure shunt atrest; independently validated withpointprocessclamps',readout='Gaussian fit(nooffset) to linear steady voltage response at75OFFcontacts, regress10umpositionbinmean RFcenters versus meanlocation',sampling='oneEM reconstruction;75synapticcontactsites, notimagingpixelorretinasampling',targets='approximate Fig2groupmean slopes0.877/0.162; descriptive comparison, not matched in the30%glycine primary illustration',sigma_um=25,bar_size_um=[50,80]),limitations=['The inferred glycine fraction is an effective co-located shunt fraction, not transmitter synapse count or measured receptor conductance.','Primary illustration uses paperlike leak,3x inherited inhibition and30% residualglycine; it improves spatialisolation/mixing but does notmatch bothFig2slopes. An optional exactvoltage-surrogatecalibration requires only1.22%residualglycine and is notselected asprimary.','The surrogate uses small-signal voltage, fixed shunts and OFFcontact sampling, whereas experimental RF centers are derived from flashed calcium signals; complete protocol matching remains outstanding.','Primary100umspot simulations use the samefixed pure-shunt background; phasicECl-60mV comparisons are separate diagnostic protocols.','Paper functional calcium lambda16um is contextual; this voltage transfer profile is neither its estimate nor its fittedtarget.','One source and one reconstruction do not estimate biological variability.'],elapsed_seconds=time.time()-start,hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),ROOT/'scripts/fig3/calibrate_spatial_v3.py',ROOT/'data/cid3.swc',ROOT/'data/sm2nrn_cid3.mat']})
 (OUT/(prefix+'selected_spatial_summary.json')).write_text(json.dumps(result,indent=2));(OUT/(prefix+'validation.json')).write_text(json.dumps(checks,indent=2));print(json.dumps(dict(summary=summary,checks=checks),indent=2),flush=True)
 if not all(c['passed'] for c in checks):raise SystemExit('validation failed')
if __name__=='__main__':main()
