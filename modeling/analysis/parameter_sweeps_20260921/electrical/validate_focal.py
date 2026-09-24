"""Independent SEClamp validation of impedance-derived finite focal maps."""
from pathlib import Path
import sys,json
import numpy as np
from neuron import h
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];MODEL=ROOT
import argparse
ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument('--outdir',type=Path,default=HERE/'reproduction')
args=ap.parse_args()
HERE=args.outdir.resolve();HERE.mkdir(parents=True,exist_ok=True)
sys.path.insert(0,str(MODEL/'scripts/fig3'))
from run_fig3_pilot import build,load_fits,load_named_config,inh_eff_weight,gaussian_spot_weights
cfg,_=load_named_config('sep_g22');fits=load_fits();checks=[]
cases=[(4,4,0,1),(4,4,2,1),(0,0,1,0),(8,8,0,2),(8,0,3,2)]
for ii,ll,ff,ee in cases:
 path=HERE/'states'/f'i{ii}_l{ll}_f{ff}.npz'
 if not path.exists():raise SystemExit(f'Waiting for {path.name}')
 a=np.load(path);mi=float(a['mi']);ml=float(a['ml']);f=float(a['retained_fraction']);e=[.5,1,2][ee]
 m=build(ml*1e-4);segs=m['segments'];inh=m['inh_nodes'];xy=m['node_coords'];iw=inh_eff_weight(fits,'OFF',100)*cfg['G_inh']*3*mi*f*gaussian_spot_weights(xy[inh,:2],m['soma_xy'],100,25)
 clamps=[]
 for idx,w in zip(inh,iw):
  c=h.SEClamp(segs[int(idx)]);c.dur1=151;c.amp1=-45;c.rs=1/float(w);clamps.append(c)
 c=h.SEClamp(segs[4705]);c.dur1=151;c.amp1=0;c.rs=1/(.012785038085113628*e);clamps.append(c)
 h.CVode().active(0);h.secondorder=0;h.dt=.1;h.finitialize(-45);h.continuerun(150)
 observed=np.array([s.v+45 for s in segs]);expected=a['focal_normalized_voltage_map']*a['focal_source_dv_mV'][ee]
 error=float(np.max(np.abs(observed-expected)));row=dict(inhibition_multiplier=mi,leak_multiplier=ml,retained_fraction=f,excitation_multiplier=e,max_all_arbor_voltage_error_mV=error,passed=error<1e-7);checks.append(row);print(json.dumps(row),flush=True)
 del c,clamps
(HERE/'focal_validation.json').write_text(json.dumps(checks,indent=2)+'\n')
if not all(r['passed'] for r in checks):raise AssertionError('focal validation failed')
