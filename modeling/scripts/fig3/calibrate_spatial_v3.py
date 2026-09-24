#!/usr/bin/env python3
"""Real NEURON spatial-mixing diagnostic. Fixed shunts, infinitesimal bar drive.
Compared with Figure2 RF-center slope summaries, not a fitted calcium assay.
"""
from pathlib import Path
import argparse,json,time
import numpy as np
from scipy.special import ndtr
from scipy.optimize import curve_fit
from scipy.stats import linregress
from neuron import h
from run_fig3_pilot import ROOT,build,load_fits,load_named_config,inh_eff_weight,gaussian_spot_weights,path_distance
OUT=ROOT/'outputs/fig3_v3_spatial';OUT.mkdir(parents=True,exist_ok=True)
def bar_weights(xy,center,positions,sigma=25.):
 x=xy[:,0,None]-center[0]-positions[None,:];y=xy[:,1,None]-center[1]
 return (ndtr((25-x)/sigma)-ndtr((-25-x)/sigma))*(ndtr((40-y)/sigma)-ndtr((-40-y)/sigma))
def gaussian(x,amp,mu,sigma):return amp*np.exp(-.5*((x-mu)/sigma)**2)
def rf_centers(response,positions):
 centers=[];r2=[]
 for r in response:
  r=r/max(r.max(),1e-30);center=float(np.sum(r*positions)/r.sum())
  try:
   fit,_=curve_fit(gaussian,positions,r,p0=[1,center,40],bounds=([0,-150,5],[2,150,150]),maxfev=1000)
   centers.append(float(fit[1]));r2.append(float(1-np.sum((r-gaussian(positions,*fit))**2)/np.sum((r-r.mean())**2)))
  except RuntimeError:centers.append(center);r2.append(0.)
 return np.array(centers),np.array(r2)
def regress_binned(x,y):
 bx=[];by=[];n=[]
 for lo in np.arange(np.floor(x.min()/10)*10,np.ceil(x.max()/10)*10,10):
  mask=(x>=lo)&(x<lo+10)
  if mask.any():bx.append(float(x[mask].mean()));by.append(float(y[mask].mean()));n.append(int(mask.sum()))
 reg=linregress(bx,by)
 return dict(slope=float(reg.slope),intercept_um=float(reg.intercept),r_squared=float(reg.rvalue**2),bin_x_um=bx,bin_rf_center_um=by,bin_n=n)
def apply_shunts(model,gbase,pernode,areas):
 for i,s in enumerate(model['segments']):s.g_pas=gbase+100*pernode[i]/areas[i];s.e_pas=model['rest_v']
 h.finitialize(model['rest_v'])
def spatial_transfer(model,record_nodes,source_nodes):
 segs=model['segments'];out=np.empty((len(record_nodes),len(source_nodes)));imp=h.Impedance()
 for j,idx in enumerate(source_nodes):
  s=segs[int(idx)];imp.loc(s.x,sec=s.sec);imp.compute(0.)
  out[:,j]=[imp.transfer(segs[int(i)].x,sec=segs[int(i)].sec) for i in record_nodes]
 return out
def local_map(model,source=4705):
 segs=model['segments'];imp=h.Impedance();s=segs[source];imp.loc(s.x,sec=s.sec);imp.compute(0.)
 rin=float(imp.input(s.x,sec=s.sec));transfer=np.array([imp.transfer(x.x,sec=x.sec) for x in segs])
 return transfer/rin,rin
def metrics(y,d,lengths):
 xx=[];yy=[]
 for lo in np.arange(0,d.max()+5,5):
  mask=(d>=lo)&(d<lo+5)
  if mask.any():xx.append(float(np.average(d[mask],weights=lengths[mask])));yy.append(float(np.average(y[mask],weights=lengths[mask])))
 x=np.r_[0,xx];b=np.r_[1,yy];cross=np.flatnonzero((b[:-1]>1/np.e)&(b[1:]<=1/np.e));crossing=None
 if len(cross):
  i=cross[0];crossing=float(x[i]+(1/np.e-b[i])*(x[i+1]-x[i])/(b[i+1]-b[i]))
 return dict(weighted_cable_um=float(np.sum(lengths*y)),path_1e_crossing_um=crossing,path_bins_um=xx,path_bin_transfer=yy)
def main():
 parser=argparse.ArgumentParser();parser.add_argument('--quick',action='store_true');args=parser.parse_args()
 cfg,_=load_named_config('sep_g22');model=build(.0001);fits=load_fits();segs=model['segments'];coords=model['node_coords'];off=model['off_nodes'];inh=model['inh_nodes'];areas=np.array([h.area(s.x,sec=s.sec) for s in segs]);lengths=np.array([s.sec.L/s.sec.nseg for s in segs]);d=path_distance(model,4705)
 iw=inh_eff_weight(fits,'OFF',100)*cfg['G_inh']*gaussian_spot_weights(coords[inh,:2],model['soma_xy'],100,25)
 base=np.bincount(inh,weights=iw,minlength=len(segs));positions=np.arange(-150,151,10.);bw=bar_weights(coords[off,:2],model['soma_xy'],positions)
 rows=[];maps={};start=time.time();leaks=[1e-4,3e-5,1e-5] if not args.quick else [1e-4];strengths=[.01,.03,.1,.3,1,3,10] if not args.quick else [.01,.03,.1,1,3]
 for leak in leaks:
  for strength in strengths:
   apply_shunts(model,leak,base*strength,areas);z=spatial_transfer(model,off,off);response=z@bw;centers,r2=rf_centers(response,positions)
   key=f'g{leak:g}_i{strength:g}';lm,rin=local_map(model)
   row=dict(key=key,g_pas_S_cm2=leak,inhibition_scale=strength,total_inhibition_uS=float(iw.sum()*strength),Rin_MOhm=rin,gaussian_rf_r2_min=float(r2.min()),gaussian_rf_r2_mean=float(r2.mean()),**regress_binned(coords[off,0]-model['soma_xy'][0],centers),**metrics(lm,d,lengths))
   rows.append(row);maps[key+'_local']=lm;maps[key+'_rf_centers']=centers;maps[key+'_rf_response']=response;maps[key+'_transfer_mohm']=z
   print(json.dumps({k:v for k,v in row.items() if not isinstance(v,list)}),flush=True)
   (OUT/'spatial_scan.json').write_text(json.dumps(dict(status='diagnostic constraint, not calcium imaging fit',rows=rows,elapsed_s=time.time()-start),indent=2))
 np.savez_compressed(OUT/'spatial_scan.npz',**maps,off_nodes=off,off_x_um=coords[off,0]-model['soma_xy'][0],bar_positions_um=positions,bar_weights=bw,path_distance_um=d,cable_lengths=lengths,source_segment=4705)
 print('DONE',time.time()-start,flush=True)
if __name__=='__main__':main()
