#!/usr/bin/env python3
"""Mechanistic spread sensitivity. Does not alter any production parameters.

Compare saved weak-EPSP source-peak snapshots with a DC passive transfer solution.
DC shunts are inserted as exactly equivalent segment passive conductances; an
independent SEClamp + IClamp run validates that representation. Attenuation
statistics describe this branched reconstruction and are NOT a fit to imaging.
"""
from pathlib import Path
import json
import numpy as np
from scipy.stats import linregress
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from neuron import h
from run_fig3_pilot import ROOT, build, load_fits, load_named_config, inh_eff_weight, gaussian_spot_weights, path_distance, simulate

OUT = ROOT / 'outputs/fig3_diagnostics_20260909'
OUT.mkdir(exist_ok=True, parents=True)
A = np.load(ROOT / 'outputs/fig3_anatomy/anatomy_maps.npz')
source = int(A['local_source_segment'])
d = A['path_distance_from_source_um']; lengths=A['cable_lengths']
source_xyz=A['xyz_segments'][source]
deuclid=np.linalg.norm(A['xyz_segments']-source_xyz,axis=1)


def binned(y, distance=d, binsize=5):
    xx=[];yy=[];ww=[]
    for lo in np.arange(0, distance.max()+binsize, binsize):
        use=(distance>=lo)&(distance<lo+binsize)
        if use.any():
            xx.append(float(np.average(distance[use],weights=lengths[use])))
            yy.append(float(np.average(y[use],weights=lengths[use])))
            ww.append(float(lengths[use].sum()))
    return np.array(xx),np.array(yy),np.array(ww)


def crossing(x,y,level=1/np.e):
    # Source 0,1 anchors the response; report first crossing, not a fitted lambda.
    x=np.r_[0,x];y=np.r_[1,y]
    candidates=np.where((y[:-1]>level)&(y[1:]<=level))[0]
    if len(candidates)==0:return None
    i=candidates[0]
    return float(x[i]+(level-y[i])*(x[i+1]-x[i])/(y[i+1]-y[i]))


def metrics(y):
    x,b,w=binned(y); xe,be,we=binned(y,deuclid)
    fits=[]
    for lo,hi in [(5,60),(5,100),(10,40),(20,80)]:
        use=(x>=lo)&(x<=hi)&(b>0)
        fit=linregress(x[use],np.log(b[use]))
        fits.append(dict(range_um=[lo,hi],lambda_um=float(-1/fit.slope),amplitude=float(np.exp(fit.intercept)),log_r_squared=float(fit.rvalue**2),n_bins=int(use.sum())))
    return dict(weighted_cable_um=float(np.sum(lengths*y)),path_1e_crossing_um=crossing(x,b),euclidean_1e_crossing_um=crossing(xe,be),log_exponential_fits=fits,
                path_bins_um=x.tolist(),path_bin_transfer=b.tolist(),path_bin_cable_um=w.tolist())

cfg,_=load_named_config('sep_g22');model=build(cfg['g_pas']);fits=load_fits()
segs=model['segments'];coords=model['node_coords'];inh=model['inh_nodes'];rest=model['rest_v']
areas=np.array([h.area(s.x,sec=s.sec) for s in segs])
iweights=inh_eff_weight(fits,'OFF',100)*cfg['G_inh']*gaussian_spot_weights(coords[inh,:2],model['soma_xy'],100,25)
per_node=np.bincount(inh,weights=iweights,minlength=len(segs))
assert np.allclose(path_distance(model,source),d)

saved={k:metrics(A[f'local_{k}_at_source_peak_normalized']) for k in ['ctrl','gaba']}

# Independent current step with original point-process shunts. At gpas=1e-3,
# intrinsic tau <= 1 ms, so 140-145 ms is long after source current starts at50.
t,v=simulate(model,[],[],iweights,np.arange(len(segs)),tonic=True,e_inh=rest,gaba_fraction=.7,gaba_remaining=0.,current_node=source,current_nA=.001,tstop=150.,sample_dt=1.)
step_response=np.mean(v[:,(t>=140)&(t<=145)],axis=1)-rest
step_normalized=step_response/step_response[source]
del v

rows=[];maps={};validation=[]
for gbase in [.001,.0001]:
    for retained in [1.,.5,.3,.1,0.]:
        for i,s in enumerate(segs):
            s.g_pas=gbase+100*per_node[i]*retained/areas[i]
        h.finitialize(rest)
        imp=h.Impedance();imp.loc(segs[source].x,sec=segs[source].sec);imp.compute(0.)
        input_mohm=float(imp.input(segs[source].x,sec=segs[source].sec))
        transfer=np.array([imp.transfer(s.x,sec=s.sec) for s in segs])/input_mohm
        key=f'gpas_{gbase:g}_retained_{retained:g}'
        maps[key]=transfer
        row=dict(key=key,g_pas_S_per_cm2=gbase,inhibitory_fraction_retained=retained,
                 input_resistance_MOhm=input_mohm,total_leak_uS=float(gbase*areas.sum()*.01),
                 total_inhibition_uS=float(iweights.sum()*retained),
                 infinite_unbranched_cable_lambda_um=float(1e4*np.sqrt(.5e-4/(4*100*gbase))))
        row.update(metrics(transfer));rows.append(row)
        print(key,'Rin',input_mohm,'cross',row['path_1e_crossing_um'],'cable',row['weighted_cable_um'],flush=True)
        if gbase==.001 and retained==.3:
            err=float(np.max(np.abs(step_normalized-transfer)))
            rin_err=float(abs(step_response[source]/.001-input_mohm)/input_mohm)
            validation.append(dict(check='DC distributed passive shunts match original point SEClamp current-step transfer',max_absolute_normalized_transfer_error=err,relative_Rin_error=rin_err,passed=err<1e-8 and rin_err<1e-8))
            validation.append(dict(check='DC vs saved weak EPSP source-peak snapshot',max_absolute_normalized_transfer_difference=float(np.max(np.abs(A['local_gaba_at_source_peak_normalized']-transfer))),weighted_cable_fractional_difference=float((metrics(transfer)['weighted_cable_um']/saved['gaba']['weighted_cable_um'])-1),interpretation='Protocol/time variation, not a numerical failure.'))
        if gbase==.001 and retained==1.:
            validation.append(dict(check='Control DC vs saved weak EPSP source-peak snapshot',max_absolute_normalized_transfer_difference=float(np.max(np.abs(A['local_ctrl_at_source_peak_normalized']-transfer))),weighted_cable_fractional_difference=float((metrics(transfer)['weighted_cable_um']/saved['ctrl']['weighted_cable_um'])-1),interpretation='Protocol/time variation, not a numerical failure.'))

np.savez_compressed(OUT/'spread_sensitivity_maps.npz',distance_um=d,source_xyz=source_xyz,**maps)
result=dict(status='diagnostic only; no selected/fitted changes to production',source_segment=source,n_segments=len(segs),total_membrane_area_um2=float(areas.sum()),full_inhibitory_load_uS=float(iweights.sum()),saved_weak_EPSP=saved,sensitivity=rows,validation=validation,
 caveats=["Source-normalized voltage transfer on one branched morphology. Neither first 1/e crossing nor log-slope fit is the physiological calcium length constant.","First 1/e crossing is from cable-length-weighted 5-um bins at the cable-weighted mean distance, with source (0,1) anchor; not a fitted exponential.","Fitted lambdas are descriptive, use log response versus path distance in explicit windows, equal weight per nonempty bin, free amplitude, zero asymptotic floor.","Each DC point-source stimulus is infinitesimal; excitation amplitude does not affect normalized passive DC transfer.","Inhibitory load is the inherited 100-um OFF-spot peak conductance held tonically at rest, not endogenous background inhibition.","Removal retains specified uniform fraction of every inhibitory contact. Identity/location/strength of glycinergic input are not measured. 0% retained is an all-inhibition-off diagnostic ceiling, not the GABA-specific cKO.","Lowering leak without refitting changes input resistance and all voltage responses. Sensitivity does not justify using new parameters without independent calibration."])
(OUT/'spread_diagnostic.json').write_text(json.dumps(result,indent=2))
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
fig,axs=plt.subplots(2,2,figsize=(10,7),layout='constrained')
colors={1.:'#333333',.5:'#c99462',.3:'#e47c29',.1:'#bb3d13',0.:'#8d1111'}
for ax,gbase in zip(axs[0],[.001,.0001]):
    for row in rows:
        if row['g_pas_S_per_cm2']!=gbase:continue
        r=row['inhibitory_fraction_retained'];ax.plot(row['path_bins_um'],row['path_bin_transfer'],color=colors[r],label=f'{100*r:g}% inhibition retained')
    ax.axhline(1/np.e,color='.7',ls=':',lw=.8)
    ax.set(xlim=(0,105),ylim=(0,1),xlabel='Dendritic path distance (µm)',ylabel='Remote / source ΔV',title=f'Passive DC transfer; gₚₐₛ = {gbase:g} S/cm²')
axs[0,0].legend(frameon=False,fontsize=8)
for gbase,style in [(.001,'-'),(.0001,'--')]:
    rr=[r for r in rows if r['g_pas_S_per_cm2']==gbase]
    retained=np.array([r['inhibitory_fraction_retained'] for r in rr])
    base=rr[0]['weighted_cable_um']
    axs[1,0].plot(100*(1-retained),[r['weighted_cable_um']/base for r in rr],style+'o',label=f'gₚₐₛ = {gbase:g}')
    axs[1,1].plot(100*(1-retained),[r['path_1e_crossing_um'] for r in rr],style+'o',label=f'gₚₐₛ = {gbase:g}')
axs[1,0].set(xlabel='Inhibitory conductance removed (%)',ylabel='Transfer-weighted cable / control',title='Spread ratio; fixed inhibitory load')
axs[1,1].set(xlabel='Inhibitory conductance removed (%)',ylabel='First binned 1/e crossing (µm)',title='Descriptive attenuation distance')
axs[1,1].text(.98,.04,'Low leak, 90–100% removal:\nno 1/e crossing within this arbor (100 µm)',transform=axs[1,1].transAxes,ha='right',va='bottom',fontsize=8,color='#9c5300')
for ax in axs[1]:ax.legend(frameon=False,fontsize=8);ax.axvline(70,color='.8',ls=':',lw=.8)
fig.suptitle('Diagnostic sensitivity: leak and residual inhibition\nOne local source, inherited spatial inhibitory weights; not a fit to physiological calcium spread',fontsize=12)
fig.savefig(OUT/'spread_sensitivity.png',dpi=180,facecolor='white');fig.savefig(OUT/'spread_sensitivity.svg',facecolor='white');plt.close(fig)
print(json.dumps(dict(saved={k:{x:v for x,v in s.items() if not x.startswith('path_bin')} for k,s in saved.items()},validation=validation,total_area=areas.sum(),inhibition=iweights.sum()),indent=2),flush=True)
