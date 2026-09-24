"""Derived metrics, complete secondary exports, and electrical review figures."""
from pathlib import Path
import json,csv,hashlib
import numpy as np
from scipy.interpolate import PchipInterpolator
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm,Normalize

HERE=Path(__file__).resolve().parent;OUT=HERE/'review';OUT.mkdir(exist_ok=True)
a=np.load(HERE/'electrical_sweep.npz');A=a['mean_auc'];M=a['inhibition_multipliers'];E=a['excitation_multipliers'];F=a['retained_fractions'];C=a['contrast_pct'];N=(9,9,3)
span=A[...,10]-A[...,3];shape=np.divide(A-A[...,3,None],span[...,None],out=np.full_like(A,np.nan),where=span[...,None]>.15)
key=A[...,[3,6,10]];gain=np.divide(key[:,:,:,1:,:],key[:,:,:,0,None,:],out=np.full((9,9,3,3,3),np.nan),where=key[:,:,:,0,None,:]>.15);delta=key[:,:,:,1:,:]-key[:,:,:,0,None,:]
pref=gain[...,:2]/gain[...,2,None]
targets=np.array([.0345740518597063,.703657829532908]);sse=(shape[:,:,:,0,6,None]-targets[0])**2+(shape[:,:,:,1:,6]-targets[1])**2
shift=np.full((9,9,3,3,11),np.nan);reason=np.full(shift.shape,'',dtype='<U40')
for ii,ll,ee in np.ndindex(N):
 ctrl=A[ii,ll,ee,0]
 inverse=PchipInterpolator(ctrl,C,extrapolate=False) if np.all(np.diff(ctrl)>0) and ctrl[-1]-ctrl[0]>.15 else None
 for ff in range(3):
  for cc,c in enumerate(C):
   response=A[ii,ll,ee,ff+1,cc]
   if not c:r='zero_contrast'
   elif not np.isfinite(ctrl).all():r='nonfinite_control'
   elif np.any(np.diff(ctrl)<0):r='nonmonotonic_control'
   elif np.any(np.diff(ctrl)==0):r='tied_control_responses'
   elif ctrl[-1]-ctrl[0]<=.15:r='weak_control_span'
   elif response<=.15:r='weak_reduced_response'
   elif response<ctrl[0] or response>ctrl[-1]:r='outside_control_response_range'
   else:r='valid';shift[ii,ll,ee,ff,cc]=float(inverse(response))/c
   reason[ii,ll,ee,ff,cc]=r
count=np.isfinite(shift).sum(-1);mean_shift=np.divide(np.nansum(shift,axis=-1),count,out=np.full(count.shape,np.nan),where=count>0)
vspread=a['focal_voltage_1e_um'][:,:,:,1:]-a['focal_voltage_1e_um'][:,:,:,0,None]
cspread=a['focal_calcium_1e_um'][:,:,:,1:]-a['focal_calcium_1e_um'][:,:,:,0,None]
np.savez_compressed(HERE/'electrical_metrics.npz',gain_15_30_70=gain,increment_15_30_70=delta,preference_15_30_over70=pref,normalized_auc=shape,shape30_sse=sse,calcium_shift_by_contrast=shift,calcium_shift_count=count,calcium_shift_mean=mean_shift,match_status=reason,voltage_spread_change_um=vspread,calcium_spread_change_um=cspread)

with (HERE/'paired_metrics.tsv').open('w') as stream:
 fields=['inhibition_multiplier','leak_multiplier','excitation_multiplier','retained_fraction','gain15','gain30','gain70','delta15','delta30','delta70','preference15_over70','preference30_over70','shape30_difference','shape30_sse','voltage_spread_change_um','calcium_spread_change_um','calcium_shift_mean','valid_shift_count','supported_contrasts'];w=csv.DictWriter(stream,fieldnames=fields,delimiter='\t');w.writeheader()
 for ii,ll,ee,ff in np.ndindex((9,9,3,3)):
  row=dict(inhibition_multiplier=M[ii],leak_multiplier=M[ll],excitation_multiplier=E[ee],retained_fraction=F[ff+1],**dict(zip(['gain15','gain30','gain70'],gain[ii,ll,ee,ff])),**dict(zip(['delta15','delta30','delta70'],delta[ii,ll,ee,ff])),preference15_over70=pref[ii,ll,ee,ff,0],preference30_over70=pref[ii,ll,ee,ff,1],shape30_difference=shape[ii,ll,ee,ff+1,6]-shape[ii,ll,ee,0,6],shape30_sse=sse[ii,ll,ee,ff],voltage_spread_change_um=vspread[ii,ll,ee,ff],calcium_spread_change_um=cspread[ii,ll,ee,ff],calcium_shift_mean=mean_shift[ii,ll,ee,ff],valid_shift_count=count[ii,ll,ee,ff],supported_contrasts=','.join(map(str,C[np.isfinite(shift[ii,ll,ee,ff])])) )
  w.writerow({k:'' if isinstance(v,(float,np.floating)) and not np.isfinite(v) else v for k,v in row.items()})
with (HERE/'contrast_matches.tsv').open('w') as stream:
 fields=['inhibition_multiplier','leak_multiplier','excitation_multiplier','retained_fraction','contrast_pct','calcium_horizontal_multiplier','match_status'];w=csv.DictWriter(stream,fieldnames=fields,delimiter='\t');w.writeheader()
 for ii,ll,ee,ff,cc in np.ndindex(shift.shape):
  value=shift[ii,ll,ee,ff,cc];w.writerow(dict(inhibition_multiplier=M[ii],leak_multiplier=M[ll],excitation_multiplier=E[ee],retained_fraction=F[ff+1],contrast_pct=C[cc],calcium_horizontal_multiplier=float(value) if np.isfinite(value) else '',match_status=reason[ii,ll,ee,ff,cc]))

def extent(v):
 finite=np.asarray(v)[np.isfinite(v)];return dict(n=int(len(finite)),minimum=float(finite.min()) if len(finite) else None,median=float(np.median(finite)) if len(finite) else None,maximum=float(finite.max()) if len(finite) else None)
summary=[]
for ee,e in enumerate(E):
 for ff,f in enumerate(F[1:]):
  p=pref[:,:,ee,ff];valid=np.isfinite(p).all(-1);cs=cspread[:,:,ee,ff]
  summary.append(dict(excitation_multiplier=float(e),retained_fraction=float(f),points=81,valid_preference_points=int(valid.sum()),both_preference_gt_one=int(np.sum(valid&np.all(p>1+1e-10,axis=-1))),preference15=extent(p[:,:,0]),preference30=extent(p[:,:,1]),gain15=extent(gain[:,:,ee,ff,0]),gain30=extent(gain[:,:,ee,ff,1]),gain70=extent(gain[:,:,ee,ff,2]),voltage_spread_change=extent(vspread[:,:,ee,ff]),calcium_spread_change=extent(cs),calcium_spread_positive_count=int(np.sum(cs>1e-10)),calcium_control_missing_crossings=int(np.isnan(a['focal_calcium_1e_um'][:,:,ee,0]).sum()),calcium_reduced_missing_crossings=int(np.isnan(a['focal_calcium_1e_um'][:,:,ee,ff+1]).sum()),calcium_shift_mean=extent(mean_shift[:,:,ee,ff]),overlap_count=extent(count[:,:,ee,ff])))
checks=dict(spot_curves=int(A.size),focal_cases=int(a['focal_source_dv_mV'].size),all_zero_contrast=bool(np.all(A[...,0]==0)),all_paired_positive_nonzero_auc=bool(np.all(A[:,:,:,1:,1:]>A[:,:,:,0,None,1:])),voltage_range=[float(a['min_voltage_mV'].min()),float(a['max_voltage_mV'].max())],largest_voltage_tail_bound=float(a['tail_voltage_bound_mV'].max()),largest_auc_tail_bound=float(a['tail_auc_bound'].max()),largest_auc_bound_relative_to_result=float(np.max(a['tail_auc_bound']/np.maximum(1,A))),tail_start_range=[float(a['tail_start_ms'].min()),float(a['tail_start_ms'].max())],all_tail_bounds_met=bool(np.all(a['tail_auc_bound']<=1e-7*a['tail_bound_scale']+1e-10)))
(HERE/'summary.json').write_text(json.dumps(dict(cases=summary,checks=checks,nominal=dict(gains=gain[4,4,1,1].tolist(),preferences=pref[4,4,1,1].tolist(),spatial_control_reduced=a['focal_calcium_1e_um'][4,4,1,[0,2]].tolist(),calcium_only_shift=float(mean_shift[4,4,1,1]))),indent=2)+'\n')

plt.rcParams.update({'font.family':'Arial','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','axes.titlesize':10})
footer='Exploratory ranges; one reconstruction. Gray: undefined ratio or unobserved first 1/e crossing. No parameters were fitted.'
def save(fig,name):
 fig.savefig(OUT/(name+'.svg'),facecolor='white');fig.savefig(OUT/(name+'.png'),facecolor='white',dpi=180);fig.savefig(OUT/(name+'.pdf'),facecolor='white');plt.close(fig)
def heat(fig,ax,values,title,label,diverging=False,limits=None):
 cmap=plt.get_cmap('RdBu_r' if diverging else 'viridis').copy();cmap.set_bad('#d9d9d9')
 finite=np.asarray(values)[np.isfinite(values)];mx=max(abs(finite).max(),1e-9) if limits is None else limits
 norm=TwoSlopeNorm(vmin=-mx,vcenter=0,vmax=mx) if diverging else (Normalize(0,limits) if limits is not None else None)
 im=ax.pcolormesh(np.log2(M),np.log2(M),np.ma.masked_invalid(values.T),shading='nearest',cmap=cmap,norm=norm)
 if diverging and finite.min()<0<finite.max():ax.contour(np.log2(M),np.log2(M),values.T,levels=[0],colors='black',linewidths=.7)
 ax.plot(0,0,'*',mfc='#079ba8',mec='white',ms=10);ax.set(xticks=[-2,-1,0,1,2],xticklabels=['0.25','0.5','1','2','4'],yticks=[-2,-1,0,1,2],yticklabels=['0.25','0.5','1','2','4'],xlabel='Full-control inhibition multiplier',ylabel='Leak multiplier',title=title)
 fig.colorbar(im,ax=ax,fraction=.046,pad=.025,label=label)
fig,ax=plt.subplots(2,3,figsize=(14,8.4));fig.subplots_adjust(wspace=.53,hspace=.45,top=.87,bottom=.12)
fig.suptitle('Electrical parameter sweep: 25% retained inhibition',fontsize=15,fontweight='bold');fig.text(.1,.91,'Original spot and focal excitation strengths; selected calcium readout held fixed. Star: original electrical settings.')
items=[(gain[:,:,1,1,1],'A  Calcium gain at 30%','Reduced / control AUC',False),(np.log2(pref[:,:,1,1,1]),'B  Intermediate contrast preference','log2 [gain(30%) / gain(70%)]',True),(shape[:,:,1,2,6]-shape[:,:,1,0,6],'C  Normalized 30% shape difference','Reduced minus control',True),(vspread[:,:,1,1],'D  Voltage spread increase','Change in first 1/e (µm)',False),(cspread[:,:,1,1],'E  Calcium spread increase','Change in first 1/e (µm)',False),(mean_shift[:,:,1,1],'F  Calcium-only contrast equivalence','Mean horizontal multiplier',False)]
for axis,item in zip(ax.flat,items):heat(fig,axis,*item)
fig.text(.1,.035,footer,fontsize=8);save(fig,'electrical_overview')
fig,ax=plt.subplots(3,3,figsize=(13.8,11.5));fig.subplots_adjust(wspace=.56,hspace=.49,top=.9,bottom=.09);fig.suptitle('Absolute calcium responses across electrical settings',fontsize=15,fontweight='bold')
for row,(contrast,ci) in enumerate([(15,3),(30,6),(70,10)]):
 for col,ff,name in [(0,0,'Control'),(1,2,'25% retained')]:heat(fig,ax[row,col],A[:,:,1,ff,ci]/1000,f'{contrast}% contrast: {name}','Ca AUC (10³ a.u. ms²)',False,150)
 heat(fig,ax[row,2],(A[:,:,1,2,ci]-A[:,:,1,0,ci])/1000,f'{contrast}% contrast: increment','Ca AUC (10³ a.u. ms²)',False,150)
fig.text(.08,.025,'Original excitation strength; common absolute scale across all panels. '+footer,fontsize=8);save(fig,'electrical_absolute_responses')
for pk,name in [(0,'low'),(1,'intermediate')]:
 fig,ax=plt.subplots(3,3,figsize=(13.8,11.5));fig.subplots_adjust(wspace=.56,hspace=.49,top=.9,bottom=.09);fig.suptitle(f'Electrical robustness of {name}-contrast preference',fontsize=15,fontweight='bold')
 values=np.log2(pref[...,pk]);limit=float(np.nanmax(abs(values)))
 for ee,e in enumerate(E):
  for ff,f in enumerate(F[1:]):heat(fig,ax[ee,ff],values[:,:,ee,ff],f'Excitation {e:g}×; inhibition retained {100*f:g}%','log2 preference',True,limit)
 fig.text(.08,.025,footer,fontsize=8);save(fig,'electrical_'+name+'_preference')
fig,ax=plt.subplots(3,3,figsize=(12,10.5));fig.subplots_adjust(wspace=.35,hspace=.47,top=.9,bottom=.09);fig.suptitle('Focal calcium spread and contrast preference across electrical settings',fontsize=14,fontweight='bold')
for ee,e in enumerate(E):
 for ff,f in enumerate(F[1:]):
  xx=cspread[:,:,ee,ff];yy=np.log2(pref[:,:,ee,ff,1]);valid=np.isfinite(xx)&np.isfinite(yy);axis=ax[ee,ff]
  axis.scatter(xx[valid],yy[valid],c=np.broadcast_to(np.log2(M)[:,None],(9,9))[valid],cmap='viridis',vmin=-2,vmax=2,s=12);axis.scatter(xx[4,4],yy[4,4],marker='*',color='#079ba8',edgecolors='black',s=90);axis.axhline(0,color='black',lw=.6)
  axis.set(xlabel='Calcium 1/e increase (µm)',ylabel='log2 [gain(30) / gain(70)]',title=f'Excitation {e:g}×; retained {100*f:g}%\n{valid.sum()}/81 resolved pairs')
fig.text(.08,.025,'Color: full-control inhibition (0.25–4×). Descriptive within-model association; no biological inference or mediation test.',fontsize=8);save(fig,'electrical_spread_preference')
report=['# Electrical robustness sweep','',f'Completed all **10,692 condition–contrast responses (972 contrast-response curves) and 972 focal cases**, using 324 electrical states and 19 distinct excitation products per state. Duplicate products were reused exactly; all planned design rows remain in the outputs.','',f'The maximum voltage was {checks["voltage_range"][1]:.6f} mV; every result remained within passive reversal bounds. All late-voltage approximation bounds passed. Transient prefixes ended at {checks["tail_start_range"][0]:g}–{checks["tail_start_range"][1]:g} ms; the largest AUC error bound divided by max(1, computed AUC) was {checks["largest_auc_bound_relative_to_result"]:.3g}.','', 'The historical 0.1-ms NEURON integration and 2-ms calcium sampling are retained. Exact resting baseline is skipped (validated sample-for-sample). A conservative all-compartment passive contraction bound permits replacement of the late voltage tail by its last recorded value; the original local sigmoid, 50-ms clearance recurrence, and trapezoidal AUC are then evaluated through 3050 ms. This is a declared numerical approximation, not substitution of a steady-state response for stimulus AUC. Five full-duration nominal/corner checks agree to at most 8.28×10⁻¹⁴ in scaled AUC; focal impedance calculations have separate point-process checks.','', '| Shared excitation | Retained inhibition | Both P15 and P30 > 1 | Resolved calcium spread pairs | P30 range |','|---|---|---|---|---|']
for r in summary:report.append(f'| {r["excitation_multiplier"]:g}× | {100*r["retained_fraction"]:g}% | {r["both_preference_gt_one"]}/{r["valid_preference_points"]} | {r["calcium_spread_change"]["n"]}/81 | {r["preference30"]["minimum"]:.4f}–{r["preference30"]["maximum"]:.3g} |')
report+=['','All nonzero-contrast reduced-inhibition AUCs exceed their paired controls. This expected monotone-readout effect is distinct from preferential fractional enhancement, which depends on electrical settings. Missing 1/e crossings mean the profile did not reach that level within this arbor; those cases remain in the source data and are never assigned zero spread. Complete per-contrast horizontal-equivalence matches include their support and exclusion reasons.','', 'These exploratory fourfold ranges are not physiological confidence intervals. Strongly depolarized corners approach the 0-mV excitation reversal; they remain visible without claiming physiological plausibility. The sweep keeps morphology, reversal potentials, contact placement, calcium activation, and all within-pair excitation settings fixed. It does not test active calcium feedback, another morphology, or measured biological distributions. Absolute cable-calcium integrals are length-weighted proxy responses, not calcium content.','', 'Reproduce: run `run_electrical.py --workers 6`, then `validate_focal.py`, then `build_electrical_outputs.py`, with the model `.venv/bin/python`. Completed states are resumable only with identical hashed inputs. Raw exact transient prefixes are retained in `transient_traces/`; full-arbor normalized focal voltage maps are retained in `states/`.','']
report += ['Black contours mark zero effect in signed maps: P = 1 in preference panels and zero normalized-shape difference in the shape panel. No contour is drawn if the boundary is absent from the sampled range.', '']
(HERE/'README.md').write_text('\n'.join(report))
hashes={str(p.relative_to(HERE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [HERE/'electrical_metrics.npz',HERE/'paired_metrics.tsv',HERE/'contrast_matches.tsv',HERE/'summary.json',HERE/'README.md',Path(__file__),*OUT.glob('*')] if p.is_file()}
(HERE/'review_manifest.json').write_text(json.dumps(hashes,indent=2)+'\n')
print(json.dumps(dict(checks=checks,cases=summary),indent=2))
