#!/usr/bin/env python3
"""Detailed exports and scientific review plots; no new simulation or fitting."""
from pathlib import Path
import csv,json,hashlib,sys
import numpy as np
from scipy.special import expit
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm,Normalize
from matplotlib.lines import Line2D

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
SOURCE=ROOT/'outputs/fig4_v9_retained_inhibition_20260916'
import argparse
ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument('--results',type=Path,default=HERE/'results')
ap.add_argument('--outdir',type=Path,default=HERE/'generated/calcium_review')
args=ap.parse_args()
RESULTS=args.results.resolve()
OUT=args.outdir.resolve();OUT.mkdir(parents=True,exist_ok=True)
with np.load(RESULTS/'calcium_sweep.npz') as archive:
    z={key:archive[key] for key in archive.files}
s=json.loads((RESULTS/'summary.json').read_text())
pairs=z['parameter_pairs'];n=2601;half=z['v_half_grid_mV'];slope=z['slope_grid_mV'];contrasts=z['contrast_pct']
auc=z['mean_auc'];shape=z['normalized_auc'];gain=z['gain_15_30_70'];pref=z['preference_15_30_over70'];inc=z['increment_15_30_70'];shift=z['calcium_shift_mean'];counts=z['calcium_shift_count'];cross=z['spatial_1e_ctrl25_um'];profiles=z['spatial_profiles_ctrl25'];x=z['spatial_distance_um']
star=pairs[-1];historical=np.flatnonzero((pairs[:,0]==-40)&np.isclose(pairs[:,1],4))[0]
plt.rcParams.update({'font.family':'Arial','font.size':9,'axes.titlesize':10,'axes.labelsize':9,'axes.linewidth':.6,'xtick.labelsize':8,'ytick.labelsize':8,'svg.fonttype':'none','savefig.dpi':200,'axes.spines.top':False,'axes.spines.right':False})
GRAY='#333333';ORANGE='#d55e00';TEAL='#079ba8'


def save(fig,name):
    fig.savefig(OUT/f'{name}.svg',facecolor='white')
    fig.savefig(OUT/f'{name}.png',facecolor='white',dpi=180)
    plt.close(fig)


def mark(ax):
    ax.plot(*star,'*',markersize=10,mfc=TEAL,mec='white',mew=.7,zorder=10)
    ax.plot(-40,4,'D',markersize=5,mfc='none',mec='black',mew=.9,zorder=10)


def heat(fig,ax,data,title,label,cmap='viridis',norm=None,contour=None):
    a=np.asarray(data[:n]).reshape(51,51).T
    palette=plt.get_cmap(cmap).copy();palette.set_bad('#d9d9d9')
    h=ax.pcolormesh(half,slope,np.ma.masked_invalid(a),shading='nearest',cmap=palette,norm=norm,rasterized=False)
    if contour is not None and np.nanmin(a)<contour<np.nanmax(a):
        ax.contour(half,slope,a,levels=[contour],colors='black',linewidths=.8)
    ax.set(xlabel='Activation midpoint (mV)',ylabel='Slope (mV)',title=title,xlim=(-45.25,-19.75),ylim=(.95,6.05))
    ax.set_xticks([-45,-40,-35,-30,-25,-20]);ax.set_yticks([1,2,3,4,5,6]);mark(ax)
    cb=fig.colorbar(h,ax=ax,fraction=.046,pad=.025);cb.set_label(label,fontsize=8);cb.ax.tick_params(labelsize=7)
    return h


def centered(data):
    lim=float(np.nanmax(np.abs(data)));return TwoSlopeNorm(vmin=-lim,vcenter=0,vmax=lim)


def base(rows,cols,size,title,subtitle):
    fig,axes=plt.subplots(rows,cols,figsize=size)
    fig.subplots_adjust(left=.075,right=.945,bottom=.12,top=.89,wspace=.47,hspace=.55)
    fig.suptitle(title,x=.075,ha='left',y=.98,fontsize=15,fontweight='bold')
    fig.text(.075,.943,subtitle,fontsize=9,ha='left')
    return fig,axes

footer='Star: exact original parameters. Diamond: historical (-40, 4 mV). Gray: technically unavailable or no 1/e crossing.\nOne reconstruction; parameter-grid coverage is descriptive, not a biological probability. No parameters were refitted.'

fig,ax=base(3,3,(13.5,11.2),'Calcium activation sweep','25% retained inhibition | 2,601 midpoint / slope combinations | shared fixed voltage responses and 50-ms clearance')
heat(fig,ax[0,0],np.log10(gain[:,1,1]),'A  Intermediate-contrast gain','log10 reduced / control AUC',norm=Normalize(0,np.nanmax(np.log10(gain[:n,1,1]))))
for a,k,label in [(ax[0,1],0,'B  Low-contrast preference'),(ax[0,2],1,'C  Intermediate-contrast preference')]:
    v=np.log2(pref[:,1,k]);heat(fig,a,v,label,f'log2 [gain({[15,30][k]}%) / gain(70%)]',cmap='RdBu_r',norm=centered(v[:n]),contour=0)
dshape=shape[:,2,6]-shape[:,0,6]
heat(fig,ax[1,0],dshape,'D  Normalized 30% shape separation','reduced minus control',cmap='RdBu_r',norm=centered(dshape[:n]),contour=0)
dspread=cross[:,1]-cross[:,0]
heat(fig,ax[1,1],dspread,'E  Focal calcium spread increase','change in first 1/e distance (um)',norm=Normalize(0,np.nanmax(dspread[:n])))
heat(fig,ax[1,2],shift[:,1],'F  Calcium-only contrast equivalence','mean horizontal multiplier',norm=Normalize(np.nanmin(shift[:n,1]),np.nanmax(shift[:n,1])) )
cs=ax[1,2].contour(half,slope,counts[:n,1].reshape(51,51).T,levels=[.5,1.5,2.5,3.5],colors='black',linewidths=.5,linestyles='dotted')
ax[1,2].text(.02,.98,'Dotted: transitions in match count',transform=ax[1,2].transAxes,va='top',fontsize=6.5,bbox=dict(facecolor='white',alpha=.8,edgecolor='none'))
for i,style,label in [(n,'-','Original'),(historical,'--','Historical')]:
    for j,col,name in [(0,GRAY,'Control'),(2,ORANGE,'Reduced')]:
        ax[2,0].plot(contrasts,auc[i,j]/1000,style,color=col,lw=1.5,label=f'{label}: {name}')
ax[2,0].set(title='G  Absolute contrast responses',xlabel='Contrast (%)',ylabel='Ca AUC (10³ a.u. ms²)',ylim=(0,155));ax[2,0].legend(fontsize=7,frameon=False,loc='lower right')
for i,style,label in [(n,'-','Original'),(historical,'--','Historical')]:
    for j,col in [(0,GRAY),(1,ORANGE)]:ax[2,1].plot(x,profiles[i,j],style,color=col,lw=1.5)
ax[2,1].axhline(1/np.e,color='#888888',ls=':',lw=.7)
ax[2,1].set(title='H  Source-normalized focal profiles',xlabel='Cable path distance (um)',ylabel='Relative calcium',ylim=(-.02,1.03),xlim=(0,x[-1]))
valid=np.isfinite(dspread[:n])&np.isfinite(pref[:n,1,1])
h=ax[2,2].scatter(dspread[:n][valid],np.log2(pref[:n,1,1][valid]),c=pairs[:n,0][valid],s=3,cmap='viridis',alpha=.5,rasterized=True)
ax[2,2].scatter(dspread[-1],np.log2(pref[-1,1,1]),marker='*',s=100,color=TEAL,edgecolors='white',zorder=5)
ax[2,2].scatter(dspread[historical],np.log2(pref[historical,1,1]),marker='D',s=25,facecolors='none',edgecolors='black',zorder=5)
ax[2,2].axhline(0,color='black',lw=.6)
ax[2,2].set(title='I  Spread and contrast preference',xlabel='Change in calcium 1/e distance (um)',ylabel='log2 [gain(30%) / gain(70%)]')
cb=fig.colorbar(h,ax=ax[2,2],fraction=.046,pad=.025);cb.set_label('Midpoint (mV)',fontsize=8)
fig.text(.075,.025,footer,fontsize=8)
save(fig,'calcium_sweep_overview')

fig,ax=base(3,3,(13.5,10.5),'Dependence on retained inhibition','Rows: 20%, 25%, and 30% of control inhibitory conductance | identical color scales across rows')
allp15=np.log2(pref[:n,:,0]);allp30=np.log2(pref[:n,:,1]);alls=shape[:n,1:,6]-shape[:n,0,6,None]
for row,f in enumerate([20,25,30]):
    heat(fig,ax[row,0],np.log2(pref[:,row,0]),f'{f}% retained: low-contrast preference','log2 [gain(15%) / gain(70%)]',cmap='RdBu_r',norm=centered(allp15),contour=0)
    heat(fig,ax[row,1],np.log2(pref[:,row,1]),f'{f}% retained: intermediate preference','log2 [gain(30%) / gain(70%)]',cmap='RdBu_r',norm=centered(allp30),contour=0)
    heat(fig,ax[row,2],shape[:,row+1,6]-shape[:,0,6],f'{f}% retained: 30% shape separation','reduced minus control',cmap='RdBu_r',norm=centered(alls),contour=0)
fig.text(.075,.025,footer,fontsize=8);save(fig,'calcium_sweep_residual_sensitivity')

fig,ax=base(3,4,(16,10.5),'Absolute responses and gain','25% retained inhibition | raw calcium proxy responses shown beside ratios | all contrast conditions use the same amplitude scale')
for row,(c,ci) in enumerate([(15,3),(30,6),(70,10)]):
    for col,j,label in [(0,0,'Control'),(1,2,'Reduced')]:
        heat(fig,ax[row,col],auc[:,j,ci]/1000,f'{c}% contrast: {label}','Ca AUC (10³ a.u. ms²)',norm=Normalize(0,150))
    heat(fig,ax[row,2],inc[:,1,row]/1000,f'{c}% contrast: absolute increment','Ca AUC (10³ a.u. ms²)',norm=Normalize(0,150))
    heat(fig,ax[row,3],np.log10(gain[:,1,row]),f'{c}% contrast: fractional gain','log10 reduced / control AUC',norm=Normalize(0,np.nanmax(np.log10(gain[:n,1]))))
fig.text(.075,.025,footer+'\nRatio denominators at or below 0.15 a.u. ms² are masked; all raw responses remain in the data.',fontsize=8);save(fig,'calcium_sweep_absolute_responses')

examples=[(n,'Exact original'),(historical,'Historical reference')]+[(np.flatnonzero((pairs[:,0]==mid)&(pairs[:,1]==k))[0],f'Domain corner {mid} / {k}') for mid,k in [(-45,1),(-45,6),(-20,1),(-20,6)]]
fig,axes=base(2,3,(13.5,8.5),'Original settings and all four domain corners','Complete corner set shown alongside original and historical settings | control (gray), 25% retained inhibition (orange)')
v=np.linspace(-45,0,451)
for j,(i,label) in enumerate(examples):
    ax=axes.flat[j]
    mid,k=pairs[i]
    for ci,color,condition in [(0,GRAY,'Control'),(2,ORANGE,'Reduced')]:ax.plot(contrasts,auc[i,ci]/1000,color=color,lw=1.8,label=condition)
    ax.set(title=f'{label}\nMidpoint {mid:.2f} mV; slope {k:.2f} mV',xlabel='Contrast (%)',ylabel='Ca AUC (10³ a.u. ms²)',ylim=(-2,155))
    labels=[f'{value:.3g}' if np.isfinite(value) else 'undefined' for value in pref[i,1]]
    ax.text(.03,.96,f'P15: {labels[0]}\nP30: {labels[1]}',transform=ax.transAxes,va='top',fontsize=8)
fig.text(.075,.025,'P15 and P30 compare fractional gain with gain at 70% contrast; values >1 indicate preferential enhancement.\nUndefined ratios have technically weak denominators; raw responses are retained. Responses are model proxies, not calibrated fluorescence.',fontsize=8)
save(fig,'calcium_sweep_examples')

# Full long-form exports explicitly record each equivalence exclusion.
with (OUT/'contrast_responses_and_matches.tsv').open('w') as f:
    columns=['point_type','v_half_mV','slope_mV','retained_fraction','contrast_pct','control_auc','reduced_auc','delta_auc','gain','control_shape','reduced_shape','shape_difference','calcium_horizontal_multiplier','match_status']
    w=csv.DictWriter(f,fieldnames=columns,delimiter='\t');w.writeheader()
    for i,p in enumerate(pairs):
        ctrl=auc[i,0]
        for j,r in enumerate([.2,.25,.3]):
            for ci,c in enumerate(contrasts):
                red=auc[i,j+1,ci];m=z['calcium_shift_by_contrast'][i,j,ci]
                if c==0:reason='zero_contrast'
                elif not np.isfinite(ctrl).all():reason='nonfinite_control'
                elif np.any(np.diff(ctrl)<0):reason='nonmonotonic_control'
                elif np.any(np.diff(ctrl)==0):reason='tied_control_responses'
                elif ctrl[-1]-ctrl[0]<=.15:reason='weak_control_span'
                elif red<=.15:reason='weak_reduced_response'
                elif not ctrl[0]<=red<=ctrl[-1]:reason='outside_control_response_range'
                elif np.isfinite(m):reason='valid'
                else:raise RuntimeError('Unexplained invalid match')
                row=dict(point_type='grid' if i<n else 'original_exact',v_half_mV=p[0],slope_mV=p[1],retained_fraction=r,contrast_pct=c,control_auc=ctrl[ci],reduced_auc=red,delta_auc=red-ctrl[ci],gain=red/ctrl[ci] if ctrl[ci]>.15 else np.nan,control_shape=shape[i,0,ci],reduced_shape=shape[i,j+1,ci],shape_difference=shape[i,j+1,ci]-shape[i,0,ci],calcium_horizontal_multiplier=m,match_status=reason)
                w.writerow({k:'' if isinstance(v,(float,np.floating)) and not np.isfinite(v) else v for k,v in row.items()})

maps=np.load(SOURCE/'nominal25/finite_local_maps.npz')
length=maps['cable_lengths']
# Whole-arbor absolute cable integral = source response * normalized integral.
absolute_cable=z['spatial_source_ca_ctrl25']*z['spatial_weighted_cable_ctrl25_um']
metadata=[]
voltages=np.load(SOURCE/'sensitivity_voltage_traces.npz')
max_v=max(float(voltages[k].max()) for k in voltages.files if k!='t_ms')
for i,p in enumerate(pairs):
    meta=dict(point_type='grid' if i<n else 'original_exact',v_half_mV=p[0],slope_mV=p[1],activation_at_rest=float(expit((-45-p[0])/p[1])),activation_at_archived_global_max=float(expit((max_v-p[0])/p[1])),archived_global_max_mV=max_v)
    for j,key in enumerate(['control','reduced25']):
        meta.update({key+'_source_ca':z['spatial_source_ca_ctrl25'][i,j],key+'_absolute_cable_ca':absolute_cable[i,j],key+'_first1e_um':cross[i,j],key+'_first1e_status':'observed' if np.isfinite(cross[i,j]) else 'no_crossing_in_sampled_arbor'})
    for j,key in enumerate(['residual20','residual25','residual30']):
        m=z['calcium_shift_by_contrast'][i,j];good=m[np.isfinite(m)]
        meta.update({key+'_shift_min':good.min() if len(good) else np.nan,key+'_shift_max':good.max() if len(good) else np.nan})
    metadata.append(meta)
with (OUT/'activation_spatial_metadata.tsv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(metadata[0]),delimiter='\t');w.writeheader()
    for row in metadata:w.writerow({k:'' if isinstance(v,(float,np.floating)) and not np.isfinite(v) else v for k,v in row.items()})

# Explicit precision margin for interpreting strict directional counts.
margin={}
for j,r in enumerate([.2,.25,.3]):
    finite=np.isfinite(pref[:n,j]).all(axis=1)
    values=np.abs(pref[:n,j][np.isfinite(pref[:n,j])]-1)
    margin[str(r)]=dict(closest_preference_to_one=float(values.min()),strict_both=int(np.sum(finite&np.all(pref[:n,j]>1,axis=1))),both_beyond_1e_minus10=int(np.sum(finite&np.all(pref[:n,j]>1+1e-10,axis=1))))
(OUT/'directional_precision.json').write_text(json.dumps(margin,indent=2)+'\n')
# Track a common contrast support as a transparent secondary descriptor, not extrapolation.
common={}
for j,r in enumerate([.2,.25,.3]):
    m=z['calcium_shift_by_contrast'][:n,j]
    selected=(contrasts>=5)&(contrasts<=15)
    valid=np.isfinite(m[:,selected]).all(axis=1)
    value=m[valid][:,selected].mean(axis=1)
    common[str(r)]=dict(contrasts=[5,10,15],n_valid=int(valid.sum()),minimum=float(value.min()),median=float(np.median(value)),maximum=float(value.max()))
(OUT/'common_support_shift.json').write_text(json.dumps(common,indent=2)+'\n')
hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in OUT.iterdir() if p.is_file() and p.name!='review_manifest.json'}
(OUT/'review_manifest.json').write_text(json.dumps(dict(source_archive_sha256=hashlib.sha256((RESULTS/'calcium_sweep.npz').read_bytes()).hexdigest(),builder_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),files=hashes),indent=2)+'\n')
print('Saved four review figures (SVG/PNG), long-form data, activation/spatial metadata, and audit metadata.')
