#!/usr/bin/env python3
"""Detailed analytical exports and individual calcium-sweep panels."""
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
    fig.savefig(OUT/f'{name}.svg',facecolor='white',bbox_inches='tight')
    fig.savefig(OUT/f'{name}.png',facecolor='white',dpi=180,bbox_inches='tight')
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


def panel():
    return plt.subplots(figsize=(5.4, 4.2), layout='constrained')


def heat_panel(data, title, label, **kwargs):
    fig, ax = panel()
    heat(fig, ax, data, title, label, **kwargs)
    ax.legend(handles=[
        Line2D([], [], ls='', marker='*', ms=8, mfc=TEAL, mec='white',
               label='Original parameters'),
        Line2D([], [], ls='', marker='D', ms=4, mfc='none', mec='black',
               label='Historical (-40, 4 mV)')],
        loc='upper center', bbox_to_anchor=(.5, -.23), ncol=2,
        fontsize=7, frameon=False)
    return fig, ax


# Each heatmap or curve is saved separately; color ranges remain shared.
fig, ax = heat_panel(np.log10(gain[:,1,1]),
    'Intermediate-contrast gain; 25% retained', 'log10 reduced / control AUC',
    norm=Normalize(0,np.nanmax(np.log10(gain[:n,1,1]))))
save(fig, 'calcium_gain30')
for k, label in [(0, 'low'), (1, 'intermediate')]:
    v=np.log2(pref[:,1,k])
    fig, ax = heat_panel(v, f'{label.title()}-contrast preference; 25% retained',
        f'log2 [gain({[15,30][k]}%) / gain(70%)]', cmap='RdBu_r', norm=centered(v[:n]), contour=0)
    save(fig, f'calcium_preference_{label}')
dshape=shape[:,2,6]-shape[:,0,6]
fig, ax = heat_panel(dshape, 'Normalized 30% shape separation; 25% retained',
    'reduced minus control', cmap='RdBu_r', norm=centered(dshape[:n]), contour=0)
save(fig, 'calcium_shape30_difference')
dspread=cross[:,1]-cross[:,0]
fig, ax = heat_panel(dspread, 'Focal calcium spread increase; 25% retained',
    'change in first 1/e distance (µm)', norm=Normalize(0,np.nanmax(dspread[:n])))
save(fig, 'calcium_spread_increase')
fig, ax = heat_panel(shift[:,1], 'Calcium-only contrast equivalence; 25% retained',
    'mean horizontal multiplier', norm=Normalize(np.nanmin(shift[:n,1]),np.nanmax(shift[:n,1])))
ax.contour(half,slope,counts[:n,1].reshape(51,51).T,levels=[.5,1.5,2.5,3.5],colors='black',linewidths=.5,linestyles='dotted')
ax.text(.02,.98,'Dotted: transitions in match count',transform=ax.transAxes,va='top',fontsize=6.5,bbox=dict(facecolor='white',alpha=.8,edgecolor='none'))
save(fig, 'calcium_contrast_equivalence')
fig, ax = panel()
for i,style,label in [(n,'-','Original'),(historical,'--','Historical')]:
    for j,col,name in [(0,GRAY,'Control'),(2,ORANGE,'Reduced')]:
        ax.plot(contrasts,auc[i,j]/1000,style,color=col,lw=1.5,label=f'{label}: {name}')
ax.set(title='Absolute contrast responses; 25% retained',xlabel='Contrast (%)',ylabel='Ca AUC (10³ a.u. ms²)',ylim=(0,155))
ax.legend(fontsize=7,frameon=False,loc='lower right')
save(fig, 'calcium_contrast_responses')
fig, ax = panel()
for i,style,label in [(n,'-','Original'),(historical,'--','Historical')]:
    for j,col,name in [(0,GRAY,'Control'),(1,ORANGE,'Reduced')]:
        ax.plot(x,profiles[i,j],style,color=col,lw=1.5,label=f'{label}: {name}')
ax.axhline(1/np.e,color='#888888',ls=':',lw=.7)
ax.set(title='Source-normalized focal profiles; 25% retained',xlabel='Cable path distance (µm)',ylabel='Relative calcium',ylim=(-.02,1.03),xlim=(0,x[-1]))
ax.legend(fontsize=7,frameon=False)
save(fig, 'calcium_focal_profiles')
fig, ax = panel()
valid=np.isfinite(dspread[:n])&np.isfinite(pref[:n,1,1])
h=ax.scatter(dspread[:n][valid],np.log2(pref[:n,1,1][valid]),c=pairs[:n,0][valid],s=3,cmap='viridis',alpha=.5,rasterized=True)
ax.scatter(dspread[-1],np.log2(pref[-1,1,1]),marker='*',s=100,color=TEAL,edgecolors='white',zorder=5,label='Original')
ax.scatter(dspread[historical],np.log2(pref[historical,1,1]),marker='D',s=25,facecolors='none',edgecolors='black',zorder=5,label='Historical')
ax.axhline(0,color='black',lw=.6)
ax.set(title='Spread and contrast preference; 25% retained',xlabel='Change in calcium 1/e distance (µm)',ylabel='log2 [gain(30%) / gain(70%)]')
fig.colorbar(h,ax=ax,fraction=.046,pad=.025,label='Midpoint (mV)')
ax.legend(fontsize=7,frameon=False)
save(fig, 'calcium_spread_preference')

allp15=np.log2(pref[:n,:,0]);allp30=np.log2(pref[:n,:,1]);alls=shape[:n,1:,6]-shape[:n,0,6,None]
for row,f in enumerate([20,25,30]):
    items=[('preference15',np.log2(pref[:,row,0]),'Low-contrast preference','log2 [gain(15%) / gain(70%)]',centered(allp15)),
           ('preference30',np.log2(pref[:,row,1]),'Intermediate preference','log2 [gain(30%) / gain(70%)]',centered(allp30)),
           ('shape30',shape[:,row+1,6]-shape[:,0,6],'30% shape separation','reduced minus control',centered(alls))]
    for name,values,title,label,norm in items:
        fig, ax = heat_panel(values,f'{title}; {f}% retained',label,cmap='RdBu_r',norm=norm,contour=0)
        save(fig,f'calcium_retained{f}_{name}')

for row,(c,ci) in enumerate([(15,3),(30,6),(70,10)]):
    for j,label in [(0,'control'),(2,'reduced')]:
        fig,ax=heat_panel(auc[:,j,ci]/1000,f'{c}% contrast: {label}; 25% retained','Ca AUC (10³ a.u. ms²)',norm=Normalize(0,150))
        save(fig,f'calcium_contrast{c}_{label}')
    fig,ax=heat_panel(inc[:,1,row]/1000,f'{c}% contrast: absolute increment; 25% retained','Ca AUC (10³ a.u. ms²)',norm=Normalize(0,150))
    save(fig,f'calcium_contrast{c}_increment')
    fig,ax=heat_panel(np.log10(gain[:,1,row]),f'{c}% contrast: fractional gain; 25% retained','log10 reduced / control AUC',norm=Normalize(0,np.nanmax(np.log10(gain[:n,1]))))
    save(fig,f'calcium_contrast{c}_gain')

examples=[(n,'Exact original'),(historical,'Historical reference')]+[(np.flatnonzero((pairs[:,0]==mid)&(pairs[:,1]==k))[0],f'Domain corner {mid} / {k}') for mid,k in [(-45,1),(-45,6),(-20,1),(-20,6)]]
for j,(i,label) in enumerate(examples):
    fig,ax=panel()
    mid,k=pairs[i]
    for ci,color,condition in [(0,GRAY,'Control'),(2,ORANGE,'25% retained')]:
        ax.plot(contrasts,auc[i,ci]/1000,color=color,lw=1.8,label=condition)
    ax.set(title=f'{label}\nMidpoint {mid:.2f} mV; slope {k:.2f} mV',xlabel='Contrast (%)',ylabel='Ca AUC (10³ a.u. ms²)',ylim=(-2,155))
    labels=[f'{value:.3g}' if np.isfinite(value) else 'undefined' for value in pref[i,1]]
    ax.text(.03,.96,f'P15: {labels[0]}\nP30: {labels[1]}',transform=ax.transAxes,va='top',fontsize=8)
    ax.legend(fontsize=7,frameon=False,loc='lower right')
    save(fig,f'calcium_example_{j+1}')

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
(OUT/'README.md').write_text(
    '# Calcium sweep panels\n\n'
    'These 36 SVG/PNG plots are exported individually from the archived sweep. '
    'Related heatmaps keep shared color scales. Stars identify the exact original '
    'parameters; diamonds identify the historical (-40, 4 mV) setting. Gray cells '
    'mark unavailable metrics or unobserved first 1/e crossings.\n\n'
    'P15 and P30 compare fractional calcium gain at 15% and 30% contrast with '
    'gain at 70% contrast. Values above one indicate preferential enhancement. '
    'Ratio denominators at or below 0.15 a.u. ms² are masked; all raw responses '
    'remain in the data. Horizontal equivalence uses only supported matches; '
    'dotted boundaries mark changes in match count.\n\n'
    'These are proxy responses from one reconstruction. Parameter-grid coverage '
    'is descriptive, not a biological probability. No parameters were refitted. '
    'The accompanying TSV/JSON exports retain full responses, match exclusions, '
    'activation and spatial metadata, and numerical precision checks.\n')
hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in OUT.iterdir() if p.is_file() and p.name!='review_manifest.json'}
(OUT/'review_manifest.json').write_text(json.dumps(dict(source_archive_sha256=hashlib.sha256((RESULTS/'calcium_sweep.npz').read_bytes()).hexdigest(),builder_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),files=hashes),indent=2)+'\n')
print('Saved 36 individual panels (SVG/PNG), long-form data, activation/spatial metadata, and audit metadata.')
