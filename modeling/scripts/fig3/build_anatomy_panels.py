#!/usr/bin/env python3
"""Render actual compartment geometry for the Fig. 3 anatomy panels.

All map scales are shared across conditions. No morphology smoothing, coordinate
changes, per-map normalization, or percentile clipping is used.
"""
from pathlib import Path
import argparse
import json
import hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from matplotlib.patches import Circle

ROOT = Path(__file__).resolve().parents[2]
PROJECT = ROOT.parent
DATA = ROOT/'outputs/fig3_anatomy'
OUT = PROJECT/'figures/main/assets/Fig3_v2'
GRAY = '#333333'
ORANGE = '#d55e00'
CYAN = '#079ba8'
CONDITIONS = [('ctrl','Control',GRAY), ('gaba','GABA removed',ORANGE)]


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--outdir',type=Path,default=OUT)
    args=ap.parse_args(); out=args.outdir;out.mkdir(parents=True,exist_ok=True)
    a=np.load(DATA/'anatomy_maps.npz')
    meta=json.loads((DATA/'metadata.json').read_text())
    xyz=a['xyz_segments']; edges=a['geometry_edges']; idx=a['edge_segment_index']
    # Order crossings consistently by physical depth, retaining actual XY view.
    order=np.argsort(edges[:,:,2].mean(1),kind='stable')
    edges=edges[order]; idx=idx[order]
    minxy=edges[:,:,:2].min((0,1));maxxy=edges[:,:,:2].max((0,1))
    limits=(minxy[0]-5,maxxy[0]+5,minxy[1]-9,maxxy[1]+7)
    source=int(a['local_source_segment']);center=a['stimulus_center_xy']
    plt.rcParams.update({'font.family':'Arial','font.size':7.5,'axes.labelsize':7.5,
        'xtick.labelsize':6.5,'ytick.labelsize':6.5,'axes.linewidth':.5,
        'axes.spines.top':False,'axes.spines.right':False,'xtick.major.width':.5,
        'ytick.major.width':.5,'xtick.major.size':2,'ytick.major.size':2,
        'svg.fonttype':'none','pdf.fonttype':42,'savefig.facecolor':'white'})
    cmap=plt.get_cmap('magma')
    map_extrema={}

    def anatomy(ax,values,vmax,spot=False,source_mark=False,scale=False,lw=.46):
        assert np.isfinite(values).all()
        # Peak depolarizations include the pre-stimulus baseline, so they are >=0.
        # Local signed snapshots may have roundoff at zero; retain the data.
        assert values.min()>-1e-6 and values.max()<=vmax+1e-6,(values.min(),values.max(),vmax)
        lc=LineCollection(edges[:,:,:2],array=values[idx],cmap=cmap,norm=Normalize(0,vmax),
                          linewidths=lw,capstyle='round',joinstyle='round')
        ax.add_collection(lc)
        if spot:
            ax.plot(*center,marker='+',ms=3.5,mew=.55,color='#777777',zorder=4)
        if source_mark:
            sx,sy=xyz[source,:2]
            ax.plot(sx,sy,marker='o',ms=3.6,mfc='white',mec=CYAN,mew=.9,zorder=5)
            ax.annotate('',xy=(sx,sy+1.2),xytext=(sx-3,sy+15),
                        arrowprops=dict(arrowstyle='-|>',color=CYAN,lw=.75,mutation_scale=5))
        if scale:
            x=limits[1]-27;y=limits[2]+8
            ax.plot([x,x+20],[y,y],color=GRAY,lw=1,solid_capstyle='butt')
            ax.text(x+10,y-2.4,'20 µm',ha='center',va='top',fontsize=6.5)
        ax.set(xlim=limits[:2],ylim=limits[2:],aspect='equal')
        ax.set_axis_off()
        return lc

    def label(fig,letter,title,x,y):
        fig.text(x,y,letter,fontsize=10,fontweight='bold',va='top')
        fig.text(x+.032,y-.001,title,fontsize=8.1,va='top')

    def save(fig,name):
        for ext in ['svg','pdf','png']:
            fig.savefig(out/f'{name}.{ext}',dpi=300)
        plt.close(fig)

    fig=plt.figure(figsize=(7,4.65))
    label(fig,'D','Spread of a local input',.014,.985)
    for column,(condition,name,color) in enumerate(CONDITIONS):
        x=.015+column*.285
        ax=fig.add_axes([x,.697,.267,.215])
        values=a[f'local_{condition}_at_source_peak_normalized']
        lc=anatomy(ax,values,1,source_mark=True,scale=column==0,lw=.55)
        fig.text(x+.134,.924,name,color=color,ha='center',fontsize=7.6)
        map_extrema[f'local_{condition}']=[float(values.min()),float(values.max())]
    cb=fig.colorbar(lc,cax=fig.add_axes([.173,.682,.235,.017]),orientation='horizontal',ticks=[0,.5,1])
    cb.outline.set_visible(False);cb.ax.tick_params(length=1.5,pad=1)
    cb.set_label('ΔV / local ΔV',labelpad=1.7,fontsize=7)

    label(fig,'E','Voltage transfer',.644,.985)
    ax=fig.add_axes([.723,.733,.25,.187])
    d=a['path_distance_from_source_um'];length=a['cable_lengths']
    for condition,_,color in CONDITIONS:
        values=a[f'local_{condition}_at_source_peak_normalized']
        bx=[];by=[]
        for low in np.arange(0,110,10):
            mask=(d>=low)&(d<low+10)
            if mask.any():
                bx.append(low+5);by.append(np.average(values[mask],weights=length[mask]))
        ax.plot(bx,by,color=color,lw=1.15)
    ax.set(xlim=(0,110),ylim=(0,1),xticks=[0,50,100],yticks=[0,.5,1],
           xlabel='Path distance (µm)',ylabel='ΔV / local ΔV')
    ax.xaxis.labelpad=2;ax.yaxis.labelpad=2

    label(fig,'F','100-µm OFF spot',.014,.594)
    sax=fig.add_axes([.571,.55,.06,.062])
    sax.add_patch(Circle((0,0),1,fill=False,ec='#999999',lw=.65,ls=(0,(2,2))))
    sax.plot(0,0,marker='+',color='#777777',ms=3,mew=.55)
    sax.set(xlim=(-1.2,1.2),ylim=(-1.2,1.2),aspect='equal');sax.set_axis_off()
    volts=[a[f'spot_{condition}_{contrast}_dv'] for condition,_,_ in CONDITIONS for contrast in [15,30,70]]
    vmax=float(np.ceil(max(v.max() for v in volts)/5)*5)
    for j,c in enumerate([15,30,70]):
        fig.text(.142+j*.208,.536,f'{c}% contrast',ha='center',fontsize=7.2)
    for row,(condition,name,color) in enumerate(CONDITIONS):
        ypos=.335-row*.175
        fig.text(.016,ypos+.08,name,rotation=90,ha='center',va='center',color=color,fontsize=7.2)
        for col,c in enumerate([15,30,70]):
            ax=fig.add_axes([.043+col*.208,ypos,.199,.17])
            values=a[f'spot_{condition}_{c}_dv']
            lc=anatomy(ax,values,vmax,spot=True,scale=row==1 and col==0,lw=.42)
            map_extrema[f'spot_{condition}_{c}']=[float(values.min()),float(values.max())]
    cb=fig.colorbar(lc,cax=fig.add_axes([.178,.132,.29,.017]),orientation='horizontal',ticks=np.arange(0,vmax+1,5))
    cb.outline.set_visible(False);cb.ax.tick_params(length=1.5,pad=1)
    cb.set_label('Peak ΔV (mV)',fontsize=7,labelpad=1.7)

    rows=meta['spot_summary']
    label(fig,'G','Voltage',.70,.594)
    axv=fig.add_axes([.772,.391,.205,.145])
    label(fig,'H','Calcium',.70,.32)
    axc=fig.add_axes([.772,.118,.205,.145])
    for condition,_,color in CONDITIONS:
        rr=sorted([r for r in rows if r['condition']==condition],key=lambda r:r['contrast_pct'])
        c=[0]+[r['contrast_pct'] for r in rr]
        dv=[0]+[r['mean_off_peak_dv_mV'] for r in rr]
        ca=[0]+[r['mean_off_peak_ca_au_ms'] for r in rr]
        axv.plot(c,dv,'-o',color=color,ms=2.3,lw=1)
        axc.plot(c,ca,'-o',color=color,ms=2.3,lw=1)
    axv.set(xlim=(0,75),ylim=(0,22),xticks=[0,15,30,70],yticks=[0,10,20],ylabel='Peak ΔV (mV)')
    axc.set(xlim=(0,75),ylim=(0,40),xticks=[0,15,30,70],yticks=[0,20,40],ylabel='Ca proxy (a.u.)',xlabel='Contrast (%)')
    for ax in [axv,axc]:ax.yaxis.labelpad=2;ax.xaxis.labelpad=2
    # Small explicit transfer curve below maps, instead of hiding the mechanism
    # inside an anatomy color scale. This is a parameterized model, not data.
    iax=axc.inset_axes([.57,.10,.40,.24])
    v=np.linspace(-55,-10,200);act=1/(1+np.exp(-(v+40)/4))
    iax.plot(v,act,color=GRAY,lw=.85)
    iax.set(xlim=(-55,-10),ylim=(0,1.05),xticks=[],yticks=[])
    iax.text(-57,1.2,'Activation',fontsize=5.3,ha='left')
    iax.set_xlabel('V',fontsize=5.2,labelpad=.5)
    fig.text(.019,.038,'Unfitted model: GABA removed; 30% glycinergic conductance retained.\nG, H: mean local peaks; calcium activation + linear clearance; shared scales.',fontsize=6.2,color='#555555',va='bottom',linespacing=1.35)
    save(fig,'Fig3_modeling_panels')

    # Full-arbor calcium maps as a companion, using the same scientific data
    # and physical rendering, with one common arbitrary-unit scale.
    fig=plt.figure(figsize=(7,4.05))
    fig.text(.02,.965,'100-µm spot: saturating channel activation',fontsize=10,fontweight='bold',va='top')
    cmax=float(np.ceil(max(a[f'spot_{condition}_{c}_ca'].max() for condition,_,_ in CONDITIONS for c in [15,30,70])/5)*5)
    for col,c in enumerate([15,30,70]):fig.text(.22+col*.28,.862,f'{c}% contrast',ha='center',fontsize=8)
    for row,(condition,name,color) in enumerate(CONDITIONS):
        y=.54-row*.30
        fig.text(.022,y+.12,name,rotation=90,ha='center',va='center',color=color,fontsize=8)
        for col,c in enumerate([15,30,70]):
            ax=fig.add_axes([.08+col*.28,y,.27,.26])
            lc=anatomy(ax,a[f'spot_{condition}_{c}_ca'],cmax,spot=True,scale=row==1 and col==0,lw=.6)
    cb=fig.colorbar(lc,cax=fig.add_axes([.31,.18,.4,.02]),orientation='horizontal',ticks=np.arange(0,cmax+1,10))
    cb.outline.set_visible(False);cb.set_label('Peak calcium proxy (a.u.); one shared scale',fontsize=8)
    fig.text(.5,.045,'Activation midpoint −40 mV, slope 4 mV; clearance 50 ms.\nNo driving-force or indicator saturation; illustrative parameters, not a fit to imaging.',ha='center',fontsize=7,color='#555555')
    save(fig,'Fig3_calcium_anatomy')
    (out/'render_manifest.json').write_text(json.dumps(dict(
        geometry='Actual pt3d polyline edges; planar XY projection; same camera, limits, and linewidth per map group',
        edge_count=len(edges),map_extrema=map_extrema,voltage_color_scale_mV=[0,vmax],
        local_color_scale=[0,1],calcium_color_scale_au=[0,cmax],
        calcium_readout='activation only + linear clearance, unchanged across conditions',
        source_hash=hashlib.sha256((DATA/'anatomy_maps.npz').read_bytes()).hexdigest(),
        local_protocol='Snapshot at source peak for each condition, normalized by simultaneous source response',
        axes_limits_um=limits),indent=2))
    print(out)


if __name__=='__main__':main()
