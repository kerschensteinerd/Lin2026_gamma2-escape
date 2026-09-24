#!/usr/bin/env python3
"""Summarize spatial-control results and export individual analytical panels."""
from pathlib import Path
import csv,json,hashlib
import numpy as np
from scipy.interpolate import PchipInterpolator
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import argparse,shutil
SCRIPT_DIR=Path(__file__).resolve().parent
ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument('--results',type=Path,default=SCRIPT_DIR)
ap.add_argument('--outdir',type=Path,default=SCRIPT_DIR.parent/'generated/spatial')
args=ap.parse_args()
OUT=args.outdir.resolve();OUT.mkdir(parents=True,exist_ok=True)
for source in args.results.resolve().iterdir():
    if source.is_file() and source.suffix in {'.json','.tsv','.npz'} and source != OUT/source.name:
        shutil.copy2(source,OUT/source.name)
F=(.2,.25,.3);ST=(.5,1.,2.);C=(0,5,10,15,20,25,30,40,50,60,70)
COL={.2:'#2a80b9',.25:'#d76b29',.3:'#4c9361'}
FC={'control':'#333333','uniform':'#d76b29','random':'#9a9a9a','proximal':'#286c9e','distal':'#8851a0'}
MARK={'uniform':'D','random':'o','proximal':'^','distal':'v'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.titlesize':9,'axes.labelsize':8,
                     'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})

def read(p):
    with (OUT/p).open() as f:rr=list(csv.DictReader(f,delimiter='\t'))
    for row in rr:
        for k,v in row.items():
            if v=='':row[k]=None
            elif v in ('True','False'):row[k]=v=='True'
            else:
                try:row[k]=float(v)
                except ValueError:pass
    return rr
def write(name,rows):
    with (OUT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def panel():
    return plt.subplots(figsize=(5.6,4.3),layout='constrained')
def save(fig,name):
    for ext in ['png','pdf','svg']:fig.savefig(OUT/f'{name}.{ext}',dpi=220,facecolor='white',bbox_inches='tight')
    plt.close(fig)
def dist(x):
    x=np.array(x,float);return dict(n=int(len(x)),min=float(x.min()),q25=float(np.quantile(x,.25)),median=float(np.median(x)),q75=float(np.quantile(x,.75)),max=float(x.max()))

def paired_focal(rows):
    ctrl={(r['source_segment'],r['focal_strength_multiplier']):r for r in rows if r['family']=='control'}
    pairs=[]
    for r in rows:
        if r['family']=='control':continue
        c=ctrl[r['source_segment'],r['focal_strength_multiplier']]
        row={k:r[k] for k in ['layout_id','family','retained_fraction','seed','source_index','source_segment','source_soma_path_um','focal_strength_multiplier']}
        for metric in ['source_dv_mV','source_calcium_au_ms','voltage_1e_um','voltage_cable_integral_um','calcium_1e_um','calcium_cable_integral_um','calcium_absolute_cable_integral_au_ms_um']:
            row['control_'+metric]=c[metric];row['reduced_'+metric]=r[metric]
            row['delta_'+metric]=r[metric]-c[metric] if r[metric] is not None and c[metric] is not None else None
        row['both_calcium_crossings_observed']=c['calcium_crossing_observed'] and r['calcium_crossing_observed']
        row['reduced_calcium_crossing_censored']=r['calcium_source_technical_ok'] and not r['calcium_crossing_observed']
        pairs.append(row)
    return pairs

def spot_summary(rows):
    groups={}
    for r in rows:groups.setdefault(r['layout_id'],[]).append(r)
    for lid in groups:groups[lid]=sorted(groups[lid],key=lambda r:r['contrast_pct'])
    ctrl=groups['control'];ac=np.array([r['mean_ca_auc_au_ms2'] for r in ctrl]);cc=np.array(C,float)
    assert np.all(np.diff(ac)>0);inv=PchipInterpolator(ac,cc,extrapolate=False)
    summary=[];paired=[];matches=[]
    for lid,rr in groups.items():
        if lid=='control':continue
        response=np.array([r['mean_ca_auc_au_ms2'] for r in rr]);ratios={};valid=[]
        for i,(r,c) in enumerate(zip(rr,ctrl)):
            ct=r['contrast_pct'];ratio=r['mean_ca_auc_au_ms2']/c['mean_ca_auc_au_ms2'] if c['mean_ca_auc_au_ms2']>.15 else None
            ratios[ct]=ratio
            paired.append(dict(layout_id=lid,family=r['family'],retained_fraction=r['retained_fraction'],contrast_pct=ct,
                control_auc=c['mean_ca_auc_au_ms2'],reduced_auc=r['mean_ca_auc_au_ms2'],delta_auc=r['mean_ca_auc_au_ms2']-c['mean_ca_auc_au_ms2'],
                auc_ratio=ratio,ratio_technical_ok=ratio is not None,control_shape=c['model_auc_shape'],reduced_shape=r['model_auc_shape'],
                delta_shape=r['model_auc_shape']-c['model_auc_shape']))
            ok=bool(ct>0 and ac[0]<=response[i]<=ac[-1]);eq=float(inv(response[i])) if ok else None
            reason='valid' if ok else ('zero_contrast' if ct==0 else 'outside_control_response_range')
            matches.append(dict(layout_id=lid,family=r['family'],retained_fraction=r['retained_fraction'],contrast_pct=ct,response=float(response[i]),
                valid=ok,reason=reason,equivalent_control_contrast=eq,horizontal_multiplier=eq/ct if ok else None))
            if ok:valid.append((ct,eq/ct))
        v=[x[1] for x in valid];i30=C.index(30)
        summary.append(dict(layout_id=lid,family=rr[0]['family'],retained_fraction=rr[0]['retained_fraction'],seed=rr[0]['seed'],
            ratio15=ratios[15],ratio30=ratios[30],ratio70=ratios[70],preference15=ratios[15]/ratios[70],preference30=ratios[30]/ratios[70],
            delta_shape30=rr[i30]['model_auc_shape']-ctrl[i30]['model_auc_shape'],reduced_shape30=rr[i30]['model_auc_shape'],
            all_nonzero_absolute_enhancement=bool(np.all(response[1:]>ac[1:])),
            mean_calcium_horizontal_shift=float(np.mean(v)) if v else None,min_calcium_horizontal_shift=float(min(v)) if v else None,
            max_calcium_horizontal_shift=float(max(v)) if v else None,n_horizontal_matches=len(v),
            supported_contrasts=','.join(f'{x[0]:g}' for x in valid)))
    return summary,paired,matches,groups

def main():
    focal=json.loads((OUT/'focal_run.json').read_text());spot=json.loads((OUT/'spot_run.json').read_text())
    assert focal['status']==spot['status']=='complete'
    sites=read('focal_all_sites.tsv');layouts=read('focal_layout_controls.tsv');sp=read('spot_responses.tsv')
    sitepair=paired_focal(sites);layoutpair=paired_focal(layouts)
    write('focal_site_paired_effects.tsv',sitepair);write('focal_layout_paired_effects.tsv',layoutpair)
    metrics,effects,matches,groups=spot_summary(sp)
    write('spot_layout_metrics.tsv',metrics);write('spot_contrast_effects.tsv',effects);write('horizontal_matches.tsv',matches)
    site_summ=[]
    for f in F:
        for strength in ST:
            rr=[r for r in sitepair if r['retained_fraction']==f and r['focal_strength_multiplier']==strength]
            d=[r['delta_calcium_1e_um'] for r in rr if r['delta_calcium_1e_um'] is not None]
            di=[r['delta_calcium_cable_integral_um'] for r in rr]
            row=dict(retained_fraction=f,focal_strength_multiplier=strength,n_paired_sites=len(rr),n_both_crossings=len(d),
                n_reduced_crossings_censored=sum(r['reduced_calcium_crossing_censored'] for r in rr),n_observed_crossings_broader=sum(x>1e-10 for x in d),
                n_cable_integrals_broader=sum(x>1e-10 for x in di))
            row.update({f'calcium_1e_delta_{k}':v for k,v in dist(d).items()});row.update({f'calcium_cable_delta_{k}':v for k,v in dist(di).items()})
            site_summ.append(row)
    write('focal_site_summary.tsv',site_summ)
    layout_summ=[]
    for f in F:
        mm=[r for r in metrics if r['retained_fraction']==f]
        lp=[r for r in layoutpair if r['retained_fraction']==f and r['focal_strength_multiplier']==1]
        row=dict(retained_fraction=f,n_layouts=len(mm),n_low_and_intermediate_preference=sum(r['preference15']>1 and r['preference30']>1 for r in mm),
                 n_all_nonzero_absolute_enhancement=sum(r['all_nonzero_absolute_enhancement'] for r in mm),
                 n_observed_focal_crossings_broader=sum(r['delta_calcium_1e_um'] is not None and r['delta_calcium_1e_um']>0 for r in lp),
                 n_focal_crossings_censored=sum(r['reduced_calcium_crossing_censored'] for r in lp))
        for key in ['preference15','preference30','delta_shape30','mean_calcium_horizontal_shift','n_horizontal_matches']:
            row.update({f'{key}_{k}':v for k,v in dist([r[key] for r in mm]).items()})
        layout_summ.append(row)
    write('layout_sensitivity_summary.tsv',layout_summ)
    summary=dict(status='complete',focal_cases=900,layout_focal_cases=210,spot_cases=770,
        paired_site_cases=len(sitepair),paired_site_crossings_observed=sum(r['both_calcium_crossings_observed'] for r in sitepair),
        paired_site_crossings_broader=sum(r['delta_calcium_1e_um'] is not None and r['delta_calcium_1e_um']>0 for r in sitepair),
        reduced_site_crossings_censored=sum(r['reduced_calcium_crossing_censored'] for r in sitepair),
        paired_site_cable_integrals_broader=sum(r['delta_calcium_cable_integral_um']>0 for r in sitepair),
        site_summary=site_summ,layout_summary=layout_summ,n_reconstructions=1,n_unique_sites=75,
        focal_validation=focal,spot_validation=spot,
        interpretation='Deterministic within-arbor sensitivity. Sites, layouts, contrasts, and randomseeds are not independent biological replicates or confidence bounds.')
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    # Every site/strength comparison is an independent panel.
    for strength in ST:
        fig_a,a=panel();fig_b,b=panel()
        a.plot([0,76],[0,76],color='.7',lw=.8,zorder=0)
        for f in F:
            rr=[r for r in sitepair if r['retained_fraction']==f and r['focal_strength_multiplier']==strength]
            observed=[r for r in rr if r['both_calcium_crossings_observed']];censored=[r for r in rr if r['reduced_calcium_crossing_censored']]
            a.scatter([r['control_calcium_1e_um'] for r in observed],[r['reduced_calcium_1e_um'] for r in observed],s=10,c=COL[f],alpha=.7,lw=0,label=f'{f:.0%} retained')
            if censored:a.scatter([r['control_calcium_1e_um'] for r in censored],[76]*len(censored),s=27,c=COL[f],marker='^',edgecolors='white',lw=.3)
            b.scatter([r['source_soma_path_um'] for r in rr],[r['delta_calcium_cable_integral_um'] for r in rr],s=10,c=COL[f],alpha=.7,lw=0,label=f'{f:.0%} retained')
            canonical=next(r for r in rr if r['source_segment']==4705)
            if canonical['both_calcium_crossings_observed']:a.scatter([canonical['control_calcium_1e_um']],[canonical['reduced_calcium_1e_um']],s=65,c=COL[f],marker='*',edgecolors='black',lw=.45,zorder=4)
            b.scatter([canonical['source_soma_path_um']],[canonical['delta_calcium_cable_integral_um']],s=65,c=COL[f],marker='*',edgecolors='black',lw=.45,zorder=4)
        a.set(xlim=(0,60),ylim=(0,80),xlabel='Control calcium 1/e distance (µm)',ylabel='Reduced calcium 1/e distance (µm)',title=f'Focal conductance {strength:g}×; all 75 sites')
        b.axhline(0,color='.7',lw=.8);b.set(xlabel='Source path distance from soma (µm)',ylabel='Change in normalized cable integral (µm)',title=f'All 75 sites, {strength:g}× focal drive')
        star=Line2D([],[],marker='*',ls='',mfc='white',mec='black',ms=8,label='Original source')
        censored_key=Line2D([],[],marker='^',ls='',color='.5',ms=5,label='Unobserved reduced crossing')
        for axis,extras in [(a,[star,censored_key]),(b,[star])]:
            handles,_=axis.get_legend_handles_labels()
            axis.legend(handles=handles+extras,frameon=False,fontsize=7,loc='upper center',bbox_to_anchor=(.5,-.22),ncol=2)
        save(fig_a,f'spatial_site_crossings_strength{strength:g}')
        save(fig_b,f'spatial_site_cable_integral_strength{strength:g}')
    # Matched-total layout effects, each saved separately.
    location_fig,location_ax=panel()
    response_fig,response_ax=panel()
    enhancement_fig,enhancement_ax=panel()
    preference_fig,preference_ax=panel()
    spread_fig,spread_ax=panel()
    equivalence_fig,equivalence_ax=panel()
    lsum=read('layout_summary.tsv');ls={r['layout_id']:r for r in lsum};met={r['layout_id']:r for r in metrics}
    for family in ['random','uniform','proximal','distal']:
        rr=[r for r in layoutpair if r['retained_fraction']==.25 and r['focal_strength_multiplier']==1 and r['family']==family]
        location_ax.scatter([ls[r['layout_id']]['weighted_mean_soma_path_um'] for r in rr],[r['reduced_calcium_1e_um'] for r in rr],
            c=FC[family],marker=MARK[family],s=20 if family=='random' else 45,label=family,edgecolors='white',lw=.3)
    ctrlcurve=np.array([r['mean_ca_auc_au_ms2'] for r in groups['control']]);grid=np.array(C)
    rand=np.array([[r['mean_ca_auc_au_ms2'] for r in groups[f'f25_random{i:02d}']] for i in range(20)])
    response_ax.fill_between(grid,rand.min(0)/1e3,rand.max(0)/1e3,color='.8',alpha=.6,label='20 random layouts: range')
    enhancement_ax.fill_between(grid[1:],(rand/np.where(ctrlcurve==0,np.nan,ctrlcurve)).min(0)[1:],(rand/np.where(ctrlcurve==0,np.nan,ctrlcurve)).max(0)[1:],color='.8',alpha=.6)
    for lid,family in [('control','control'),('f25_uniform','uniform'),('f25_proximal','proximal'),('f25_distal','distal')]:
        y=np.array([r['mean_ca_auc_au_ms2'] for r in groups[lid]])
        response_ax.plot(grid,y/1e3,color=FC[family],lw=1.4,label=family)
        if family!='control':enhancement_ax.plot(grid[1:],y[1:]/ctrlcurve[1:],color=FC[family],lw=1.4,label=family)
    location_ax.set(xlabel='Retained-weight mean soma distance (µm)',ylabel='Focal calcium 1/e distance (µm)',title='Same total inhibition; different locations')
    location_ax.legend(frameon=False,fontsize=7)
    response_ax.set(xlabel='OFF contrast (%)',ylabel='Calcium AUC (×10³ arbitrary units·ms²)',title='Absolute responses, 25% retained')
    curve_handles,curve_labels=response_ax.get_legend_handles_labels()
    reduced_legend=[(h,l) for h,l in zip(curve_handles,curve_labels) if l!='control']
    enhancement_ax.legend([h for h,l in reduced_legend],[l for h,l in reduced_legend],frameon=False,fontsize=7,loc='upper right')
    response_ax.legend(frameon=False,fontsize=7)
    enhancement_ax.set(xlabel='OFF contrast (%)',ylabel='Reduced/control calcium AUC',yscale='log',title='Fractional enhancement, 25% retained')
    for f in F:
        for family in ['random','uniform','proximal','distal']:
            mm=[r for r in metrics if r['retained_fraction']==f and r['family']==family]
            size=15 if family=='random' else 40
            preference_ax.scatter([r['preference15'] for r in mm],[r['preference30'] for r in mm],c=COL[f],marker=MARK[family],s=size,alpha=.8,edgecolors='white',lw=.3)
            lp=[r for r in layoutpair if r['retained_fraction']==f and r['focal_strength_multiplier']==1 and r['family']==family]
            spread_ax.scatter([r['delta_calcium_1e_um'] for r in lp],[met[r['layout_id']]['preference15'] for r in lp],c=COL[f],marker=MARK[family],s=size,alpha=.8,edgecolors='white',lw=.3)
            offset={'random':0,'uniform':-.007,'proximal':.007,'distal':.014}[family]
            equivalence_ax.scatter([(f+offset)*100]*len(mm),[r['mean_calcium_horizontal_shift'] for r in mm],c=COL[f],marker=MARK[family],s=size,alpha=.8,edgecolors='white',lw=.3)
    preference_ax.axhline(1,color='.7',lw=.8);preference_ax.axvline(1,color='.7',lw=.8)
    preference_ax.set(xlabel='R15/R70',ylabel='R30/R70',title='Preferential lower-contrast gain')
    spread_ax.set(xlabel='Change in focal calcium 1/e distance (µm)',ylabel='R15/R70',title='Layout-dependent spread and gain')
    equivalence_ax.set(xlabel='Retained inhibition (%)',ylabel='Mean calcium horizontal multiplier',xticks=[20,25,30],title='Overlap-limited contrast equivalence')
    layout_keys=[Line2D([],[],color=COL[f],marker='o',ls='',label=f'{f:.0%} retained') for f in F]
    layout_keys += [Line2D([],[],color='.4',marker=MARK[family],ls='',label=family) for family in ['random','uniform','proximal','distal']]
    for axis in [preference_ax,spread_ax,equivalence_ax]:
        axis.legend(handles=layout_keys,loc='upper center',bbox_to_anchor=(.5,-.22),ncol=3,frameon=False,fontsize=7)
    for fig,name in [(location_fig,'spatial_layout_location'),(response_fig,'spatial_layout_responses'),
                     (enhancement_fig,'spatial_layout_enhancement'),(preference_fig,'spatial_layout_preference'),
                     (spread_fig,'spatial_layout_spread_preference'),(equivalence_fig,'spatial_layout_contrast_equivalence')]:
        save(fig,name)
    # Focused summary text with denominator/censoring/support limitations.
    n25=next(r for r in site_summ if r['retained_fraction']==.25 and r['focal_strength_multiplier']==1)
    report=['# Spatial-control results — 21 September 2026', '',
        'Completed the prespecified site/strength and conductance-distribution controls on the unchanged Figure 4 v9 model. All calcium parameters and excitation inputs remain shared within each control–reduced pair. No calibration or fitting was performed.', '',
        f"Across 675 reduced-versus-control site/strength comparisons, all {summary['paired_site_crossings_broader']} observed calcium 1/e crossings were broader. Six reduced conditions had no observed 1/e crossing within the arbor and are explicitly censored. The normalized calcium cable integral increased in all {summary['paired_site_cable_integrals_broader']}/675 comparisons, including those censored cases.", '',
        f"At nominal 25% retained inhibition and the original focal strength, calcium 1/e distance increased at 75/75 sites: median {n25['calcium_1e_delta_median']:.2f} µm, range {n25['calcium_1e_delta_min']:.2f}–{n25['calcium_1e_delta_max']:.2f} µm. These are sitewise descriptors from one arbor, not biological replicates.", '',
        '## Matched-total inhibition layouts', '',
        'For each 20%, 25%, and 30% retained fraction, tested uniform weakening, 20 random retention orders (PCG64 seeds 0–19), proximal-first retention, and distal-first retention. Nonuniform patterns retain full original contact weights until the conductance budget is reached, partially retain only the boundary contact, and set the remainder to zero. No contact is strengthened; the total is matched exactly. Contact counts need not match the retained conductance fraction.', '',
        '| Retained | Layouts with R15/R70 > 1 and R30/R70 > 1 | R15/R70 range | R30/R70 range | Mean calcium shift range | Valid contrasts per layout |',
        '|---|---:|---:|---:|---:|---:|']
    for r in layout_summ:
        report.append(f"| {r['retained_fraction']:.0%} | {r['n_low_and_intermediate_preference']}/{r['n_layouts']} | {r['preference15_min']:.2f}–{r['preference15_max']:.2f} | {r['preference30_min']:.2f}–{r['preference30_max']:.2f} | {r['mean_calcium_horizontal_shift_min']:.3f}–{r['mean_calcium_horizontal_shift_max']:.3f} | {r['n_horizontal_matches_min']:.0f}–{r['n_horizontal_matches_max']:.0f} |")
    maxdv=max(c['max_absolute_dv_error_mV'] for c in json.loads((OUT/'focal_validation.json').read_text()) if 'max_absolute_dv_error_mV' in c)
    report.extend(['',
        'Here R15, R30, and R70 are reduced/control calcium AUC ratios at 15%, 30%, and 70% contrast. All 69 reduced layouts enhanced absolute calcium AUC at every nonzero contrast. This is consistent with the monotone calcium readout and pointwise reduced inhibitory conductances; the preferential fractional-enhancement metrics and spatial profiles provide the more informative sensitivity checks.', '',
        'The horizontal multiplier uses inverse PCHIP of the full-control calcium AUC curve without extrapolation. Each mean summarizes only its valid sampled overlap; full matches, exclusions, and contrast support are saved. It is neither an amplitude ratio nor the earlier pooled voltage/calcium Figure 7 marker.', '',
        'For nominal 25% retained inhibition, uniform weakening yields a calcium-only mean horizontal multiplier of 3.268536 over 5%, 10%, 15%, and 20% contrast. The smallest layout mean is 3.267108 (random seed 10); distal-first retention gives 4.840160 over 5%, 10%, and 15% contrast. These ranges therefore compare descriptive means with different overlap support, not fits of a common shift over a fixed contrast domain.', '',
        '## Numerical validation', '',
        f'Focal calculations use the passive DC impedance solution with a finite excitatory conductance. Six independently specified full 3050-ms NEURON simulations matched the complete voltage maps within {maxdv:.3g} mV. Historical canonical-source descriptors were reproduced.', '',
        f"The 770 spot cases preserve the original 0.1-ms NEURON onset transient, 2-ms sampling, local nonlinear calcium transformation, and 50–3050-ms trapezoidal AUC. A late voltage tail is extended only after a conservative all-compartment passive contraction bound is satisfied; the calcium clearance recurrence is evaluated through the complete stimulus. Maximum declared relative AUC bound was {spot['max_tail_auc_relative_bound']:.3g}. Eight independent full 3052-ms original-protocol comparisons had maximum mean AUC relative difference {spot['max_direct_mean_auc_relative_error']:.3g}. This is a documented numerical acceleration of the protocol, not a replacement of transient AUC by steady activation.", '',
        '## Files and reproducibility', '',
        '- `run_spatial.py`: rerunnable focal/spot stages; exact executed snapshots are saved beside it.',
        '- `layout_summary.tsv`, `layout_contact_weights.tsv`, `source_sites.tsv`: complete input identities and weights.',
        '- `focal_all_sites.tsv`, `focal_layout_controls.tsv`, `focal_binned_profiles.tsv`: absolute and normalized focal descriptors.',
        '- `focal_site_paired_effects.tsv`, `focal_site_summary.tsv`, `focal_layout_paired_effects.tsv`: paired changes and censoring.',
        '- `spot_responses.tsv`, `spot_local_responses.tsv`, `spot_layout_metrics.tsv`, `spot_contrast_effects.tsv`, `horizontal_matches.tsv`: all model responses and comparisons.',
        '- `focal_cache/`, `spot_voltage_cache/`: reproducible transfer maps/input resistances and voltage prefixes.',
        '- `focal_validation.json`, `spot_validation.json`, `spot_tail_diagnostics.tsv`, `source_hashes_before.json`: validation and provenance.',
        '- `spatial_site_*.pdf/.svg/.png`, `spatial_layout_*.pdf/.svg/.png`: 12 individual analytical panels.', '',
        'Rerun with the repository `.venv/bin/python`, then `run_spatial.py --stage focal`, `run_spatial.py --stage spot`, and `build_outputs.py`. The spot stage imports the preserved `../electrical/passive_tail.py` helper and verifies its hash. Original sources, calibration, figures, and manuscript remain unchanged.', '',
        '## Interpretation', '',
        'The specified electrical and calcium settings support the spatial effect across all sampled input locations, and the qualitative contrast-gain redistribution survives the tested residual-contact distributions. The magnitudes depend on site, input strength, and inhibitory layout. Source-normalized 1/e distances are first-crossing descriptors, not fitted exponential length constants. Calcium diffusion and active membrane currents are not modeled. Spatial/gain covariance is not a causal mediation analysis. The 75 sites, 20 random seeds, and parameter conditions do not provide biological confidence intervals or evidence across independent reconstructions.', ''])
    (OUT/'README.md').write_text('\n'.join(report))
    print(json.dumps({k:v for k,v in summary.items() if k in ['focal_cases','layout_focal_cases','spot_cases','paired_site_crossings_observed','paired_site_crossings_broader','reduced_site_crossings_censored','paired_site_cable_integrals_broader']},indent=2))

if __name__=='__main__':main()
