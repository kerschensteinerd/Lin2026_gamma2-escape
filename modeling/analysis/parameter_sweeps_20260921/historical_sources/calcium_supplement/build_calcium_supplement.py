#!/usr/bin/env python3
"""Publication-scale plots of validated calcium sweep; no simulation or fitting."""
from pathlib import Path
import csv
import hashlib
import json
import platform
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / 'results/calcium_sweep.npz'
OUT = HERE
DATA = OUT / 'source_data'
W, H = 180.0, 225.0
MM = 25.4
GRAY, ORANGE, CYAN = '#333333', '#d55e00', '#079ba8'
MISSING = '#d8d8d8'
with np.load(SOURCE) as archive:
    z = {key: archive[key] for key in archive.files}
source_hash = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
pairs = z['parameter_pairs']
n = len(z['v_half_grid_mV']) * len(z['slope_grid_mV'])
assert n == 2601 and len(pairs) == 2602
ix_selected = n
def point(half, slope):
    found = np.flatnonzero(np.all(pairs == [half, slope], axis=1))
    assert len(found) == 1
    return int(found[0])
ix_corner = point(-20, 6)
params = pairs[ix_selected]
half, slopes, contrast = (z[k] for k in ['v_half_grid_mV', 'slope_grid_mV', 'contrast_pct'])
auc, norm_auc = z['mean_auc'], z['normalized_auc']
pref = np.log2(z['preference_15_30_over70'])
shape_delta = norm_auc[:, 1:, 6] - norm_auc[:, 0, 6, None]
crossing = z['spatial_1e_ctrl25_um']
spread = crossing[:, 1] - crossing[:, 0]
shift = z['calcium_shift_mean'][:, 1]
count = z['calcium_shift_count'][:, 1]
PREF_NORM = TwoSlopeNorm(vmin=-7, vcenter=0, vmax=7)
SHAPE_NORM = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)

plt.rcParams.update({
    'font.family': 'Arial', 'font.size': 7.3, 'axes.labelsize': 7.3,
    'axes.titlesize': 7.8, 'xtick.labelsize': 6.8, 'ytick.labelsize': 6.8,
    'axes.linewidth': .5, 'axes.spines.top': False, 'axes.spines.right': False,
    'xtick.major.width': .5, 'ytick.major.width': .5,
    'xtick.major.size': 2, 'ytick.major.size': 2,
    'svg.fonttype': 'none', 'pdf.fonttype': 42,
    'savefig.facecolor': 'white', 'savefig.dpi': 300,
    'axes.unicode_minus': True,
})

checks = []
plot_manifest = {}
figures = []
def record(name, actual, expected):
    ok = np.array_equal(np.asarray(actual), np.asarray(expected), equal_nan=True)
    checks.append({'check': name, 'passed': bool(ok)})
    assert ok, name

def text(fig, x, top, string, **kw):
    return fig.text(x/W, 1-top/fig._height_mm, string, va='top', **kw)

def axes(fig, x, top, width, height):
    return fig.add_axes([x/W, 1-(top+height)/fig._height_mm, width/W, height/fig._height_mm])

def panel(fig, letter, title, x, top, subtitle=None):
    text(fig, x, top, letter, fontsize=10, fontweight='bold')
    text(fig, x+5.0, top+.3, title, fontsize=7.8)
    if subtitle:
        text(fig, x+5.0, top+4.8, subtitle, fontsize=6.8, color='#555555')

def base(height, title, detail):
    fig = plt.figure(figsize=(W/MM, height/MM))
    fig._height_mm = height
    text(fig, 1, 1.1, title, fontsize=9.3, fontweight='bold')
    text(fig, 179, 1.5, 'Related to Figure 4', fontsize=7.3, ha='right')
    text(fig, 1, 6.2, detail, fontsize=7.0, color='#444444')
    return fig

def symbols(ax):
    ax.plot(*params, '*', ms=6.8, mfc=CYAN, mec='white', mew=.55,
            clip_on=False, zorder=8)
    ax.plot(-20, 6, 's', ms=3.2, mfc='none', mec='black', mew=.6,
            clip_on=False, zorder=8)

def map_axes(fig, x, top, vals, cmap, norm, label, name, zero=False,
             cbar=True, width=59.5, height=29.5, cbticks=None, markers=True):
    ax = axes(fig, x, top, width, height)
    arr = np.asarray(vals[:n]).reshape(len(half), len(slopes)).T
    palette = plt.get_cmap(cmap).copy()
    palette.set_bad(MISSING)
    mesh = ax.pcolormesh(half, slopes, np.ma.masked_invalid(arr),
                        shading='nearest', cmap=palette, norm=norm,
                        edgecolors='none', linewidth=0, rasterized=True,
                        antialiased=False)
    if zero and np.nanmin(arr) < 0 < np.nanmax(arr):
        ax.contour(half, slopes, arr, levels=[0], colors='black', linewidths=.65)
    ax.set(xlim=(-45, -20), ylim=(1, 6), xlabel='Activation midpoint (mV)',
           ylabel='Slope (mV)', xticks=[-45, -35, -25, -20], yticks=[1, 3, 6])
    ax.xaxis.labelpad = 2
    ax.yaxis.labelpad = 2
    if markers:
        symbols(ax)
    if cbar:
        cax = axes(fig, x+width+2.2, top, 2.1, height)
        cb = fig.colorbar(mesh, cax=cax, ticks=cbticks)
        cb.outline.set_linewidth(.4)
        cb.ax.tick_params(width=.4, length=1.6, pad=1.4, labelsize=6.3)
        cb.set_label(label, fontsize=6.7, labelpad=2)
    record(name, mesh.get_array().filled(np.nan), arr)
    plot_manifest[name] = {'source': 'source_data/plotted_arrays.npz',
                           'array': name, 'reshape': [51, 51], 'transpose': True,
                           'x_limits': [-45, -20], 'y_limits': [1, 6],
                           'cmap': cmap, 'color_min': norm.vmin, 'color_max': norm.vmax,
                           'masked_count': int(np.sum(~np.isfinite(vals[:n])))}
    return ax, mesh

def response_ax(fig, x, top, point_ids, styles, ylim):
    ax = axes(fig, x, top, 64, 29.5)
    for i, style in zip(point_ids, styles):
        for c, color in [(0, GRAY), (2, ORANGE)]:
            line, = ax.plot(contrast, auc[i, c]/1000, style, color=color, lw=1.0)
            record(f'contrast_curve_{i}_{c}', line.get_ydata(), auc[i,c]/1000)
            markers = np.isin(contrast, [15, 30, 70])
            ax.plot(contrast[markers], auc[i,c,markers]/1000, 'o', ms=2.1,
                    color=color, mec=color, mew=.25)
    ax.set(xlim=(0, 72), ylim=ylim, xticks=[0, 15, 30, 70],
           xlabel='Contrast (%)', ylabel='Ca AUC (10³ a.u. ms²)')
    ax.xaxis.labelpad = 2
    ax.yaxis.labelpad = 2
    return ax

fig = base(H, 'Figure Sx. Sensitivity of the calcium transformation',
           '25% retained inhibition; 2,601 parameter pairs; shared 50-ms calcium clearance')
row = [14, 64, 114, 164]
panel(fig, 'A', 'Absolute calcium responses', 1, row[0],
      'Selected calcium parameters')
ax = response_ax(fig, 14, row[0]+11, [ix_selected], ['-'], (0, 155))
ax.set_yticks([0, 50, 100, 150])
ax.legend(handles=[Line2D([],[],color=GRAY,lw=1,label='Control'),
                   Line2D([],[],color=ORANGE,lw=1,label='25% retained')],
          loc='upper left', frameon=False, fontsize=6.5, labelspacing=.3,
          borderpad=0, handlelength=1.4)
panel(fig, 'B', 'A regime with reversed preference', 92, row[0],
      'Midpoint -20 mV; slope 6 mV (domain corner)')
ax = response_ax(fig, 105, row[0]+11, [ix_corner], ['-'], (0, 45))
ax.set_yticks([0, 15, 30, 45])
ax.text(.035, .94, 'P15 = 0.72\nP30 = 0.87', transform=ax.transAxes,
        va='top', fontsize=7.0, linespacing=1.3)

panel(fig, 'C', 'Preference at low contrast', 1, row[1],
      'P15 = gain(15%) / gain(70%)')
map_axes(fig, 14, row[1]+11, pref[:,1,0], 'RdBu_r', PREF_NORM,
         'log2 P15', 'main_C_log2_preference15', zero=True, cbticks=[-6,0,6])
panel(fig, 'D', 'Preference at intermediate contrast', 92, row[1],
      'P30 = gain(30%) / gain(70%)')
map_axes(fig, 105, row[1]+11, pref[:,1,1], 'RdBu_r', PREF_NORM,
         'log2 P30', 'main_D_log2_preference30', zero=True, cbticks=[-6,0,6])

panel(fig, 'E', 'Normalized contrast-response shape', 1, row[2],
      '30% response after separate endpoint scaling')
map_axes(fig, 14, row[2]+11, shape_delta[:,1], 'RdBu_r', SHAPE_NORM,
         'Reduced - control', 'main_E_shape30_difference', zero=True, cbticks=[-1,0,1])
panel(fig, 'F', 'Spatial calcium spread', 92, row[2],
      'Original focal site; 27 reduced crossings censored')
map_axes(fig, 105, row[2]+11, spread, 'viridis', Normalize(0,80),
         'Δ 1/e distance (µm)', 'main_F_spread_difference_um', cbticks=[0,40,80])
assert np.nanmax(spread[:n]) < 80

panel(fig, 'G', 'Focal calcium profiles', 1, row[3],
      'Selected calcium parameters')
ax = axes(fig, 14, row[3]+11, 64, 29.5)
distance = z['spatial_distance_um']
for i, style in [(ix_selected,'-')]:
    for j, color in [(0,GRAY),(1,ORANGE)]:
        line, = ax.plot(distance, z['spatial_profiles_ctrl25'][i,j], style,
                        color=color, lw=1.0)
        record(f'spatial_profile_{i}_{j}', line.get_ydata(), z['spatial_profiles_ctrl25'][i,j])
ax.axhline(1/np.e, ls=':', color='#777777', lw=.6)
ax.text(104, 1/np.e+.02, '1/e', ha='right', va='bottom', fontsize=6.4, color='#555555')
ax.set(xlim=(0,105), ylim=(0,1.025), xticks=[0,50,100], yticks=[0,.5,1],
       xlabel='Cable path distance (µm)', ylabel='Remote / source Ca')
ax.xaxis.labelpad=2
ax.yaxis.labelpad=2
panel(fig, 'H', 'Calcium-only contrast equivalence', 92, row[3],
      'Dotted lines: boundaries in valid match count')
ax, mesh = map_axes(fig, 105, row[3]+11, shift, 'viridis', Normalize(3.1,3.5),
         'Horizontal multiplier', 'main_H_calcium_shift', cbticks=[3.1,3.3,3.5])
ax.contour(half, slopes, count[:n].reshape(51,51).T,
           levels=[.5,1.5,2.5,3.5], colors='black', linewidths=.45, linestyles=':')

legend_ax = axes(fig, 1, 215.0, 178, 4)
legend_ax.axis('off')
handles = [Line2D([],[],ls='',marker='*',ms=7,mfc=CYAN,mec='white',mew=.5,label='Selected parameters'),
           Line2D([],[],ls='',marker='s',ms=3.5,mfc='none',mec='black',mew=.6,label='Domain corner (-20, 6 mV)'),
           Patch(facecolor=MISSING,edgecolor='none',label='Unavailable metric')]
legend_ax.legend(handles=handles,loc='center left',ncol=3,fontsize=6.4,frameon=False,
                 handlelength=1.05,handletextpad=.4,columnspacing=1.1,borderpad=0)
text(fig,1,220.3,'Black contours in C-E: zero effect. All maps retain the full prescribed domain; no refitting.',
     fontsize=6.8,color='#444444')
figures.append((fig, 'Figure_Sx_calcium_sweep'))

# A companion page is needed to show all retained fractions without overcrowding
# the publication-sized main figure. Every column uses the same color range.
HS = 181.0
fig2 = base(HS, 'Figure Sx - companion. Retained-inhibition sensitivity',
            'Same midpoint/slope grid and calcium transformation in all conditions; no electrical refitting')
xs = [12,70,128]
tops = [17,66,115]
column_titles = ['Low-contrast preference','Intermediate preference','Normalized 30% shape']
map_width, map_height = 40.5, 28.0
for rr,fraction in enumerate([20,25,30]):
    for cc in range(3):
        letter = chr(ord('A')+rr*3+cc)
        top,x = tops[rr],xs[cc]
        panel(fig2,letter,f'{fraction}% retained',x-10,top)
        text(fig2,x-5,top+5.0,column_titles[cc],fontsize=6.6,color='#555555')
        data = pref[:,rr,cc] if cc<2 else shape_delta[:,rr]
        normalizer = PREF_NORM if cc<2 else SHAPE_NORM
        key = f'companion_{letter}_{fraction}_'+(['log2_preference15','log2_preference30','shape30_difference'][cc])
        ax, mesh = map_axes(fig2,x,top+11,data,'RdBu_r',normalizer,'',key,zero=True,
                            cbar=False,width=map_width,height=map_height)
        ax.set_xticks([-45,-35,-20])
        ax.xaxis.label.set_fontsize(6.6)
        ax.yaxis.label.set_fontsize(6.6)
        ax.tick_params(labelsize=6.4)
        if rr==2:
            cax=axes(fig2,x,162.3,map_width,2.2)
            ticks=[-6,0,6] if cc<2 else [-1,0,1]
            cb=fig2.colorbar(mesh,cax=cax,orientation='horizontal',ticks=ticks)
            cb.outline.set_linewidth(.4)
            cb.ax.tick_params(width=.4,length=1.5,pad=1.3,labelsize=6.4)
            cb.set_label(['log2 P15','log2 P30','Reduced - control'][cc],fontsize=7.0,labelpad=1.7)
text(fig2,1,175.7,'Symbols and technical masks as in Figure Sx. Black contours: zero effect; color ranges shared by column.',
     fontsize=6.7,color='#444444')
figures.append((fig2, 'Figure_Sx_calcium_retained_inhibition'))

# Source data use the exact arrays given to artists, plus the full raw responses
# and all individual shift matches so each plotted summary is reconstructable.
export = {
    'v_half_grid_mV':half, 'slope_grid_mV':slopes, 'parameter_pairs':pairs,
    'example_indices':np.array([ix_selected,ix_corner]),
    'contrast_pct':contrast, 'mean_auc_all_conditions':auc,
    'main_C_log2_preference15':pref[:,1,0],
    'main_D_log2_preference30':pref[:,1,1],
    'main_E_shape30_difference':shape_delta[:,1],
    'main_F_spread_difference_um':spread,
    'main_H_calcium_shift':shift,
    'main_H_shift_count':count,
    'all_horizontal_matches':z['calcium_shift_by_contrast'],
    'all_fractional_gain_15_30_70':z['gain_15_30_70'],
    'spatial_distance_um':distance,
    'spatial_profiles_ctrl25':z['spatial_profiles_ctrl25'],
    'spatial_crossing_ctrl25_um':crossing,
}
for rr,fraction in enumerate([20,25,30]):
    for cc in range(3):
        letter=chr(ord('A')+rr*3+cc)
        key=f'companion_{letter}_{fraction}_'+(['log2_preference15','log2_preference30','shape30_difference'][cc])
        export[key]=pref[:,rr,cc] if cc<2 else shape_delta[:,rr]
np.savez_compressed(DATA/'plotted_arrays.npz',**export)
map_columns=['main_C_log2_preference15','main_D_log2_preference30',
             'main_E_shape30_difference','main_F_spread_difference_um',
             'main_H_calcium_shift','main_H_shift_count']
map_columns += [key for key in export if key.startswith('companion_')]
def safe(value):
    return '' if not np.isfinite(value) else float(value)
with (DATA/'parameter_maps.tsv').open('w') as stream:
    writer=csv.writer(stream,delimiter='\t')
    writer.writerow(['point_type','midpoint_mV','slope_mV']+map_columns)
    for i,(hh,kk) in enumerate(pairs):
        writer.writerow(['grid' if i<n else 'selected_exact',hh,kk]+[safe(export[key][i]) for key in map_columns])
with (DATA/'example_contrast_curves.tsv').open('w') as stream:
    writer=csv.writer(stream,delimiter='\t')
    writer.writerow(['parameter_example','midpoint_mV','slope_mV','retained_fraction','contrast_pct','mean_ca_auc','plotted_ca_auc_div1000'])
    for name,i in [('selected',ix_selected),('domain_corner',ix_corner)]:
        for c,fraction in [(0,1),(2,.25)]:
            for ci,cc in enumerate(contrast):
                writer.writerow([name,*pairs[i],fraction,cc,auc[i,c,ci],auc[i,c,ci]/1000])
with (DATA/'example_spatial_profiles.tsv').open('w') as stream:
    writer=csv.writer(stream,delimiter='\t')
    writer.writerow(['parameter_example','midpoint_mV','slope_mV','retained_fraction','distance_um','normalized_ca'])
    for name,i in [('selected',ix_selected)]:
        for c,fraction in [(0,1),(1,.25)]:
            for ci,dd in enumerate(distance):
                writer.writerow([name,*pairs[i],fraction,dd,z['spatial_profiles_ctrl25'][i,c,ci]])

# Geometric bounds check for figure text; visual review remains required.
bounds=[]
for fig,name in figures:
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    canvas=fig.bbox
    for t in fig.findobj(matplotlib.text.Text):
        if not t.get_visible() or not t.get_text(): continue
        b=t.get_window_extent(renderer)
        if b.width and b.height:
            outside=b.x0<canvas.x0-.6 or b.y0<canvas.y0-.6 or b.x1>canvas.x1+.6 or b.y1>canvas.y1+.6
            if outside: bounds.append({'figure':name,'text':t.get_text(),'bounds':list(b.bounds)})
    for ext in ['svg','png','pdf']:
        fig.savefig(OUT/f'{name}.{ext}',dpi=300 if ext=='png' else 600,
                    metadata={'Creator':Path(__file__).name} if ext=='pdf' else None)
    plt.close(fig)
checks.append({'check':'all_text_inside_page','passed':not bool(bounds),'outside':bounds})
checks.append({'check':'source_archive_unchanged','passed':hashlib.sha256(SOURCE.read_bytes()).hexdigest()==source_hash})
assert all(c['passed'] for c in checks),checks
(OUT/'qa/plot_validation.json').write_text(json.dumps(checks,indent=2)+'\n')
(DATA/'plot_mapping.json').write_text(json.dumps({
    'source_archive':str(SOURCE),'source_archive_sha256':source_hash,
    'parameter_grid_order':'midpoint outer, slope inner; reshape(51,51).T for pcolormesh',
    'condition_order':[1.0,.20,.25,.30],
    'example_indices':{'selected':ix_selected,'domain_corner':ix_corner},
    'published_exact_point_is_not_grid_member':True,
    'main_A':{'point_indices':[ix_selected],'conditions':[0,2],'array':'mean_auc_all_conditions','scale':.001},
    'main_B':{'point_indices':[ix_corner],'conditions':[0,2],'array':'mean_auc_all_conditions','scale':.001},
    'main_G':{'point_indices':[ix_selected],'array':'spatial_profiles_ctrl25'},
    'maps':plot_manifest,
    'raw_missing_values':'NaN in NPZ; empty TSV cells. Spatial NaN means censored crossing.'
},indent=2)+'\n')
manifest={'source_archive_sha256':source_hash,'builder_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
          'figure_dimensions_mm':{'main':[W,H],'companion':[W,HS]},
          'font':'Arial','base_font_pt':7.3,'panel_letter_pt':10,
          'heatmap_export':'Exact parameter cells embedded at 600 dpi in PDF/SVG to eliminate vector seam artifacts; text, axes, contours and curves remain vector.',
          'python':platform.python_version(),'numpy':np.__version__,'matplotlib':matplotlib.__version__,
          'artifact_hashes':{str(p.relative_to(OUT)):hashlib.sha256(p.read_bytes()).hexdigest()
                            for p in OUT.rglob('*') if p.is_file() and p.name!='manifest.json'},
          'all_array_checks_passed':all(c['passed'] for c in checks),
          'visual_QA':'See qa/visual_review.md; PDF rendering required after build.'}
(OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps({'figures':[str(OUT/f'{name}.pdf') for _,name in figures],
                  'checks_passed':len(checks),'source_sha256':source_hash},indent=2))
