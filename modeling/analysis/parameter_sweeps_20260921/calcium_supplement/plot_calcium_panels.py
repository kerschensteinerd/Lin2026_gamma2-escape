#!/usr/bin/env python3
"""Plot individual calcium-analysis panels from archived results; no simulation or fitting."""
from pathlib import Path
import argparse
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
ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument('--results', type=Path, default=HERE.parent / 'results')
ap.add_argument('--outdir', type=Path, default=HERE.parent / 'generated/calcium_panels')
args = ap.parse_args()
SOURCE = args.results.resolve() / 'calcium_sweep.npz'
OUT = args.outdir.resolve()
DATA = OUT / 'source_data'
DATA.mkdir(parents=True, exist_ok=True)
(OUT / 'qa').mkdir(parents=True, exist_ok=True)
GRAY, ORANGE, CYAN = '#333333', '#d55e00', '#079ba8'
MISSING = '#d8d8d8'
with np.load(SOURCE) as archive:
    z = {key: archive[key] for key in archive.files}
source_hash = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
pairs = z['parameter_pairs']
half, slopes, contrast = (z[k] for k in ['v_half_grid_mV', 'slope_grid_mV', 'contrast_pct'])
n = len(half) * len(slopes)
assert n == 2601 and len(pairs) == 2602
ix_selected = n
corner_indices = np.flatnonzero(np.all(pairs == [-20, 6], axis=1))
assert len(corner_indices) == 1
ix_corner = int(corner_indices[0])
params = pairs[ix_selected]
auc, norm_auc = z['mean_auc'], z['normalized_auc']
pref = np.log2(z['preference_15_30_over70'])
shape_delta = norm_auc[:, 1:, 6] - norm_auc[:, 0, 6, None]
crossing = z['spatial_1e_ctrl25_um']
spread = crossing[:, 1] - crossing[:, 0]
shift = z['calcium_shift_mean'][:, 1]
count = z['calcium_shift_count'][:, 1]
distance = z['spatial_distance_um']
PREF_NORM = TwoSlopeNorm(vmin=-7, vcenter=0, vmax=7)
SHAPE_NORM = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)

plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 10, 'axes.labelsize': 10,
    'axes.titlesize': 11, 'xtick.labelsize': 9, 'ytick.labelsize': 9,
    'axes.linewidth': .7, 'axes.spines.top': False, 'axes.spines.right': False,
    'xtick.major.width': .7, 'ytick.major.width': .7,
    'svg.fonttype': 'none', 'savefig.facecolor': 'white', 'savefig.dpi': 200,
    'axes.unicode_minus': True,
})

checks = []
plot_manifest = {}
panels = []


def record(name, actual, expected):
    ok = np.array_equal(np.asarray(actual), np.asarray(expected), equal_nan=True)
    checks.append({'check': name, 'passed': bool(ok)})
    assert ok, name


def panel_axes(title):
    fig, ax = plt.subplots(figsize=(7.2, 4.8), layout='constrained')
    ax.set_title(title, pad=10)
    return fig, ax


def condition_legend(ax, location='upper left'):
    ax.legend(handles=[
        Line2D([], [], color=GRAY, lw=1.5, label='Control'),
        Line2D([], [], color=ORANGE, lw=1.5, label='25% retained inhibition'),
    ], loc=location, frameon=False, fontsize=9)


def map_panel(vals, cmap, norm, label, name, title, zero=False, cbticks=None,
              missing_label='Unavailable metric'):
    fig, ax = panel_axes(title)
    arr = np.asarray(vals[:n]).reshape(len(half), len(slopes)).T
    palette = plt.get_cmap(cmap).copy()
    palette.set_bad(MISSING)
    mesh = ax.pcolormesh(half, slopes, np.ma.masked_invalid(arr),
                        shading='nearest', cmap=palette, norm=norm,
                        edgecolors='none', linewidth=0, rasterized=True,
                        antialiased=False)
    if zero and np.nanmin(arr) < 0 < np.nanmax(arr):
        ax.contour(half, slopes, arr, levels=[0], colors='black', linewidths=.8)
    ax.set(xlim=(-45, -20), ylim=(1, 6), xlabel='Activation midpoint (mV)',
           ylabel='Slope (mV)', xticks=[-45, -35, -25, -20], yticks=[1, 3, 6])
    ax.plot(*params, '*', ms=10, mfc=CYAN, mec='white', mew=.7,
            clip_on=False, zorder=8)
    ax.plot(-20, 6, 's', ms=5, mfc='none', mec='black', mew=.8,
            clip_on=False, zorder=8)
    cb = fig.colorbar(mesh, ax=ax, ticks=cbticks, pad=.04)
    cb.outline.set_linewidth(.6)
    cb.set_label(label)
    ax.legend(handles=[
        Line2D([], [], ls='', marker='*', ms=10, mfc=CYAN, mec='white',
               mew=.7, label='Selected parameters'),
        Line2D([], [], ls='', marker='s', ms=5, mfc='none', mec='black',
               mew=.8, label='Domain corner (-20, 6 mV)'),
        Patch(facecolor=MISSING, edgecolor='none', label=missing_label),
    ], loc='upper center', bbox_to_anchor=(.5, -.21), ncol=2,
       frameon=False, fontsize=8, columnspacing=1)
    record(name, mesh.get_array().filled(np.nan), arr)
    plot_manifest[name] = {
        'source': 'source_data/plotted_arrays.npz', 'array': name,
        'reshape': [len(half), len(slopes)], 'transpose': True,
        'x_limits': [-45, -20], 'y_limits': [1, 6],
        'cmap': cmap, 'color_min': norm.vmin, 'color_max': norm.vmax,
        'zero_contour': zero,
        'masked_count': int(np.sum(~np.isfinite(vals[:n]))),
        'missing_label': missing_label,
    }
    panels.append((fig, name))
    return ax


for example, i, ylim, yticks, title in [
    ('selected', ix_selected, (0, 155), [0, 50, 100, 150],
     'Absolute calcium responses\nSelected calcium parameters'),
    ('domain_corner', ix_corner, (0, 45), [0, 15, 30, 45],
     'Absolute calcium responses with reversed preference\nMidpoint −20 mV; slope 6 mV'),
]:
    fig, ax = panel_axes(title)
    for c, color in [(0, GRAY), (2, ORANGE)]:
        line, = ax.plot(contrast, auc[i, c] / 1000, '-', color=color, lw=1.5)
        record(f'contrast_curve_{i}_{c}', line.get_ydata(), auc[i, c] / 1000)
        markers = np.isin(contrast, [15, 30, 70])
        ax.plot(contrast[markers], auc[i, c, markers] / 1000, 'o', ms=4,
                color=color, mec=color, mew=.3)
    ax.set(xlim=(0, 72), ylim=ylim, xticks=[0, 15, 30, 70], yticks=yticks,
           xlabel='Contrast (%)', ylabel='Ca AUC (10³ a.u. ms²)')
    condition_legend(ax)
    if example == 'domain_corner':
        preferences = z['preference_15_30_over70'][i, 1]
        ax.text(.03, .76, f'P15 = {preferences[0]:.2f}\nP30 = {preferences[1]:.2f}',
                transform=ax.transAxes, ha='left', va='top', linespacing=1.4)
    panels.append((fig, f'calcium_response_{example}'))

metric_names = ['log2_preference15', 'log2_preference30', 'shape30_difference']
for rr, fraction in enumerate([20, 25, 30]):
    for cc, metric in enumerate(metric_names):
        vals = pref[:, rr, cc] if cc < 2 else shape_delta[:, rr]
        title = [
            'Low-contrast preference: P15 = gain(15%) / gain(70%)',
            'Intermediate preference: P30 = gain(30%) / gain(70%)',
            'Normalized contrast-response shape at 30%',
        ][cc]
        map_panel(vals, 'RdBu_r', PREF_NORM if cc < 2 else SHAPE_NORM,
                  ['log2 P15', 'log2 P30', 'Reduced − control'][cc],
                  f'retained_{fraction}_{metric}',
                  f'{title}\n{fraction}% retained inhibition; black contour: zero',
                  zero=True, cbticks=[-6, 0, 6] if cc < 2 else [-1, 0, 1])

map_panel(spread, 'viridis', Normalize(0, 80), 'Δ 1/e distance (µm)',
          'spatial_spread_difference_um',
          'Spatial calcium spread at the original focal site\n25% retained inhibition',
          cbticks=[0, 40, 80], missing_label='Censored crossing')
assert np.nanmax(spread[:n]) < 80

fig, ax = panel_axes('Focal calcium profiles\nSelected calcium parameters')
for j, color in [(0, GRAY), (1, ORANGE)]:
    line, = ax.plot(distance, z['spatial_profiles_ctrl25'][ix_selected, j], '-',
                    color=color, lw=1.5)
    record(f'spatial_profile_{ix_selected}_{j}', line.get_ydata(),
           z['spatial_profiles_ctrl25'][ix_selected, j])
ax.axhline(1 / np.e, ls=':', color='#777777', lw=.8)
ax.text(104, 1 / np.e + .02, '1/e', ha='right', va='bottom', color='#555555')
ax.set(xlim=(0, 105), ylim=(0, 1.025), xticks=[0, 50, 100], yticks=[0, .5, 1],
       xlabel='Cable path distance (µm)', ylabel='Remote / source Ca')
condition_legend(ax, 'upper right')
panels.append((fig, 'calcium_spatial_profile_selected'))

ax = map_panel(shift, 'viridis', Normalize(3.1, 3.5), 'Horizontal multiplier',
               'calcium_contrast_equivalence',
               'Calcium-only contrast equivalence; 25% retained inhibition\n'
               'Dotted contours: valid-match count boundaries', cbticks=[3.1, 3.3, 3.5])
ax.contour(half, slopes, count[:n].reshape(len(half), len(slopes)).T,
           levels=[.5, 1.5, 2.5, 3.5], colors='black', linewidths=.7, linestyles=':')
plot_manifest['calcium_contrast_equivalence']['count_contour_levels'] = [.5, 1.5, 2.5, 3.5]

# Exact plotted arrays, raw responses, and individual matches permit independent
# reconstruction of every plotted summary without rerunning the simulations.
export = {
    'v_half_grid_mV': half, 'slope_grid_mV': slopes, 'parameter_pairs': pairs,
    'example_indices': np.array([ix_selected, ix_corner]),
    'contrast_pct': contrast, 'mean_auc_all_conditions': auc,
    'spatial_spread_difference_um': spread,
    'calcium_contrast_equivalence': shift, 'calcium_shift_count': count,
    'all_horizontal_matches': z['calcium_shift_by_contrast'],
    'all_fractional_gain_15_30_70': z['gain_15_30_70'],
    'spatial_distance_um': distance,
    'spatial_profiles_ctrl25': z['spatial_profiles_ctrl25'],
    'spatial_crossing_ctrl25_um': crossing,
}
for rr, fraction in enumerate([20, 25, 30]):
    for cc, metric in enumerate(metric_names):
        export[f'retained_{fraction}_{metric}'] = pref[:, rr, cc] if cc < 2 else shape_delta[:, rr]
np.savez_compressed(DATA / 'plotted_arrays.npz', **export)
map_columns = ['spatial_spread_difference_um', 'calcium_contrast_equivalence', 'calcium_shift_count']
map_columns += [key for key in export if key.startswith('retained_')]


def safe(value):
    return '' if not np.isfinite(value) else float(value)


with (DATA / 'parameter_maps.tsv').open('w') as stream:
    writer = csv.writer(stream, delimiter='\t')
    writer.writerow(['point_type', 'midpoint_mV', 'slope_mV'] + map_columns)
    for i, (hh, kk) in enumerate(pairs):
        writer.writerow(['grid' if i < n else 'selected_exact', hh, kk] +
                        [safe(export[key][i]) for key in map_columns])
with (DATA / 'example_contrast_curves.tsv').open('w') as stream:
    writer = csv.writer(stream, delimiter='\t')
    writer.writerow(['parameter_example', 'midpoint_mV', 'slope_mV', 'retained_fraction',
                     'contrast_pct', 'mean_ca_auc', 'plotted_ca_auc_div1000'])
    for name, i in [('selected', ix_selected), ('domain_corner', ix_corner)]:
        for c, fraction in [(0, 1), (2, .25)]:
            for ci, cc in enumerate(contrast):
                writer.writerow([name, *pairs[i], fraction, cc, auc[i, c, ci], auc[i, c, ci] / 1000])
with (DATA / 'example_spatial_profiles.tsv').open('w') as stream:
    writer = csv.writer(stream, delimiter='\t')
    writer.writerow(['parameter_example', 'midpoint_mV', 'slope_mV', 'retained_fraction',
                     'distance_um', 'normalized_ca'])
    for c, fraction in [(0, 1), (1, .25)]:
        for ci, dd in enumerate(distance):
            writer.writerow(['selected', *pairs[ix_selected], fraction, dd,
                             z['spatial_profiles_ctrl25'][ix_selected, c, ci]])

# Check the rendered labels and exact artist arrays before saving each panel.
bounds = []
for fig, name in panels:
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    canvas = fig.bbox
    for t in fig.findobj(matplotlib.text.Text):
        if not t.get_visible() or not t.get_text():
            continue
        b = t.get_window_extent(renderer)
        if b.width and b.height:
            outside = (b.x0 < canvas.x0 - .6 or b.y0 < canvas.y0 - .6 or
                       b.x1 > canvas.x1 + .6 or b.y1 > canvas.y1 + .6)
            if outside:
                bounds.append({'panel': name, 'text': t.get_text(), 'bounds': list(b.bounds)})
    for ext in ['svg', 'png']:
        fig.savefig(OUT / f'{name}.{ext}', dpi=200)
    plt.close(fig)
checks.append({'check': 'all_text_inside_panel', 'passed': not bool(bounds), 'outside': bounds})
checks.append({'check': 'source_archive_unchanged',
               'passed': hashlib.sha256(SOURCE.read_bytes()).hexdigest() == source_hash})
assert all(c['passed'] for c in checks), checks
(OUT / 'qa/plot_validation.json').write_text(json.dumps(checks, indent=2) + '\n')
(DATA / 'plot_mapping.json').write_text(json.dumps({
    'source_archive': str(SOURCE), 'source_archive_sha256': source_hash,
    'parameter_grid_order': 'midpoint outer, slope inner; reshape(51,51).T for pcolormesh',
    'condition_order': [1.0, .20, .25, .30],
    'example_indices': {'selected': ix_selected, 'domain_corner': ix_corner},
    'selected_exact_point_is_not_grid_member': True,
    'calcium_response_selected': {'point_indices': [ix_selected], 'conditions': [0, 2],
                                  'array': 'mean_auc_all_conditions', 'scale': .001},
    'calcium_response_domain_corner': {'point_indices': [ix_corner], 'conditions': [0, 2],
                                       'array': 'mean_auc_all_conditions', 'scale': .001},
    'calcium_spatial_profile_selected': {'point_indices': [ix_selected],
                                         'array': 'spatial_profiles_ctrl25'},
    'maps': plot_manifest,
    'raw_missing_values': 'NaN in NPZ; empty TSV cells. Spatial NaN means censored crossing.',
}, indent=2) + '\n')
manifest = {
    'source_archive_sha256': source_hash,
    'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'panel_names': [name for _, name in panels], 'formats': ['svg', 'png'],
    'font': 'DejaVu Sans',
    'heatmaps': 'Exact parameter cells; text, axes, contours, and curves are vector in SVG.',
    'python': platform.python_version(), 'numpy': np.__version__, 'matplotlib': matplotlib.__version__,
    'artifact_hashes': {str(p.relative_to(OUT)): hashlib.sha256(p.read_bytes()).hexdigest()
                       for p in OUT.rglob('*') if p.is_file() and p.name != 'manifest.json'},
    'all_checks_passed': all(c['passed'] for c in checks),
}
(OUT / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
print(json.dumps({'panels': [str(OUT / f'{name}.png') for _, name in panels],
                  'checks_passed': len(checks), 'source_sha256': source_hash}, indent=2))
