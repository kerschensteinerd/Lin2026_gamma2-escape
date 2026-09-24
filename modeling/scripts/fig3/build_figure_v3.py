#!/usr/bin/env python3
"""Draw calibrated Fig. 3 v3 results with shared anatomy scales and fit labels."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize, ListedColormap

ROOT = Path(__file__).resolve().parents[2]
PROJECT = ROOT.parent
DATA = ROOT / 'outputs/fig3_v3'
OUT = PROJECT / 'figures/main/assets/Fig3_v3'
GRAY, ORANGE, CYAN = '#333333', '#d55e00', '#079ba8'
CONDITIONS = [('ctrl', 'Control', GRAY), ('gaba', 'GABA removed', ORANGE)]
# Match final-size summary circles in Figures 5 v4 and 7 v2 (1.65 pt radius).
SUMMARY_MARKER_SIZE_PT = 3.3
SUMMARY_MARKER_EDGE_PT = .2
SUMMARY_ERROR_WIDTH_PT = .6


def main(model_first=False):
    out = OUT.with_name('Fig3_v4') if model_first else OUT
    out.mkdir(parents=True, exist_ok=True)
    validation = json.loads((DATA / 'validation.json').read_text())
    checks = validation['checks'] if isinstance(validation, dict) else validation
    assert all(check['passed'] for check in checks), 'Simulation validation failed; refusing to render.'
    with np.load(DATA / 'anatomy_maps.npz') as archive:
        a = {key: archive[key] for key in archive.files}
    with np.load(ROOT / 'outputs/fig3_v3_spatial/selected_spatial_maps.npz') as spatial:
        assert np.allclose(a['geometry_edges'], spatial['geometry_edges']), 'Spatial and contrast geometry differ.'
        for key in ['local_source_segment', 'local_ctrl_normalized', 'local_gaba_normalized', 'path_distance_from_source_um']:
            a[key] = spatial[key]
    with np.load(ROOT / 'outputs/fig3_v3_spatial/finite_local_maps.npz') as spatial:
        for key in ('local_ctrl_ca', 'local_gaba_ca'):
            a[key] = spatial[key]
    meta = json.loads((DATA / 'metadata.json').read_text())
    display = json.loads((DATA / 'figure_data.json').read_text())
    experiment = json.loads((PROJECT / 'analysis/fig3_model/experimental_fig3.json').read_text())
    edges, index = a['geometry_edges'], a['edge_segment_index']
    order = np.argsort(edges[:, :, 2].mean(1), kind='stable')
    edges, index = edges[order], index[order]
    xymin, xymax = edges[:, :, :2].min((0, 1)), edges[:, :, :2].max((0, 1))
    limits = (xymin[0]-4, xymax[0]+4, xymin[1]-9, xymax[1]+6)
    source = int(a['local_source_segment'])
    center = a['stimulus_center_xy']
    plt.rcParams.update({'font.family': 'Arial', 'font.size': 7.3, 'axes.labelsize': 7.3,
        'xtick.labelsize': 6.4, 'ytick.labelsize': 6.4, 'axes.linewidth': .5,
        'axes.spines.top': False, 'axes.spines.right': False, 'xtick.major.width': .5,
        'ytick.major.width': .5, 'xtick.major.size': 2, 'ytick.major.size': 2,
        'svg.fonttype': 'none', 'pdf.fonttype': 42, 'savefig.facecolor': 'white'})
    # Keep high responses visible as thin colored cables on white paper.
    # The complete data range is mapped to this palette; no values are clipped.
    palette = ListedColormap(plt.get_cmap('magma')(np.linspace(0, .82, 256)), name='magma_print')

    def anatomy(ax, values, maximum, mark_source=False, spot=False, scale=False, lw=.5, scale_raise_um=0):
        assert np.isfinite(values).all()
        assert values.min() > -1e-5 and values.max() <= maximum + 1e-5
        lc = LineCollection(edges[:, :, :2], array=values[index], cmap=palette,
                            norm=Normalize(0, maximum), linewidths=lw,
                            capstyle='round', joinstyle='round')
        ax.add_collection(lc)
        if mark_source:
            sx, sy = a['xyz_segments'][source, :2]
            ax.plot(sx, sy, 'o', ms=3.5, mfc='white', mec=CYAN, mew=.9, zorder=5)
            ax.annotate('', xy=(sx, sy+1.1), xytext=(sx-4, sy+15),
                        arrowprops=dict(arrowstyle='-|>', color=CYAN, lw=.7, mutation_scale=5))
        if spot:
            ax.plot(*center, '+', ms=3.5, mew=.6, color='#777777', zorder=5)
        if scale:
            x, y = limits[1]-27, limits[2]+7+scale_raise_um
            ax.plot([x, x+20], [y, y], color=GRAY, lw=1, solid_capstyle='butt')
            ax.text(x+10, y-2.4, '20 µm', ha='center', va='top', fontsize=6.2)
        ax.set(xlim=limits[:2], ylim=limits[2:], aspect='equal')
        ax.set_axis_off()
        return lc

    def label(fig, letter, title, x, y):
        fig.text(x, y, letter, fontsize=10, fontweight='bold', va='top')
        fig.text(x+.032, y-.001, title, fontsize=8.0, va='top')

    def colorbar(fig, lc, box, ticks, text):
        cb = fig.colorbar(lc, cax=fig.add_axes(box), orientation='horizontal', ticks=ticks)
        cb.outline.set_visible(False)
        cb.ax.tick_params(length=1.5, pad=1)
        cb.set_label(text, fontsize=6.7, labelpad=1.5)
        return cb

    def save(fig, name):
        for extension in ('svg', 'pdf', 'png'):
            fig.savefig(out / f'{name}.{extension}', dpi=300)
        plt.close(fig)

    fig = plt.figure(figsize=(7, 5.6))
    label(fig, 'A' if model_first else 'D', 'Spread of local calcium activation', .014, .989)
    for column, (key, name, color) in enumerate(CONDITIONS):
        x = .015 + column*.292
        fig.text(x+.137, .947, name, color=color, ha='center', fontsize=7.5)
        lc = anatomy(fig.add_axes([x, .757, .275, .178]),
                     a[f'local_{key}_ca'], 50, mark_source=True, scale=column == 0, scale_raise_um=5)
    colorbar(fig, lc, [.183, .744, .24, .014], [0, 25, 50], 'Ca proxy (a.u.)')
    label(fig, 'B' if model_first else 'E', 'Spatial attenuation', .658, .989)
    ax = fig.add_axes([.74, .790, .235, .145])
    for key, name, color in CONDITIONS:
        curve = display['spatial'][key]
        ax.plot(curve['distance_um'], curve['transfer_ca'], color=color, lw=1.15)
        ax.plot(curve['distance_um'], curve['transfer_v'], color=color, lw=.8, ls=(0, (3, 2)))
    ax.axhline(1/np.e, color='.8', lw=.6, ls=':')
    ax.set(xlim=(0, 105), ylim=(0, 1), xticks=[0, 50, 100], yticks=[0, .5, 1],
           xlabel='Path distance (µm)', ylabel='Remote / source response')
    ax.xaxis.labelpad = 2
    ax.yaxis.labelpad = 2
    ax.plot([], [], color=GRAY, lw=1.1, label='Ca')
    ax.plot([], [], color=GRAY, lw=.8, ls=(0, (3, 2)), label='V')
    ax.legend(loc='upper right', frameon=False, fontsize=5.9, handlelength=1.6, borderpad=0, labelspacing=.2)

    label(fig, 'C' if model_first else 'F', '100-µm OFF spot: local calcium responses', .014, .683)
    calcium_max = max(a[f'spot_{key}_{c}_ca'].max() for key, _, _ in CONDITIONS for c in (15, 30, 70))
    calcium_max = float(np.ceil(calcium_max/5)*5)
    for col, contrast in enumerate((15, 30, 70)):
        fig.text(.208+col*.303, .636, f'{contrast}% contrast', ha='center', fontsize=7.4)
    for row, (key, name, color) in enumerate(CONDITIONS):
        y = .466 - row*.177
        fig.text(.023, y+.083, name, rotation=90, color=color, ha='center', va='center', fontsize=7.2)
        for col, contrast in enumerate((15, 30, 70)):
            lc = anatomy(fig.add_axes([.065+col*.303, y, .287, .165]),
                         a[f'spot_{key}_{contrast}_ca'], calcium_max,
                         spot=True, scale=row == 1 and col == 0, lw=.47)
    colorbar(fig, lc, [.369, .278, .27, .013], np.linspace(0, calcium_max, 3), 'Peak Ca proxy (a.u.)')

    axes = []
    summary_panels = [('D', 'Voltage', .014), ('E', 'Calcium', .513)] if model_first else [
        ('G', 'Voltage', .014), ('H', 'Calcium', .351), ('I', 'OFF imaging comparison', .685)]
    summary_width = .397 if model_first else .244
    for letter, title, x in summary_panels:
        label(fig, letter, title, x, .231)
        axes.append(fig.add_axes([x+.065, .077, summary_width, .113]))
    axv, axc = axes[:2]
    axfit = None if model_first else axes[2]
    summary_checks = []
    for key, name, color in CONDITIONS:
        r = display['contrast'][key]
        contrast = np.asarray(r['contrast_pct'])
        for ax, field in [(axv, 'mean_peak_dv_mV'), (axc, 'mean_peak_ca_au')]:
            values = np.asarray(r[field])
            line, = ax.plot(contrast, values, color=color, lw=1.05)
            assert np.array_equal(line.get_xdata(), contrast)
            assert np.array_equal(line.get_ydata(), values)
            ix = np.isin(contrast, [15, 30, 70])
            points, = ax.plot(contrast[ix], values[ix], 'o', ms=SUMMARY_MARKER_SIZE_PT,
                             mew=SUMMARY_MARKER_EDGE_PT, color=color, clip_on=False)
            assert np.array_equal(points.get_xdata(), contrast[ix])
            assert np.array_equal(points.get_ydata(), values[ix])
            summary_checks.append(dict(condition=key, field=field, line_points=len(contrast),
                                       marker_points=int(ix.sum())))
        if axfit is not None:
            axfit.plot(contrast, r['auc_shape'], color=color, lw=1.05)
            group = experiment['groups'][('Ctrl' if key == 'ctrl' else 'cKO')+'_OFF']['contrasts']
            exp = [group[str(c)]['normalized_mean'] for c in (.15, .3, .7)]
            sem = [group[str(c)]['normalized_sem'] for c in (.15, .3, .7)]
            axfit.errorbar([15, 30, 70], exp, yerr=sem, fmt='o', color=color, mfc=color,
                           ms=SUMMARY_MARKER_SIZE_PT, mew=SUMMARY_MARKER_EDGE_PT,
                           elinewidth=SUMMARY_ERROR_WIDTH_PT, capsize=0, zorder=5, clip_on=False)
    for ax in axes:
        ax.set(xlim=(0, 73), xticks=[0, 30, 70], xlabel='Contrast (%)')
        ax.xaxis.labelpad = 1.5
        ax.yaxis.labelpad = 1.5
    axv.set(ylim=(0, max(max(display['contrast'][key]['mean_peak_dv_mV']) for key, _, _ in CONDITIONS)*1.08), ylabel='Peak ΔV (mV)')
    axc.set(ylim=(0, calcium_max), ylabel='Ca proxy (a.u.)')
    if axfit is not None:
        axfit.set(xlim=(12, 73), xticks=[15, 30, 70], ylim=(-.03, 1.2), yticks=[0, .5, 1], ylabel='Normalized AUC')
        fig.text(.736, .207, 'Lines: model; points: imaging', fontsize=5.7, va='top', color='#555555')
    footer = ('Illustrative calibration; fixed inhibition, 30% glycine retained. Shared calcium scales in A, C and E.'
              if model_first else display['footer'])
    fig.text(.020, .004, footer, fontsize=6.05, color='#555555', va='bottom', linespacing=1.35)
    save(fig, 'Fig3_modeling_panels')

    # Companion voltage maps display the same runs before calcium transformation.
    fig = plt.figure(figsize=(7, 4.2))
    fig.text(.02, .965, '100-µm OFF spot: voltage before calcium activation', fontsize=10, fontweight='bold', va='top')
    vmax = float(np.ceil(max(a[f'spot_{key}_{c}_dv'].max() for key, _, _ in CONDITIONS for c in (15, 30, 70))/5)*5)
    for col, contrast in enumerate((15, 30, 70)):
        fig.text(.218+col*.28, .854, f'{contrast}% contrast', ha='center', fontsize=8)
    for row, (key, name, color) in enumerate(CONDITIONS):
        y = .535 - row*.295
        fig.text(.025, y+.13, name, rotation=90, color=color, ha='center', va='center', fontsize=8)
        for col, contrast in enumerate((15, 30, 70)):
            lc = anatomy(fig.add_axes([.08+col*.28, y, .27, .26]), a[f'spot_{key}_{contrast}_dv'], vmax,
                         spot=True, scale=row == 1 and col == 0, lw=.6)
    colorbar(fig, lc, [.32, .183, .38, .022], np.arange(0, vmax+1, 10), 'Peak ΔV (mV); shared scale')
    companion_footer = display['companion_footer'].replace('Figure 3 v3', 'Figure 3 v4') if model_first else display['companion_footer']
    fig.text(.5, .041, companion_footer, ha='center', fontsize=7, color='#555555', linespacing=1.4)
    save(fig, 'Fig3_voltage_anatomy')
    (out/'render_manifest.json').write_text(json.dumps(dict(
        panels='A-E' if model_first else 'D-I',
        anatomy_sha256=hashlib.sha256((DATA/'anatomy_maps.npz').read_bytes()).hexdigest(),
        figure_data_sha256=hashlib.sha256((DATA/'figure_data.json').read_bytes()).hexdigest(),
        edges=len(edges), geometry='Exact pt3d XY projection, fixed geometry and shared scales; no spatial smoothing',
        limits_um=limits, calcium_scale=[0, calcium_max], voltage_scale=[0, vmax],
        local_calcium_scale=[0, 50], local_profiles='Finite conductance input, source-normalized V and Ca',
        color_palette='magma 0 to 0.82, covering complete data range; high-value cables remain visible on white',
        plotted_summary_values_exactly_preserved=True, summary_checks=summary_checks,
        summary_style=dict(panels='D-E' if model_first else 'G-I',
                           axes_width_pt=summary_width*504,
                           marker_diameter_pt=SUMMARY_MARKER_SIZE_PT,
                           marker_edge_width_pt=SUMMARY_MARKER_EDGE_PT,
                           error_bar_width_pt=SUMMARY_ERROR_WIDTH_PT, error_bar_caps=False,
                           references=['Fig5_v4', 'Fig7_v2'],
                           marker_clipping=False,
                           marker_fill='Filled model points in D-E' if model_first else 'Filled model points in G-H and filled imaging points in I')), indent=2)+'\n')
    print(out)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model-first', action='store_true',
                        help='Render A–E for the model-first Figure 3 v4 layout, omitting the imaging comparison.')
    main(model_first=parser.parse_args().model_first)
