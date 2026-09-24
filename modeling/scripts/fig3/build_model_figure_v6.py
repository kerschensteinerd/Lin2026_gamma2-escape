#!/usr/bin/env python3
"""Draw the model-first Fig. 3 A–E block at native 1.5-column typography."""
from pathlib import Path
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize, ListedColormap
from build_figure_v3 import (ROOT, PROJECT, DATA, GRAY, ORANGE, CYAN,
                            CONDITIONS, SUMMARY_MARKER_SIZE_PT,
                            SUMMARY_MARKER_EDGE_PT, SUMMARY_ERROR_WIDTH_PT)

OUT = PROJECT / 'figures/main/assets/Fig3_v6'
WIDTH_PT = 135 * 72 / 25.4 - 16
HEIGHT_PT = 384.


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    validation = json.loads((DATA / 'validation.json').read_text())
    checks = validation['checks'] if isinstance(validation, dict) else validation
    assert all(check['passed'] for check in checks), 'Simulation validation failed; refusing to render.'
    with np.load(DATA / 'anatomy_maps.npz') as archive:
        a = {key: archive[key] for key in archive.files}
    with np.load(ROOT / 'outputs/fig3_v3_spatial/selected_spatial_maps.npz') as spatial:
        assert np.allclose(a['geometry_edges'], spatial['geometry_edges'])
        for key in ['local_source_segment', 'local_ctrl_normalized', 'local_gaba_normalized', 'path_distance_from_source_um']:
            a[key] = spatial[key]
    with np.load(ROOT / 'outputs/fig3_v3_spatial/finite_local_maps.npz') as spatial:
        for key in ('local_ctrl_ca', 'local_gaba_ca'):
            a[key] = spatial[key]
    display = json.loads((DATA / 'figure_data.json').read_text())
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
    palette = ListedColormap(plt.get_cmap('magma')(np.linspace(0, .82, 256)), name='magma_print')
    fig = plt.figure(figsize=(WIDTH_PT/72, HEIGHT_PT/72))
    plot_boxes = {}

    def text(x, top, label, **kwargs):
        return fig.text(x/WIDTH_PT, 1-top/HEIGHT_PT, label, va='top', **kwargs)

    def axes(name, x, top, width, height):
        plot_boxes[name] = dict(x=x, top=top, width=width, height=height)
        return fig.add_axes([x/WIDTH_PT, 1-(top+height)/HEIGHT_PT, width/WIDTH_PT, height/HEIGHT_PT])

    def panel(letter, title, x, top):
        text(x, top, letter, fontsize=10, fontweight='bold')
        text(x+14, top+.4, title, fontsize=8.)

    # Geometry, depth order, color scale, markers and scale bars are inherited
    # exactly from the v3/v4 anatomy helper. Only physical placement changes.
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

    panel('A', 'Spread of local calcium activation', 0, 0)
    for col, (key, name, color) in enumerate(CONDITIONS):
        x = col * 115.
        text(x+56.5, 17, name, color=color, ha='center', fontsize=7.5)
        anatomy(axes(f'A_{key}', x, 26, 113, 75), a[f'local_{key}_ca'], 50,
                mark_source=True, scale=col == 0, scale_raise_um=5)

    panel('B', 'Spatial attenuation', 245, 0)
    ax = axes('B', 279, 21, 86, 80)
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
    ax.legend(loc='upper right', frameon=False, fontsize=5.9, handlelength=1.6,
              borderpad=0, labelspacing=.2)

    panel('C', '100-µm OFF spot: local calcium responses', 0, 115)
    calcium_max = float(np.ceil(max(a[f'spot_{key}_{c}_ca'].max()
                        for key, _, _ in CONDITIONS for c in (15, 30, 70))/5)*5)
    assert calcium_max == 50, 'A/C share a single scale; re-evaluate if data range changes.'
    for col, contrast in enumerate((15, 30, 70)):
        text(71+col*115, 129, f'{contrast}% contrast', ha='center', fontsize=7.4)
    for row, (key, name, color) in enumerate(CONDITIONS):
        top = 140+row*52
        fig.text(5/WIDTH_PT, 1-(top+24.5)/HEIGHT_PT, name, rotation=90, color=color,
                 ha='center', va='center', fontsize=7.2)
        for col, contrast in enumerate((15, 30, 70)):
            lc = anatomy(axes(f'C_{key}_{contrast}', 16+col*115, top, 110, 49),
                         a[f'spot_{key}_{contrast}_ca'], calcium_max,
                         spot=True, scale=row == 1 and col == 0, lw=.47)
    cb = fig.colorbar(lc, cax=axes('A_C_shared_colorbar', 130, 246, 112, 4),
                      orientation='horizontal', ticks=[0, 25, 50])
    cb.outline.set_visible(False)
    cb.ax.tick_params(length=1.5, pad=1)
    cb.set_label('Ca proxy (a.u.)', fontsize=6.7, labelpad=1.5)

    panel('D', 'Voltage', 0, 274)
    panel('E', 'Calcium', 184, 274)
    axv = axes('D', 43, 290, 105, 70)
    axc = axes('E', 227, 290, 105, 70)
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
                             mew=SUMMARY_MARKER_EDGE_PT, color=color, mfc=color,
                             clip_on=False)
            assert np.array_equal(points.get_xdata(), contrast[ix])
            assert np.array_equal(points.get_ydata(), values[ix])
            summary_checks.append(dict(condition=key, field=field, line_points=len(contrast),
                                       marker_points=int(ix.sum())))
    for ax in (axv, axc):
        ax.set(xlim=(0, 73), xticks=[0, 30, 70], xlabel='Contrast (%)')
        ax.xaxis.labelpad = 1.5
        ax.yaxis.labelpad = 1.5
    axv.set(ylim=(0, max(max(display['contrast'][key]['mean_peak_dv_mV'])
                        for key, _, _ in CONDITIONS)*1.08), ylabel='Peak ΔV (mV)')
    axc.set(ylim=(0, calcium_max), yticks=[0, 25, 50], ylabel='Ca proxy (a.u.)')
    fig.canvas.draw()
    for name, box in plot_boxes.items():
        # Record the dimensions after the equal-aspect anatomy adjustment too.
        ax = fig.axes[list(plot_boxes).index(name)]
        bbox = ax.get_position()
        box['rendered_width'] = float(bbox.width*WIDTH_PT)
        box['rendered_height'] = float(bbox.height*HEIGHT_PT)
    for extension in ('svg', 'pdf', 'png'):
        fig.savefig(OUT / f'Fig3_modeling_panels.{extension}', dpi=300)
    plt.close(fig)
    manifest = dict(
        panels='A-E', canvas_points=[WIDTH_PT, HEIGHT_PT], plot_boxes_points=plot_boxes,
        anatomy_sha256=hashlib.sha256((DATA/'anatomy_maps.npz').read_bytes()).hexdigest(),
        figure_data_sha256=hashlib.sha256((DATA/'figure_data.json').read_bytes()).hexdigest(),
        local_calcium_sha256=hashlib.sha256((ROOT/'outputs/fig3_v3_spatial/finite_local_maps.npz').read_bytes()).hexdigest(),
        edges=len(edges), geometry='Exact pt3d XY projection, equal aspect, fixed geometry and shared scales; no spatial smoothing',
        limits_um=limits, calcium_scale=[0, calcium_max], shared_colorbar_panels=['A', 'C'],
        local_profiles='Finite conductance input, source-normalized V and Ca',
        color_palette='magma 0 to 0.82, covering complete data range; high-value cables remain visible on white',
        plotted_summary_values_exactly_preserved=True, summary_checks=summary_checks,
        summary_style=dict(panels='D-E', axes_width_pt=105, axes_height_pt=70,
                           marker_diameter_pt=SUMMARY_MARKER_SIZE_PT,
                           marker_edge_width_pt=SUMMARY_MARKER_EDGE_PT,
                           error_bar_width_pt=SUMMARY_ERROR_WIDTH_PT, error_bar_caps=False,
                           references=['Fig5_v4', 'Fig7_v2'], marker_clipping=False,
                           marker_fill='Filled model points in D-E'),
        artifact_changes='135-mm layout; D/E reduced from 115x90 to 105x70 pt to balance the middle row; 7 pt clearance below the colorbar label; shared A/C calcium scale and native font/marker sizes preserved')
    (OUT/'model_render_manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(OUT)


if __name__ == '__main__':
    main()
