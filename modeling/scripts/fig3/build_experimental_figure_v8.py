#!/usr/bin/env python3
"""Extend v6 imaging trace heights so the two intact figure blocks align."""
from pathlib import Path
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from build_figure_v3 import GRAY, ORANGE, SUMMARY_MARKER_SIZE_PT, SUMMARY_MARKER_EDGE_PT, SUMMARY_ERROR_WIDTH_PT

PROJECT = Path(__file__).resolve().parents[3]
SOURCE = PROJECT / 'analysis/fig3_model/experimental_fig3.json'
OUT = PROJECT / 'figures/main/assets/Fig3_v8'
W, H = 135*72/25.4 - 16, 384.


def main():
    data = json.loads(SOURCE.read_text())
    plt.rcParams.update({'font.family': 'Arial', 'font.size': 7.3,
        'axes.labelsize': 7.3, 'xtick.labelsize': 6.4, 'ytick.labelsize': 6.4,
        'axes.linewidth': .5, 'axes.spines.top': False, 'axes.spines.right': False,
        'xtick.major.width': .5, 'ytick.major.width': .5,
        'xtick.major.size': 2, 'ytick.major.size': 2,
        'svg.fonttype': 'none', 'pdf.fonttype': 42, 'savefig.facecolor': 'white'})
    fig = plt.figure(figsize=(W/72, H/72))
    plot_boxes = {}

    def text(x, top, s, **kwargs):
        return fig.text(x/W, 1-top/H, s, va='top', **kwargs)

    def axes(name, x, top, width, height):
        plot_boxes[name] = dict(x=x, top=top, width=width, height=height)
        return fig.add_axes([x/W, 1-(top+height)/H, width/W, height/H])

    text(0, 0, 'F', fontsize=10, fontweight='bold')
    text(14, .4, '100-µm spot: dendritic calcium', fontsize=8)
    handles = [Line2D([], [], color=color, lw=1.05, label=name)
               for name, color in [('Ctrl', GRAY), ('VG3-Gabrg2cKO', ORANGE)]]
    fig.legend(handles=handles, loc='upper left', bbox_to_anchor=(.47, 1),
               bbox_transform=fig.transFigure, frameon=False, ncol=2, fontsize=7.3,
               handlelength=1.4, columnspacing=.8, labelspacing=.3, borderpad=0)
    columns = [(21.+i*W/3, c) for i,c in enumerate((15,30,70))]
    trace_width, trace_height = 99., 40.45
    trace_artists, summary_artists = [], []
    for x, contrast in columns:
        text(x+trace_width/2, 19, f'{contrast}% contrast', ha='center', fontsize=7.4)

    for layer, heading_top, row_tops in [('OFF', 33, [47., 91.75]), ('ON', 143.5, [157.5, 202.25])]:
        text(21., heading_top, f'{layer} sublamina', fontsize=8)
        for x, contrast in columns:
            band = Rectangle(((x+trace_width/2)/W, 1-(row_tops[1]+trace_height)/H),
                             trace_width/(2*W), (row_tops[1]+trace_height-row_tops[0])/H,
                             transform=fig.transFigure, fc='#e5e6e7', ec='none', zorder=-1)
            fig.patches.append(band)
            for row, (genotype, color) in enumerate([('Ctrl', GRAY), ('cKO', ORANGE)]):
                source = data['groups'][genotype+'_'+layer]['average_traces']
                time = np.asarray(source['time_s']['values'], dtype=float)
                values = np.asarray(source[str(contrast/100)]['values'], dtype=float)
                ax = axes(f'F_{genotype}_{layer}_{contrast}', x, row_tops[row], trace_width, trace_height)
                ax.patch.set_alpha(0)
                ax.axhline(0, color='#777777', lw=.5, ls=(0, (3, 3)))
                line, = ax.plot(time, values, color=color, lw=1.05)
                ax.set(xlim=(0, 6), ylim=(-.04, 1.05))
                ax.set_axis_off()
                assert np.array_equal(line.get_xdata(), time)
                assert np.array_equal(line.get_ydata(), values)
                trace_artists.append(dict(group=genotype+'_'+layer, contrast=contrast,
                                          count=len(values), source_range=source[str(contrast/100)]['range']))

    fig.add_artist(Line2D([154/W, (154+trace_width/2)/W], [1-34/H]*2,
                          transform=fig.transFigure, color=GRAY, lw=.7, solid_capstyle='butt'))
    text(154+trace_width/4, 36, '3 s', ha='center', fontsize=6.4)
    scale_height = trace_height*.3/1.09
    fig.add_artist(Line2D([252/W]*2, [1-43/H, 1-(43-scale_height)/H],
                          transform=fig.transFigure, color=GRAY, lw=.7, solid_capstyle='butt'))
    text(256, 35, '0.3 Norm. F', fontsize=6.4)

    for letter, layer, left in [('G', 'OFF', 0.), ('H', 'ON', (W/3)/W)]:
        text(left*W, 256, letter, fontsize=10, fontweight='bold')
        text(left*W+14, 256.4, f'{layer} sublamina', fontsize=8)
        ax = axes(letter, left*W+33, 280, 82, 80)
        for genotype, color in [('Ctrl', GRAY), ('cKO', ORANGE)]:
            group = data['groups'][genotype+'_'+layer]['contrasts']
            means = np.array([group[str(c)]['normalized_mean'] for c in (.15, .3, .7)])
            sems = np.array([group[str(c)]['normalized_sem'] for c in (.15, .3, .7)])
            error = ax.errorbar([15, 30, 70], means, yerr=sems, fmt='o-', color=color,
                mfc=color, ms=SUMMARY_MARKER_SIZE_PT, mew=SUMMARY_MARKER_EDGE_PT,
                lw=1.05, elinewidth=SUMMARY_ERROR_WIDTH_PT, capsize=0, clip_on=False)
            assert np.array_equal(error.lines[0].get_ydata(), means)
            assert not error.lines[1]
            summary_artists.append(dict(group=genotype+'_'+layer, means=means.tolist(), sems=sems.tolist()))
        cko = data['groups']['cKO_'+layer]['contrasts']['0.3']
        ax.text(30, cko['normalized_mean']+cko['normalized_sem']+.025, '***',
                ha='center', va='bottom', fontsize=7.3)
        ax.set(xlim=(12, 73), ylim=(-.03, 1.2), xticks=[15, 30, 70], yticks=[0, .5, 1],
               xlabel='Contrast (%)', ylabel='Normalized AUC')
        ax.xaxis.labelpad = ax.yaxis.labelpad = 1.5

    text(2*W/3, 256, 'I', fontsize=10, fontweight='bold')
    text(14+2*W/3, 256.4, 'OFF comparison', fontsize=8)
    text(33+2*W/3, 266, 'Lines: calibrated model', fontsize=5.5, color='#555555')
    text(33+2*W/3, 272, 'Points: imaging', fontsize=5.5, color='#555555')
    ax = axes('I', 33+2*W/3, 280, 82, 80)
    model = json.loads((PROJECT/'vg3_cid3_rf_simulation-main/outputs/fig3_v3/figure_data.json').read_text())
    for key, genotype, color in [('ctrl','Ctrl',GRAY), ('gaba','cKO',ORANGE)]:
        curve = model['contrast'][key]
        line, = ax.plot(curve['contrast_pct'], curve['auc_shape'], color=color, lw=1.05)
        assert np.array_equal(line.get_xdata(), curve['contrast_pct'])
        assert np.array_equal(line.get_ydata(), curve['auc_shape'])
        group = data['groups'][genotype+'_OFF']['contrasts']
        means = [group[str(c)]['normalized_mean'] for c in (.15,.3,.7)]
        sems = [group[str(c)]['normalized_sem'] for c in (.15,.3,.7)]
        error = ax.errorbar([15,30,70], means, yerr=sems, fmt='o', color=color, mfc=color,
                    ms=SUMMARY_MARKER_SIZE_PT, mew=SUMMARY_MARKER_EDGE_PT,
                    elinewidth=SUMMARY_ERROR_WIDTH_PT, capsize=0, clip_on=False, zorder=5)
        assert np.array_equal(error.lines[0].get_ydata(), means)
        assert not error.lines[1]
    ax.set(xlim=(12,73), ylim=(-.03,1.2), xticks=[15,30,70], yticks=[0,.5,1],
           xlabel='Contrast (%)', ylabel='Normalized AUC')
    ax.xaxis.labelpad = ax.yaxis.labelpad = 1.5
    fig.canvas.draw()
    for name, box in plot_boxes.items():
        bbox = fig.axes[list(plot_boxes).index(name)].get_position()
        box['rendered_width'] = float(bbox.width*W)
        box['rendered_height'] = float(bbox.height*H)
    OUT.mkdir(parents=True, exist_ok=True)
    for extension in ('pdf', 'svg', 'png'):
        fig.savefig(OUT/f'Fig3_experimental_panels.{extension}', dpi=300)
    plt.close(fig)
    (OUT/'experimental_render_manifest.json').write_text(json.dumps(dict(
        source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        plot_boxes_points=plot_boxes,
        traces=trace_artists, summaries=summary_artists, trace_data_exactly_preserved=True,
        normalized_summary_values_exactly_preserved=True,
        source='Existing source-workbook export; no refit, resampling, or renormalization',
        layout=dict(page_points=[W,H], trace_width_pt=trace_width, trace_height_pt=trace_height, stimulus_cartoon='Omitted in narrow format; spot size in F title, contrast in column headings, OFF phase shaded',
                    former_trace_width_pt=56.693*504/393,
                    sublamina_labels='Headings above the paired genotype trace rows; same full names in G/H',
                    summary_axes_width_pt=82, summary_axes_height_pt=80, panel_order='F full-width traces; G/H imaging; I calibrated comparison last',
                    source_layout='v6, with all horizontal placements unchanged; only trace heights and downstream vertical positions adjusted',
                    prior_page_height_pt=354, prior_trace_height_pt=32.7,
                    trace_height_increase_pt=7.75, summary_vertical_shift_pt=31,
                    final_summary_axis_bottom_pt=360, aligns_with_model_panels='D/E',
                    fluorescence_scale_height_pt=scale_height),
        style=dict(colors=[GRAY,ORANGE], marker_diameter_pt=SUMMARY_MARKER_SIZE_PT,
                   marker_edge_width_pt=SUMMARY_MARKER_EDGE_PT, error_bar_caps=False,
                   font='Arial', headings_pt=8, axis_pt=7.3, ticks_pt=6.4)), indent=2)+'\n')
    print(OUT/'Fig3_experimental_panels.pdf')


if __name__ == '__main__':
    main()
