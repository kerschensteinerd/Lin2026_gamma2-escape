#!/usr/bin/env python3
"""Reflow the unchanged exported Figure 3 imaging values into a common grid."""
from pathlib import Path
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle
from matplotlib.lines import Line2D
from build_figure_v3 import GRAY, ORANGE, SUMMARY_MARKER_SIZE_PT, SUMMARY_MARKER_EDGE_PT, SUMMARY_ERROR_WIDTH_PT

PROJECT = Path(__file__).resolve().parents[3]
SOURCE = PROJECT / 'analysis/fig3_model/experimental_fig3.json'
OUT = PROJECT / 'figures/main/assets/Fig3_v3'
W, H = 504., 306.


def main():
    data = json.loads(SOURCE.read_text())
    plt.rcParams.update({'font.family': 'Arial', 'font.size': 7.3,
        'axes.labelsize': 7.3, 'xtick.labelsize': 6.4, 'ytick.labelsize': 6.4,
        'axes.linewidth': .5, 'axes.spines.top': False, 'axes.spines.right': False,
        'xtick.major.width': .5, 'ytick.major.width': .5,
        'xtick.major.size': 2, 'ytick.major.size': 2,
        'svg.fonttype': 'none', 'pdf.fonttype': 42, 'savefig.facecolor': 'white'})
    fig = plt.figure(figsize=(W/72, H/72))

    def text(x, top, s, **kwargs):
        return fig.text(x/W, 1-top/H, s, va='top', **kwargs)

    def axes(x, top, width, height):
        return fig.add_axes([x/W, 1-(top+height)/H, width/W, height/H])

    text(7.056, 4.5, 'A', fontsize=10, fontweight='bold')
    handles = [Line2D([], [], color=color, lw=1.05, label=name)
               for name, color in [('Ctrl', GRAY), ('VG3-Gabrg2cKO', ORANGE)]]
    fig.legend(handles=handles, loc='upper left', bbox_to_anchor=(.14, 1),
               bbox_transform=fig.transFigure, frameon=False, ncol=2, fontsize=7.3,
               handlelength=1.8, columnspacing=2, labelspacing=.3, borderpad=0)
    columns = [(7.5, 15), (112.5, 30), (217.5, 70)]
    trace_width, trace_height = 89., 43.6
    trace_artists, summary_artists = [], []
    for x, contrast in columns:
        text(x+trace_width/2, 24, f'{contrast}% contrast', ha='center', fontsize=7.4)
        stim = axes(x, 36, trace_width, 23)
        stim.set(xlim=(0, trace_width), ylim=(0, 23), aspect='equal')
        stim.add_patch(Rectangle((0, 0), trace_width, 23, color=(.5, .5, .5), lw=0))
        for cx, sign in [(trace_width/4, 1), (3*trace_width/4, -1)]:
            gray = .5*(1+sign*contrast/100)
            stim.add_patch(Circle((cx, 11.5), 4.2, fc=(gray, gray, gray), lw=0))
        stim.axvline(trace_width/2, color='white', lw=.45)
        stim.set_axis_off()

    for layer, heading_top, row_tops in [('OFF', 69, [82., 131.]), ('ON', 186, [199., 248.])]:
        text(7.5, heading_top, f'{layer} sublamina', fontsize=8)
        for x, contrast in columns:
            band = Rectangle(((x+trace_width/2)/W, 1-(row_tops[1]+trace_height)/H),
                             trace_width/(2*W), (row_tops[1]+trace_height-row_tops[0])/H,
                             transform=fig.transFigure, fc='#e5e6e7', ec='none', zorder=-1)
            fig.patches.append(band)
            for row, (genotype, color) in enumerate([('Ctrl', GRAY), ('cKO', ORANGE)]):
                source = data['groups'][genotype+'_'+layer]['average_traces']
                time = np.asarray(source['time_s']['values'], dtype=float)
                values = np.asarray(source[str(contrast/100)]['values'], dtype=float)
                ax = axes(x, row_tops[row], trace_width, trace_height)
                ax.patch.set_alpha(0)
                ax.axhline(0, color='#777777', lw=.5, ls=(0, (3, 3)))
                line, = ax.plot(time, values, color=color, lw=1.05)
                ax.set(xlim=(0, 6), ylim=(-.04, 1.05))
                ax.set_axis_off()
                assert np.array_equal(line.get_xdata(), time)
                assert np.array_equal(line.get_ydata(), values)
                trace_artists.append(dict(group=genotype+'_'+layer, contrast=contrast,
                                          count=len(values), source_range=source[str(contrast/100)]['range']))

    fig.add_artist(Line2D([143/W, (143+trace_width/2)/W], [1-68/H]*2,
                          transform=fig.transFigure, color=GRAY, lw=.7, solid_capstyle='butt'))
    text(143+trace_width/4, 70, '3 s', ha='center', fontsize=6.4)
    scale_height = trace_height*.3/1.09
    fig.add_artist(Line2D([314/W]*2, [1-119/H, 1-(119-scale_height)/H],
                          transform=fig.transFigure, color=GRAY, lw=.7, solid_capstyle='butt'))
    text(317, 109, '0.3', fontsize=6.4)
    fig.text(330/W, 1-112/H, 'Norm. F', rotation=90, ha='center', va='center', fontsize=6.4)

    for letter, layer, title_top, plot_top in [('B', 'OFF', 69, 84), ('C', 'ON', 186, 201)]:
        text(.658*W, title_top, letter, fontsize=10, fontweight='bold')
        text(.690*W, title_top+.4, f'{layer} sublamina', fontsize=8)
        ax = axes(.74*W, plot_top, .235*W, 85)
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

    OUT.mkdir(parents=True, exist_ok=True)
    for extension in ('pdf', 'svg', 'png'):
        fig.savefig(OUT/f'Fig3_experimental_panels.{extension}', dpi=300)
    plt.close(fig)
    (OUT/'experimental_render_manifest.json').write_text(json.dumps(dict(
        source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        traces=trace_artists, summaries=summary_artists, trace_data_exactly_preserved=True,
        normalized_summary_values_exactly_preserved=True,
        source='Existing source-workbook export; no refit, resampling, or renormalization',
        layout=dict(page_points=[W,H], trace_width_pt=trace_width,
                    former_trace_width_pt=56.693*504/393,
                    sublamina_labels='Headings above the paired genotype trace rows; same full names in B/C',
                    summary_axes_x_pt=.74*W, summary_axes_width_pt=.235*W),
        style=dict(colors=[GRAY,ORANGE], marker_diameter_pt=SUMMARY_MARKER_SIZE_PT,
                   marker_edge_width_pt=SUMMARY_MARKER_EDGE_PT, error_bar_caps=False,
                   font='Arial', headings_pt=8, axis_pt=7.3, ticks_pt=6.4)), indent=2)+'\n')
    print(OUT/'Fig3_experimental_panels.pdf')


if __name__ == '__main__':
    main()
