"""
04_plot_calcium_sweep.py - Plot the Figure S6 and S7 panels

This script plots the results of the calcium transformation sweep as separate panels:
- Absolute peak calcium responses across contrast at the selected parameters and at a domain corner
- Maps over the midpoint and slope grid of the low-contrast preference (P15), the
  intermediate-contrast preference (P30) and the normalized 30% response difference,
  for 20%, 25% and 30% retained inhibition
- Change in calcium spread from the focal input, and the source-normalized profiles
- The calcium-only horizontal contrast multiplier

P15 and P30 are the gain (reduced / control peak calcium) at 15% or 30% contrast divided by the
gain at 70% contrast. The star marks the selected parameters and the open square marks
the corner example. Gray cells are undefined. Each panel is saved as SVG and PNG.

Prerequisites:
    - output/calcium_sweep/ from 03_run_calcium_sweep.py

Usage:
    python 04_plot_calcium_sweep.py

Output:
    - Panels in figures/calcium_sweep_panels/
    - Plotted arrays and tables in figures/calcium_sweep_panels/source_data/

Author: Pratyush Ramakrishna
Date: 2026
"""

import os
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle

# =============================================================================
# CONFIGURATION
# =============================================================================

# Directory paths
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(SCRIPT_DIR, 'output', 'calcium_sweep')
FIGURES_DIR = os.path.join(SCRIPT_DIR, 'figures', 'calcium_sweep_panels')
DATA_DIR = os.path.join(FIGURES_DIR, 'source_data')

# Colors
GRAY, ORANGE, CYAN = '#333333', '#d55e00', '#079ba8'
MISSING = '#d8d8d8'

# Color limits for the preference and shape maps
PREF_NORM = TwoSlopeNorm(vmin=-7, vcenter=0, vmax=7)
SHAPE_NORM = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)

RETAINED_PERCENT = [20, 25, 30]
METRIC_NAMES = ['log2_preference15', 'log2_preference30', 'shape30_difference']

plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 10, 'axes.labelsize': 10,
    'axes.titlesize': 11, 'xtick.labelsize': 9, 'ytick.labelsize': 9,
    'axes.linewidth': .7, 'axes.spines.top': False, 'axes.spines.right': False,
    'xtick.major.width': .7, 'ytick.major.width': .7,
    'svg.fonttype': 'none', 'savefig.facecolor': 'white', 'savefig.dpi': 200,
    'axes.unicode_minus': True,
})

# =============================================================================
# FUNCTIONS
# =============================================================================

def save_panel(fig, name):
    """Save a figure as SVG and PNG in the figures folder, then close it."""
    # Draw once first so the constrained layout settles before saving
    fig.canvas.draw()
    for ext in ['svg', 'png']:
        fig.savefig(os.path.join(FIGURES_DIR, f'{name}.{ext}'), dpi=200)
    plt.close(fig)
    print(f"      Saved: {name}")


def panel_axes(title):
    """Create an empty panel with a title."""
    fig, ax = plt.subplots(figsize=(7.2, 4.8), layout='constrained')
    ax.set_title(title, pad=10)
    return fig, ax


def condition_legend(ax, location='upper left'):
    """Add the control and reduced-inhibition legend."""
    ax.legend(handles=[
        Line2D([], [], color=GRAY, lw=1.5, label='Control'),
        Line2D([], [], color=ORANGE, lw=1.5, label='25% retained inhibition'),
    ], loc=location, frameon=False, fontsize=9)


def map_panel(vals, half, slopes, n, params, cmap, norm, label, title,
              zero=False, cbticks=None, missing_label='Unavailable metric', unbounded=None):
    """
    Plot a map over the midpoint and slope grid.

    Parameters
    ----------
    vals : ndarray
        One value per parameter pair (the first n are the grid)
    half, slopes : ndarray
        Midpoint and slope grid axes (mV)
    n : int
        Number of grid points
    params : ndarray
        Selected (midpoint, slope), drawn as a star
    cmap : str
        Color map name
    norm : Normalize
        Color limits
    label : str
        Color bar label
    title : str
        Panel title
    zero : bool
        Draw a black contour at zero (only if the map changes sign; the title says so)
    cbticks : list, optional
        Color bar ticks
    missing_label : str
        Legend text for gray (undefined) cells
    unbounded : ndarray of bool, optional
        Cells where control is silent and the ratio is unbounded; hatched

    Returns
    -------
    fig, ax
    """
    arr = np.asarray(vals[:n]).reshape(len(half), len(slopes)).T
    crosses_zero = zero and np.nanmin(arr) < 0 < np.nanmax(arr)
    if zero:
        title += '; black contour: zero' if crosses_zero else '; no sign change in the domain'
    fig, ax = panel_axes(title)
    palette = plt.get_cmap(cmap).copy()
    palette.set_bad(MISSING)
    mesh = ax.pcolormesh(half, slopes, np.ma.masked_invalid(arr),
                         shading='nearest', cmap=palette, norm=norm,
                         edgecolors='none', linewidth=0, rasterized=True,
                         antialiased=False)
    if crosses_zero:
        ax.contour(half, slopes, arr, levels=[0], colors='black', linewidths=.8)
    handles = []
    if unbounded is not None and np.any(unbounded[:n]):
        dx, dy = half[1] - half[0], slopes[1] - slopes[0]
        mask = np.asarray(unbounded[:n]).reshape(len(half), len(slopes)).T
        for j, i in zip(*np.nonzero(mask)):
            ax.add_patch(Rectangle((half[i] - dx / 2, slopes[j] - dy / 2), dx, dy, fill=False,
                                   hatch='////', edgecolor='white', linewidth=0, zorder=5))
        handles.append(Patch(facecolor=palette(norm(norm.vmax)), hatch='////', edgecolor='white',
                             label='Control silent (ratio unbounded)'))
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
    ] + handles, loc='upper center', bbox_to_anchor=(.5, -.21), ncol=2,
       frameon=False, fontsize=8, columnspacing=1)
    return fig, ax


def nice_max(value):
    """Round a positive value up to 1, 1.5, 2, 3 or 5 times a power of ten."""
    exponent = np.floor(np.log10(value))
    for m in (1, 1.5, 2, 3, 5, 10):
        if value <= m * 10 ** exponent:
            return float(m * 10 ** exponent)


def safe(value):
    """Blank for undefined values in tables."""
    return '' if not np.isfinite(value) else float(value)


# =============================================================================
# MAIN SCRIPT
# =============================================================================

def main():
    print("=" * 60)
    print("VG3 Amacrine Cell Model - Calcium Sweep Panels")
    print("=" * 60)
    os.makedirs(DATA_DIR, exist_ok=True)

    # Load the sweep results
    print("\n[1/5] Loading sweep results")
    with np.load(os.path.join(RESULTS_DIR, 'calcium_sweep.npz')) as archive:
        z = {key: archive[key] for key in archive.files}
    pairs = z['parameter_pairs']
    half, slopes, contrast = (z[k] for k in ['v_half_grid_mV', 'slope_grid_mV', 'contrast_pct'])
    n = len(half) * len(slopes)
    ix_selected = n
    ix_corner = int(np.flatnonzero(np.all(pairs == [-20, 6], axis=1))[0])
    params = pairs[ix_selected]

    auc, norm_auc = z['mean_peak_ca'], z['normalized_peak_ca']
    # Control silent, reduced inhibition responding (unbounded preference): drawn at the color
    # limit and hatched
    pref = np.log2(z['preference_15_30_over70'])
    unbounded = np.isposinf(pref)
    pref[unbounded] = PREF_NORM.vmax
    shape_delta = norm_auc[:, 1:, 6] - norm_auc[:, 0, 6, None]
    crossing = z['spatial_1e_ctrl25_um']
    spread = crossing[:, 1] - crossing[:, 0]
    shift = z['calcium_shift_mean'][:, 1]
    count = z['calcium_shift_count'][:, 1]
    distance = z['spatial_distance_um']
    print(f"      {n} grid points plus the selected point")

    # Absolute responses at the selected parameters and at the corner
    print("\n[2/5] Plotting absolute calcium responses")
    for example, i, title in [
        ('selected', ix_selected, 'Absolute calcium responses\nSelected calcium parameters'),
        ('domain_corner', ix_corner,
         'Absolute calcium responses at the domain corner\nMidpoint −20 mV; slope 6 mV'),
    ]:
        fig, ax = panel_axes(title)
        for c, color in [(0, GRAY), (2, ORANGE)]:
            ax.plot(contrast, auc[i, c], '-', color=color, lw=1.5)
            markers = np.isin(contrast, [15, 30, 70])
            ax.plot(contrast[markers], auc[i, c, markers], 'o', ms=4,
                    color=color, mec=color, mew=.3)
        top = nice_max(1.05 * max(auc[i, 0].max(), auc[i, 2].max()))
        ax.set(xlim=(0, 72), ylim=(0, top), xticks=[0, 15, 30, 70], yticks=[0, top / 2, top],
               xlabel='Contrast (%)', ylabel='Peak Ca proxy (a.u.)')
        condition_legend(ax)
        if example == 'domain_corner':
            preferences = z['preference_15_30_over70'][i, 1]
            ax.text(.03, .76, f'P15 = {preferences[0]:.2f}\nP30 = {preferences[1]:.2f}',
                    transform=ax.transAxes, ha='left', va='top', linespacing=1.4)
        save_panel(fig, f'calcium_response_{example}')

    # Preference and shape maps for each retained fraction
    print("\n[3/5] Plotting preference and shape maps")
    for rr, fraction in enumerate(RETAINED_PERCENT):
        for cc, metric in enumerate(METRIC_NAMES):
            vals = pref[:, rr, cc] if cc < 2 else shape_delta[:, rr]
            title = [
                'Low-contrast preference: P15 = gain(15%) / gain(70%)',
                'Intermediate preference: P30 = gain(30%) / gain(70%)',
                'Normalized contrast-response shape at 30%',
            ][cc]
            fig, ax = map_panel(vals, half, slopes, n, params, 'RdBu_r',
                                PREF_NORM if cc < 2 else SHAPE_NORM,
                                ['log2 P15', 'log2 P30', 'Reduced − control'][cc],
                                f'{title}\n{fraction}% retained inhibition',
                                zero=True, cbticks=[-6, 0, 6] if cc < 2 else [-1, 0, 1],
                                unbounded=unbounded[:, rr, cc] if cc < 2 else None)
            save_panel(fig, f'retained_{fraction}_{metric}')

    # Calcium spread from the focal input
    print("\n[4/5] Plotting calcium spread and contrast equivalence")
    spread_top = float(np.ceil(np.nanmax(spread[:n]) / 10) * 10)   # next 10 um above the maximum, no clipping
    fig, ax = map_panel(spread, half, slopes, n, params, 'viridis', Normalize(0, spread_top),
                        'Δ 1/e distance (µm)',
                        'Spatial calcium spread at the original focal site\n25% retained inhibition',
                        cbticks=[0, spread_top / 2, spread_top], missing_label='Censored crossing')
    save_panel(fig, 'spatial_spread_difference_um')

    fig, ax = panel_axes('Focal calcium profiles\nSelected calcium parameters')
    for j, color in [(0, GRAY), (1, ORANGE)]:
        ax.plot(distance, z['spatial_profiles_ctrl25'][ix_selected, j], '-',
                color=color, lw=1.5)
    ax.axhline(1 / np.e, ls=':', color='#777777', lw=.8)
    ax.text(104, 1 / np.e + .02, '1/e', ha='right', va='bottom', color='#555555')
    ax.set(xlim=(0, 105), ylim=(0, 1.025), xticks=[0, 50, 100], yticks=[0, .5, 1],
           xlabel='Cable path distance (µm)', ylabel='Remote / source Ca')
    condition_legend(ax, 'upper right')
    save_panel(fig, 'calcium_spatial_profile_selected')

    # Horizontal contrast multiplier, with the number of matched contrasts contoured
    low, high = np.nanpercentile(shift[:n], [2, 98])
    low, high = np.floor(low * 10) / 10, np.ceil(high * 10) / 10
    fig, ax = map_panel(shift, half, slopes, n, params, 'viridis', Normalize(low, high),
                        'Horizontal multiplier',
                        'Calcium-only contrast equivalence; 25% retained inhibition\n'
                        'Dotted contours: valid-match count boundaries',
                        cbticks=[low, (low + high) / 2, high])
    ax.contour(half, slopes, count[:n].reshape(len(half), len(slopes)).T,
               levels=[.5, 1.5, 2.5, 3.5], colors='black', linewidths=.7, linestyles=':')
    save_panel(fig, 'calcium_contrast_equivalence')

    # Save the plotted arrays and example tables
    print("\n[5/5] Saving the plotted data")
    export = {
        'v_half_grid_mV': half, 'slope_grid_mV': slopes, 'parameter_pairs': pairs,
        'example_indices': np.array([ix_selected, ix_corner]),
        'contrast_pct': contrast, 'mean_peak_ca_all_conditions': auc,
        'spatial_spread_difference_um': spread,
        'calcium_contrast_equivalence': shift, 'calcium_shift_count': count,
        'all_horizontal_matches': z['calcium_shift_by_contrast'],
        'all_fractional_gain_15_30_70': z['gain_15_30_70'],
        'spatial_distance_um': distance,
        'spatial_profiles_ctrl25': z['spatial_profiles_ctrl25'],
        'spatial_crossing_ctrl25_um': crossing,
    }
    for rr, fraction in enumerate(RETAINED_PERCENT):
        for cc, metric in enumerate(METRIC_NAMES):
            export[f'retained_{fraction}_{metric}'] = pref[:, rr, cc] if cc < 2 else shape_delta[:, rr]
    np.savez_compressed(os.path.join(DATA_DIR, 'plotted_arrays.npz'), **export)

    map_columns = ['spatial_spread_difference_um', 'calcium_contrast_equivalence', 'calcium_shift_count']
    map_columns += [key for key in export if key.startswith('retained_')]
    with open(os.path.join(DATA_DIR, 'parameter_maps.tsv'), 'w') as stream:
        writer = csv.writer(stream, delimiter='\t')
        writer.writerow(['point_type', 'midpoint_mV', 'slope_mV'] + map_columns)
        for i, (hh, kk) in enumerate(pairs):
            writer.writerow(['grid' if i < n else 'selected_exact', hh, kk] +
                            [safe(export[key][i]) for key in map_columns])

    with open(os.path.join(DATA_DIR, 'example_contrast_curves.tsv'), 'w') as stream:
        writer = csv.writer(stream, delimiter='\t')
        writer.writerow(['parameter_example', 'midpoint_mV', 'slope_mV', 'retained_fraction',
                         'contrast_pct', 'mean_peak_ca'])
        for name, i in [('selected', ix_selected), ('domain_corner', ix_corner)]:
            for c, fraction in [(0, 1), (2, .25)]:
                for ci, cc in enumerate(contrast):
                    writer.writerow([name, *pairs[i], fraction, cc, auc[i, c, ci]])

    with open(os.path.join(DATA_DIR, 'example_spatial_profiles.tsv'), 'w') as stream:
        writer = csv.writer(stream, delimiter='\t')
        writer.writerow(['parameter_example', 'midpoint_mV', 'slope_mV', 'retained_fraction',
                         'distance_um', 'normalized_ca'])
        for c, fraction in [(0, 1), (1, .25)]:
            for ci, dd in enumerate(distance):
                writer.writerow(['selected', *pairs[ix_selected], fraction, dd,
                                 z['spatial_profiles_ctrl25'][ix_selected, c, ci]])

    print("\n" + "=" * 60)
    print("PLOTTING COMPLETE")
    print("=" * 60)
    print(f"\n  Figures: {FIGURES_DIR}")
    print(f"  Data:    {DATA_DIR}")


if __name__ == '__main__':
    main()
