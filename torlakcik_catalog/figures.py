"""
Figure generation module for the Torlakcık Catalog.

Generates all figures for the paper — PUBLICATION-READY VERSION.

"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
from matplotlib.gridspec import GridSpec
from astropy.coordinates import SkyCoord
import astropy.units as u
from typing import Optional, Dict, List

# - Font setup (cross-platform)
_FONT_CANDIDATES = [
    '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf',
    '/System/Library/Fonts/Helvetica.ttc',
    'C:/Windows/Fonts/arial.ttf',
]
_font_loaded = False
for _fpath in _FONT_CANDIDATES:
    try:
        fm.fontManager.addfont(_fpath)
        _font_loaded = True
        break
    except Exception:
        continue

if _font_loaded:
    plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Helvetica']
else:
    plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'sans-serif']
plt.rcParams['axes.unicode_minus'] = False

# - Publication-quality defaults
plt.rcParams.update({
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'font.size': 9,
    'font.family': 'sans-serif',
    'axes.labelsize': 10,
    'axes.titlesize': 11,
    'axes.titleweight': 'bold',
    'legend.fontsize': 8,
    'legend.frameon': True,
    'legend.edgecolor': '0.8',
    'legend.fancybox': True,
    'figure.facecolor': 'white',
    'axes.facecolor': 'white',
    'axes.edgecolor': '0.3',
    'axes.linewidth': 0.6,
    'xtick.major.width': 0.6,
    'ytick.major.width': 0.6,
    'xtick.minor.width': 0.4,
    'ytick.minor.width': 0.4,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.top': False,
    'ytick.right': False,
    'grid.alpha': 0.3,
    'grid.linewidth': 0.5,
    'lines.markersize': 4,
    'lines.linewidth': 1.2,
})

# - Journal column width
COL_WIDTH = 3.5
COL2_WIDTH = 7.0

# - Color palette (colorblind-friendly)─
C_EMPIRICAL  = '#2176AE'
C_SYNTHETIC  = '#D64933'
C_RETAINED   = '#2A9D8F'
C_EXCLUDED   = '#E76F51'
C_AGE        = '#E63946'
C_MASS       = '#457B9D'
C_METALL     = '#2A9D8F'
C_VAR        = '#7B2D8E'
C_GRID       = '#CCCCCC'
C_ISAACSON   = '#E63946'   # red for Isaacson
C_MEERKAT   = '#2176AE'   # blue for MeerKAT

CRITERIA_COLORS = ['#E63946', '#F4A261', '#E9C46A', '#2A9D8F', '#2176AE', '#7B2D8E', '#264653', '#6B705C']

# - Output directory
FIGURES_DIR = os.path.dirname(os.path.abspath(__file__))
if not os.path.exists(FIGURES_DIR):
    os.makedirs(FIGURES_DIR)


# ---
#  Figure 2: Exclusion by Criterion
# ---
def figure2_exclusion_by_criterion(catalog: pd.DataFrame,
                                    save_path: Optional[str] = None):
    excluded = catalog[catalog['decision'] == 'EXCLUDE']
    total = len(catalog)

    criteria = {
        'R1': r'Mass $>$ 1.5 M$_\odot$',
        'R2': 'Age upper $<$ 3 Gyr',
        'R3': 'Spectral O/B/A/F0–F4',
        'R4': '[Fe/H] $<$ −0.4',
        'R5': 'Multiplicity',
        'R6': 'Phot. variability',
        'R7': 'M-dwarf activity',
    }

    counts = {}
    for code in criteria:
        counts[code] = len(excluded[excluded['reason_code'].str.contains(code, na=False)])

    codes = list(criteria.keys())
    values = [counts[c] for c in codes]
    fractions = [100 * c / total for c in values]

    fig, ax = plt.subplots(figsize=(COL2_WIDTH, 3.5))

    bars = ax.bar(range(len(codes)), values, color=CRITERIA_COLORS,
                  edgecolor='white', linewidth=0.5, width=0.7)

    for bar, val, frac in zip(bars, values, fractions):
        height = bar.get_height()
        ax.annotate(f'{val:,}\n({frac:.1f}%)',
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 4), textcoords="offset points",
                    ha='center', va='bottom', fontsize=7, fontweight='bold')

    ax.set_xticks(range(len(codes)))
    ax.set_xticklabels([f"{c}" for c in codes], fontsize=8)
    ax.set_ylabel('Stars Excluded', fontsize=9)
    ax.set_ylim(0, max(values) * 1.22)

    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(range(len(codes)))
    ax2.set_xticklabels(list(criteria.values()), fontsize=5.5, style='italic', color='0.4')
    ax2.tick_params(axis='x', length=0, pad=12)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    ax2.spines['left'].set_visible(False)

    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f'{x/1000:.0f}k' if x >= 1000 else f'{x:.0f}'))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', color=C_GRID, linewidth=0.4)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, bbox_inches='tight', facecolor='white')
        print(f"Saved: {save_path}")
    plt.close(fig)
    return fig, ax


# ---
#  Figure 3: Spectral Type Breakdown
# ---
def figure3_spectral_type_breakdown(catalog: pd.DataFrame,
                                      save_path: Optional[str] = None):
    sp_order = ['O', 'B', 'A', 'F0-F4', 'F5-F9', 'G', 'K', 'M']
    excluded_counts, retained_counts, exclusion_rates = [], [], []

    for sp in sp_order:
        sp_df = catalog[catalog['spectral_type'] == sp]
        n_total = len(sp_df)
        n_excl = len(sp_df[sp_df['decision'] == 'EXCLUDE'])
        excluded_counts.append(n_excl)
        retained_counts.append(n_total - n_excl)
        exclusion_rates.append(100 * n_excl / n_total if n_total > 0 else 0)

    fig, ax = plt.subplots(figsize=(COL2_WIDTH, 3.5))
    x = np.arange(len(sp_order))
    width = 0.65

    bars1 = ax.bar(x, retained_counts, width, label='Retained',
                   color=C_RETAINED, edgecolor='white', linewidth=0.5)
    bars2 = ax.bar(x, excluded_counts, width, bottom=retained_counts,
                   label='Excluded', color=C_EXCLUDED, edgecolor='white', linewidth=0.5)

    for i, rate in enumerate(exclusion_rates):
        total_h = excluded_counts[i] + retained_counts[i]
        ax.annotate(f'{rate:.0f}%',
                    xy=(i, total_h), xytext=(0, 3), textcoords="offset points",
                    ha='center', va='bottom', fontsize=7.5, fontweight='bold', color='0.2')

    ax.set_xticks(x)
    ax.set_xticklabels(sp_order, fontsize=8)
    ax.set_xlabel('Spectral Type', fontsize=9)
    ax.set_ylabel('Number of Stars', fontsize=9)
    ax.legend(loc='upper right', fontsize=7)
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f'{x/1000:.0f}k' if x >= 1000 else f'{x:.0f}'))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', color=C_GRID, linewidth=0.4)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, bbox_inches='tight', facecolor='white')
        print(f"Saved: {save_path}")
    plt.close(fig)
    return fig, ax


# ---
#  HR / Color-Magnitude Diagram
# ---
def figure_hr_diagram(catalog: pd.DataFrame,
                       save_path: Optional[str] = None):
    from matplotlib.colors import LogNorm
    from matplotlib.patches import Rectangle

    cols = set(catalog.columns)
    has_bp_rp = 'bp_rp' in cols
    has_g_mag = 'phot_g_mean_mag' in cols
    has_teff = 'teff_gspphot' in cols
    has_logg = 'logg_gspphot' in cols
    has_parallax = 'parallax' in cols

    can_compute_Mg = has_g_mag and has_parallax
    use_mode = None

    if has_bp_rp and can_compute_Mg:
        use_mode = 'bp_rp_Mg'
    elif has_teff and can_compute_Mg:
        use_mode = 'teff_Mg'
    elif has_teff and has_logg:
        use_mode = 'teff_logg'

    if use_mode is None:
        print("WARNING: Insufficient columns for HR diagram.")
        return None, None

    print(f"  HR diagram mode: {use_mode}")
    df = catalog.copy()

    if can_compute_Mg:
        valid_plx = df['parallax'] > 0
        df = df[valid_plx].copy()
        df['abs_g'] = df['phot_g_mean_mag'] + 5 * np.log10(df['parallax'] / 100)

    fig, ax = plt.subplots(figsize=(COL2_WIDTH, 6.0))

    if use_mode == 'bp_rp_Mg':
        x_col = 'bp_rp'
        xlabel = r'$G_\mathrm{BP} - G_\mathrm{RP}$'
        ylabel = r'$M_G$ (from $d = 1/\varpi$)'

        valid = df.dropna(subset=[x_col, 'abs_g']).copy()
        valid = valid[(valid[x_col] > -0.3) & (valid[x_col] < 3.8)]
        valid = valid[(valid['abs_g'] > -1) & (valid['abs_g'] < 12)]
        ret_valid = valid[valid['decision'] == 'RETAIN']

        x_bins = np.linspace(-0.3, 3.8, 250)
        y_bins = np.linspace(-1, 12, 250)

        H_all, xedges, yedges = np.histogram2d(
            valid[x_col], valid['abs_g'], bins=[x_bins, y_bins])
        H_all_plot = H_all.T.copy()
        H_all_plot[H_all_plot == 0] = np.nan
        ax.pcolormesh(xedges, yedges, H_all_plot,
                      cmap='Greys', norm=LogNorm(vmin=1, vmax=H_all.max()),
                      rasterized=True, alpha=0.5)

        H_ret, _, _ = np.histogram2d(
            ret_valid[x_col], ret_valid['abs_g'], bins=[x_bins, y_bins])
        H_ret_plot = H_ret.T.copy()
        H_ret_plot[H_ret_plot == 0] = np.nan
        ax.pcolormesh(xedges, yedges, H_ret_plot,
                      cmap='GnBu', norm=LogNorm(vmin=1, vmax=H_ret.max()),
                      rasterized=True, alpha=0.7)

        ax.annotate('Upper MS', xy=(0.35, 1.5), fontsize=7, color='#555', ha='center', style='italic')
        ax.annotate('Giants', xy=(1.0, 0.5), fontsize=7, color='#555', ha='center', style='italic')
        ax.annotate('M dwarfs', xy=(2.8, 9.5), fontsize=7, color='#555', ha='center', style='italic')
        ax.set_xlim(-0.3, 3.5)
        ax.set_ylim(11, -1)

    elif use_mode == 'teff_Mg':
        x_col = 'teff_gspphot'
        xlabel = r'$T_\mathrm{eff}$ (K)'
        ylabel = r'$M_G$ (from $d = 1/\varpi$)'

        valid = df.dropna(subset=[x_col, 'abs_g']).copy()
        ret_valid = valid[valid['decision'] == 'RETAIN']

        x_bins = np.linspace(2800, 10500, 250)
        y_bins = np.linspace(-1, 11.5, 250)

        H_all, xedges, yedges = np.histogram2d(
            valid[x_col], valid['abs_g'], bins=[x_bins, y_bins])
        H_all_plot = H_all.T.copy()
        H_all_plot[H_all_plot == 0] = np.nan
        ax.pcolormesh(xedges, yedges, H_all_plot,
                      cmap='Greys', norm=LogNorm(vmin=1, vmax=H_all.max()),
                      rasterized=True, alpha=0.5)

        H_ret, _, _ = np.histogram2d(
            ret_valid[x_col], ret_valid['abs_g'], bins=[x_bins, y_bins])
        H_ret_plot = H_ret.T.copy()
        H_ret_plot[H_ret_plot == 0] = np.nan
        ax.pcolormesh(xedges, yedges, H_ret_plot,
                      cmap='GnBu', norm=LogNorm(vmin=1, vmax=H_ret.max()),
                      rasterized=True, alpha=0.7)

        # Structural annotations
        ax.annotate('Upper MS\n(excluded)', xy=(8500, 1), fontsize=7, color='#555',
                    ha='center', style='italic')
        ax.annotate('Main Sequence\n(retained)', xy=(5500, 7), fontsize=7, color=C_RETAINED,
                    ha='center', style='italic', fontweight='bold')
        ax.annotate('Red Giants', xy=(4200, 1.5), fontsize=7, color='#555',
                    ha='center', style='italic')
        ax.annotate('M dwarfs', xy=(3400, 10), fontsize=7, color='#555',
                    ha='center', style='italic')

        ax.set_xlim(10000, 3000)
        ax.set_ylim(11, -1)

    elif use_mode == 'teff_logg':
        x_col = 'teff_gspphot'
        y_col = 'logg_gspphot'
        xlabel = r'$T_\mathrm{eff}$ (K)'
        ylabel = r'$\log\,g$'

        valid = df.dropna(subset=[x_col, y_col]).copy()
        valid = valid[(valid[x_col] > 2500) & (valid[x_col] < 10500)]
        valid = valid[(valid[y_col] > 0.5) & (valid[y_col] < 5.5)]
        ret_valid = valid[valid['decision'] == 'RETAIN']

        x_bins = np.linspace(2800, 10500, 200)
        y_bins = np.linspace(0.5, 5.5, 200)

        H_all, xedges, yedges = np.histogram2d(
            valid[x_col], valid[y_col], bins=[x_bins, y_bins])
        H_all_plot = H_all.T.copy()
        H_all_plot[H_all_plot == 0] = np.nan
        ax.pcolormesh(xedges, yedges, H_all_plot,
                      cmap='Greys', norm=LogNorm(vmin=1, vmax=H_all.max()),
                      rasterized=True, alpha=0.4)

        from scipy.ndimage import gaussian_filter
        H_ret, _, _ = np.histogram2d(
            ret_valid[x_col], ret_valid[y_col], bins=[x_bins, y_bins])
        H_ret_smooth = gaussian_filter(H_ret.T, sigma=1.5)

        xcenter = 0.5 * (xedges[:-1] + xedges[1:])
        ycenter = 0.5 * (yedges[:-1] + yedges[1:])

        n_levels = 8
        levels_ret = np.linspace(H_ret_smooth.max()*0.02, H_ret_smooth.max()*0.95, n_levels)
        ax.contourf(xcenter, ycenter, H_ret_smooth, levels=levels_ret,
                    cmap='GnBu', alpha=0.35, extend='max')
        ax.contour(xcenter, ycenter, H_ret_smooth,
                   levels=levels_ret[3:], colors='#1a7a6d', linewidths=0.6, alpha=0.8)

        ax.text(7000, 4.3, 'Main Sequence', fontsize=7, color='#444', ha='center', style='italic', alpha=0.8)
        ax.text(4200, 1.8, 'Red Giants', fontsize=7, color='#444', ha='center', style='italic', alpha=0.8)
        ax.text(3200, 4.9, 'M dwarfs', fontsize=7, color='#444', ha='center', style='italic', alpha=0.8)
        ax.axvline(x=8500, color='0.6', linestyle=':', linewidth=0.5, alpha=0.5)

        ax.set_xlim(10000, 2800)
        ax.set_ylim(5.5, 0.5)

    # - Common formatting
    ax.set_xlabel(xlabel, fontsize=10)
    ax.set_ylabel(ylabel, fontsize=10)

    if use_mode == 'teff_logg':
        legend_elements = [
            mpatches.Patch(facecolor='grey', alpha=0.4, label='Full sample (density)'),
            mpatches.Patch(facecolor=C_RETAINED, alpha=0.35, label='Retained (density)'),
            plt.Line2D([0], [0], color='#1a7a6d', linewidth=0.8, label='Retained (contour)'),
        ]
    else:
        legend_elements = [
            mpatches.Patch(facecolor='grey', alpha=0.5, label='All catalog (density)'),
            mpatches.Patch(facecolor=C_RETAINED, alpha=0.7, label='Retained (density)'),
        ]

    ax.legend(handles=legend_elements, loc='upper right', fontsize=6.5,
              framealpha=0.9, edgecolor='0.7')
    ax.grid(True, color=C_GRID, linewidth=0.3, alpha=0.5)

    if can_compute_Mg:
        fig.text(0.12, 0.01,
                 r'Note: $M_G$ computed via $d = 1/\varpi$; justified by $\varpi/\sigma_\varpi > 5$ cut.',
                 fontsize=6.5, color='0.4')

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, bbox_inches='tight', facecolor='white')
        print(f"Saved: {save_path}")
    plt.close(fig)
    return fig, ax


# ---
#  Figure 4: Empirical vs Synthetic Proxy Comparison
# ---
def figure4_empirical_vs_synthetic(catalog: pd.DataFrame,
                                    save_path: Optional[str] = None):
    from analysis import compare_empirical_vs_synthetic

    proxy_comp = compare_empirical_vs_synthetic(catalog)

    categories = ['Multiplicity\n(R5)', 'Variability\n(R6)', 'Overall\nExclusion']

    empirical_values = [
        proxy_comp['emp_r5_pct'],
        proxy_comp['emp_r6_pct'],
        proxy_comp['emp_overall_pct']
    ]
    synthetic_values = [
        proxy_comp['syn_r5_pct'],
        proxy_comp['syn_r6_pct'],
        proxy_comp['syn_overall_pct']
    ]

    diff = proxy_comp['syn_overall_pct'] - proxy_comp['emp_overall_pct']

    fig, ax = plt.subplots(figsize=(COL2_WIDTH, 3.6))
    y = np.arange(len(categories))
    height = 0.28

    bars1 = ax.barh(y + height/2, empirical_values, height,
                    label='Empirical (Gaia DR3 flags)', color=C_EMPIRICAL,
                    edgecolor='white', linewidth=0.5)
    bars2 = ax.barh(y - height/2, synthetic_values, height,
                    label='Synthetic proxy', color=C_SYNTHETIC,
                    edgecolor='white', linewidth=0.5)

    max_bar_val = max(max(empirical_values), max(synthetic_values))
    label_offset = max_bar_val * 0.02

    for bar, val in zip(bars1, empirical_values):
        w = bar.get_width()
        ax.text(w + label_offset, bar.get_y() + bar.get_height()/2,
                f'{val:.1f}%', ha='left', va='center',
                fontsize=8, fontweight='bold', color=C_EMPIRICAL)

    for bar, val in zip(bars2, synthetic_values):
        w = bar.get_width()
        ax.text(w + label_offset, bar.get_y() + bar.get_height()/2,
                f'{val:.1f}%', ha='left', va='center',
                fontsize=8, fontweight='bold', color=C_SYNTHETIC)

    max_label_right = max_bar_val + max_bar_val * 0.13
    ratio_x = max_label_right + max_bar_val * 0.10

    ax.text(ratio_x, y[-1] - 0.48, 'Syn/Emp\nratio', ha='center', va='center',
            fontsize=6, color='0.3', fontweight='bold', linespacing=1.2)

    for i, (emp, syn) in enumerate(zip(empirical_values, synthetic_values)):
        if emp > 0 and syn > 0:
            ratio = syn / emp
            ratio_text = f'{ratio:.1f}x' if ratio >= 1 else f'{ratio:.2f}x'
        elif syn > 0:
            ratio_text = '—'
        else:
            ratio_text = '0x'
        ax.text(ratio_x, y[i], ratio_text,
                ha='center', va='center', fontsize=7.5, color='0.4', style='italic')

    ax.set_yticks(y)
    ax.set_yticklabels(categories, fontsize=9)
    ax.set_xlabel('Exclusion Rate (%)', fontsize=9)
    ax.set_xlim(0, ratio_x + max_bar_val * 0.08)
    ax.invert_yaxis()

    ax.legend(loc='lower center', ncol=2, fontsize=7.5,
              framealpha=0.95, edgecolor='0.8', borderpad=0.4,
              handletextpad=0.5, columnspacing=1.0, bbox_to_anchor=(0.5, 1.02))

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='x', color=C_GRID, linewidth=0.4)

    direction = 'overestimate' if diff > 0 else 'underestimate'
    fig.text(0.12, 0.01,
             f'Note: Synthetic proxies {direction} overall exclusion by {abs(diff):.1f} percentage points.',
             fontsize=7, color='0.35')

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, bbox_inches='tight', facecolor='white')
        print(f"Saved: {save_path}")
    plt.close(fig)
    return fig, ax


# ---
#  Figure 5: Sky Distribution (Mollweide) — DENSITY PLOT
# ---
def figure5_sky_distribution(catalog: pd.DataFrame,
                              save_path: Optional[str] = None):
    fig = plt.figure(figsize=(COL2_WIDTH, 4.0))
    ax = fig.add_subplot(111, projection='mollweide')

    if 'ra' not in catalog.columns or 'dec' not in catalog.columns:
        ax.text(0.5, 0.5, 'RA/Dec data not available',
                transform=ax.transAxes, ha='center', va='center', fontsize=10)
        if save_path:
            plt.savefig(save_path, bbox_inches='tight', facecolor='white')
            print(f"Saved: {save_path}")
        plt.close(fig)
        return fig, ax

    retained = catalog[catalog['decision'] == 'RETAIN']
    excluded = catalog[catalog['decision'] == 'EXCLUDE']

    n_bins_ra = 120
    n_bins_dec = 60

    ra_excl = np.radians(180 - excluded['ra'].values)
    dec_excl = np.radians(excluded['dec'].values)
    H_excl, xedges, yedges = np.histogram2d(
        ra_excl, dec_excl, bins=[n_bins_ra, n_bins_dec],
        range=[[-np.pi, np.pi], [-np.pi/2, np.pi/2]])
    from scipy.ndimage import gaussian_filter
    H_excl_smooth = gaussian_filter(H_excl.T, sigma=1.2)

    ra_ret = np.radians(180 - retained['ra'].values)
    dec_ret = np.radians(retained['dec'].values)
    H_ret, _, _ = np.histogram2d(
        ra_ret, dec_ret, bins=[n_bins_ra, n_bins_dec],
        range=[[-np.pi, np.pi], [-np.pi/2, np.pi/2]])
    H_ret_smooth = gaussian_filter(H_ret.T, sigma=1.2)

    xcenter = 0.5 * (xedges[:-1] + xedges[1:])
    ycenter = 0.5 * (yedges[:-1] + yedges[1:])

    levels_excl = np.linspace(H_excl_smooth.max()*0.05, H_excl_smooth.max()*0.9, 6)
    ax.contourf(xcenter, ycenter, H_excl_smooth, levels=levels_excl,
                cmap='Oranges', alpha=0.4, extend='max')

    levels_ret = np.linspace(H_ret_smooth.max()*0.05, H_ret_smooth.max()*0.9, 6)
    ax.contourf(xcenter, ycenter, H_ret_smooth, levels=levels_ret,
                cmap='GnBu', alpha=0.6, extend='max')
    ax.contour(xcenter, ycenter, H_ret_smooth, levels=levels_ret[2:],
               colors=C_RETAINED, linewidths=0.4, alpha=0.6)

    tick_positions = np.radians(np.arange(-150, 180, 30))
    tick_labels = ['22h', '20h', '18h', '16h', '14h', '12h', '10h', '8h', '6h', '4h', '2h']
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels, fontsize=7)
    ax.set_xlabel('Right Ascension', fontsize=8)
    ax.set_ylabel('Declination', fontsize=8)
    ax.tick_params(axis='y', labelsize=7)
    ax.grid(True, color=C_GRID, linewidth=0.3)

    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=C_EXCLUDED, alpha=0.4, label='Excluded (density)'),
        Patch(facecolor=C_RETAINED, alpha=0.6, label='Retained (density)'),
    ]
    ax.legend(handles=legend_elements, loc='lower left', fontsize=7,
              framealpha=0.9, edgecolor='0.7')

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, bbox_inches='tight', facecolor='white')
        print(f"Saved: {save_path}")
    plt.close(fig)
    return fig, ax


# ---
#  Figure 6: Sensitivity Analysis
# ---
def figure6_sensitivity_analysis(sensitivity_results: Dict[str, pd.DataFrame],
                                  baseline_exclusion: float,
                                  save_path: Optional[str] = None):
    fig, axes = plt.subplots(2, 2, figsize=(COL2_WIDTH, COL2_WIDTH * 0.85))

    configs = [
        ('age',         C_AGE,   'Age Threshold (Gyr)',           3.0,  r'Age'),
        ('mass',        C_MASS,  r'Mass Threshold (M$_\odot$)',   1.5,  r'Mass'),
        ('metallicity', C_METALL,'[Fe/H] Threshold',             -0.4,  'Metallicity'),
        ('variability', C_VAR,   'Variability Threshold (mag)',   0.01,  'Variability'),
    ]

    for ax, (param, color, xlabel, baseline, title) in zip(axes.flat, configs):
        df = sensitivity_results[param]
        y_vals = 100 * df['excluded_fraction']

        ax.plot(df['threshold'], y_vals, 'o-', color=color, markersize=3, linewidth=1.2)
        ax.axhline(y=100 * baseline_exclusion, color='0.5', linestyle=':',
                   linewidth=0.8, label=f'Baseline ({100*baseline_exclusion:.1f}%)')
        ax.axvline(x=baseline, color='0.5', linestyle='--', linewidth=0.8,
                   alpha=0.7, label=f'Threshold = {baseline}')

        ax.set_xlabel(xlabel, fontsize=8)
        ax.set_ylabel('Exclusion Rate (%)', fontsize=8)
        ax.set_title(title, fontsize=9, fontweight='bold')
        ax.legend(fontsize=6, loc='best')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(color=C_GRID, linewidth=0.3)
        ax.tick_params(labelsize=7)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, bbox_inches='tight', facecolor='white')
        print(f"Saved: {save_path}")
    plt.close(fig)
    return fig, axes


# ---
#  FOV Analysis Figure
# ---
def figure_fov_analysis(catalog: pd.DataFrame,
                         save_path: Optional[str] = None):
    from analysis import analyze_all_telescopes

    fig, axes = plt.subplots(1, 2, figsize=(COL2_WIDTH, 3.2))

    ax = axes[0]
    telescope_df = analyze_all_telescopes(catalog)
    telescopes = telescope_df['Telescope'].values
    expected = telescope_df['Expected Stars'].values
    p_values = telescope_df['P(>=1 star)'].values

    colors = [C_EXCLUDED, C_EMPIRICAL, C_RETAINED, C_VAR]
    bars = ax.bar(telescopes, 100 * p_values, color=colors,
                  edgecolor='white', linewidth=0.5, width=0.6)

    for bar, exp, p in zip(bars, expected, p_values):
        ax.annotate(f'{100*p:.1f}%\n($\\lambda$={exp:.2f})',
                    xy=(bar.get_x() + bar.get_width()/2, bar.get_height()),
                    xytext=(0, 3), textcoords="offset points",
                    ha='center', va='bottom', fontsize=7, fontweight='bold')

    ax.set_ylabel('P($\\geq$1 retained star)', fontsize=8)
    ax.set_ylim(0, 115)
    ax.tick_params(labelsize=7)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', color=C_GRID, linewidth=0.3)

    ax = axes[1]
    if 'parallax' in catalog.columns:
        retained = catalog[catalog['decision'] == 'RETAIN'].copy()
        retained['distance_pc'] = 1000.0 / retained['parallax']
        valid_dist = retained[retained['distance_pc'] > 0]['distance_pc']

        # Limit x-axis to ~3000 pc for cleaner view
        upper_limit = 3000
        ax.hist(valid_dist[valid_dist <= upper_limit], bins=60,
                color=C_EMPIRICAL, edgecolor='none', alpha=0.8)
        ax.axvline(valid_dist.median(), color=C_SYNTHETIC, linestyle='--',
                   linewidth=1.5, label=f'Median: {valid_dist.median():.0f} pc')
        ax.set_xlim(0, upper_limit)
        ax.legend(loc='best', fontsize=7)

    ax.set_xlabel('Distance (pc)', fontsize=8)
    ax.set_ylabel('Retained Stars', fontsize=8)
    ax.tick_params(labelsize=7)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', color=C_GRID, linewidth=0.3)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, bbox_inches='tight', facecolor='white')
        print(f"Saved: {save_path}")
    plt.close(fig)
    return fig, axes


# ---
#  Figure 7: GCNS Comparison
# ---
def figure7_gcns_comparison(catalog: pd.DataFrame,
                             save_path: Optional[str] = None):
    """GCNS retention comparison figure (Howard feedback #1)."""
    if 'in_gcns' not in catalog.columns or not catalog['in_gcns'].any():
        print("  Skipping GCNS comparison figure (no GCNS data)")
        return None, None

    fig, axes = plt.subplots(1, 2, figsize=(COL2_WIDTH, 3.2))

    # Panel 1: Retention rate comparison
    ax = axes[0]
    gcns_cat = catalog[catalog['in_gcns']]
    full_cat = catalog

    gcns_ret = len(gcns_cat[gcns_cat['decision'] == 'RETAIN'])
    gcns_total = len(gcns_cat)
    full_ret = len(full_cat[full_cat['decision'] == 'RETAIN'])
    full_total = len(full_cat)

    labels = ['Full Catalog', 'GCNS\n(<=100 pc)']
    retention_pcts = [100*full_ret/full_total, 100*gcns_ret/gcns_total]
    colors = [C_EMPIRICAL, C_RETAINED]

    bars = ax.bar(labels, retention_pcts, color=colors, edgecolor='white', linewidth=0.5, width=0.55)
    for bar, pct, n_ret, n_tot in zip(bars, retention_pcts,
                                        [full_ret, gcns_ret], [full_total, gcns_total]):
        ax.annotate(f'{pct:.1f}%\n({n_ret:,}/{n_tot:,})',
                    xy=(bar.get_x() + bar.get_width()/2, bar.get_height()),
                    xytext=(0, 3), textcoords="offset points",
                    ha='center', va='bottom', fontsize=7.5, fontweight='bold')

    ax.set_ylabel('Retention Rate (%)', fontsize=8)
    ax.set_ylim(0, max(retention_pcts) * 1.3)
    ax.tick_params(labelsize=7)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', color=C_GRID, linewidth=0.3)

    # Panel 2: GCNS exclusion reasons
    ax = axes[1]
    gcns_excl = gcns_cat[gcns_cat['decision'] == 'EXCLUDE']
    codes = ['R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7']
    counts = [len(gcns_excl[gcns_excl['reason_code'].str.contains(c, na=False)]) for c in codes]

    bars = ax.bar(codes, counts, color=CRITERIA_COLORS, edgecolor='white', linewidth=0.5, width=0.6)
    for bar, val in zip(bars, counts):
        if val > 0:
            ax.annotate(f'{val}',
                        xy=(bar.get_x() + bar.get_width()/2, bar.get_height()),
                        xytext=(0, 3), textcoords="offset points",
                        ha='center', va='bottom', fontsize=7)

    ax.set_xlabel('Exclusion Criterion', fontsize=8)
    ax.set_ylabel('GCNS Stars Excluded', fontsize=8)
    ax.tick_params(labelsize=7)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', color=C_GRID, linewidth=0.3)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, bbox_inches='tight', facecolor='white')
        print(f"Saved: {save_path}")
    plt.close(fig)
    return fig, axes


# ---
#  Figure 8: Parallax Quality
# ---
def figure8_parallax_quality(catalog: pd.DataFrame,
                              save_path: Optional[str] = None):
    """Parallax quality figure (Howard feedback #3)."""
    from analysis import analyze_parallax_quality

    poe_stats = analyze_parallax_quality(catalog)
    if 'error' in poe_stats:
        print("  Skipping parallax quality figure")
        return None, None

    fig, axes = plt.subplots(1, 2, figsize=(COL2_WIDTH, 3.2))

    # Panel 1: CDF of parallax_over_error
    ax = axes[0]
    poe = catalog['parallax_over_error'].dropna()
    sorted_poe = np.sort(poe.values)
    cdf = np.arange(1, len(sorted_poe)+1) / len(sorted_poe)

    ax.plot(sorted_poe, cdf, color=C_EMPIRICAL, linewidth=1.2,
               label=r'CDF of $\varpi/\sigma_\varpi$')
    ax.axvline(x=5, color=C_SYNTHETIC, linestyle='--', linewidth=1.0,
               label=r'$\varpi/\sigma_\varpi = 5$')
    ax.axvline(x=10, color='#F4A261', linestyle=':', linewidth=1.0,
               label=r'$\varpi/\sigma_\varpi = 10$')

    # Shade region where naive 1/pi is reliable
    ax.axvspan(5, sorted_poe.max(), alpha=0.1, color=C_RETAINED,
               label=r'Reliable $d = 1/\varpi$')

    ax.set_xlabel(r'$\varpi/\sigma_\varpi$', fontsize=9)
    ax.set_ylabel('Cumulative Fraction', fontsize=9)
    ax.set_xlim(0, min(500, sorted_poe.max()))
    ax.legend(loc='lower right', fontsize=6.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(color=C_GRID, linewidth=0.3)

    # Panel 2: Summary statistics bar chart
    ax = axes[1]
    thresholds = [5, 10, 20, 50]
    pcts = [poe_stats['pct_gt_5'], poe_stats['pct_gt_10'],
            poe_stats['pct_gt_20'], poe_stats['pct_gt_50']]

    bars = ax.bar([str(t) for t in thresholds], pcts,
                  color=[C_RETAINED, C_EMPIRICAL, '#F4A261', C_EXCLUDED],
                  edgecolor='white', linewidth=0.5, width=0.6)
    for bar, pct in zip(bars, pcts):
        ax.annotate(f'{pct:.1f}%',
                    xy=(bar.get_x() + bar.get_width()/2, bar.get_height()),
                    xytext=(0, 3), textcoords="offset points",
                    ha='center', va='bottom', fontsize=8, fontweight='bold')

    ax.set_xlabel(r'Threshold $\varpi/\sigma_\varpi$', fontsize=9)
    ax.set_ylabel('% of Stars Above Threshold', fontsize=9)
    ax.set_ylim(0, 115)
    ax.tick_params(labelsize=7)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', color=C_GRID, linewidth=0.3)

    fig.text(0.12, 0.01,
             r'With $\varpi/\sigma_\varpi > 5$, naive $d = 1/\varpi$ and Bayesian distances agree within 5% (Bailer-Jones+ 2021).',
             fontsize=6.5, color='0.4')

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, bbox_inches='tight', facecolor='white')
        print(f"Saved: {save_path}")
    plt.close(fig)
    return fig, axes


# ---
#  BL Cross-Match Overview Figure
# ---
def figure_bl_crossmatch(bl_results: dict,
                          save_path: Optional[str] = None):
    """
    Breakthrough Listen cross-match summary figure.

    Shows matched/excluded/retained counts for each catalog
    in separate panels (vastly different sample sizes).
    """
    isaacson = bl_results.get('isaacson')
    meerkat = bl_results.get('meerkat')

    if not isaacson and not meerkat:
        print("  Skipping BL cross-match figure (no BL data)")
        return None, None

    has_both = isaacson is not None and meerkat is not None

    if has_both:
        fig, axes = plt.subplots(1, 2, figsize=(COL2_WIDTH, 3.5))

        # - Panel 1: Isaacson et al. (2017) summary
        ax = axes[0]
        isa_match = isaacson['matched_count']
        isa_excl = isaacson['excluded_count']
        isa_ret  = isaacson['retained_count']

        bars = ax.bar(['Matched', 'Excluded', 'Retained'],
                     [isa_match, isa_excl, isa_ret],
                     color=['#95a5a6', C_EXCLUDED, C_RETAINED],
                     edgecolor='white', linewidth=0.5, width=0.55)
        for bar, val in zip(bars, [isa_match, isa_excl, isa_ret]):
            pct = 100 * val / isa_match
            ax.annotate(f'{val:,}\n({pct:.1f}%)',
                        xy=(bar.get_x() + bar.get_width()/2, bar.get_height()),
                        xytext=(0, 3), textcoords="offset points",
                        ha='center', va='bottom', fontsize=8, fontweight='bold')
        ax.set_title('Isaacson et al. (2017)\nPASP 129, 054501', fontsize=8.5, fontweight='bold')
        ax.set_ylabel('Number of Stars', fontsize=8)
        ax.tick_params(labelsize=7)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', color=C_GRID, linewidth=0.3)

        # - Panel 2: MeerKAT summary
        ax = axes[1]
        mkat_match = meerkat['matched_count']
        mkat_excl = meerkat['excluded_count']
        mkat_ret  = meerkat['retained_count']

        bars = ax.bar(['Matched', 'Excluded', 'Retained'],
                       [mkat_match, mkat_excl, mkat_ret],
                       color=['#95a5a6', C_EXCLUDED, C_RETAINED],
                       edgecolor='white', linewidth=0.5, width=0.55)
        for bar, val in zip(bars, [mkat_match, mkat_excl, mkat_ret]):
            pct = 100 * val / mkat_match
            val_str = f'{val:,}' if val < 1e6 else f'{val/1e6:.1f}M'
            ax.annotate(f'{val_str}\n({pct:.1f}%)',
                        xy=(bar.get_x() + bar.get_width()/2, bar.get_height()),
                        xytext=(0, 3), textcoords="offset points",
                        ha='center', va='bottom', fontsize=8, fontweight='bold')
        ax.set_title('MeerKAT 1M\nCzech et al. (2021)', fontsize=8.5, fontweight='bold')
        ax.tick_params(labelsize=7)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', color=C_GRID, linewidth=0.3)

    else:
        # Single catalog: 1-panel layout
        fig, ax = plt.subplots(figsize=(COL_WIDTH, 3.5))
        overlap = isaacson if isaacson else meerkat

        labels = ['Matched\nin Gaia', 'Excluded\nby Model', 'Retained\nby Model']
        values = [overlap['matched_count'], overlap['excluded_count'], overlap['retained_count']]
        colors = ['#95a5a6', C_EXCLUDED, C_RETAINED]

        bars = ax.bar(labels, values, color=colors, edgecolor='white', linewidth=0.5, width=0.55)
        for bar, val in zip(bars, values):
            val_str = f'{val:,}' if val < 1e6 else f'{val/1e6:.1f}M'
            ax.annotate(val_str,
                        xy=(bar.get_x() + bar.get_width()/2, bar.get_height()),
                        xytext=(0, 3), textcoords="offset points",
                        ha='center', va='bottom', fontsize=8, fontweight='bold')

        ax.set_ylabel('Number of Stars', fontsize=8)
        ax.tick_params(labelsize=7)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', color=C_GRID, linewidth=0.3)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, bbox_inches='tight', facecolor='white')
        print(f"Saved: {save_path}")
    plt.close(fig)
    return fig, axes if has_both else (fig, ax)


# ---
#  BL Cross-Match Per-Criterion Figures (2 separate figures)
# ---
def figure_bl_criteria_isaacson(bl_results: dict,
                                 save_path: Optional[str] = None):
    """
    Per-criterion exclusion breakdown for Isaacson et al. (2017).
    Separate figure so the small sample (N~1700) is fully visible.
    """
    isaacson = bl_results.get('isaacson')
    if not isaacson:
        return None, None

    codes = ['R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7']
    code_labels = {
        'R1': 'Mass',
        'R2': 'Age',
        'R3': 'Spectral',
        'R4': '[Fe/H]',
        'R5': 'Multipl.',
        'R6': 'Variab.',
        'R7': 'M-dwarf',
    }
    counts = [isaacson['reason_counts'].get(c, 0) for c in codes]
    total_excl = isaacson['excluded_count']

    fig, ax = plt.subplots(figsize=(COL2_WIDTH, 3.5))

    bars = ax.bar(codes, counts, color=CRITERIA_COLORS[:7],
                  edgecolor='white', linewidth=0.5, width=0.65)
    for bar, val in zip(bars, counts):
        if val > 0:
            pct = 100 * val / total_excl if total_excl > 0 else 0
            ax.annotate(f'{val}\n({pct:.1f}%)',
                        xy=(bar.get_x() + bar.get_width()/2, bar.get_height()),
                        xytext=(0, 4), textcoords="offset points",
                        ha='center', va='bottom', fontsize=7.5, fontweight='bold')

    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(range(len(codes)))
    ax2.set_xticklabels([code_labels[c] for c in codes],
                        fontsize=7, style='italic', color='0.4')
    ax2.tick_params(axis='x', length=0, pad=10)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    ax2.spines['left'].set_visible(False)

    ax.set_ylabel('Stars Excluded', fontsize=9)
    ax.set_xlabel('Exclusion Criterion', fontsize=8)
    ax.set_ylim(0, max(counts) * 1.28)
    ax.set_title('Isaacson et al. (2017) — Per-Criterion Breakdown',
                 fontsize=9.5, fontweight='bold', pad=18)
    ax.tick_params(labelsize=8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', color=C_GRID, linewidth=0.4)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, bbox_inches='tight', facecolor='white')
        print(f"Saved: {save_path}")
    plt.close(fig)
    return fig, ax


def figure_bl_criteria_meerkat(bl_results: dict,
                                save_path: Optional[str] = None):
    """
    Per-criterion exclusion breakdown for MeerKAT 1M catalog.
    Separate figure so the large sample (N~98k) is fully visible.
    """
    meerkat = bl_results.get('meerkat')
    if not meerkat:
        return None, None

    codes = ['R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7']
    code_labels = {
        'R1': 'Mass',
        'R2': 'Age',
        'R3': 'Spectral',
        'R4': '[Fe/H]',
        'R5': 'Multipl.',
        'R6': 'Variab.',
        'R7': 'M-dwarf',
    }
    counts = [meerkat['reason_counts'].get(c, 0) for c in codes]
    total_excl = meerkat['excluded_count']

    fig, ax = plt.subplots(figsize=(COL2_WIDTH, 3.5))

    bars = ax.bar(codes, counts, color=CRITERIA_COLORS[:7],
                  edgecolor='white', linewidth=0.5, width=0.65)
    for bar, val in zip(bars, counts):
        if val > 0:
            pct = 100 * val / total_excl if total_excl > 0 else 0
            val_str = f'{val:,}' if val < 1e6 else f'{val/1e6:.1f}M'
            ax.annotate(f'{val_str}\n({pct:.1f}%)',
                        xy=(bar.get_x() + bar.get_width()/2, bar.get_height()),
                        xytext=(0, 4), textcoords="offset points",
                        ha='center', va='bottom', fontsize=7.5, fontweight='bold')

    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(range(len(codes)))
    ax2.set_xticklabels([code_labels[c] for c in codes],
                        fontsize=7, style='italic', color='0.4')
    ax2.tick_params(axis='x', length=0, pad=10)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    ax2.spines['left'].set_visible(False)

    ax.yaxis.set_major_formatter(mticker.FuncFormatter(
        lambda x, _: f'{x/1000:.0f}k' if x >= 1000 else f'{x:.0f}'))
    ax.set_ylabel('Stars Excluded', fontsize=9)
    ax.set_xlabel('Exclusion Criterion', fontsize=8)
    ax.set_ylim(0, max(counts) * 1.28)
    ax.set_title('MeerKAT 1M (Czech et al. 2021) — Per-Criterion Breakdown',
                 fontsize=9.5, fontweight='bold', pad=18)
    ax.tick_params(labelsize=8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', color=C_GRID, linewidth=0.4)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, bbox_inches='tight', facecolor='white')
        print(f"Saved: {save_path}")
    plt.close(fig)
    return fig, ax


# ---
#  Age Comparison Figure
# ---
def figure_age_comparison(catalog: pd.DataFrame,
                           save_path: Optional[str] = None):
    from analysis import compare_age_criteria

    comparison = compare_age_criteria(catalog)

    fig, axes = plt.subplots(1, 2, figsize=(COL2_WIDTH, 3.2))

    ax = axes[0]
    valid_ages = catalog[catalog['age_flame_spec'] > 0]['age_flame_spec']

    ax.hist(valid_ages, bins=100, color=C_EMPIRICAL, edgecolor='none', alpha=0.75)
    ax.axvline(x=3.0, color=C_SYNTHETIC, linestyle='--', linewidth=1.5,
               label='Threshold (3 Gyr)')
    ax.axvspan(2.0, 4.0, alpha=0.15, color='#F4A261', label='Boundary zone (2–4 Gyr)')

    ax.set_xlabel('Stellar Age (Gyr)', fontsize=8)
    ax.set_ylabel('Number of Stars', fontsize=8)
    ax.legend(loc='best', fontsize=6.5)
    ax.set_xlim(0, 15)
    ax.tick_params(labelsize=7)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', color=C_GRID, linewidth=0.3)

    ax = axes[1]
    labels = ['Hard Cut\n(age < 3 Gyr)', 'Uncertainty-Aware\n(age$_{upper}$ < 3 Gyr)']
    values = [comparison['old_r2_count'], comparison['new_r2_count']]
    colors = [C_EXCLUDED, C_RETAINED]

    bars = ax.bar(labels, values, color=colors, edgecolor='white', linewidth=0.5, width=0.55)
    for bar, val in zip(bars, values):
        ax.annotate(f'{val:,}',
                    xy=(bar.get_x() + bar.get_width()/2, bar.get_height()),
                    xytext=(0, 3), textcoords="offset points",
                    ha='center', va='bottom', fontsize=9, fontweight='bold')

    diff = comparison['difference']
    delta_pct = 100 * diff / comparison['old_r2_count'] if comparison['old_r2_count'] > 0 else 0
    ax.annotate(f'$\\Delta$ = +{diff:,} ({delta_pct:.1f}% fewer excluded)',
                xy=(0.5, 0.88), xycoords='axes fraction',
                ha='center', fontsize=7.5, fontweight='bold', color=C_RETAINED,
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor=C_RETAINED, alpha=0.85))

    ax.set_ylabel('Stars Excluded by R2', fontsize=8)
    ax.tick_params(labelsize=7)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', color=C_GRID, linewidth=0.3)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, bbox_inches='tight', facecolor='white')
        print(f"Saved: {save_path}")
    plt.close(fig)
    return fig, axes


# ---
#  Generate All Figures
# ---
def generate_all_figures(catalog: pd.DataFrame,
                          sensitivity_results: dict,
                          bl_results: dict = None):
    """Generate all publication-ready figures.

    Parameters
    ----------
    catalog : pd.DataFrame
        Torlakcık catalog with decision/reason_code columns
    sensitivity_results : dict
        Sensitivity analysis results
    bl_results : dict, optional
        BL cross-match results: {'isaacson': dict|None, 'meerkat': dict|None}
    """
    if bl_results is None:
        bl_results = {}

    print("\n" + "=" * 60)
    print("GENERATING ALL FIGURES (Publication-Ready)")
    print("=" * 60)

    baseline_exclusion = len(catalog[catalog['decision'] == 'EXCLUDE']) / len(catalog)

    print("\n  Figure 2: Exclusion by criterion...")
    figure2_exclusion_by_criterion(catalog,
                                    os.path.join(FIGURES_DIR, 'figure2_exclusion_by_criterion.png'))

    print("  Figure 3: Spectral type breakdown...")
    figure3_spectral_type_breakdown(catalog,
                                     os.path.join(FIGURES_DIR, 'figure3_spectral_breakdown.png'))

    print("  HR Diagram...")
    figure_hr_diagram(catalog,
                       os.path.join(FIGURES_DIR, 'figure_hr_diagram.png'))

    print("  Figure 4: Empirical vs Synthetic...")
    figure4_empirical_vs_synthetic(catalog,
                                    os.path.join(FIGURES_DIR, 'figure4_empirical_vs_synthetic.png'))

    print("  Figure 5: Sky distribution...")
    figure5_sky_distribution(catalog,
                              os.path.join(FIGURES_DIR, 'figure5_sky_distribution.png'))

    # Figure 6: Sensitivity (only if data available)
    if sensitivity_results and len(sensitivity_results) == 4:
        print("  Figure 6: Sensitivity analysis...")
        figure6_sensitivity_analysis(sensitivity_results, baseline_exclusion,
                                      os.path.join(FIGURES_DIR, 'figure6_sensitivity_analysis.png'))

    print("  Figure 7: GCNS comparison...")
    figure7_gcns_comparison(catalog,
                             os.path.join(FIGURES_DIR, 'figure7_gcns_comparison.png'))

    print("  Figure 8: Parallax quality...")
    figure8_parallax_quality(catalog,
                              os.path.join(FIGURES_DIR, 'figure8_parallax_quality.png'))

    # BL Cross-match figures
    if bl_results.get('isaacson') or bl_results.get('meerkat'):
        print("  BL Cross-match overview...")
        figure_bl_crossmatch(bl_results,
                              os.path.join(FIGURES_DIR, 'figure_bl_crossmatch.png'))
        if bl_results.get('isaacson'):
            print("  BL Isaacson per-criterion breakdown...")
            figure_bl_criteria_isaacson(bl_results,
                                         os.path.join(FIGURES_DIR, 'figure_bl_criteria_isaacson.png'))
        if bl_results.get('meerkat'):
            print("  BL MeerKAT per-criterion breakdown...")
            figure_bl_criteria_meerkat(bl_results,
                                        os.path.join(FIGURES_DIR, 'figure_bl_criteria_meerkat.png'))

    print("  Age comparison...")
    figure_age_comparison(catalog,
                           os.path.join(FIGURES_DIR, 'figure_age_comparison.png'))

    print("  FOV analysis...")
    figure_fov_analysis(catalog,
                         os.path.join(FIGURES_DIR, 'figure_fov_analysis.png'))

    print(f"\n  All figures saved to: {FIGURES_DIR}")
