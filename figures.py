# SETI Parametric Avoidance Model — Figure Generation
# Input: df (processed DataFrame with exclusion flags)
# Output: PNG figures saved to Desktop/seti_project/

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
import os

save_dir = os.path.join(os.path.expanduser('~'), 'Desktop', 'seti_project')
os.makedirs(save_dir, exist_ok=True)


def plot_criteria_breakdown(df):
    fig, ax = plt.subplots(figsize=(10, 6))
    crit_labels = ['Mass\n>1.5M⊙', 'Age\n<3Gyr', 'Spectral\nO/B/A/F0-4',
                   '[Fe/H]\n<-0.4', 'Multiplicity\n≥1', 'Photometric\nVariability', 'Active\nM Dwarf']
    crit_values = [
        df['exclude_mass'].sum(), df['exclude_age'].sum(),
        df['exclude_spec'].sum(), df['exclude_feh'].sum(),
        df['exclude_mult'].sum(), df['exclude_var'].sum(),
        df['exclude_act'].sum()
    ]
    colors = ['#d62728','#ff7f0e','#ffd700','#2ca02c','#1f77b4','#9467bd','#8c564b']
    bars = ax.bar(crit_labels, crit_values, color=colors, edgecolor='none')
    for bar, val in zip(bars, crit_values):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 5000,
                f'{val:,}', ha='center', va='bottom', fontsize=9, fontweight='bold')
    ax.set_ylabel('Number of Stars Excluded')
    ax.set_title(f'Stars Excluded per Avoidance Criterion\n(Gaia DR3, N = {len(df):,})')
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{int(x):,}'))
    ax.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, 'figure1_criteria_breakdown.png'), dpi=200, bbox_inches='tight')
    plt.close()
    print("Saved: figure1_criteria_breakdown.png")


def plot_spectype_breakdown(df):
    spec_order = ['O','B','A','F0-F4','F5-F9','G','K','M']
    fig, ax = plt.subplots(figsize=(10, 6))
    x = np.arange(len(spec_order))
    width = 0.35
    retained_counts = []
    excluded_counts = []
    for sp in spec_order:
        mask = df['spec_type'] == sp
        retained_counts.append((mask & ~df['excluded']).sum())
        excluded_counts.append((mask & df['excluded']).sum())
    ax.bar(x - width/2, retained_counts, width, label='Retained', color='#2ca02c', alpha=0.85)
    bars_e = ax.bar(x + width/2, excluded_counts, width, label='Excluded', color='#d62728', alpha=0.85)
    for bar, exc, ret in zip(bars_e, excluded_counts, retained_counts):
        total_sp = exc + ret
        if total_sp > 0:
            pct = exc / total_sp * 100
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1000,
                    f'{pct:.0f}%', ha='center', va='bottom', fontsize=8,
                    color='#d62728', fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(spec_order)
    ax.set_ylabel('Number of Stars')
    ax.set_title(f'Retained vs Excluded Stars by Spectral Type\n(Gaia DR3, N = {len(df):,})')
    ax.legend()
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{int(x):,}'))
    ax.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, 'figure2_spectype_breakdown.png'), dpi=200, bbox_inches='tight')
    plt.close()
    print("Saved: figure2_spectype_breakdown.png")


def plot_synthetic_vs_empirical():
    fig, ax = plt.subplots(figsize=(9, 6))
    labels = ['Synthetic\n(N=10,000)', 'Synthetic\n(N=1,754,135)', 'Empirical Gaia DR3\n(N=1,754,135)']
    values = [65.7, 71.1, 66.2]
    colors = ['#aec7e8', '#7fbc8c', '#1f77b4']
    bars = ax.bar(labels, values, color=colors, edgecolor='none', width=0.5)
    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.2,
                f'{val}%', ha='center', va='bottom', fontsize=13, fontweight='bold')
    ax.axhline(66.2, color='#1f77b4', linestyle='--', linewidth=1.2, alpha=0.6)
    ax.set_ylabel('Exclusion Rate (%)')
    ax.set_ylim(60, 74)
    ax.set_title('Exclusion Rate: Synthetic vs Empirical Parameters\nfor Multiplicity and Photometric Variability')
    ax.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, 'figure3_synthetic_vs_empirical.png'), dpi=200, bbox_inches='tight')
    plt.close()
    print("Saved: figure3_synthetic_vs_empirical.png")


def plot_skymap(df):
    fig = plt.figure(figsize=(18, 9), facecolor='white')
    ax = fig.add_subplot(111, projection='mollweide')
    ax.set_facecolor('white')
    ra_rad  = np.deg2rad(df['ra'].values - 180)
    dec_rad = np.deg2rad(df['dec'].values)
    excluded_mask = df['excluded'].values
    retained_mask = ~excluded_mask
    excl_count = excluded_mask.sum()
    ret_count  = retained_mask.sum()
    total = len(df)
    ax.scatter(ra_rad[excluded_mask], dec_rad[excluded_mask],
               s=0.08, c='#d62728', alpha=0.2, linewidths=0, rasterized=True)
    ax.scatter(ra_rad[retained_mask], dec_rad[retained_mask],
               s=0.12, c='#1f77b4', alpha=0.4, linewidths=0, rasterized=True)
    ax.grid(True, color='gray', alpha=0.25, linewidth=0.5)
    ax.set_xticklabels(['210°','240°','270°','300°','330°','0°','30°','60°','90°','120°','150°'], fontsize=9)
    legend_elements = [
        Line2D([0],[0], marker='o', color='w', markerfacecolor='#d62728',
               markersize=10,
               label=f'Excluded (N = {excl_count:,} | {excl_count/total*100:.1f}%)',
               linewidth=0),
        Line2D([0],[0], marker='o', color='w', markerfacecolor='#1f77b4',
               markersize=10,
               label=f'Retained (N = {ret_count:,} | {ret_count/total*100:.1f}%)',
               linewidth=0),
    ]
    ax.legend(handles=legend_elements, loc='lower right',
              facecolor='white', edgecolor='black', fontsize=11)
    ax.set_title(f'Sky Distribution of SETI Candidate Stars — Gaia DR3\n'
                 f'Parametric Avoidance Model  |  N = {total:,}',
                 fontsize=14, fontweight='bold', pad=20)
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, 'skymap.png'), dpi=200, bbox_inches='tight', facecolor='white')
    plt.close()
    print("Saved: skymap.png")


if __name__ == '__main__':
    print("Run analysis.py first to generate df, then call these functions.")
