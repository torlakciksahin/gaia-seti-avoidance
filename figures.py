# SETI Parametric Avoidance Model — Figure Generation
# Run after analysis.py or after loading df from CSV
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import numpy as np
from matplotlib.lines import Line2D
import os

save_dir = r'C:\Users\PC\OneDrive\Desktop\avoidance model'
os.makedirs(save_dir, exist_ok=True)


# FIGURE 1 — Stars Excluded per Avoidance Criterion

def plot_criteria_breakdown(df):
    fig, ax = plt.subplots(figsize=(11, 6))
    crit_labels = [
        'Mass\n>1.5M⊙', 'Age\n<3Gyr', 'Spectral\nO/B/A/F0-4',
        '[Fe/H]\n<−0.4', 'Multiplicity\n≥1 (NSS)', 'Photometric\nVariability', 'Active\nM Dwarf'
    ]
    crit_values = [
        df['exclude_mass'].sum(), df['exclude_age'].sum(),
        df['exclude_spec'].sum(), df['exclude_feh'].sum(),
        df['exclude_mult'].sum(), df['exclude_var'].sum(),
        df['exclude_act'].sum()
    ]
    colors = ['#d62728', '#ff7f0e', '#ffd700', '#2ca02c', '#1f77b4', '#9467bd', '#8c564b']
    bars = ax.bar(crit_labels, crit_values, color=colors, edgecolor='none', width=0.6)
    for bar, val in zip(bars, crit_values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 8000,
                f'{val:,}', ha='center', va='bottom', fontsize=9, fontweight='bold')
    ax.set_ylabel('Number of Stars Excluded', fontsize=11)
    ax.set_title(f'Stars Excluded per Avoidance Criterion\n(Gaia DR3, N = {len(df):,})', fontsize=12, fontweight='bold')
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{int(x):,}'))
    ax.grid(axis='y', alpha=0.3)
    ax.set_ylim(0, max(crit_values) * 1.15)
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, 'figure1_criteria_breakdown.png'), dpi=200, bbox_inches='tight')
    plt.close()
    print("Saved: figure1_criteria_breakdown.png")



# FIGURE 2 — Retained vs Excluded by Spectral Type

def plot_spectype_breakdown(df):
    spec_order = ['O', 'B', 'A', 'F0-F4', 'F5-F9', 'G', 'K', 'M']
    fig, ax = plt.subplots(figsize=(11, 6))
    x = np.arange(len(spec_order))
    width = 0.35
    retained_counts, excluded_counts = [], []
    for sp in spec_order:
        mask = df['spec_type'] == sp
        retained_counts.append((mask & ~df['excluded']).sum())
        excluded_counts.append((mask & df['excluded']).sum())
    ax.bar(x - width / 2, retained_counts, width, label='Retained', color='#2ca02c', alpha=0.88)
    bars_e = ax.bar(x + width / 2, excluded_counts, width, label='Excluded', color='#d62728', alpha=0.88)
    for bar, exc, ret in zip(bars_e, excluded_counts, retained_counts):
        total_sp = exc + ret
        if total_sp > 0:
            pct = exc / total_sp * 100
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1000,
                    f'{pct:.0f}%', ha='center', va='bottom', fontsize=8,
                    color='#d62728', fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(spec_order, fontsize=10)
    ax.set_ylabel('Number of Stars', fontsize=11)
    ax.set_title(f'Retained vs Excluded Stars by Spectral Type\n(Gaia DR3, N = {len(df):,})', fontsize=12, fontweight='bold')
    ax.legend(fontsize=10)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{int(x):,}'))
    ax.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, 'figure2_spectype_breakdown.png'), dpi=200, bbox_inches='tight')
    plt.close()
    print("Saved: figure2_spectype_breakdown.png")



# FIGURE 3 — Synthetic vs Empirical Comparison

def plot_synthetic_vs_empirical(df):
    total = len(df)
    excl = df['excluded'].sum()
    empirical_rate = excl / total * 100

    fig, ax = plt.subplots(figsize=(9, 6))
    labels = ['Synthetic\n(N=10,000)', 'Synthetic\n(N=1,754,135)', f'Empirical Gaia DR3\n(N=1,754,135)']
    values = [65.7, 71.1, round(empirical_rate, 1)]
    colors = ['#aec7e8', '#7fbc8c', '#1f77b4']
    bars = ax.bar(labels, values, color=colors, edgecolor='none', width=0.5)
    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.2,
                f'{val}%', ha='center', va='bottom', fontsize=13, fontweight='bold')
    ax.axhline(empirical_rate, color='#1f77b4', linestyle='--', linewidth=1.3, alpha=0.7,
               label=f'Empirical baseline ({empirical_rate:.1f}%)')
    ax.set_ylabel('Exclusion Rate (%)', fontsize=11)
    ax.set_ylim(60, 75)
    ax.set_title('Exclusion Rate: Synthetic vs Empirical Parameters\nfor Multiplicity and Photometric Variability',
                 fontsize=12, fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, 'figure3_synthetic_vs_empirical.png'), dpi=200, bbox_inches='tight')
    plt.close()
    print("Saved: figure3_synthetic_vs_empirical.png")



# FIGURE 4 — Decision Flow Diagram

def plot_decision_flow():
    fig, ax = plt.subplots(figsize=(10, 14))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 14)
    ax.axis('off')
    ax.set_facecolor('white')
    fig.patch.set_facecolor('white')

    box_w, box_h = 5.5, 0.6
    cx = 5.0
    exc_x = 8.5
    exc_w = 2.6

    criteria = [
        ('INPUT: Stellar Parameters', '#2c3e50', 'white', 13.0, False),
        ('Mass > 1.5 M⊙?',           '#2980b9', 'white', 11.5, True),
        ('Age < 3 Gyr?',             '#2980b9', 'white', 10.0, True),
        ('Spectral Type O/B/A/F0-F4?','#2980b9', 'white',  8.5, True),
        ('[Fe/H] < −0.4?',           '#2980b9', 'white',  7.0, True),
        ('non_single_star ≥ 1?',     '#2980b9', 'white',  5.5, True),
        ('Variability > 0.01 mag\nor VARIABLE flag?', '#2980b9', 'white', 4.0, True),
        ('Active M Dwarf?\n(rotation/short-timescale flag)', '#2980b9', 'white', 2.5, True),
    ]

    reason_codes = ['R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7']

    for i, (label, color, tcolor, y, has_exclude) in enumerate(criteria):
        fancy = mpatches.FancyBboxPatch(
            (cx - box_w / 2, y - box_h / 2), box_w, box_h,
            boxstyle='round,pad=0.1', facecolor=color, edgecolor='none', zorder=3
        )
        ax.add_patch(fancy)
        ax.text(cx, y, label, ha='center', va='center',
                color=tcolor, fontsize=9, fontweight='bold', zorder=4,
                multialignment='center')

        if has_exclude:
            rc = reason_codes[i - 1]
            exc_box = mpatches.FancyBboxPatch(
                (exc_x - exc_w / 2, y - box_h / 2), exc_w, box_h,
                boxstyle='round,pad=0.1', facecolor='#e74c3c', edgecolor='none', zorder=3
            )
            ax.add_patch(exc_box)
            ax.text(exc_x, y, f'EXCLUDE\n+ {rc}', ha='center', va='center',
                    color='white', fontsize=8, fontweight='bold', zorder=4,
                    multialignment='center')
            ax.annotate('', xy=(exc_x - exc_w / 2, y), xytext=(cx + box_w / 2, y),
                        arrowprops=dict(arrowstyle='->', color='#e74c3c', lw=1.5), zorder=2)
            ax.text(cx + box_w / 2 + 0.15, y + 0.08, 'YES', fontsize=7,
                    color='#e74c3c', fontweight='bold')

        if i < len(criteria) - 1:
            next_y = criteria[i + 1][3]
            ax.annotate('', xy=(cx, next_y + box_h / 2), xytext=(cx, y - box_h / 2),
                        arrowprops=dict(arrowstyle='->', color='#555555', lw=1.5), zorder=2)

    retain_y = 1.0
    retain_box = mpatches.FancyBboxPatch(
        (cx - box_w / 2, retain_y - box_h / 2), box_w, box_h,
        boxstyle='round,pad=0.1', facecolor='#27ae60', edgecolor='none', zorder=3
    )
    ax.add_patch(retain_box)
    ax.text(cx, retain_y, 'RETAIN: Suitable for SETI Search',
            ha='center', va='center', color='white', fontsize=10, fontweight='bold', zorder=4)
    ax.annotate('', xy=(cx, retain_y + box_h / 2), xytext=(cx, criteria[-1][3] - box_h / 2),
                arrowprops=dict(arrowstyle='->', color='#555555', lw=1.5), zorder=2)

    ax.set_title('Parametric Avoidance Model — Decision Flow', fontsize=13, fontweight='bold', pad=10)
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, 'figure4_decision_flow.png'), dpi=200, bbox_inches='tight',
                facecolor='white')
    plt.close()
    print("Saved: figure4_decision_flow.png")



# FIGURE 5 — Sensitivity Analysis

def plot_sensitivity(df):
    baseline_rate = df['excluded'].mean() * 100

    mass_thresholds = np.arange(1.0, 2.5, 0.1)
    mass_excl = [((df['mass_flame_spec'] > t) | df['exclude_age'] | df['exclude_spec'] |
                  df['exclude_feh'] | df['exclude_mult'] | df['exclude_var'] | df['exclude_act']
                  ).mean() * 100 for t in mass_thresholds]

    age_thresholds = np.arange(1.0, 6.0, 0.25)
    age_excl = [(df['exclude_mass'] | (df['age_flame_spec'] < t) | df['exclude_spec'] |
                 df['exclude_feh'] | df['exclude_mult'] | df['exclude_var'] | df['exclude_act']
                 ).mean() * 100 for t in age_thresholds]

    feh_thresholds = np.arange(-0.8, 0.0, 0.05)
    feh_excl = [(df['exclude_mass'] | df['exclude_age'] | df['exclude_spec'] |
                 (df['mh_gspphot'] < t) | df['exclude_mult'] | df['exclude_var'] | df['exclude_act']
                 ).mean() * 100 for t in feh_thresholds]

    var_thresholds = np.arange(0.005, 0.1, 0.005)
    var_excl = [(df['exclude_mass'] | df['exclude_age'] | df['exclude_spec'] |
                 df['exclude_feh'] | df['exclude_mult'] |
                 ((df['range_mag_g_fov'].fillna(0) > t) | (df['phot_variable_flag'] == 'VARIABLE')) |
                 df['exclude_act']).mean() * 100 for t in var_thresholds]

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    fig.suptitle('Sensitivity Analysis — Exclusion Rate vs Threshold Value',
                 fontsize=13, fontweight='bold')

    configs = [
        (axes[0, 0], mass_thresholds, mass_excl, 1.5, 'Mass Threshold (M⊙)', 'Mass Threshold Sensitivity', 'Baseline (1.5 M⊙)'),
        (axes[0, 1], age_thresholds,  age_excl,  3.0, 'Age Threshold (Gyr)',  'Age Threshold Sensitivity',  'Baseline (3 Gyr)'),
        (axes[1, 0], feh_thresholds,  feh_excl, -0.4, '[Fe/H] Threshold',     'Metallicity Threshold Sensitivity', 'Baseline ([Fe/H] < −0.4)'),
        (axes[1, 1], var_thresholds,  var_excl,  0.01,'Variability Threshold (mag)', 'Variability Threshold Sensitivity', 'Baseline (0.01 mag)'),
    ]

    for ax, xdata, ydata, baseline_x, xlabel, title, blabel in configs:
        ax.plot(xdata, ydata, color='steelblue', linewidth=2)
        ax.axvline(baseline_x, color='red', linestyle='--', linewidth=1.5, label=blabel)
        ax.axhline(baseline_rate, color='gray', linestyle=':', alpha=0.7, linewidth=1.2)
        ax.set_xlabel(xlabel, fontsize=10)
        ax.set_ylabel('Exclusion Rate (%)', fontsize=10)
        ax.set_title(title, fontsize=11)
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, 'figure5_sensitivity_analysis.png'), dpi=200, bbox_inches='tight')
    plt.close()
    print("Saved: figure5_sensitivity_analysis.png")


# FIGURE 6 — Sky Distribution (Mollweide)

def plot_skymap(df):
    fig = plt.figure(figsize=(18, 9), facecolor='white')
    ax = fig.add_subplot(111, projection='mollweide')
    ax.set_facecolor('white')

    ra_rad  = np.deg2rad(df['ra'].values - 180)
    dec_rad = np.deg2rad(df['dec'].values)
    excluded_mask = df['excluded'].values.astype(bool)
    retained_mask = ~excluded_mask
    excl_count = excluded_mask.sum()
    ret_count  = retained_mask.sum()
    total = len(df)

    ax.scatter(ra_rad[excluded_mask], dec_rad[excluded_mask],
               s=0.08, c='#d62728', alpha=0.22, linewidths=0, rasterized=True)
    ax.scatter(ra_rad[retained_mask], dec_rad[retained_mask],
               s=0.12, c='#1f77b4', alpha=0.42, linewidths=0, rasterized=True)

    ax.grid(True, color='gray', alpha=0.25, linewidth=0.5)
    ax.set_xticklabels(['210°','240°','270°','300°','330°','0°','30°','60°','90°','120°','150°'], fontsize=9)

    legend_elements = [
        Line2D([0],[0], marker='o', color='w', markerfacecolor='#d62728', markersize=11,
               label=f'Excluded  (N = {excl_count:,}  |  {excl_count/total*100:.1f}%)', linewidth=0),
        Line2D([0],[0], marker='o', color='w', markerfacecolor='#1f77b4', markersize=11,
               label=f'Retained  (N = {ret_count:,}  |  {ret_count/total*100:.1f}%)', linewidth=0),
    ]
    ax.legend(handles=legend_elements, loc='lower right',
              facecolor='white', edgecolor='black', fontsize=11, framealpha=1)
    ax.set_title(f'Sky Distribution of SETI Candidate Stars — Gaia DR3\n'
                 f'Parametric Avoidance Model  |  N = {total:,}',
                 fontsize=14, fontweight='bold', pad=20)

    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, 'figure6_skymap.png'), dpi=200, bbox_inches='tight', facecolor='white')
    plt.close()
    print("Saved: figure6_skymap.png")

# RUN ALL

if __name__ == '__main__':
    import pandas as pd
    df = pd.read_csv(os.path.join(save_dir, 'gaia_seti_avoidance_robust_corrected.csv'), low_memory=False)
    print(f"Loaded: {df.shape}")

    plot_criteria_breakdown(df)
    plot_spectype_breakdown(df)
    plot_synthetic_vs_empirical(df)
    plot_decision_flow()
    plot_sensitivity(df)
    plot_skymap(df)

    print("\nAll 6 figures saved.")
