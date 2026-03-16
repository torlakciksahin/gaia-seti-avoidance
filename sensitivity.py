# SETI Parametric Avoidance Model — Sensitivity Analysis
# Input: df (processed DataFrame with exclusion
 # Output: sensitivity_analysis.png

import matplotlib.pyplot as plt
import numpy as np
import os

save_dir = r'C:\Users\PC\OneDrive\Desktop\avoidance model'
os.makedirs(save_dir, exist_ok=True)


def run_sensitivity(df):
    mass_thresholds = np.arange(1.0, 2.5, 0.1)
    mass_excl = []
    for t in mass_thresholds:
        excl = (
            (df['mass_flame_spec'] > t) | df['exclude_age'] | df['exclude_spec'] |
            df['exclude_feh'] | df['exclude_mult'] | df['exclude_var'] | df['exclude_act']
        ).mean() * 100
        mass_excl.append(excl)

    age_thresholds = np.arange(1.0, 6.0, 0.25)
    age_excl = []
    for t in age_thresholds:
        excl = (
            df['exclude_mass'] | (df['age_flame_spec'] < t) | df['exclude_spec'] |
            df['exclude_feh'] | df['exclude_mult'] | df['exclude_var'] | df['exclude_act']
        ).mean() * 100
        age_excl.append(excl)

    feh_thresholds = np.arange(-0.8, 0.0, 0.05)
    feh_excl = []
    for t in feh_thresholds:
        excl = (
            df['exclude_mass'] | df['exclude_age'] | df['exclude_spec'] |
            (df['mh_gspphot'] < t) | df['exclude_mult'] | df['exclude_var'] | df['exclude_act']
        ).mean() * 100
        feh_excl.append(excl)

    var_thresholds = np.arange(0.005, 0.1, 0.005)
    var_excl = []
    for t in var_thresholds:
        excl = (
            df['exclude_mass'] | df['exclude_age'] | df['exclude_spec'] |
            df['exclude_feh'] | df['exclude_mult'] |
            (df['range_mag_g_fov'].fillna(0) > t) | df['exclude_act']
        ).mean() * 100
        var_excl.append(excl)

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle('Sensitivity Analysis — Exclusion Rate vs Threshold Value',
                 fontsize=13, fontweight='bold')

    axes[0,0].plot(mass_thresholds, mass_excl, 'steelblue', linewidth=2)
    axes[0,0].axvline(1.5, color='red', linestyle='--', label='Baseline (1.5 M⊙)')
    axes[0,0].axhline('excluded', color='gray', linestyle=':', alpha=0.7)
    axes[0,0].set_xlabel('Mass Threshold (M⊙)')
    axes[0,0].set_ylabel('Exclusion Rate (%)')
    axes[0,0].set_title('Mass Threshold Sensitivity')
    axes[0,0].legend()
    axes[0,0].grid(alpha=0.3)

    axes[0,1].plot(age_thresholds, age_excl, 'steelblue', linewidth=2)
    axes[0,1].axvline(3.0, color='red', linestyle='--', label='Baseline (3 Gyr)')
    axes[0,1].axhline(65.3, color='gray', linestyle=':', alpha=0.7)
    axes[0,1].set_xlabel('Age Threshold (Gyr)')
    axes[0,1].set_ylabel('Exclusion Rate (%)')
    axes[0,1].set_title('Age Threshold Sensitivity')
    axes[0,1].legend()
    axes[0,1].grid(alpha=0.3)

    axes[1,0].plot(feh_thresholds, feh_excl, 'steelblue', linewidth=2)
    axes[1,0].axvline(-0.4, color='red', linestyle='--', label='Baseline ([Fe/H] < -0.4)')
    axes[1,0].axhline(65.3, color='gray', linestyle=':', alpha=0.7)
    axes[1,0].set_xlabel('[Fe/H] Threshold')
    axes[1,0].set_ylabel('Exclusion Rate (%)')
    axes[1,0].set_title('Metallicity Threshold Sensitivity')
    axes[1,0].legend()
    axes[1,0].grid(alpha=0.3)

    axes[1,1].plot(var_thresholds, var_excl, 'steelblue', linewidth=2)
    axes[1,1].axvline(0.01, color='red', linestyle='--', label='Baseline (0.01 mag)')
    axes[1,1].axhline(65.3, color='gray', linestyle=':', alpha=0.7)
    axes[1,1].set_xlabel('Variability Threshold (mag)')
    axes[1,1].set_ylabel('Exclusion Rate (%)')
    axes[1,1].set_title('Variability Threshold Sensitivity')
    axes[1,1].legend()
    axes[1,1].grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, 'figure5_sensitivity_analysis.png'),
                dpi=200, bbox_inches='tight')
    plt.close()
    print("Saved: figure5_sensitivity_analysis.png")


if __name__ == '__main__':
    print("Run analysis.py first to generate df, then call run_sensitivity(df).")
