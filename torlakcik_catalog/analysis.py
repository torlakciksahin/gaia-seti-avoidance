"""
Statistical analysis module for the Torlakcık Catalog.

"""

import pandas as pd
import numpy as np
from typing import Tuple
from scipy import stats
from config import TELESCOPES, THRESHOLDS


def calculate_distance_statistics(df):
    """Calculate distance statistics for retained stars."""
    retained = df[df['decision'] == 'RETAIN'].copy()

    # Check if parallax column exists - use list to avoid KeyError
    cols = list(df.columns)
    if 'parallax' not in cols:
        return {
            'median_distance_pc': None,
            'mean_distance_pc': None,
            'std_distance_pc': None,
            'min_distance_pc': None,
            'max_distance_pc': None,
            'q25_distance_pc': None,
            'q75_distance_pc': None,
            'n_valid': 0,
            'error': 'Parallax column not found in catalog'
        }

    retained['distance_pc'] = 1000.0 / retained['parallax']
    valid_distances = retained[retained['distance_pc'] > 0]['distance_pc']

    return {
        'median_distance_pc': valid_distances.median(),
        'mean_distance_pc': valid_distances.mean(),
        'std_distance_pc': valid_distances.std(),
        'min_distance_pc': valid_distances.min(),
        'max_distance_pc': valid_distances.max(),
        'q25_distance_pc': valid_distances.quantile(0.25),
        'q75_distance_pc': valid_distances.quantile(0.75),
        'n_valid': len(valid_distances),
    }


def calculate_sky_density(df):
    """Calculate sky density of retained stars."""
    retained = df[df['decision'] == 'RETAIN']
    n_retained = len(retained)
    total_sky_area_sqdeg = 4 * np.pi * (180 / np.pi) ** 2
    density = n_retained / total_sky_area_sqdeg
    avg_separation_arcmin = np.sqrt(1 / density) * 60
    fov_needed_sqdeg = 1 / density
    fov_radius_arcmin = np.sqrt(fov_needed_sqdeg / np.pi) * 60

    return {
        'n_retained': n_retained,
        'total_sky_area_sqdeg': total_sky_area_sqdeg,
        'density_per_sqdeg': density,
        'avg_separation_arcmin': avg_separation_arcmin,
        'fov_for_1_star_sqarcmin': fov_needed_sqdeg * 3600,
        'fov_radius_arcmin': fov_radius_arcmin,
    }


def calculate_fov_coverage(df, beam_diameter_arcmin):
    """Calculate probability of finding retained stars in a telescope beam."""
    sky_stats = calculate_sky_density(df)
    density = sky_stats['density_per_sqdeg']
    beam_radius_deg = (beam_diameter_arcmin / 2) / 60
    beam_area_sqdeg = np.pi * beam_radius_deg ** 2
    expected_stars = density * beam_area_sqdeg
    p_at_least_1 = 1 - np.exp(-expected_stars)
    p_at_least_n = {}
    for n in [1, 2, 3, 5, 10]:
        if n == 1:
            p_at_least_n[n] = p_at_least_1
        else:
            p_at_least_n[n] = 1 - stats.poisson.cdf(n - 1, expected_stars)

    return {
        'beam_diameter_arcmin': beam_diameter_arcmin,
        'beam_area_sqdeg': beam_area_sqdeg,
        'expected_stars': expected_stars,
        'p_at_least_1': p_at_least_1,
        'p_at_least_n': p_at_least_n,
    }


def analyze_all_telescopes(df):
    """Calculate FOV coverage for all configured telescopes."""
    results = []
    for name, params in TELESCOPES.items():
        coverage = calculate_fov_coverage(df, params['beam_fwhm_arcmin'])
        results.append({
            'Telescope': name,
            'Beam FWHM (arcmin)': params['beam_fwhm_arcmin'],
            'Expected Stars': coverage['expected_stars'],
            'P(>=1 star)': coverage['p_at_least_1'],
            'P(>=2 stars)': coverage['p_at_least_n'][2],
        })
    return pd.DataFrame(results)


def compare_age_criteria(df):
    """Compare old vs new age criterion approaches."""
    old_young = df[df['age_flame_spec'] < 3.0]
    new_young = df[df['age_flame_spec_upper'] < 3.0]
    boundary_zone = df[(df['age_flame_spec'] >= 2.0) & (df['age_flame_spec'] <= 4.0)]

    return {
        'old_r2_count': len(old_young),
        'new_r2_count': len(new_young),
        'difference': len(old_young) - len(new_young),
        'boundary_zone_count': len(boundary_zone),
        'boundary_zone_fraction': len(boundary_zone) / len(df),
    }


def analyze_exclusion_by_spectral_type(df):
    """Create detailed breakdown by spectral type."""
    results = []
    for sp_type in ['O', 'B', 'A', 'F0-F4', 'F5-F9', 'G', 'K', 'M']:
        sp_df = df[df['spectral_type'] == sp_type]
        n_total = len(sp_df)
        n_excluded = len(sp_df[sp_df['decision'] == 'EXCLUDE'])
        n_retained = n_total - n_excluded

        results.append({
            'Spectral Type': sp_type,
            'Total': n_total,
            'Excluded': n_excluded,
            'Retained': n_retained,
            'Exclusion Rate (%)': 100 * n_excluded / n_total if n_total > 0 else 0,
        })

    return pd.DataFrame(results)


def compare_empirical_vs_synthetic(df):
    """
    Compare empirical Gaia DR3 flags against synthetic proxy methods.

    Empirical proxies use Gaia's published classification flags:
      - R5: non_single_star >= 1
      - R6: range_mag_g_fov > variability_threshold OR phot_variable_flag == 'VARIABLE'

    Synthetic proxies use derived quantities:
      - R5: RUWE > 1.4 (Renormalized Unit Weight Error, indicates multiplicity)
      - R6: sigma_G/G > synthetic_variability_threshold
            (fractional flux error as variability indicator)

    The sigma_G/G threshold is set to 0.001 (0.1%) because Gaia's mean
    photometry is extremely precise (median fractional error ~0.02%).
    The original 5% threshold was inappropriate for Gaia DR3 mean flux
    errors, which measure photometric precision rather than variability
    amplitude. Even at 0.1%, this proxy captures far fewer variables
    than empirical flags, reflecting the fundamental difference between
    precision metrics and variability classification.

    Parameters
    ----------
    df : pd.DataFrame
        Catalog with decision, reason_code, and Gaia columns

    Returns
    -------
    dict
        Dictionary with empirical and synthetic percentages
    """
    total = len(df)

    # --- Empirical R5: Gaia non_single_star flag ---
    emp_r5_mask = df['non_single_star'] >= 1

    # --- Synthetic R5: RUWE > 1.4 ---
    if 'ruwe' in df.columns:
        syn_r5_mask = df['ruwe'] > 1.4
    else:
        print("WARNING: 'ruwe' column missing. Synthetic R5 skipped.")
        syn_r5_mask = pd.Series(False, index=df.index)

    # --- Empirical R6: Gaia variability flags ---
    emp_r6_mask = (
        (df['range_mag_g_fov'] > THRESHOLDS['variability']) |
        (df['phot_variable_flag'] == 'VARIABLE')
    )

    # --- Synthetic R6: fractional flux error (sigma_G / G) ---
    if 'phot_g_mean_flux' in df.columns and 'phot_g_mean_flux_error' in df.columns:
        fractional_error = df['phot_g_mean_flux_error'].astype(np.float64) / df['phot_g_mean_flux'].astype(np.float64)
        syn_r6_mask = fractional_error > THRESHOLDS['synthetic_variability']
        # Also compute stats for the paper
        fe_mean = fractional_error.mean()
        fe_median = fractional_error.median()
        fe_gt_001 = (fractional_error > 0.001).sum()
        fe_gt_005 = (fractional_error > 0.005).sum()
        fe_gt_01 = (fractional_error > 0.01).sum()
    else:
        print("WARNING: Flux columns missing. Synthetic R6 skipped.")
        syn_r6_mask = pd.Series(False, index=df.index)
        fe_mean = fe_median = 0.0
        fe_gt_001 = fe_gt_005 = fe_gt_01 = 0

    # --- Overall exclusion ---
    # Empirical: from the actual pipeline decision
    emp_excluded = len(df[df['decision'] == 'EXCLUDE'])

    # Synthetic: use base exclusion (R1,R2,R3,R4,R7) + synthetic R5,R6
    base_exclusion_mask = df['reason_code'].str.contains('R1|R2|R3|R4|R7', regex=True, na=False)
    syn_overall_mask = base_exclusion_mask | syn_r5_mask | syn_r6_mask
    syn_excluded = syn_overall_mask.sum()

    result = {
        'total': total,
        'emp_r5_pct': (emp_r5_mask.sum() / total) * 100,
        'syn_r5_pct': (syn_r5_mask.sum() / total) * 100,
        'emp_r6_pct': (emp_r6_mask.sum() / total) * 100,
        'syn_r6_pct': (syn_r6_mask.sum() / total) * 100,
        'emp_overall_pct': (emp_excluded / total) * 100,
        'syn_overall_pct': (syn_excluded / total) * 100,
        # Additional stats for the paper
        'fe_mean': fe_mean,
        'fe_median': fe_median,
        'fe_gt_001_pct': (fe_gt_001 / total) * 100,
        'fe_gt_005_pct': (fe_gt_005 / total) * 100,
        'fe_gt_01_pct': (fe_gt_01 / total) * 100,
    }

    return result


def analyze_nearby_proxy(df, distance_limit_pc=100.0):
    """
    Analyze stars within a distance limit as a proxy for GCNS comparison.
    
    Per H. Isaacson feedback: compare with Gaia Catalog of Nearby Stars.
    Since GCNS is not directly available in our data, we use a distance
    cut (d <= 100 pc) as a volume-complete proxy. This approximates the
    GCNS sample (Gaia Collaboration et al. 2021, 331,312 stars within 100 pc).
    
    Parameters
    ----------
    df : pd.DataFrame
        Catalog with 'decision', 'parallax' columns
    distance_limit_pc : float
        Distance limit in parsecs (default: 100 pc for GCNS proxy)
    
    Returns
    -------
    dict
        Nearby sample analysis results
    """
    if 'parallax' not in df.columns:
        return {'error': 'Parallax column not found'}
    
    # Compute distance
    df_copy = df.copy()
    valid_plx = df_copy['parallax'] > 0
    df_copy.loc[valid_plx, 'distance_pc'] = 1000.0 / df_copy.loc[valid_plx, 'parallax']
    
    nearby = df_copy[df_copy['distance_pc'] <= distance_limit_pc]
    n_nearby = len(nearby)
    
    if n_nearby == 0:
        return {'error': f'No stars found within {distance_limit_pc} pc'}
    
    nearby_retained = nearby[nearby['decision'] == 'RETAIN']
    nearby_excluded = nearby[nearby['decision'] == 'EXCLUDE']
    
    # Exclusion reasons
    reason_counts = {}
    for code in ['R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7']:
        count = len(nearby_excluded[nearby_excluded['reason_code'].str.contains(code, na=False)])
        reason_counts[code] = count
    
    # Spectral type breakdown
    sp_breakdown = {}
    if 'spectral_type' in nearby.columns:
        for sp in ['O', 'B', 'A', 'F0-F4', 'F5-F9', 'G', 'K', 'M']:
            sp_df = nearby[nearby['spectral_type'] == sp]
            if len(sp_df) > 0:
                sp_excl = len(sp_df[sp_df['decision'] == 'EXCLUDE'])
                sp_breakdown[sp] = {
                    'total': len(sp_df),
                    'retained': len(sp_df) - sp_excl,
                    'excluded': sp_excl,
                    'exclusion_rate': sp_excl / len(sp_df)
                }
    
    # Full catalog stats for comparison
    full_retained = len(df_copy[df_copy['decision'] == 'RETAIN'])
    
    return {
        'nearby_total': n_nearby,
        'nearby_retained': len(nearby_retained),
        'nearby_excluded': len(nearby_excluded),
        'nearby_retention_rate': len(nearby_retained) / n_nearby,
        'nearby_exclusion_rate': len(nearby_excluded) / n_nearby,
        'full_retention_rate': full_retained / len(df_copy) if len(df_copy) > 0 else 0,
        'distance_limit_pc': distance_limit_pc,
        'reason_counts': reason_counts,
        'spectral_breakdown': sp_breakdown,
    }


def analyze_parallax_quality(df):
    """
    Analyze parallax quality to justify naive distance methodology.
    
    Per H. Isaacson feedback: clarify why naive 1/π is sufficient.
    Key argument: with ϖ/σ_ϖ > 5, the fractional distance uncertainty
    is < 20%, and Bailer-Jones et al. (2021) show that for such
    high-SNR parallaxes, naive and Bayesian distances agree within ~5%.
    """
    if 'parallax_over_error' not in df.columns:
        return {'error': 'parallax_over_error not available'}
    
    poe = df['parallax_over_error'].dropna()
    
    return {
        'n_valid': len(poe),
        'median_poe': poe.median(),
        'mean_poe': poe.mean(),
        'min_poe': poe.min(),
        'pct_gt_5': (poe > 5).mean() * 100,
        'pct_gt_10': (poe > 10).mean() * 100,
        'pct_gt_20': (poe > 20).mean() * 100,
        'pct_gt_50': (poe > 50).mean() * 100,
        'p25': poe.quantile(0.25),
        'p75': poe.quantile(0.75),
    }


def print_analysis_report(df):
    """Print comprehensive analysis report."""
    print("\n" + "=" * 70)
    print("COMPREHENSIVE ANALYSIS REPORT")
    print("=" * 70)

    # Distance statistics
    print("\n--- Distance Statistics (Retained Stars) ---")
    dist_stats = calculate_distance_statistics(df)
    if dist_stats.get('error'):
        print(f"Note: {dist_stats['error']}")
    else:
        print(f"Median distance: {dist_stats['median_distance_pc']:.0f} pc")
        print(f"Mean distance: {dist_stats['mean_distance_pc']:.0f} pc")
        print(f"Valid distances: {dist_stats['n_valid']:,}")

    # Parallax quality (justifies naive 1/π methodology)
    print("\n--- Parallax Quality (Distance Methodology Justification) ---")
    poe_stats = analyze_parallax_quality(df)
    if 'error' not in poe_stats:
        print(f"  Median ϖ/σ_ϖ: {poe_stats['median_poe']:.1f}")
        print(f"  Mean ϖ/σ_ϖ: {poe_stats['mean_poe']:.1f}")
        print(f"  > 10: {poe_stats['pct_gt_10']:.1f}%  |  > 20: {poe_stats['pct_gt_20']:.1f}%  |  > 50: {poe_stats['pct_gt_50']:.1f}%")
        print(f"  → With ϖ/σ_ϖ > 5 cut, naive 1/π distances are reliable.")
        print(f"    Bailer-Jones et al. (2021) show <5% bias for such high-SNR parallaxes.")

    # Nearby sample (GCNS proxy)
    print("\n--- Nearby Sample (≤100 pc, GCNS Proxy) ---")
    nearby_stats = analyze_nearby_proxy(df)
    if 'error' not in nearby_stats:
        print(f"  Stars within 100 pc: {nearby_stats['nearby_total']:,}")
        print(f"  Retained: {nearby_stats['nearby_retained']:,} ({100*nearby_stats['nearby_retention_rate']:.1f}%)")
        print(f"  Excluded: {nearby_stats['nearby_excluded']:,} ({100*nearby_stats['nearby_exclusion_rate']:.1f}%)")
        print(f"  Full catalog retention: {100*nearby_stats['full_retention_rate']:.1f}%")
        print(f"  Primary exclusion reasons for nearby stars:")
        for code, count in nearby_stats['reason_counts'].items():
            if count > 0:
                print(f"    {code}: {count}")

    # Sky density
    print("\n--- Sky Density Analysis ---")
    density_stats = calculate_sky_density(df)
    print(f"Retained stars: {density_stats['n_retained']:,}")
    print(f"Sky density: {density_stats['density_per_sqdeg']:.2f} stars/sq deg")
    print(f"Average separation: {density_stats['avg_separation_arcmin']:.1f} arcmin")

    # Telescope coverage
    print("\n--- Telescope FOV Coverage ---")
    telescope_df = analyze_all_telescopes(df)
    print(telescope_df.to_string(index=False))

    # Age criteria comparison
    print("\n--- Age Criteria Comparison ---")
    age_comparison = compare_age_criteria(df)
    print(f"Old criterion (age < 3 Gyr): {age_comparison['old_r2_count']:,}")
    print(f"New criterion (age_upper < 3 Gyr): {age_comparison['new_r2_count']:,}")
    print(f"Difference: {age_comparison['difference']:,} stars now retained")
    print(f"Boundary zone (2-4 Gyr): {age_comparison['boundary_zone_count']:,} ({100*age_comparison['boundary_zone_fraction']:.1f}%)")

    # Spectral type breakdown
    print("\n--- Spectral Type Breakdown ---")
    sp_df = analyze_exclusion_by_spectral_type(df)
    print(sp_df.to_string(index=False))

    # Empirical vs Synthetic proxy comparison
    print("\n--- Empirical vs. Synthetic Proxy Comparison ---")
    proxy_comp = compare_empirical_vs_synthetic(df)
    print(f"R5 (Multiplicity): Empirical = {proxy_comp['emp_r5_pct']:.1f}% | Synthetic (RUWE>1.4) = {proxy_comp['syn_r5_pct']:.1f}%")
    print(f"R6 (Variability):  Empirical = {proxy_comp['emp_r6_pct']:.1f}% | Synthetic (sigma_G/G>{THRESHOLDS['synthetic_variability']}) = {proxy_comp['syn_r6_pct']:.1f}%")
    print(f"Overall Exclusion: Empirical = {proxy_comp['emp_overall_pct']:.1f}% | Synthetic = {proxy_comp['syn_overall_pct']:.1f}%")
    diff = proxy_comp['syn_overall_pct'] - proxy_comp['emp_overall_pct']
    print(f"Difference: Synthetic proxies {'overestimate' if diff > 0 else 'underestimate'} exclusion by {abs(diff):.1f} percentage points.")
    # Print fractional error stats for interpretation
    print(f"\n  Fractional flux error (sigma_G/G) distribution:")
    print(f"    Mean: {proxy_comp['fe_mean']:.6f} ({proxy_comp['fe_mean']*100:.4f}%)")
    print(f"    Median: {proxy_comp['fe_median']:.6f} ({proxy_comp['fe_median']*100:.4f}%)")
    print(f"    > 0.1%: {proxy_comp['fe_gt_001_pct']:.2f}% of stars")
    print(f"    > 0.5%: {proxy_comp['fe_gt_005_pct']:.2f}% of stars")
    print(f"    > 1.0%: {proxy_comp['fe_gt_01_pct']:.2f}% of stars")
