"""
Torlakcık Catalog: SETI Target Avoidance Model
Main Execution Script

This script runs the complete analysis pipeline including:
- Data loading & derived columns (flux→mag, distance)
- GCNS cross-match by source_id
- Avoidance model application with age uncertainty
- Statistical analysis
- Empirical vs. Synthetic proxy comparison
- Sensitivity analysis (fast mode)
- Cross-match with Breakthrough Listen (Isaacson 2017 + MeerKAT)
- Figure generation

"""

import os
import sys
import pandas as pd
import numpy as np
from astropy.table import Table

# Import local modules
from config import (
    FITS_FILE, OUTPUT_CATALOG, OUTPUT_RETAINED, BASE_DIR,
    GCNS_FILE, BL_ISAACSON_FILE, BL_MEERKAT_1M_FILE,
    G_ZEROPOINT, OUTPUT_FIGURES_DIR, OUTPUT_DIR
)
from avoidance_model import (
    apply_avoidance_model,
    calculate_exclusion_statistics,
    print_statistics,
)
from analysis import print_analysis_report, compare_age_criteria, compare_empirical_vs_synthetic
from sensitivity import run_full_sensitivity_analysis
from cross_match import (
    load_bl_catalog_isaacson,
    load_bl_catalog_meerkat,
    analyze_bl_overlap,
    print_bl_comparison_report,
)
from figures import generate_all_figures


def load_data(filepath: str) -> pd.DataFrame:
    """Load stellar catalog from FITS file."""
    print(f"Loading data from {filepath}...")
    t = Table.read(filepath)
    df = t.to_pandas()
    print(f"Loaded {len(df):,} stars")
    print(f"  Columns: {list(df.columns)}")
    return df


def compute_derived_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute derived columns: phot_g_mean_mag, abs_g, distance_pc, bp_rp.
    Also normalizes Gaia DR3 column names (e.g., age_flame_spec_upp → upper).
    """
    # - Column name normalization (Gaia DR3 uses short names)
    rename_map = {}
    if 'age_flame_spec_upp' in df.columns and 'age_flame_spec_upper' not in df.columns:
        rename_map['age_flame_spec_upp'] = 'age_flame_spec_upper'
    if 'age_flame_spec_low' in df.columns and 'age_flame_spec_lower' not in df.columns:
        rename_map['age_flame_spec_low'] = 'age_flame_spec_lower'
    if rename_map:
        df = df.rename(columns=rename_map)
        print(f"  Renamed columns: {rename_map}")

    # - phot_g_mean_mag from flux (if not already present)
    if 'phot_g_mean_mag' not in df.columns and 'phot_g_mean_flux' in df.columns:
        valid_flux = df['phot_g_mean_flux'] > 0
        df.loc[valid_flux, 'phot_g_mean_mag'] = (
            -2.5 * np.log10(df.loc[valid_flux, 'phot_g_mean_flux']) + G_ZEROPOINT
        )
        print(f"  Computed phot_g_mean_mag from flux ({valid_flux.sum():,} stars)")

    # - bp_rp from BP/RP magnitudes (if not already present)
    if 'bp_rp' not in df.columns:
        has_bp = 'phot_bp_mean_mag' in df.columns
        has_rp = 'phot_rp_mean_mag' in df.columns
        if has_bp and has_rp:
            df['bp_rp'] = df['phot_bp_mean_mag'] - df['phot_rp_mean_mag']
            valid_bp_rp = df['bp_rp'].notna().sum()
            print(f"  Computed bp_rp from BP/RP magnitudes ({valid_bp_rp:,} valid)")
        else:
            print(f"  WARNING: Cannot compute bp_rp — missing phot_bp_mean_mag or phot_rp_mean_mag")
            print(f"  Available columns: {list(df.columns)[:20]}...")

    # - abs_g and distance_pc from parallax
    if 'parallax' in df.columns and 'phot_g_mean_mag' in df.columns:
        valid_plx = df['parallax'] > 0
        df.loc[valid_plx, 'distance_pc'] = 1000.0 / df.loc[valid_plx, 'parallax']
        df.loc[valid_plx, 'abs_g'] = (
            df.loc[valid_plx, 'phot_g_mean_mag']
            + 5 * np.log10(df.loc[valid_plx, 'parallax'] / 100)
        )
        print(f"  Computed abs_g and distance_pc ({valid_plx.sum():,} stars)")
    else:
        missing = []
        if 'parallax' not in df.columns:
            missing.append('parallax')
        if 'phot_g_mean_mag' not in df.columns:
            missing.append('phot_g_mean_mag')
        print(f"  WARNING: Cannot compute abs_g — missing: {missing}")

    return df


def load_gcns(filepath: str) -> set:
    """
    Load GCNS (Gaia Catalog of Nearby Stars) and return source_ids.
    Supports VOTable format.
    """
    if not os.path.exists(filepath):
        print(f"  GCNS file not found: {filepath}")
        return set()

    try:
        t = Table.read(filepath)
        gcns_df = t.to_pandas()

        # Find source_id column
        sid_col = None
        for cand in ['source_id', 'Source', 'SOURCE_ID', 'sid']:
            if cand in gcns_df.columns:
                sid_col = cand
                break

        if sid_col is None:
            # Try case-insensitive
            for col in gcns_df.columns:
                if 'source' in col.lower() and 'id' in col.lower():
                    sid_col = col
                    break

        if sid_col is None:
            print(f"  WARNING: Could not find source_id column in GCNS")
            print(f"  GCNS columns: {list(gcns_df.columns)[:15]}")
            return set()

        source_ids = set(gcns_df[sid_col].astype(np.int64))
        print(f"  Loaded GCNS: {len(source_ids):,} source_ids (column: {sid_col})")
        return source_ids

    except Exception as e:
        print(f"  ERROR loading GCNS: {e}")
        return set()


def cross_match_gcns(df: pd.DataFrame, gcns_ids: set) -> pd.DataFrame:
    """Add GCNS membership flag to catalog."""
    if not gcns_ids:
        df['in_gcns'] = False
        return df

    df['in_gcns'] = df['source_id'].isin(gcns_ids)
    gcns_members = df['in_gcns'].sum()
    print(f"  GCNS members in catalog: {gcns_members:,}")
    return df


def save_results(catalog: pd.DataFrame):
    """Save catalog in CSV and FITS formats."""
    from astropy.table import Table as AstropyTable

    print("\nSaving results...")
    retained = catalog[catalog['decision'] == 'RETAIN']

    catalog_dir = os.path.join(OUTPUT_DIR, 'catalogs')
    os.makedirs(catalog_dir, exist_ok=True)

    base_full = os.path.join(catalog_dir, 'torlakcik_catalog')
    base_ret  = os.path.join(catalog_dir, 'retained_stars')

    def _sanitize_for_fits(df):
        """Convert problematic dtypes for FITS compatibility."""
        df_out = df.copy()
        for col in df_out.columns:
            if pd.api.types.is_integer_dtype(df_out[col]) and df_out[col].isna().any():
                df_out[col] = df_out[col].astype('float64')
            elif pd.api.types.is_bool_dtype(df_out[col]) and df_out[col].isna().any():
                df_out[col] = df_out[col].astype('object')
            elif df_out[col].dtype == object:
                try:
                    df_out[col] = df_out[col].fillna('')
                except Exception:
                    pass
        return df_out

    for subset_name, df, base in [('full', catalog, base_full), ('retained', retained, base_ret)]:
        label = 'Full catalog' if subset_name == 'full' else 'Retained stars'
        n = len(df)
        print(f"\n  -- {label} (N={n:,})")

        csv_path = base + '.csv'
        df.to_csv(csv_path, index=False)
        print(f"    CSV:  {csv_path}")

        fits_path = base + '.fits'
        try:
            df_clean = _sanitize_for_fits(df)
            t = AstropyTable.from_pandas(df_clean)
            t.write(fits_path, format='fits', overwrite=True)
            print(f"    FITS: {fits_path}")
        except Exception as e:
            print(f"    FITS: FAILED ({e})")

    print(f"\n  Catalogs saved to: {catalog_dir}/")


def main():
    """Main execution function."""
    print("=" * 70)
    print("TORLAKCIK CATALOG: SETI TARGET AVOIDANCE MODEL")
    print("With Age Uncertainty-Aware Filtering + Howard Isaacson Feedback")
    print("=" * 70)

    # - Step 1: Load data
    print("\n[Step 1] Loading data...")
    try:
        catalog = load_data(FITS_FILE)
    except FileNotFoundError:
        print(f"Error: FITS file not found at {FITS_FILE}")
        sys.exit(1)

    # - Step 2: Compute derived columns
    print("\n[Step 2] Computing derived columns...")
    catalog = compute_derived_columns(catalog)

    # - Step 3: Load GCNS
    print("\n[Step 3] Loading GCNS...")
    gcns_ids = load_gcns(GCNS_FILE)
    catalog = cross_match_gcns(catalog, gcns_ids)

    # - Step 4: Apply avoidance model
    print("\n[Step 4] Applying avoidance model (age_upper < 3 Gyr)...")
    catalog = apply_avoidance_model(catalog, use_age_upper=True)

    # - Step 5: Calculate statistics
    print("\n[Step 5] Calculating statistics...")
    stats = calculate_exclusion_statistics(catalog)
    print_statistics(stats)

    # - Step 6: Comprehensive analysis
    print("\n[Step 6] Running comprehensive analysis...")
    print_analysis_report(catalog)

    # - Step 7: Age criteria comparison
    print("\n[Step 7] Comparing age criteria approaches...")
    age_comparison = compare_age_criteria(catalog)
    print(f"  Old approach (age < 3 Gyr): {age_comparison['old_r2_count']:,} excluded")
    print(f"  New approach (age_upper < 3 Gyr): {age_comparison['new_r2_count']:,} excluded")
    print(f"  Difference: {age_comparison['difference']:,} stars now retained")

    # - Step 8: Empirical vs. Synthetic proxy comparison
    print("\n[Step 8] Empirical vs Synthetic proxy comparison...")
    proxy_comp = compare_empirical_vs_synthetic(catalog)
    print(f"  R5 (Multiplicity): Empirical = {proxy_comp['emp_r5_pct']:.1f}% | Synthetic (RUWE>1.4) = {proxy_comp['syn_r5_pct']:.1f}%")
    print(f"  R6 (Variability):  Empirical = {proxy_comp['emp_r6_pct']:.1f}% | Synthetic = {proxy_comp['syn_r6_pct']:.1f}%")
    print(f"  Overall Exclusion: Empirical = {proxy_comp['emp_overall_pct']:.1f}% | Synthetic = {proxy_comp['syn_overall_pct']:.1f}%")
    diff = proxy_comp['syn_overall_pct'] - proxy_comp['emp_overall_pct']
    print(f"  Difference: Synthetic proxies {'overestimate' if diff > 0 else 'underestimate'} by {abs(diff):.1f} pp")

    # - Step 9: Parallax quality
    print("\n[Step 9] Parallax quality analysis (distance methodology)...")
    from analysis import analyze_parallax_quality
    poe_stats = analyze_parallax_quality(catalog)
    if 'error' not in poe_stats:
        print(f"  Median ϖ/σ_ϖ: {poe_stats['median_poe']:.1f}")
        print(f"  > 5: {poe_stats['pct_gt_5']:.1f}%  |  > 10: {poe_stats['pct_gt_10']:.1f}%  |  > 20: {poe_stats['pct_gt_20']:.1f}%")
        print(f"  → Naive 1/π distances justified by high ϖ/σ_ϖ.")
    else:
        print(f"  {poe_stats['error']}")

    # - Step 10: GCNS cross-match analysis
    print("\n[Step 10] GCNS cross-match analysis...")
    if catalog['in_gcns'].any():
        gcns_cat = catalog[catalog['in_gcns']]
        gcns_retained = len(gcns_cat[gcns_cat['decision'] == 'RETAIN'])
        gcns_excluded = len(gcns_cat[gcns_cat['decision'] == 'EXCLUDE'])
        full_retained_pct = len(catalog[catalog['decision'] == 'RETAIN']) / len(catalog) * 100
        print(f"  Method: source_id cross-match")
        print(f"  GCNS members: {len(gcns_cat):,}")
        print(f"  Retained: {gcns_retained:,} ({100*gcns_retained/len(gcns_cat):.1f}%)")
        print(f"  Excluded: {gcns_excluded:,} ({100*gcns_excluded/len(gcns_cat):.1f}%)")
        print(f"  Full catalog retention: {full_retained_pct:.1f}%")

        # GCNS exclusion reasons
        gcns_excl_df = gcns_cat[gcns_cat['decision'] == 'EXCLUDE']
        print(f"  Primary GCNS exclusion reasons:")
        for code in ['R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7']:
            count = len(gcns_excl_df[gcns_excl_df['reason_code'].str.contains(code, na=False)])
            if count > 0:
                print(f"    {code}: {count}")
    else:
        print("  No GCNS members found in catalog.")

    # - Step 11: Sensitivity analysis
    SKIP_SENSITIVITY = True  # ← Change to False to run

    if SKIP_SENSITIVITY:
        print("\n[Step 11] Sensitivity analysis SKIPPED (set SKIP_SENSITIVITY=False to run)")
        sensitivity_results = {}
        for param in ['age', 'mass', 'metallicity', 'variability']:
            csv_path = os.path.join(BASE_DIR, f'sensitivity_{param}.csv')
            if os.path.exists(csv_path):
                sensitivity_results[param] = pd.read_csv(csv_path)
                print(f"  Loaded saved: {csv_path}")
        if len(sensitivity_results) < 4:
            print("  Note: Some/all sensitivity CSVs missing. Figures will use empty results.")
    else:
        print("\n[Step 11] Running sensitivity analysis (FAST mode)...")
        sensitivity_results = run_full_sensitivity_analysis(catalog)

        # Print results
        print("\n  Sensitivity Inflection Points:")
        from sensitivity import find_inflection_points
        inflection = find_inflection_points(sensitivity_results)
        for param, info in inflection.items():
            print(f"    {param}: threshold = {info['inflection_threshold']:.2f}, "
                  f"exclusion = {100*info['exclusion_at_inflection']:.1f}%")

        # Save
        for param in ['age', 'mass', 'metallicity', 'variability']:
            csv_path = os.path.join(BASE_DIR, f'sensitivity_{param}.csv')
            sensitivity_results[param].to_csv(csv_path, index=False)
            print(f"  Saved: {csv_path}")

    # - Step 12: Cross-match with Breakthrough Listen
    print("\n[Step 12] Cross-match with Breakthrough Listen...")

    bl_results = {}

    # 12a: Isaacson et al. 2017 (~1,700 nearby stars)
    print("\n  --- Isaacson et al. 2017 (PASP 129, 054501) ---")
    isaacson_df = load_bl_catalog_isaacson(BL_ISAACSON_FILE)
    if isaacson_df is not None:
        isaacson_overlap = analyze_bl_overlap(
            catalog, isaacson_df,
            tolerance_arcsec=5.0,
            catalog_name="Isaacson et al. 2017"
        )
        print_bl_comparison_report(isaacson_overlap)
        bl_results['isaacson'] = isaacson_overlap
    else:
        print("  Skipping Isaacson cross-match.")
        bl_results['isaacson'] = None

    # 12b: MeerKAT 1-million-star sample
    print("\n  --- MeerKAT 1-Million-Star Sample ---")
    meerkat_df = load_bl_catalog_meerkat(BL_MEERKAT_1M_FILE)
    if meerkat_df is not None:
        # MeerKAT has source_id — use direct cross-match (much faster + exact)
        has_sid = 'source_id' in meerkat_df.columns and 'source_id' in catalog.columns
        if has_sid:
            print(f"  Using source_id cross-match for {len(meerkat_df):,} MeerKAT targets...")
        else:
            print(f"  Cross-matching {len(meerkat_df):,} MeerKAT targets by coordinates (this may take a minute)...")
        meerkat_overlap = analyze_bl_overlap(
            catalog, meerkat_df,
            tolerance_arcsec=5.0,
            catalog_name="MeerKAT 1M",
            use_source_id=has_sid
        )
        print_bl_comparison_report(meerkat_overlap)
        bl_results['meerkat'] = meerkat_overlap
    else:
        print("  Skipping MeerKAT cross-match.")
        bl_results['meerkat'] = None

    if not bl_results.get('isaacson') and not bl_results.get('meerkat'):
        print("  No BL catalogs found. Place files in the working directory:")
        print(f"    - {BL_ISAACSON_FILE}")
        print(f"    - {BL_MEERKAT_1M_FILE}")

    # - Step 13: Generate figures
    print("\n[Step 13] Generating all figures...")
    generate_all_figures(catalog, sensitivity_results, bl_results)

    # - Step 14: Save results
    print("\n[Step 14] Saving results...")
    save_results(catalog)

    # - Final summary
    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)
    print(f"Total stars: {len(catalog):,}")
    print(f"Excluded: {stats['excluded_count']:,} ({100*stats['excluded_fraction']:.1f}%)")
    print(f"Retained: {stats['retained_count']:,} ({100*stats['retained_fraction']:.1f}%)")
    print(f"\nKey improvement: {age_comparison['difference']:,} additional stars retained")
    print("due to uncertainty-aware age filtering.")

    # BL summary
    if bl_results.get('isaacson'):
        isa = bl_results['isaacson']
        print(f"\nIsaacson et al. 2017: {isa['retained_count']:,} / {isa['matched_count']:,} retained")
    if bl_results.get('meerkat'):
        mkat = bl_results['meerkat']
        print(f"MeerKAT 1M: {mkat['retained_count']:,} / {mkat['matched_count']:,} retained")

    # GCNS summary
    if catalog['in_gcns'].any():
        gcns_ret = len(catalog[(catalog['in_gcns']) & (catalog['decision'] == 'RETAIN')])
        gcns_tot = catalog['in_gcns'].sum()
        print(f"\nGCNS: {gcns_ret:,} / {gcns_tot:,} retained")

    print(f"\nFigures:  {OUTPUT_FIGURES_DIR}")
    print(f"Catalogs: {os.path.join(OUTPUT_DIR, 'catalogs')}/")
    print(f"  → CSV | FITS")
    print("=" * 70)

    return catalog, stats, sensitivity_results, bl_results


if __name__ == "__main__":
    catalog, stats, sensitivity_results, bl_results = main()
