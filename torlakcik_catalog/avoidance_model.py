"""
Torlakcık Catalog: SETI Target Avoidance Model

This module implements the parametric avoidance function for SETI target
selection based on seven stellar criteria.

VECTORIZED VERSION: Uses numpy/pandas vectorized operations instead of
row-by-row iteration for ~100x speedup on large catalogs.

"""

import pandas as pd
import numpy as np
from typing import Tuple, List, Optional
from config import THRESHOLDS, SPECTRAL_BOUNDARIES


def classify_spectral_type(teff: float) -> str:
    """
    Classify stellar spectral type based on effective temperature.

    Parameters
    ----------
    teff : float
        Effective temperature in Kelvin

    Returns
    -------
    str
        Spectral type classification
    """
    if pd.isna(teff):
        return 'Unknown'

    if teff >= SPECTRAL_BOUNDARIES['O']:
        return 'O'
    elif teff >= SPECTRAL_BOUNDARIES['B']:
        return 'B'
    elif teff >= SPECTRAL_BOUNDARIES['A']:
        return 'A'
    elif teff >= SPECTRAL_BOUNDARIES['F0_F4']:
        return 'F0-F4'
    elif teff >= SPECTRAL_BOUNDARIES['F5_F9']:
        return 'F5-F9'
    elif teff >= SPECTRAL_BOUNDARIES['G']:
        return 'G'
    elif teff >= SPECTRAL_BOUNDARIES['K']:
        return 'K'
    else:
        return 'M'


def classify_spectral_types_vectorized(teff_series: pd.Series) -> pd.Series:
    """
    Vectorized spectral type classification for entire column.

    Parameters
    ----------
    teff_series : pd.Series
        Series of effective temperatures in Kelvin

    Returns
    -------
    pd.Series
        Series of spectral type strings
    """
    result = pd.Series('Unknown', index=teff_series.index, dtype='object')
    valid = teff_series.notna()

    t = teff_series[valid]
    result[valid] = np.select(
        [
            t >= SPECTRAL_BOUNDARIES['O'],
            t >= SPECTRAL_BOUNDARIES['B'],
            t >= SPECTRAL_BOUNDARIES['A'],
            t >= SPECTRAL_BOUNDARIES['F0_F4'],
            t >= SPECTRAL_BOUNDARIES['F5_F9'],
            t >= SPECTRAL_BOUNDARIES['G'],
            t >= SPECTRAL_BOUNDARIES['K'],
        ],
        ['O', 'B', 'A', 'F0-F4', 'F5-F9', 'G', 'K'],
        default='M'
    )

    return result


def evaluate_star(row: pd.Series,
                  use_age_upper: bool = True,
                  custom_thresholds: Optional[dict] = None) -> Tuple[str, str]:
    """
    Evaluate a single star against all avoidance criteria.
    Kept for backward compatibility, but apply_avoidance_model_vectorized
    is recommended for large catalogs.

    Parameters
    ----------
    row : pd.Series
        Star data row containing all required parameters
    use_age_upper : bool
        If True, use age upper bound for R2 criterion
    custom_thresholds : dict, optional
        Custom threshold values for sensitivity analysis

    Returns
    -------
    tuple
        (decision, reason_codes) where decision is 'EXCLUDE' or 'RETAIN'
    """
    thresholds = custom_thresholds if custom_thresholds else THRESHOLDS
    reasons = []

    # R1 - Mass criterion
    mass = row.get('mass_flame_spec')
    if pd.notna(mass) and mass > thresholds['mass']:
        reasons.append('R1')

    # R2 - Age criterion (using upper bound for uncertainty-aware filtering)
    if use_age_upper:
        age_upper = row.get('age_flame_spec_upper')
        if pd.notna(age_upper) and age_upper < thresholds['age']:
            reasons.append('R2')
    else:
        age = row.get('age_flame_spec')
        if pd.notna(age) and age < thresholds['age']:
            reasons.append('R2')

    # R3 - Spectral type criterion
    spectral_type = row.get('spectral_type', 'Unknown')
    if spectral_type in ['O', 'B', 'A', 'F0-F4']:
        reasons.append('R3')

    # R4 - Metallicity criterion
    metallicity = row.get('mh_gspphot')
    if pd.notna(metallicity) and metallicity < thresholds['metallicity']:
        reasons.append('R4')

    # R5 - Multiplicity criterion
    nss = row.get('non_single_star')
    if pd.notna(nss) and nss >= 1:
        reasons.append('R5')

    # R6 - Photometric variability criterion
    range_mag = row.get('range_mag_g_fov')
    var_flag = row.get('phot_variable_flag')
    if (pd.notna(range_mag) and range_mag > thresholds['variability']) or var_flag == 'VARIABLE':
        reasons.append('R6')

    # R7 - M-dwarf activity criterion
    if spectral_type == 'M':
        rot_mod = row.get('in_vari_rotation_modulation')
        short_ts = row.get('in_vari_short_timescale')
        if (pd.notna(rot_mod) and rot_mod == True) or (pd.notna(short_ts) and short_ts == True):
            reasons.append('R7')

    if reasons:
        return 'EXCLUDE', ';'.join(reasons)
    else:
        return 'RETAIN', ''


def apply_avoidance_model(df: pd.DataFrame,
                          use_age_upper: bool = True,
                          custom_thresholds: Optional[dict] = None) -> pd.DataFrame:
    """
    Apply the avoidance model to a stellar catalog.
    Uses VECTORIZED operations for performance.

    Parameters
    ----------
    df : pd.DataFrame
        Stellar catalog with required columns
    use_age_upper : bool
        If True, use age upper bound for R2 criterion
    custom_thresholds : dict, optional
        Custom threshold values for sensitivity analysis

    Returns
    -------
    pd.DataFrame
        Catalog with added 'spectral_type', 'decision', and 'reason_code' columns
    """
    df = df.copy()
    thresholds = custom_thresholds if custom_thresholds else THRESHOLDS

    # Classify spectral types (vectorized)
    print("Classifying spectral types...")
    df['spectral_type'] = classify_spectral_types_vectorized(df['teff_gspphot'])

    # Apply avoidance criteria (vectorized)
    print("Applying avoidance criteria...")
    reason_masks = {}

    # R1 - Mass criterion: mass > threshold
    mass = df['mass_flame_spec']
    reason_masks['R1'] = mass.notna() & (mass > thresholds['mass'])

    # R2 - Age criterion
    if use_age_upper:
        age_col = df['age_flame_spec_upper']
        reason_masks['R2'] = age_col.notna() & (age_col < thresholds['age'])
    else:
        age_col = df['age_flame_spec']
        reason_masks['R2'] = age_col.notna() & (age_col < thresholds['age'])

    # R3 - Spectral type criterion
    reason_masks['R3'] = df['spectral_type'].isin(['O', 'B', 'A', 'F0-F4'])

    # R4 - Metallicity criterion
    mh = df['mh_gspphot']
    reason_masks['R4'] = mh.notna() & (mh < thresholds['metallicity'])

    # R5 - Multiplicity criterion
    nss = df['non_single_star']
    reason_masks['R5'] = nss.notna() & (nss >= 1)

    # R6 - Photometric variability criterion
    range_mag = df['range_mag_g_fov']
    var_flag = df['phot_variable_flag']
    reason_masks['R6'] = (range_mag.notna() & (range_mag > thresholds['variability'])) | (var_flag == 'VARIABLE')

    # R7 - M-dwarf activity criterion
    is_m_dwarf = df['spectral_type'] == 'M'
    rot_mod = df.get('in_vari_rotation_modulation', pd.Series(False, index=df.index))
    short_ts = df.get('in_vari_short_timescale', pd.Series(False, index=df.index))
    # Handle boolean-like values
    rot_mod_bool = rot_mod.fillna(False).astype(bool)
    short_ts_bool = short_ts.fillna(False).astype(bool)
    reason_masks['R7'] = is_m_dwarf & (rot_mod_bool | short_ts_bool)

    # Build reason_code string and decision
    reason_code_parts = []
    for code in ['R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7']:
        mask = reason_masks[code]
        part = np.where(mask, code, '')
        reason_code_parts.append(part)

    # Combine reason codes
    stacked = np.stack(reason_code_parts, axis=0)
    # Join non-empty codes with ';'
    reason_codes = []
    n_criteria = stacked.shape[0]
    for i in range(len(df)):
        codes = [stacked[j, i] for j in range(n_criteria) if stacked[j, i] != '']
        reason_codes.append(';'.join(codes))

    df['reason_code'] = reason_codes
    df['decision'] = np.where(df['reason_code'] == '', 'RETAIN', 'EXCLUDE')

    return df


def calculate_exclusion_statistics(df: pd.DataFrame) -> dict:
    """
    Calculate comprehensive exclusion statistics.

    Parameters
    ----------
    df : pd.DataFrame
        Catalog with 'decision' and 'reason_code' columns

    Returns
    -------
    dict
        Dictionary containing all statistics
    """
    total = len(df)
    excluded = df[df['decision'] == 'EXCLUDE']
    retained = df[df['decision'] == 'RETAIN']

    stats = {
        'total_stars': total,
        'excluded_count': len(excluded),
        'retained_count': len(retained),
        'excluded_fraction': len(excluded) / total,
        'retained_fraction': len(retained) / total,
    }

    # Per-criterion counts
    for code in ['R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7']:
        count = excluded[excluded['reason_code'].str.contains(code, na=False)].shape[0]
        stats[f'{code}_count'] = count
        stats[f'{code}_fraction'] = count / total

    # Per-spectral-type statistics
    for sp_type in df['spectral_type'].unique():
        sp_total = len(df[df['spectral_type'] == sp_type])
        sp_excluded = len(df[(df['spectral_type'] == sp_type) & (df['decision'] == 'EXCLUDE')])
        stats[f'{sp_type}_total'] = sp_total
        stats[f'{sp_type}_excluded'] = sp_excluded
        stats[f'{sp_type}_exclusion_rate'] = sp_excluded / sp_total if sp_total > 0 else 0

    return stats


def print_statistics(stats: dict):
    """
    Print formatted exclusion statistics.
    """
    print("\n" + "=" * 60)
    print("EXCLUSION STATISTICS")
    print("=" * 60)
    print(f"Total stars: {stats['total_stars']:,}")
    print(f"Excluded: {stats['excluded_count']:,} ({100*stats['excluded_fraction']:.1f}%)")
    print(f"Retained: {stats['retained_count']:,} ({100*stats['retained_fraction']:.1f}%)")

    print("\n" + "-" * 60)
    print("Per-Criterion Exclusion Counts")
    print("-" * 60)
    criterion_names = {
        'R1': 'Mass > 1.5 M_sun',
        'R2': 'Age upper < 3 Gyr',
        'R3': 'Spectral O/B/A/F0-F4',
        'R4': '[Fe/H] < -0.4',
        'R5': 'Multiplicity',
        'R6': 'Phot. variability',
        'R7': 'M-dwarf activity',
    }
    for code, name in criterion_names.items():
        count = stats.get(f'{code}_count', 0)
        fraction = stats.get(f'{code}_fraction', 0)
        print(f"{code}: {count:,} ({100*fraction:.1f}%) - {name}")

    print("\n" + "-" * 60)
    print("Per-Spectral-Type Exclusion Rates")
    print("-" * 60)
    for sp_type in ['O', 'B', 'A', 'F0-F4', 'F5-F9', 'G', 'K', 'M', 'Unknown']:
        if f'{sp_type}_total' in stats:
            total = stats[f'{sp_type}_total']
            rate = stats[f'{sp_type}_exclusion_rate']
            print(f"{sp_type}: {100*rate:.1f}% excluded (N={total:,})")
