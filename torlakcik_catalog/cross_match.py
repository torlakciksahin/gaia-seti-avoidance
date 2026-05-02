"""
Cross-match module for comparing with Breakthrough Listen targets.

Supports two BL target catalogs:
  1. Isaacson et al. 2017 (PASP 129, 054501) — ~1,700 nearby stars (FITS)
  2. MeerKAT 1-million-star sample (CSV)

Handles special column formats:
  - Isaacson: sexagesimal component columns (RAh, RAm, RAs, DE-, DEd, DEm, DEs)
  - MeerKAT: 'decl' for Dec, 'ra' for RA; also supports source_id cross-match

"""

import os
import pandas as pd
import numpy as np
from astropy.table import Table
from astropy.coordinates import SkyCoord, Angle
import astropy.units as u
from typing import Tuple, Optional, Dict


# ---
#  Column name auto-detection
# ---
_RA_CANDIDATES = ['ra', 'RA', 'ra_deg', 'ra_degree', 'raj2000',
                  'RAJ2000', 'ra_j2000', 'right_ascension',
                  'RA_deg', 'ra_hour', 'ra_hms']
_DEC_CANDIDATES = ['dec', 'DEC', 'decl', 'dec_deg', 'dec_degree',
                   'dej2000', 'DEJ2000', 'dec_j2000', 'declination',
                   'DECL', 'DEC_deg', 'dec_dms']


def _find_column(columns: list, candidates: list) -> Optional[str]:
    """Find first matching column name from a candidate list."""
    cols_lower = {c.lower(): c for c in columns}
    for cand in candidates:
        if cand.lower() in cols_lower:
            return cols_lower[cand.lower()]
    for cand in candidates:
        for col in columns:
            if col.lower() == cand.lower():
                return col
    return None


def _has_sexagesimal_components(columns: list) -> bool:
    """Check if columns contain sexagesimal RA/Dec components (Isaacson format)."""
    cols_lower = {c.lower() for c in columns}
    return ('rah' in cols_lower or 'ra_h' in cols_lower) and \
           ('ded' in cols_lower or 'de_d' in cols_lower)


def _build_coords_from_components(df: pd.DataFrame) -> pd.DataFrame:
    """
    Build ra_deg and dec_deg from sexagesimal component columns.

    Handles Isaacson et al. 2017 format:
      RAh, RAm, RAs  → RA in hours, minutes, seconds
      DE-, DEd, DEm, DEs → Dec sign, degrees, arcmin, arcsec

    Also handles CDS/VizieR variants:
      RAh, RAm, RAs, DE+, DE-, DEd, DEm, DEs
    """
    # Find RA component columns
    rah_col = None
    for c in ['RAh', 'ra_h', 'RA_H']:
        if c in df.columns:
            rah_col = c
            break
    ram_col = None
    for c in ['RAm', 'ra_m', 'RA_M']:
        if c in df.columns:
            ram_col = c
            break
    ras_col = None
    for c in ['RAs', 'ra_s', 'RA_S']:
        if c in df.columns:
            ras_col = c
            break

    # Find Dec component columns
    ded_col = None
    for c in ['DEd', 'de_d', 'DE_D', 'DE+', 'DE-']:
        if c in df.columns:
            ded_col = c
            break
    dem_col = None
    for c in ['DEm', 'de_m', 'DE_M']:
        if c in df.columns:
            dem_col = c
            break
    des_col = None
    for c in ['DEs', 'de_s', 'DE_S']:
        if c in df.columns:
            des_col = c
            break

    # Find Dec sign column
    sign_col = None
    for c in ['DE-', 'DE+', 'DE-']:
        if c in df.columns:
            sign_col = c
            break

    if rah_col is None or ded_col is None:
        raise ValueError("Could not find sexagesimal component columns")

    print(f"  Sexagesimal components: RA=({rah_col}, {ram_col}, {ras_col}), "
          f"Dec sign={sign_col}, Dec=({ded_col}, {dem_col}, {des_col})")

    # Build RA in decimal hours → degrees
    # NOTE: .copy() is essential — FITS → astropy Table → pandas preserves
    #       memory-mapped read-only numpy arrays; in-place += would crash.
    ra_hours = df[rah_col].astype(float).values.copy()
    if ram_col is not None:
        ra_hours = ra_hours + df[ram_col].astype(float).values / 60.0
    if ras_col is not None:
        ra_hours = ra_hours + df[ras_col].astype(float).values / 3600.0

    df['ra_deg'] = ra_hours * 15.0  # hours → degrees

    # Build Dec in decimal degrees — always compute as positive first,
    # then apply sign at the end.
    dec_abs = df[ded_col].astype(float).values.copy()
    if dem_col is not None:
        dec_abs = dec_abs + df[dem_col].astype(float).values / 60.0
    if des_col is not None:
        dec_abs = dec_abs + df[des_col].astype(float).values / 3600.0

    # Determine sign from DE-/DE+ column
    # The column may contain string characters ('-'/'+'/'nan') or float
    # values (NaN, -1, 1).  We must handle both cases safely.
    is_negative = np.zeros(len(df), dtype=bool)

    if sign_col is not None:
        raw = df[sign_col].values
        for i, val in enumerate(raw):
            s = str(val)                       # guaranteed string
            if '-' in s:                       # catches '-', '-1.0', etc.
                is_negative[i] = True
            elif isinstance(val, (int, float, np.floating)) and val < 0:
                is_negative[i] = True          # catches numeric negatives
    else:
        # Fallback: scan for any DE-/DE+ column
        for col in df.columns:
            if col in ['DE-', 'DE+']:
                raw = df[col].values
                for i, val in enumerate(raw):
                    s = str(val)
                    if '-' in s or (isinstance(val, (int, float, np.floating)) and val < 0):
                        is_negative[i] = True
                break

    # Apply sign
    dec_deg = dec_abs.copy()
    dec_deg[is_negative] = -dec_abs[is_negative]

    df['dec_deg'] = dec_deg

    return df


# ---
#  Coordinate format detection & conversion
# ---

def _detect_coordinate_format(series: pd.Series, sample_size: int = 100) -> str:
    """Detect whether a coordinate column is in decimal degrees or sexagesimal."""
    sample = series.dropna().head(sample_size)
    if len(sample) == 0:
        return 'decimal_degrees'

    first_val = str(sample.iloc[0]).strip()

    if any(c in first_val for c in [':', 'h', 'm', 's']):
        return 'sexagesimal'

    try:
        val = float(first_val)
        if 0 <= val <= 360:
            return 'decimal_degrees'
        elif val > 24:
            return 'decimal_degrees'
    except (ValueError, TypeError):
        return 'sexagesimal'

    try:
        vals = sample.astype(float)
        if vals.max() <= 24 and vals.min() >= 0:
            return 'sexagesimal_hours'
        return 'decimal_degrees'
    except (ValueError, TypeError):
        return 'sexagesimal'


def _convert_ra_to_deg(series: pd.Series, fmt: str) -> pd.Series:
    """Convert RA column to decimal degrees."""
    result = pd.Series(np.nan, index=series.index, dtype=float)

    if fmt == 'decimal_degrees':
        try:
            result = series.astype(float)
            if result.max() <= 24:
                result = result * 15.0
        except (ValueError, TypeError):
            pass
    elif fmt == 'sexagesimal_hours':
        for i, val in series.items():
            try:
                result.iloc[i] = Angle(str(val), unit='hourangle').deg
            except Exception:
                try:
                    result.iloc[i] = float(val) * 15.0
                except Exception:
                    pass
    else:
        for i, val in series.items():
            try:
                result.iloc[i] = Angle(str(val), unit='hourangle').deg
            except Exception:
                try:
                    result.iloc[i] = float(val)
                except Exception:
                    pass
    return result


def _convert_dec_to_deg(series: pd.Series, fmt: str) -> pd.Series:
    """Convert Dec column to decimal degrees."""
    result = pd.Series(np.nan, index=series.index, dtype=float)

    if fmt == 'decimal_degrees':
        try:
            result = series.astype(float)
        except (ValueError, TypeError):
            pass
    else:
        for i, val in series.items():
            try:
                result.iloc[i] = Angle(str(val), unit='deg').deg
            except Exception:
                try:
                    result.iloc[i] = float(val)
                except Exception:
                    pass
    return result


# ---
#  Load BL target catalogs
# ---

def _parse_pm_unit(unit_str: Optional[str], col_name: str) -> str:
    """
    Parse a proper motion unit string into a canonical form.

    Recognized canonical forms: 'mas/yr', 'arcsec/yr', 's/yr', 'deg/yr'

    In VizieR Hipparcos-based catalogs, pmRA is often in 's/yr' (RA time-seconds
    per year), while pmDE is in 'arcsec/yr'.  This distinction is critical
    because 1 s of RA = 15 arcsec at the equator.
    """
    if unit_str is None or unit_str == '':
        return 'unknown'

    u = str(unit_str).strip().lower()

    # mas/yr variants
    if 'mas' in u or 'milliarcsec' in u:
        return 'mas/yr'

    # arcsec/yr variants
    if ('arcsec' in u or 'arc-sec' in u or u == '"/yr' or
            u == "''/yr" or ('/yr' in u and 'sec' in u and 'mas' not in u and 'milli' not in u)):
        return 'arcsec/yr'

    # s/yr (RA time-seconds per year) — VizieR Hipparcos convention
    if u in ('s/yr', 's yr-1', 'sec/yr', 'second/yr', 'time-s/yr'):
        return 's/yr'

    # deg/yr
    if 'deg' in u or 'degree' in u:
        return 'deg/yr'

    # Cannot determine → fall back to heuristic
    return 'unknown'


def _apply_proper_motion_correction(df: pd.DataFrame,
                                    target_epoch: float = 2016.0,
                                    pm_units: Optional[Dict[str, str]] = None) -> pd.DataFrame:
    """
    Apply proper motion correction to propagate coordinates to a target epoch.

    Isaacson et al. 2017 coordinates are on the Hipparcos/Tycho epoch (~J2000).
    Gaia DR3 reference epoch is J2016.0.  For nearby stars with high proper
    motion, 16+ years of drift can shift positions by 10-100+ arcseconds,
    causing coordinate cross-match to fail at small tolerances.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain 'ra_deg', 'dec_deg', and proper motion columns.
    target_epoch : float
        Target epoch in Julian years (default 2016.0 for Gaia DR3).
    pm_units : dict, optional
        Mapping of column name → unit string (e.g. {'pmRA': 'arcsec/yr'}).
        If provided, overrides heuristic unit detection.

    Returns
    -------
    pd.DataFrame
        DataFrame with added 'ra_deg_original', 'dec_deg_original' (pre-PM)
        and updated 'ra_deg', 'dec_deg' (propagated to target_epoch).
    """
    if pm_units is None:
        pm_units = {}

    # Save original coordinates
    df['ra_deg_original'] = df['ra_deg'].copy()
    df['dec_deg_original'] = df['dec_deg'].copy()

    # Detect proper motion columns
    pmra_col = None
    for c in ['pmRA', 'pmra', 'pm_ra', 'pmRA_J2000']:
        if c in df.columns:
            pmra_col = c
            break
    pmde_col = None
    for c in ['pmDE', 'pmde', 'pm_dec', 'pmDE_J2000', 'pmdec']:
        if c in df.columns:
            pmde_col = c
            break
    ep_col = None
    for c in ['Ep', 'epoch', 'Epoch', 'ref_epoch']:
        if c in df.columns:
            ep_col = c
            break

    if pmra_col is None or pmde_col is None:
        print("  No proper motion columns found — skipping PM correction")
        return df

    # Source epoch: use per-star epoch if available, else J2000.0
    if ep_col is not None:
        source_epoch = df[ep_col].astype(float).values
        n_unique = len(np.unique(source_epoch))
        print(f"  PM correction: epoch '{ep_col}' ({n_unique} unique), "
              f"target = J{target_epoch}")
    else:
        source_epoch = np.full(len(df), 2000.0)
        print(f"  PM correction: assuming J2000.0, target = J{target_epoch}")

    # Time difference in years
    dt = target_epoch - source_epoch

    # Read proper motion values
    pmra = pd.to_numeric(df[pmra_col], errors='coerce').values
    pmde = pd.to_numeric(df[pmde_col], errors='coerce').values

    # - Determine proper motion units
    # Priority: (1) explicit pm_units dict  →  (2) heuristic from values
    pmra_unit = _parse_pm_unit(pm_units.get(pmra_col), pmra_col)
    pmde_unit = _parse_pm_unit(pm_units.get(pmde_col), pmde_col)

    # Heuristic fallback: if either unit is unknown, guess from value range
    pm_max = max(np.nanmax(np.abs(pmra)), np.nanmax(np.abs(pmde)))
    if pmra_unit == 'unknown' or pmde_unit == 'unknown':
        if pm_max > 50:
            guessed = 'mas/yr'
        elif pm_max > 0.5:
            guessed = 'arcsec/yr'
        else:
            guessed = 's/yr'   # very small values → likely RA time-seconds
        if pmra_unit == 'unknown':
            pmra_unit = guessed
        if pmde_unit == 'unknown':
            pmde_unit = guessed
        print(f"  PM units (heuristic): pmRA={pmra_unit}, pmDE={pmde_unit} "
              f"(max |pm| = {pm_max:.1f})")
    else:
        print(f"  PM units (FITS header): pmRA={pmra_unit}, pmDE={pmde_unit}")

    # - Convert to deg/yr
    # pmRA conversion (RA direction)
    if pmra_unit == 'mas/yr':
        pmra_deg_per_yr = pmra / 3.6e6
        pmra_is_angular = True       # μ_α* — already projected angular velocity
    elif pmra_unit == 'arcsec/yr':
        pmra_deg_per_yr = pmra / 3600.0
        pmra_is_angular = True       # μ_α* — Hipparcos convention
    elif pmra_unit == 's/yr':
        # RA time-seconds/yr → angular: 1 s RA = 15 arcsec at equator
        # This gives μ_α* directly (projected angular velocity)
        pmra_deg_per_yr = pmra * 15.0 / 3600.0     # s/yr × 15"/s × (1°/3600") = deg/yr
        pmra_is_angular = True
    elif pmra_unit == 'deg/yr':
        pmra_deg_per_yr = pmra
        pmra_is_angular = False      # raw RA coordinate rate → needs ÷ cos(Dec)
    else:
        pmra_deg_per_yr = pmra / 3600.0
        pmra_is_angular = True

    # pmDE conversion (Dec direction) — always angular
    if pmde_unit == 'mas/yr':
        pmde_deg_per_yr = pmde / 3.6e6
    elif pmde_unit == 'arcsec/yr':
        pmde_deg_per_yr = pmde / 3600.0
    elif pmde_unit == 'deg/yr':
        pmde_deg_per_yr = pmde
    else:
        pmde_deg_per_yr = pmde / 3600.0

    # - Apply correction
    # If pmRA is μ_α* (projected angular PM), the RA coordinate shift is:
    #   Δα* = μ_α* × Δt           (angular shift on the sky)
    #   Δα  = μ_α* / cos(δ) × Δt  (RA coordinate change)
    dec_rad = np.deg2rad(df['dec_deg'].values)
    cos_dec = np.cos(dec_rad)
    cos_dec_safe = np.where(np.abs(cos_dec) < 1e-10, 1e-10, cos_dec)

    if pmra_is_angular:
        # pmRA is μ_α* — need to divide by cos(Dec) to get RA coordinate shift
        ra_shift_deg = pmra_deg_per_yr * dt / cos_dec_safe
    else:
        # pmRA is already RA coordinate rate — no cos(Dec) correction needed
        ra_shift_deg = pmra_deg_per_yr * dt

    dec_shift_deg = pmde_deg_per_yr * dt

    df['ra_deg'] = df['ra_deg'].values + ra_shift_deg
    df['dec_deg'] = df['dec_deg'].values + dec_shift_deg

    # - Statistics
    # Total angular shift on the sky (great-circle distance)
    total_shift_arcsec = np.sqrt((pmra_deg_per_yr * dt)**2 +
                                 (pmde_deg_per_yr * dt)**2) * 3600.0

    median_shift = np.nanmedian(total_shift_arcsec)
    max_shift = np.nanmax(total_shift_arcsec)
    n_large = np.nansum(total_shift_arcsec > 5)

    print(f"  PM correction: median={median_shift:.1f}\", "
          f"max={max_shift:.1f}\", "
          f"{n_large} stars > 5\"")

    return df


def load_bl_catalog_isaacson(filepath: str) -> Optional[pd.DataFrame]:
    """
    Load Isaacson et al. 2017 (PASP 129, 054501) target list.

    Handles sexagesimal component columns (RAh, RAm, RAs, DE-, DEd, DEm, DEs)
    as well as standard single-column RA/Dec formats.

    Applies proper motion correction from catalog epoch to Gaia DR3 epoch
    (J2016.0) to ensure accurate cross-matching of nearby high-PM stars.

    Proper motion units are read from the FITS header TUNIT keywords when
    available.  Falls back to heuristic detection if TUNIT is absent.

    Reference: Isaacson, H., Siemion, A.P.V., Marcy, G.W., et al. 2017,
    PASP, 129, 054501.
    """
    if not os.path.exists(filepath):
        print(f"  Isaacson catalog not found: {filepath}")
        return None

    try:
        t = Table.read(filepath)

        # - Extract proper motion units from FITS header (TUNIT)
        pm_units = {}
        for col in t.colnames:
            if col in ('pmRA', 'pmra', 'pm_ra', 'pmRA_J2000',
                       'pmDE', 'pmde', 'pm_dec', 'pmDE_J2000', 'pmdec'):
                if t[col].unit is not None:
                    pm_units[col] = str(t[col].unit)
                    print(f"  FITS TUNIT: {col} = {t[col].unit}")

        df = t.to_pandas()
        print(f"  Loaded Isaacson et al. 2017: {len(df)} targets")
        print(f"  Columns: {list(df.columns)}")

        # - Method 1: Sexagesimal components (Isaacson format)
        if _has_sexagesimal_components(df.columns):
            print("  Detected sexagesimal component columns (RAh/RAm/RAs, DE-/DEd/DEm/DEs)")
            df = _build_coords_from_components(df)
        else:
            # - Method 2: Standard single RA/Dec columns
            ra_col = _find_column(df.columns, _RA_CANDIDATES)
            dec_col = _find_column(df.columns, _DEC_CANDIDATES)

            if ra_col is None or dec_col is None:
                print(f"  WARNING: Could not find RA/Dec columns in Isaacson catalog")
                print(f"  Available columns: {list(df.columns)}")
                return None

            print(f"  Detected RA column: '{ra_col}', Dec column: '{dec_col}'")
            ra_format = _detect_coordinate_format(df[ra_col])
            dec_format = _detect_coordinate_format(df[dec_col])
            print(f"  RA format: {ra_format}, Dec format: {dec_format}")
            df['ra_deg'] = _convert_ra_to_deg(df[ra_col], ra_format)
            df['dec_deg'] = _convert_dec_to_deg(df[dec_col], dec_format)

        # - Proper motion correction to Gaia DR3 epoch (J2016.0)
        # Isaacson stars are nearby with high proper motion; without this
        # correction, cross-match at 5" tolerance misses many targets.
        df = _apply_proper_motion_correction(df, target_epoch=2016.0,
                                             pm_units=pm_units)

        # Validate
        valid = df.dropna(subset=['ra_deg', 'dec_deg'])
        n_valid = len(valid)
        print(f"  Valid coordinates (after PM correction): {n_valid} / {len(df)}")

        if n_valid == 0:
            print("  WARNING: No valid coordinates found after conversion")
            return None

        # Print sample for verification
        if n_valid > 0:
            sample = valid.head(3)
            for _, row in sample.iterrows():
                ra_orig = row.get('ra_deg_original', row['ra_deg'])
                dec_orig = row.get('dec_deg_original', row['dec_deg'])
                print(f"    Sample: RA={row['ra_deg']:.4f}° (was {ra_orig:.4f}°), "
                      f"Dec={row['dec_deg']:.4f}° (was {dec_orig:.4f}°)"
                      f"  (Star: {row.get('Star', 'N/A')})")

        return df

    except Exception as e:
        print(f"  ERROR loading Isaacson catalog: {e}")
        import traceback
        traceback.print_exc()
        return None


def load_bl_catalog_meerkat(filepath: str) -> Optional[pd.DataFrame]:
    """
    Load MeerKAT 1-million-star sample (CSV).

    Handles 'decl' column name for Dec.
    Also detects 'source_id' column for direct cross-match option.
    """
    if not os.path.exists(filepath):
        print(f"  MeerKAT 1M catalog not found: {filepath}")
        return None

    try:
        df = pd.read_csv(filepath)
        print(f"  Loaded MeerKAT 1M sample: {len(df):,} targets")
        print(f"  Columns: {list(df.columns)[:20]}{'...' if len(df.columns) > 20 else ''}")

        # Check for source_id (allows direct cross-match)
        has_source_id = 'source_id' in df.columns
        if has_source_id:
            print(f"  Detected 'source_id' column — direct cross-match available")

        # Find RA/Dec columns
        ra_col = _find_column(df.columns, _RA_CANDIDATES)
        dec_col = _find_column(df.columns, _DEC_CANDIDATES)

        if ra_col is None or dec_col is None:
            print(f"  WARNING: Could not find RA/Dec columns in MeerKAT catalog")
            print(f"  Available columns: {list(df.columns)}")
            return None

        print(f"  Detected RA column: '{ra_col}', Dec column: '{dec_col}'")

        # Convert to decimal degrees
        ra_format = _detect_coordinate_format(df[ra_col])
        dec_format = _detect_coordinate_format(df[dec_col])
        print(f"  RA format: {ra_format}, Dec format: {dec_format}")

        df['ra_deg'] = _convert_ra_to_deg(df[ra_col], ra_format)
        df['dec_deg'] = _convert_dec_to_deg(df[dec_col], dec_format)

        # Validate
        valid = df.dropna(subset=['ra_deg', 'dec_deg'])
        print(f"  Valid coordinates: {len(valid):,} / {len(df):,}")

        if len(valid) == 0:
            print("  WARNING: No valid coordinates found after conversion")
            return None

        # Print sample for verification
        sample = valid.head(3)
        for _, row in sample.iterrows():
            print(f"    Sample: RA={row['ra_deg']:.4f}°, Dec={row['dec_deg']:.4f}°"
                  f"  (source_id: {row.get('source_id', 'N/A')})")

        return df

    except Exception as e:
        print(f"  ERROR loading MeerKAT catalog: {e}")
        return None


# ---
#  Cross-match by source_id (fast, exact)
# ---

def cross_match_by_source_id(catalog: pd.DataFrame,
                              targets: pd.DataFrame,
                              sid_col_cat: str = 'source_id',
                              sid_col_target: str = 'source_id') -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Cross-match by Gaia source_id (exact match, very fast).

    Returns
    -------
    tuple
        (matched_catalog_rows, matched_target_rows)
    """
    cat_ids = set(catalog[sid_col_cat].dropna().astype(np.int64))
    tgt_ids = targets[sid_col_target].dropna().astype(np.int64)

    matched_mask = tgt_ids.isin(cat_ids)
    matched_targets = targets[matched_mask].copy()

    # Get corresponding catalog rows
    matched_target_ids = set(matched_targets[sid_col_target].astype(np.int64))
    matched_catalog = catalog[catalog[sid_col_cat].astype(np.int64).isin(matched_target_ids)].copy()

    print(f"  source_id cross-match: {len(matched_targets):,} / {len(targets):,} "
          f"({100*len(matched_targets)/len(targets):.1f}%) matched")

    return matched_catalog, matched_targets


# ---
#  Cross-match by coordinates
# ---

def cross_match_coordinates(catalog: pd.DataFrame,
                            targets: pd.DataFrame,
                            ra_col_cat: str = 'ra',
                            dec_col_cat: str = 'dec',
                            ra_col_target: str = 'ra_deg',
                            dec_col_target: str = 'dec_deg',
                            tolerance_arcsec: float = 5.0) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Cross-match two catalogs by coordinates using astropy."""
    cat_valid = catalog.dropna(subset=[ra_col_cat, dec_col_cat]).copy()
    tgt_valid = targets.dropna(subset=[ra_col_target, dec_col_target]).copy()

    if len(cat_valid) == 0 or len(tgt_valid) == 0:
        return pd.DataFrame(), pd.DataFrame(), targets.copy()

    print(f"  Cross-matching {len(tgt_valid):,} targets against {len(cat_valid):,} catalog stars...")

    cat_coords = SkyCoord(ra=cat_valid[ra_col_cat].values * u.degree,
                          dec=cat_valid[dec_col_cat].values * u.degree)
    tgt_coords = SkyCoord(ra=tgt_valid[ra_col_target].values * u.degree,
                          dec=tgt_valid[dec_col_target].values * u.degree)

    idx, d2d, _ = tgt_coords.match_to_catalog_sky(cat_coords)
    match_mask = d2d < tolerance_arcsec * u.arcsec

    matched_targets = tgt_valid[match_mask].copy()
    matched_catalog = cat_valid.iloc[idx[match_mask]].copy()
    matched_targets['separation_arcsec'] = d2d[match_mask].arcsec
    unmatched_targets = tgt_valid[~match_mask].copy()

    print(f"  Matched: {len(matched_targets):,} / {len(tgt_valid):,} "
          f"({100*len(matched_targets)/len(tgt_valid):.1f}%) within {tolerance_arcsec}\"")

    return matched_catalog.reset_index(drop=True), \
           matched_targets.reset_index(drop=True), \
           unmatched_targets.reset_index(drop=True)


# ---
#  Analyze overlap with Torlakcık Catalog
# ---

def analyze_bl_overlap(catalog: pd.DataFrame,
                       bl_targets: pd.DataFrame,
                       tolerance_arcsec: float = 5.0,
                       catalog_name: str = "BL",
                       use_source_id: bool = False) -> dict:
    """
    Analyze overlap between Torlakcık Catalog and BL targets.

    Parameters
    ----------
    use_source_id : bool
        If True and both catalogs have source_id, use direct ID match
        instead of coordinate matching.
    """
    # Try source_id match first if available
    if use_source_id and 'source_id' in catalog.columns and 'source_id' in bl_targets.columns:
        print(f"  Using source_id cross-match for {catalog_name}...")
        matched_cat, matched_bl = cross_match_by_source_id(catalog, bl_targets)
        unmatched_bl = bl_targets[~bl_targets['source_id'].astype(np.int64).isin(
            set(matched_bl['source_id'].astype(np.int64)) if len(matched_bl) > 0 else set()
        )].copy()
    else:
        # Coordinate-based cross-match
        matched_cat, matched_bl, unmatched_bl = cross_match_coordinates(
            catalog, bl_targets, tolerance_arcsec=tolerance_arcsec)

    # Analyze decisions for matched stars
    excluded = matched_cat[matched_cat['decision'] == 'EXCLUDE'] if len(matched_cat) > 0 else pd.DataFrame()
    retained = matched_cat[matched_cat['decision'] == 'RETAIN'] if len(matched_cat) > 0 else pd.DataFrame()

    # Analyze exclusion reasons
    reason_counts = {}
    for code in ['R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7']:
        count = len(excluded[excluded['reason_code'].str.contains(code, na=False)]) if len(excluded) > 0 else 0
        reason_counts[code] = count

    # Spectral type breakdown
    sp_breakdown = {}
    if len(matched_cat) > 0 and 'spectral_type' in matched_cat.columns:
        for sp in matched_cat['spectral_type'].unique():
            sp_df = matched_cat[matched_cat['spectral_type'] == sp]
            sp_excl = len(sp_df[sp_df['decision'] == 'EXCLUDE'])
            sp_breakdown[sp] = {
                'total': len(sp_df),
                'excluded': sp_excl,
                'retained': len(sp_df) - sp_excl,
            }

    return {
        'catalog_name': catalog_name,
        'bl_total': len(bl_targets),
        'matched_count': len(matched_cat),
        'unmatched_count': len(unmatched_bl),
        'match_rate': len(matched_cat) / len(bl_targets) if len(bl_targets) > 0 else 0,
        'excluded_count': len(excluded),
        'retained_count': len(retained),
        'excluded_fraction': len(excluded) / len(matched_cat) if len(matched_cat) > 0 else 0,
        'reason_counts': reason_counts,
        'spectral_breakdown': sp_breakdown,
        'matched_catalog': matched_cat,
        'matched_bl': matched_bl,
        'unmatched_bl': unmatched_bl,
    }


def print_bl_comparison_report(overlap_results: dict):
    """Print formatted BL comparison report."""
    name = overlap_results.get('catalog_name', 'BL')
    print(f"\n  {'='*58}")
    print(f"  BREAKTHROUGH LISTEN CROSS-MATCH: {name}")
    print(f"  {'='*58}")
    print(f"  Total BL targets:        {overlap_results['bl_total']:,}")
    print(f"  Matched in Gaia DR3:     {overlap_results['matched_count']:,} "
          f"({100*overlap_results['match_rate']:.1f}%)")
    print(f"  Unmatched (no Gaia data): {overlap_results['unmatched_count']:,}")

    if overlap_results['matched_count'] > 0:
        print(f"\n  Of matched {name} targets:")
        print(f"    Excluded by our model: {overlap_results['excluded_count']:,} "
              f"({100*overlap_results['excluded_fraction']:.1f}%)")
        print(f"    Retained by our model: {overlap_results['retained_count']:,} "
              f"({100*(1-overlap_results['excluded_fraction']):.1f}%)")

        print(f"\n  Primary exclusion reasons for {name} targets:")
        reason_names = {
            'R1': 'Mass', 'R2': 'Age', 'R3': 'Spectral type',
            'R4': 'Metallicity', 'R5': 'Multiplicity',
            'R6': 'Variability', 'R7': 'M-dwarf activity',
        }
        for code, count in overlap_results['reason_counts'].items():
            if count > 0:
                print(f"    {code} ({reason_names[code]}): {count}")

        if overlap_results.get('spectral_breakdown'):
            print(f"\n  Spectral type breakdown for matched {name} targets:")
            for sp in ['O', 'B', 'A', 'F0-F4', 'F5-F9', 'G', 'K', 'M']:
                if sp in overlap_results['spectral_breakdown']:
                    info = overlap_results['spectral_breakdown'][sp]
                    rate = 100 * info['excluded'] / info['total'] if info['total'] > 0 else 0
                    print(f"    {sp}: {info['total']} matched, {info['excluded']} excluded ({rate:.0f}%)")
