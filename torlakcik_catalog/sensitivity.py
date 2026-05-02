"""
Sensitivity analysis module for the Torlakcık Catalog.

Performs one-at-a-time sensitivity analysis on all continuous thresholds.
Uses vectorized avoidance model for performance.

"""

import pandas as pd
import numpy as np
from typing import List, Tuple
from avoidance_model import apply_avoidance_model, calculate_exclusion_statistics
from config import THRESHOLDS, SENSITIVITY_RANGES


def sensitivity_analysis_age(df: pd.DataFrame,
                              age_range: Tuple[float, float, float] = (1.0, 6.0, 0.5)
                              ) -> pd.DataFrame:
    """
    Perform sensitivity analysis on the age threshold.

    Parameters
    ----------
    df : pd.DataFrame
        Full stellar catalog
    age_range : tuple
        (start, stop, step) in Gyr

    Returns
    -------
    pd.DataFrame
        Sensitivity analysis results
    """
    results = []
    ages = np.arange(age_range[0], age_range[1] + age_range[2], age_range[2])

    for age_threshold in ages:
        custom_thresholds = THRESHOLDS.copy()
        custom_thresholds['age'] = age_threshold

        result_df = apply_avoidance_model(df, use_age_upper=True, custom_thresholds=custom_thresholds)
        stats = calculate_exclusion_statistics(result_df)

        results.append({
            'threshold': age_threshold,
            'excluded_count': stats['excluded_count'],
            'excluded_fraction': stats['excluded_fraction'],
            'retained_count': stats['retained_count'],
        })

    return pd.DataFrame(results)


def sensitivity_analysis_mass(df: pd.DataFrame,
                               mass_range: Tuple[float, float, float] = (0.8, 3.0, 0.1)
                               ) -> pd.DataFrame:
    """
    Perform sensitivity analysis on the mass threshold.
    """
    results = []
    masses = np.arange(mass_range[0], mass_range[1] + mass_range[2], mass_range[2])

    for mass_threshold in masses:
        custom_thresholds = THRESHOLDS.copy()
        custom_thresholds['mass'] = mass_threshold

        result_df = apply_avoidance_model(df, use_age_upper=True, custom_thresholds=custom_thresholds)
        stats = calculate_exclusion_statistics(result_df)

        results.append({
            'threshold': mass_threshold,
            'excluded_count': stats['excluded_count'],
            'excluded_fraction': stats['excluded_fraction'],
            'retained_count': stats['retained_count'],
        })

    return pd.DataFrame(results)


def sensitivity_analysis_metallicity(df: pd.DataFrame,
                                      feh_range: Tuple[float, float, float] = (-1.0, 0.0, 0.05)
                                      ) -> pd.DataFrame:
    """
    Perform sensitivity analysis on the metallicity threshold.
    """
    results = []
    metallicities = np.arange(feh_range[0], feh_range[1] + feh_range[2], feh_range[2])

    for feh_threshold in metallicities:
        custom_thresholds = THRESHOLDS.copy()
        custom_thresholds['metallicity'] = feh_threshold

        result_df = apply_avoidance_model(df, use_age_upper=True, custom_thresholds=custom_thresholds)
        stats = calculate_exclusion_statistics(result_df)

        results.append({
            'threshold': feh_threshold,
            'excluded_count': stats['excluded_count'],
            'excluded_fraction': stats['excluded_fraction'],
            'retained_count': stats['retained_count'],
        })

    return pd.DataFrame(results)


def sensitivity_analysis_variability(df: pd.DataFrame,
                                      var_range: Tuple[float, float, float] = (0.005, 0.05, 0.005)
                                      ) -> pd.DataFrame:
    """
    Perform sensitivity analysis on the variability threshold.
    """
    results = []
    var_thresholds = np.arange(var_range[0], var_range[1] + var_range[2], var_range[2])

    for var_threshold in var_thresholds:
        custom_thresholds = THRESHOLDS.copy()
        custom_thresholds['variability'] = var_threshold

        result_df = apply_avoidance_model(df, use_age_upper=True, custom_thresholds=custom_thresholds)
        stats = calculate_exclusion_statistics(result_df)

        results.append({
            'threshold': var_threshold,
            'excluded_count': stats['excluded_count'],
            'excluded_fraction': stats['excluded_fraction'],
            'retained_count': stats['retained_count'],
        })

    return pd.DataFrame(results)


def run_full_sensitivity_analysis(df: pd.DataFrame) -> dict:
    """
    Run complete sensitivity analysis for all parameters.

    Parameters
    ----------
    df : pd.DataFrame
        Full stellar catalog

    Returns
    -------
    dict
        Dictionary of DataFrames for each parameter
    """
    print("Running sensitivity analysis...")
    print("  - Age threshold...")
    age_results = sensitivity_analysis_age(df)

    print("  - Mass threshold...")
    mass_results = sensitivity_analysis_mass(df)

    print("  - Metallicity threshold...")
    feh_results = sensitivity_analysis_metallicity(df)

    print("  - Variability threshold...")
    var_results = sensitivity_analysis_variability(df)

    return {
        'age': age_results,
        'mass': mass_results,
        'metallicity': feh_results,
        'variability': var_results,
    }


def find_inflection_points(results: dict) -> dict:
    """
    Identify inflection points in sensitivity curves.

    Parameters
    ----------
    results : dict
        Output from run_full_sensitivity_analysis

    Returns
    -------
    dict
        Inflection point analysis
    """
    inflection_info = {}

    for param, df in results.items():
        # Find the point of maximum slope change
        fractions = df['excluded_fraction'].values
        thresholds = df['threshold'].values

        if len(fractions) > 2:
            # Calculate second derivative
            second_deriv = np.diff(fractions, n=2)
            inflection_idx = np.argmax(np.abs(second_deriv)) + 1
            inflection_threshold = thresholds[inflection_idx]

            inflection_info[param] = {
                'inflection_threshold': inflection_threshold,
                'exclusion_at_inflection': fractions[inflection_idx],
            }

    return inflection_info
