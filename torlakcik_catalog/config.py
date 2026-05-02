"""
Configuration file for the Torlakcık Catalog.

"""

import os

THRESHOLDS = {
    'mass': 1.5,
    'age': 3.0,
    'metallicity': -0.4,
    'variability': 0.01,
    'synthetic_variability': 0.001,  # sigma_G/G threshold for synthetic R6 proxy
}

SPECTRAL_BOUNDARIES = {
    'O': 30000, 'B': 10000, 'A': 7500, 'F0_F4': 6600,
    'F5_F9': 6000, 'G': 5200, 'K': 3700, 'M': 0,
}

REASON_CODES = {
    'R1': 'Mass > 1.5 M_sun',
    'R2': 'Age upper bound < 3 Gyr',
    'R3': 'Spectral type O/B/A/F0-F4',
    'R4': '[Fe/H] < -0.4',
    'R5': 'Multiplicity (NSS >= 1)',
    'R6': 'Photometric variability',
    'R7': 'M-dwarf activity',
}

# G-band zero-point (Riello et al. 2021)
G_ZEROPOINT = 25.6883657251

## paths
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
FITS_FILE = os.path.join(BASE_DIR, 'gaia_result_v2.fits')
GCNS_FILE = os.path.join(BASE_DIR, 'gcns.vot')  # VOTable format

OUTPUT_DIR = os.path.join(BASE_DIR, 'output')
OUTPUT_CATALOG = os.path.join(OUTPUT_DIR, 'catalogs', 'torlakcik_catalog.csv')
OUTPUT_RETAINED = os.path.join(OUTPUT_DIR, 'catalogs', 'retained_stars.csv')
OUTPUT_FIGURES_DIR = os.path.join(OUTPUT_DIR, 'figures')

# Breakthrough Listen target catalogs
BL_ISAACSON_FILE = os.path.join(BASE_DIR, 'J_PASP_129_E4501_table1.dat.fits')
BL_MEERKAT_1M_FILE = os.path.join(BASE_DIR, '1_million_sample_complete.csv')

os.makedirs(OUTPUT_FIGURES_DIR, exist_ok=True)
os.makedirs(os.path.join(OUTPUT_DIR, 'catalogs'), exist_ok=True)

TELESCOPES = {
    'GBT': {'beam_fwhm_arcmin': 9.0, 'frequency_ghz': 1.4},
    'Parkes': {'beam_fwhm_arcmin': 14.0, 'frequency_ghz': 1.4},
    'MeerKAT': {'beam_fwhm_arcmin': 58.0, 'frequency_ghz': 1.4},
    'VLA': {'beam_fwhm_arcmin': 30.0, 'frequency_ghz': 1.4},
}

SENSITIVITY_RANGES = {
    'age': (1.0, 6.0, 0.5),
    'mass': (0.8, 3.0, 0.1),
    'metallicity': (-1.0, 0.0, 0.05),
    'variability': (0.005, 0.05, 0.005),
}
