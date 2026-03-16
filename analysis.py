from astroquery.gaia import Gaia
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# Recording Folder
save_dir = os.path.join(os.path.expanduser('~'), 'Desktop', 'seti_project')
os.makedirs(save_dir, exist_ok=True)

# Gaia DR3 Data
query = """
SELECT
    gs.source_id, gs.ra, gs.dec,
    ap.teff_gspphot, ap.mh_gspphot, ap.logg_gspphot,
    fl.mass_flame_spec, fl.age_flame_spec,
    gs.non_single_star, gs.phot_variable_flag,
    vs.range_mag_g_fov,
    vs.in_vari_rotation_modulation,
    vs.in_vari_short_timescale
FROM gaiadr3.gaia_source AS gs
JOIN gaiadr3.astrophysical_parameters AS ap ON gs.source_id = ap.source_id
JOIN gaiadr3.astrophysical_parameters_supp AS fl ON gs.source_id = fl.source_id
LEFT JOIN gaiadr3.vari_summary AS vs ON gs.source_id = vs.source_id
WHERE ap.teff_gspphot IS NOT NULL
AND ap.mh_gspphot IS NOT NULL
AND fl.mass_flame_spec IS NOT NULL
AND fl.age_flame_spec IS NOT NULL
"""
print("Data is loading...")
job = Gaia.launch_job_async(query)
df = job.get_results().to_pandas()
print(f"Veri geldi: {df.shape}")

# Spectral Type
def teff_to_spectype(teff):
    if teff >= 30000: return 'O'
    elif teff >= 10000: return 'B'
    elif teff >= 7500: return 'A'
    elif teff >= 6600: return 'F0-F4'
    elif teff >= 6000: return 'F5-F9'
    elif teff >= 5200: return 'G'
    elif teff >= 3700: return 'K'
    else: return 'M'

df['spec_type'] = df['teff_gspphot'].apply(teff_to_spectype)

# Avoidance 
df['exclude_mass'] = df['mass_flame_spec'] > 1.5
df['exclude_age']  = df['age_flame_spec'] < 3.0
df['exclude_spec'] = df['spec_type'].isin(['O','B','A','F0-F4'])
df['exclude_feh']  = df['mh_gspphot'] < -0.4
df['exclude_mult'] = df['non_single_star'].fillna(0) >= 2
df['exclude_var']  = (
    (df['range_mag_g_fov'].fillna(0) > 0.01) |
    (df['phot_variable_flag'] == 'VARIABLE')
)
mask_M = df['spec_type'] == 'M'
df['exclude_act'] = (
    mask_M & (
        df['in_vari_rotation_modulation'].fillna(False) |
        df['in_vari_short_timescale'].fillna(False)
    )
)
df['excluded'] = df[['exclude_mass','exclude_age','exclude_spec',
                      'exclude_feh','exclude_mult','exclude_var','exclude_act']].any(axis=1)

# Reason codes
def reason_codes(row):
    codes = []
    if row['exclude_mass']: codes.append('R1')
    if row['exclude_age']:  codes.append('R2')
    if row['exclude_spec']: codes.append('R3')
    if row['exclude_feh']:  codes.append('R4')
    if row['exclude_mult']: codes.append('R5')
    if row['exclude_var']:  codes.append('R6')
    if row['exclude_act']:  codes.append('R7')
    return ';'.join(codes) if codes else 'RETAIN'

df['reason_codes'] = df.apply(reason_codes, axis=1)

# Results
total    = len(df)
excluded = df['excluded'].sum()
retained = total - excluded
print(f"\nToplam:   {total:,}")
print(f"Dışlanan: {excluded:,} ({excluded/total*100:.1f}%)")
print(f"Korunan:  {retained:,} ({retained/total*100:.1f}%)")

print("\nComparison:")
for col, code, label in [
    ('exclude_mass','R1','Mass > 1.5 M⊙'),
    ('exclude_age', 'R2','Age < 3 Gyr'),
    ('exclude_spec','R3','Spectral O/B/A/F0-F4'),
    ('exclude_feh', 'R4','[Fe/H] < -0.4'),
    ('exclude_mult','R5','Multiplicity ≥ 3'),
    ('exclude_var', 'R6','Photometric variability'),
    ('exclude_act', 'R7','Active M dwarf'),
]:
    n = df[col].sum()
    print(f"  {code} {label}: {n:,} ({n/total*100:.1f}%)")

# Save
csv_path = os.path.join(save_dir, 'gaia_seti_avoidance_robust.csv')
df.to_csv(csv_path, index=False)
print(f"\nCSV saved: {csv_path}")
