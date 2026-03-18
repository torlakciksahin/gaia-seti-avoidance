https://doi.org/10.5281/zenodo.19054136

# Torlakcık Catalog
### A Parametric Avoidance Model for SETI Target Selection

---

## Overview

**Torlakcık Catalog** is an open-source stellar exclusion catalog and analysis pipeline for operational SETI target selection, applied to the full Gaia DR3 dataset.

This repository contains the analysis pipeline, figure generation scripts,
and sensitivity analysis used in the study.

The full exclusion catalog is archived on Zenodo (see Data section).

> **"Where Not to Look: A Parametric Avoidance Model for SETI Target Selection"**
> Torlakcık, Ş. (2026) — *manuscript in preparation*

The model applies seven astrophysical thresholds to **1,754,135 Gaia DR3 stars** and excludes **66.2%** of candidates while retaining high-priority targets such as stable K-type and chromospherically quiet M dwarfs.

---

## Results

| | N | Fraction |
|---|---|---|
| Total sample | 1,754,135 | 100% |
| Excluded | 1,161,623 | 66.2% |
| Retained | 592,512 | 33.8% |

---

## Exclusion Criteria

| Code | Criterion | Threshold | Note |
|------|-----------|-----------|------|
| R1 | Stellar mass | > 1.5 M⊙ | |
| R2 | Stellar age | < 3 Gyr | |
| R3 | Spectral type | O, B, A, F0–F4 | |
| R4 | Metallicity | [Fe/H] < −0.4 | |
| R5 | Multiplicity | `non_single_star` ≥ 1 | All NSS-flagged sources excluded; Gaia flag does not distinguish binaries from higher-order multiples |
| R6 | Photometric variability | range_mag_g > 0.01 mag or VARIABLE flag | |
| R7 | M dwarf activity | rotation modulation or short-timescale variability flag | |

---

## Repository Structure
```
Torlakcık Catalog/
├── analysis.py           # Main Gaia DR3 pipeline
├── figures.py            # Figure generation (all 6 figures)
├── sensitivity.py        # Sensitivity analysis
├── figures/              # Publication figures
│   ├── figure1_criteria_breakdown.png
│   ├── figure2_spectype_breakdown.png
│   ├── figure3_synthetic_vs_empirical.png
│   ├── figure4_decision_flow.png
│   ├── figure5_sensitivity_analysis.png
│   └── figure6_skymap.png
└── README.md
```

---

## Data


The **Torlakcık Catalog**, a Gaia DR3–based exclusion catalog containing
1,754,135 stars, is publicly available on Zenodo:
https://doi.org/10.5281/zenodo.19054136
The Torlakcık Catalog csv fie contains ~1.75 million rows and exceeds the row limits of common spreadsheet software such as Microsoft Excel (1,048,576 rows).
The dataset is intended for programmatic analysis using tools such as Python (pandas), R, or dedicated astronomical catalog software.
For large-table inspection, tools such as ExtremeCSV or TOPCAT are recommended.

# TAP
```sql
Gaia DR3 data accessed via ESA TAP service:
SELECT
gs.source_id, gs.ra, gs.dec,
ap.teff_gspphot, ap.mh_gspphot,
ap.logg_gspphot,
fl.mass_flame_spec, fl.age_flame_spec,
gs.non_single_star,
gs.phot_variable_flag,
vs.range_mag_g_fov,
vs.in_vari_rotation_modulation,
vs.in_vari_short_timescale
FROM gaiadr3.gaia_source AS gs
JOIN gaiadr3.astrophysical_parameters AS ap
ON gs.source_id = ap.source_id
JOIN gaiadr3.astrophysical_parameters_supp AS fl
ON gs.source_id = fl.source_id
LEFT JOIN gaiadr3.vari_summary AS vs
ON gs.source_id = vs.source_id
WHERE ap.teff_gspphot IS NOT NULL
AND ap.mh_gspphot IS NOT NULL
AND fl.mass_flame_spec IS NOT NULL
AND fl.age_flame_spec IS NOT NULL
```


---

## Requirements
```
pip install astroquery pandas numpy matplotlib
```

---

## Citation

*Citation will be updated upon publication.*

---

## License

MIT License — see LICENSE file.
