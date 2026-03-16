# Torlakcık Catalog
### A Parametric Avoidance Model for SETI Target Selection

---

## Overview

**Torlakcık Catalog** is an open-source stellar exclusion catalog and analysis pipeline for operational SETI target selection, applied to the full Gaia DR3 dataset.

This repository contains the analysis code, exclusion catalog, and publication figures associated with:

> **"Where Not to Look: A Parametric Avoidance Model for SETI Target Selection"**
> Torlakcık, Ş. (2026) — *manuscript in preparation*

The model applies seven astrophysical thresholds to **1,754,135 Gaia DR3 stars** and excludes **65.3%** of candidates while retaining high-priority targets such as stable K-type and chromospherically quiet M dwarfs.

---

## Results

| | N | Fraction |
|---|---|---|
| Total sample | 1,754,135 | 100% |
| Excluded | 1,145,218 | 65.3% |
| Retained | 608,917 | 34.7% |

---

## Exclusion Criteria

| Code | Criterion | Threshold |
|------|-----------|-----------|
| R1 | Stellar mass | > 1.5 M⊙ |
| R2 | Stellar age | < 3 Gyr |
| R3 | Spectral type | O, B, A, F0–F4 |
| R4 | Metallicity | [Fe/H] < −0.4 |
| R5 | Multiplicity | ≥ 3 components |
| R6 | Photometric variability | range_mag_g > 0.01 mag or VARIABLE flag |
| R7 | M dwarf activity | rotation modulation or short-timescale variability |

---

## Repository Structure
```
Torlakcık Catalog/
├── analysis.py          # Main Gaia DR3 pipeline
├── figures.py           # Figure generation
├── sensitivity.py       # Sensitivity analysis
├── figures/             # Publication figures
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

Full exclusion catalog (N = 1,754,135) will be made available on Zenodo upon publication.

Gaia DR3 data accessed via ESA TAP service:
`https://gea.esac.esa.int/tap-server/tap`

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
