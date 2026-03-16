# gaia-seti-avoidance
Parametric avoidance model for SETI target selection — applied to Gaia DR3 (N=1,754,135)


# gaia-seti-avoidance

A parametric avoidance model for operational SETI target selection, applied to Gaia DR3.

## Overview

This repository contains the analysis code, exclusion catalog, and figures for the paper:

> **"Where Not to Look: A Parametric Avoidance Model for SETI Target Selection"**
> S. Torlakcık (2026)

The model applies seven astrophysical thresholds to 1,754,135 Gaia DR3 stars and excludes 65.3% of candidates while retaining high-priority targets such as stable K-type and chromospherically quiet M dwarfs.

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

## Results

| | N | Fraction |
|---|---|---|
| Total sample | 1,754,135 | 100% |
| Excluded | 1,145,218 | 65.3% |
| Retained | 608,917 | 34.7% |

## Repository Structure
```
gaia-seti-avoidance/
├── analysis.py                        # Main pipeline
├── figures.py                         # Figure generation
├── sensitivity.py                     # Sensitivity analysis
├── gaia_seti_avoidance_robust.csv     # Full exclusion catalog
├── figures/
│   ├── figure1_criteria_breakdown.png
│   ├── figure2_spectype_breakdown.png
│   ├── figure3_synthetic_vs_empirical.png
│   ├── figure4_decision_flow.png
│   ├── figure5_sensitivity_analysis.png
│   └── figure6_skymap.png
└── README.md
```

## Data
 Gaia DR3 (N=1,754,135)
Full TAP query is included in `analysis.py`.

## Requirements
```
astroquery
pandas
numpy
matplotlib
```

Install with:
```
pip install astroquery pandas numpy matplotlib
```

## License

MIT License — see LICENSE file.

## Citation
Arriving in the near future
