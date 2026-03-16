# gaia-seti-avoidance

A parametric avoidance model for operational SETI target selection, applied to Gaia DR3.

## Overview

This repository contains the analysis code, exclusion catalog, and figures associated with:

> **"Where Not to Look: A Parametric Avoidance Model for SETI Target Selection"**
> S. Torlakcık (2026) — *manuscript in preparation*

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

## Requirements
```
pip install astroquery pandas numpy matplotlib
```

## License

MIT License — see LICENSE file.
