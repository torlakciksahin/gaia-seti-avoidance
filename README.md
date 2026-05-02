# Torlakcik Catalog

### *Where Not to Look:* A Parametric Avoidance Model for SETI Target Selection

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19956677.svg)](https://doi.org/10.5281/zenodo.19956677)



---

Imagine you have a radio telescope and a few thousand hours of observing time. The sky contains over a billion stars. Where do you point?

Most SETI surveys pick targets based on proximity — the closest stars, the brightest stars. But if you're looking for technological signatures around stars that could *actually* host habitable worlds, you need a different approach. You need to know **where not to look**.

The **Torlakcik Catalog** is a stellar exclusion catalog that does exactly that. It applies seven interpretable astrophysical criteria to 1.74 million Gaia DR3 stars and tells you which ones to skip — and *why* — so you can spend your telescope time on the stars that matter.

> Torlakcik, S. (2026) — *manuscript in preparation*

---

## What It Does

The model starts with **1,742,306 Gaia DR3 stars** that have reliable parallax measurements and complete astrophysical parameter coverage. It then applies a set of transparent, tweakable rules to filter out stars that are unlikely hosts for complex life:

| Rule | What it checks | Threshold | Why |
|------|---------------|-----------|-----|
| **R1** | Stellar mass | > 1.5 M<sub>&odot;</sub> | Massive stars burn out fast — not enough time for complex life to emerge |
| **R2** | Stellar age (upper bound) | < 3 Gyr | Young stars haven't had enough time either, and we use the *upper bound* to be fair to uncertain ages |
| **R3** | Spectral type | O, B, A, F0–F4 | These stars blast their planets with UV radiation |
| **R4** | Metallicity | [Fe/H] < −0.4 | Low metallicity means fewer building blocks for rocky planets |
| **R5** | Multiplicity | `non_single_star` &ge; 1 | Binary and multi-star systems destabilize planetary orbits |
| **R6** | Photometric variability | > 0.01 mag or `VARIABLE` flag | Wild brightness swings suggest a hostile environment for climate stability |
| **R7** | M dwarf activity | rotation/short-timescale flag | Active M dwarfs flare violently, which can strip planetary atmospheres |

After applying all seven rules, the model **excludes 55.4%** of the sample and **retains 777,597 stars** — stable K-type stars, quiet M dwarfs, and other promising candidates. The median distance to a retained star is just 382 pc, well within reach of current radio SETI facilities.

### A note on age (R2)

This is where the model does something a bit different. Rather than using the Gaia point-estimate age and applying a hard cut (which would unfairly penalize stars with large uncertainties), R2 uses the **age upper bound** from Gaia FLAME. This uncertainty-aware approach retains **355,086 additional stars** that a naive cut would throw away. If the upper bound on a star's age is at least 3 Gyr, we give it the benefit of the doubt.

---

## Key Numbers

| | N | Fraction |
|---|---|---|
| Total sample | 1,742,306 | 100% |
| Excluded | 964,709 | 55.4% |
| **Retained** | **777,597** | **44.6%** |
| Median distance (retained) | 382 pc | — |

### Breakthrough Listen Cross-Match

How does our habitability-driven approach compare to the proximity-driven Breakthrough Listen survey?

| Sample | Matched | Excluded | Retained |
|--------|---------|----------|---------|
| BL Primary (Isaacson et al. 2017) | 405 |229 (56.5%) | 176 (43.5%) |
| MeerKAT 1M (Czech et al. 2021) | 98,716 | 42,531 (43.1%) | 56,185 (56.9%) |

Over half of BL's primary targets would be excluded by our criteria — mostly due to low metallicity. This isn't a criticism of BL; it's a reflection of two fundamentally different survey philosophies. BL casts a wide net based on proximity. We ask a different question: *which nearby stars are actually worth listening to?*

---

## Getting the Data

The full catalog is on Zenodo:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19956677.svg)](https://doi.org/10.5281/zenodo.19956677)



Heads up — the CSV has ~1.74 million rows, which is more than Excel can handle (its limit is 1,048,576). Use Python, R, or [TOPCAT](https://www.star.bris.ac.uk/~mbt/topcat/) to work with it:

```python
import pandas as pd

df = pd.read_csv("torlakcik_catalog.csv")
retained = df[df["excluded"] == False]
print(f"Retained: {len(retained):,} stars")
```



---

## Quickstart

```bash
git clone https://github.com/torlakciksahin/gaia-seti-avoidance.git
cd gaia-seti-avoidance
pip install -r requirements.txt
```

Then download the catalog from [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19956677.svg)](https://doi.org/10.5281/zenodo.19956677)
, place it in `data/`, and start exploring:

```python
import pandas as pd

df = pd.read_csv("data/torlakcik_catalog.csv")
retained = df[df["excluded"] == False]
print(f"Retained: {len(retained):,} stars")
```

---

## Gaia DR3 TAP Query

Want to build the catalog from scratch? Here's the ADQL query that pulls the required columns from the Gaia archive:

```sql
SELECT
gs.source_id, gs.ra, gs.dec,
gs.parallax, gs.parallax_over_error,
gs.phot_g_mean_mag, gs.bp_rp,
ap.teff_gspphot, ap.mh_gspphot,
ap.logg_gspphot,
fl.mass_flame_spec,
fl.age_flame_spec,
fl.age_flame_spec_upper,
gs.non_single_star,
gs.phot_variable_flag,
vs.range_mag_g_fov,
vs.in_vari_rotation_modulation,
vs.in_vari_short_timescale,
gs.ruwe,
gs.phot_g_mean_flux,
gs.phot_g_mean_flux_error
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
AND gs.parallax_over_error > 5

```

---

## Requirements

```
pip install astroquery pandas numpy matplotlib scipy astropy
```

---

## Related

- [stellar-avoidance](https://github.com/torlakciksahin/stellar-avoidance) -- The catalog-agnostic, reusable version of this pipeline. Same R1-R7 logic, but configurable via YAML for any stellar survey. If you want to apply these (or your own) criteria to a different catalog, use that one.

## Contributing

Found a bug? Have an idea for a new criterion or a different threshold? Issues and pull requests are welcome.

If you use this catalog in your work, please cite the Zenodo record and the associated manuscript.

---

## Citation

```bibtex
@misc{torlakcik2026,
  author       = {Torlakc{\i}k, {\c{S}}ahin},
  title        = {Where Not to Look: A Parametric Avoidance Model for SETI Target Selection},
  year         = {2026},
  howpublished = {Manuscript in preparation},
  doi          = {10.5281/zenodo.19644176},
  url          = {[https://doi.org/10.5281/zenodo.19644176](https://zenodo.org/records/19956677)}
}
```

---

## Acknowledgments

This work made use of data from the European Space Agency (ESA) mission [Gaia](https://www.cosmos.esa.int/gaia), processed by the Gaia Data Processing and Analysis Consortium ([DPAC](https://www.cosmos.esa.int/web/gaia/dpac/consortium)).

---

## License

MIT License — see [LICENSE](LICENSE) file.
