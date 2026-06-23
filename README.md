# The Gambia 2074 — Bayesian Mortality & Population Projection

**An independent, fully-reproducible probabilistic projection of The Gambia's
population to 2074, built on the 2024 Census and four decades of demographic
surveillance, using the Bayesian methodology adopted by the UN Population
Division.**

> Final-year research project · Statistics, KNUST · Abdoulie Balisa

---

## The one-line pitch

> *"I produced the first independent, uncertainty-quantified population
> projection for The Gambia off the new 2024 Census — fitting a Bayesian
> Lee–Carter mortality model with PyMC, extending it with the Li–Lee coherent
> method, and benchmarking against the UN's own Bayesian hierarchical approach.
> The result is a 50-year projection of population, age structure, and
> dependency ratios with full credible intervals — the numbers a country needs
> to plan schools, jobs, pensions and clinics."*

## Why this project

The Gambia has **no complete civil/vital registration system**, so it is absent
from the [Human Mortality Database](https://www.mortality.org/). Naively fitting
a mortality model to such a country is indefensible. This project's contribution
is to confront that limitation honestly and turn it into the research question:
**how do you build a credible, uncertainty-aware population projection for a
data-sparse country, and how much do the answers depend on the method?**

What makes The Gambia uniquely tractable:

- **Health & Demographic Surveillance Systems (HDSS):** Farafenni (since **1981**)
  and Basse (since **2007**) — rare multi-decade empirical mortality data with
  cause of death by verbal autopsy.
- **A brand-new 2024 Census** (2.4M people; first fully digital census) giving a
  fresh base population.
- **DHS surveys** (2013, 2019–20) for child and adult mortality.
- **UN WPP 2024** reconstructed single-age series and an official projection to
  benchmark against.

## Methodology in one diagram

```
Data triangulation                 Mortality model              Projection
(census, DHS, HDSS, WPP)      (Bayesian, with uncertainty)   (cohort-component)
        │                              │                            │
        ▼                              ▼                            ▼
 reconstructed kx,t  ──►  Bayesian Lee–Carter (PyMC, MCMC)  ──►  Leslie-matrix
 age-specific rates       + Li–Lee coherent (W. Africa group)    projection to 2074
        │                              │                            │
        │                  validated by: backtest (≤2010→2023),     ▼
        │                  HDSS life tables, GBD cross-check    population by age/sex
        └──────────────────────────────┴──────────────►  + 90% credible intervals
                                                          → dependency ratios,
                                                            school-age, working-age,
                                                            elderly populations
        benchmark ↑ : UN WPP official Bayesian projection (bayesPop/bayesLife)
```

## Repository structure

```
gambia-population-projection/
├── README.md                 ← you are here
├── docs/
│   ├── 00-research-protocol.md   ← full proposal: methodology, data, lit, plan
│   ├── 01-data-sources.md        ← exact sources, URLs, access notes, licences
│   └── 02-literature.md          ← annotated bibliography
├── data/
│   ├── raw/        ← downloaded source files (never edited)
│   ├── processed/  ← cleaned, analysis-ready tables
│   └── external/   ← reference data (e.g. WPP comparison series)
├── src/            ← reusable Python modules (data, models, projection, viz)
├── notebooks/      ← exploratory + reproducible analysis notebooks
├── references/     ← PDFs and the .bib bibliography
├── figures/        ← generated charts
└── reports/        ← the final research report / paper drafts
```

## Status

🟢 **Phase 0 — Scoping & protocol.** Project framed, methodology decided, data
landscape verified. See [`docs/00-research-protocol.md`](docs/00-research-protocol.md).

Next: data acquisition (Phase 1).

## Reproducibility

Every figure and number in the final report will be regenerable from raw public
data by running the pipeline end-to-end. Modeling uses a pinned conda
environment (PyMC builds cleanly on Windows via conda-forge).

## Licence & data governance

Code: MIT. Data: each source retains its own licence (see
[`docs/01-data-sources.md`](docs/01-data-sources.md)). No restricted microdata is
committed to this repository.
