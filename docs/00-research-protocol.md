# Research Protocol

**Working title:** *Probabilistic Mortality and Population Projections for The
Gambia to 2074: A Bayesian Lee–Carter Approach Anchored on the 2024 Census and
Demographic Surveillance Data*

**Author:** Abdoulie Balisa · BSc Statistics, KNUST
**Status:** Draft v0.1 (Phase 0) · Last updated: 2026-06-23

---

## Abstract (working)

The Gambia lacks a complete civil registration and vital statistics (CRVS)
system, so its mortality and population dynamics must be inferred from a
patchwork of censuses, sample surveys, demographic surveillance, and
model-based reconstructions. This project produces an **independent,
fully-reproducible, uncertainty-quantified projection of The Gambia's
population to 2074**, taking the 2024 Population and Housing Census as the base.
We model age-specific mortality with a **Bayesian Lee–Carter** specification
estimated by MCMC, propagate parameter and forecast uncertainty into
**probabilistic life tables**, and extend the single-population model with the
**Li–Lee coherent method** (forecasting The Gambia jointly with a West-African
reference group to borrow strength and prevent implausible long-run divergence).
We project the population forward with a **cohort-component (Leslie-matrix)**
method and report population by age and sex, the total/old-age/child dependency
ratios, and the school-age and working-age populations, each with 80%/95%
credible intervals. Results are **validated** by out-of-sample backtesting,
comparison with Farafenni/Basse HDSS empirical life tables, cross-checked
against the Global Burden of Disease estimates, and **benchmarked** against the
official UN World Population Prospects 2024 projection, which itself uses a
Bayesian hierarchical methodology. The central methodological finding is a
quantification of how much The Gambia's demographic future depends on the
modelling choices an analyst must make under data scarcity.

---

## 1. Background and motivation

Population projections are the quiet input behind almost every long-range public
decision: how many classrooms and teachers to fund, how large the future labour
force and pension liability will be, how many clinics and how much vaccine to
buy. For The Gambia — a young, fast-growing country of 2.4 million (2024
Census) — these numbers matter acutely, yet the data foundation to produce them
is thin:

- **No complete CRVS.** Deaths are substantially under-registered, so there is
  no national time series of age-specific death rates. The Gambia is therefore
  absent from the Human Mortality Database, the standard input for classical
  mortality forecasting.
- **Sparse, heterogeneous evidence.** Mortality is known mainly through (i)
  periodic censuses, (ii) DHS surveys (child mortality directly; adult mortality
  via sibling-survival histories), (iii) two long-running HDSS sites, and (iv)
  internationally produced reconstructions (UN WPP, IHME/GBD) that themselves
  ingest the above.

This is the normal condition of most of sub-Saharan Africa, and it is precisely
why a careful, uncertainty-honest treatment is valuable rather than a
point-estimate "answer". The project does not pretend the data are better than
they are; it makes the uncertainty explicit and shows which conclusions are
robust to it.

### Why The Gambia is unusually tractable

- **Farafenni HDSS** — continuous demographic surveillance since **1981**
  (~55,000 people by 2015), with cause of death assigned by verbal autopsy
  (InterVA-4 / ICD-10). One of the longest-running such platforms in Africa.
- **Basse HDSS** — since **2007**, following ~180,000 people, among the largest
  in the region, reporting falling adult and child mortality.
- **2024 Census** — the country's first fully digital census, providing a fresh,
  high-quality base population by age and sex.

These give rare *empirical* mortality data to anchor and validate a model that
must otherwise lean on reconstructions.

## 2. Problem statement and research questions

**RQ1 (Estimation).** What is the most defensible reconstruction of The Gambia's
age-specific mortality schedule and its trend over time, given only sparse and
heterogeneous data, and with honest uncertainty?

**RQ2 (Forecasting).** Under a Bayesian Lee–Carter model, what is the projected
mortality schedule and life expectancy of The Gambia to 2074, with credible
intervals — and how does a coherent (Li–Lee, West-African reference group)
specification change the forecast?

**RQ3 (Projection).** Taking the 2024 Census as the base, what is the projected
population by age and sex to 2074, and the resulting total, child and old-age
dependency ratios, school-age (3–17) and working-age (15–64) populations, each
with credible intervals?

**RQ4 (Robustness).** How sensitive are the headline conclusions to the
modelling choices forced by data scarcity (single vs. coherent model; choice of
mortality input series; jump-off adjustment), and how do the projections compare
with the official UN WPP 2024 projection?

## 3. Significance and intended contribution

1. **Substantive.** An independent, transparent, reproducible projection for The
   Gambia off the new 2024 Census — usable by the Gambia Bureau of Statistics
   (GBoS), ministries (Education, Health, Finance), and development partners.
2. **Methodological.** A worked, honest template for probabilistic mortality
   projection in a *no-CRVS* setting: how to triangulate sources, where the
   uncertainty really comes from, and how much method choice matters.
3. **Communication.** Findings framed for both a technical audience (the report)
   and the public (a plain-language brief suitable for national media), with
   uncertainty communicated faithfully rather than hidden behind a single line.

## 4. Data sources (summary — full detail in `01-data-sources.md`)

| Source | What it provides | Role | Access |
|---|---|---|---|
| **UN WPP 2024** | Reconstructed single-age population & life tables 1950–2023; official projection to 2100 | Primary mortality input + benchmark | Open (CSV bulk + API) |
| **2024 PHC (GBoS)** | Base population by age & sex; recent change | Projection base year | Open (preliminary out; detailed tables as released) |
| **DHS 2013, 2019–20** | Under-5 mortality (direct); adult mortality (sibling survival); fertility | Empirical anchors & inputs | Open (DHS Program) |
| **Farafenni & Basse HDSS** | Empirical age-specific mortality, multi-decade | Validation / reality-check | Published life tables open; microdata via iSHARE request (stretch) |
| **IHME / GBD** | Independent age-specific mortality estimates w/ uncertainty | Cross-check | Open (GHDx) |
| **World Bank / WHO** | Headline indicators (e0, CDR, IMR) | Sanity checks | Open |

**Data honesty note.** UN WPP and GBD series are *model-based reconstructions*,
not raw observations. Where the Bayesian Lee–Carter is fit to a reconstructed
series, this is stated explicitly, the reconstruction's own assumptions are
summarised, and the HDSS empirical life tables are used as an external check
that the reconstruction (and our forecast) is not detached from observed
Gambian mortality.

## 5. Methodology

### 5.0 Pipeline overview

```
raw sources → harmonised age-specific death rates m(x,t) and exposures
            → Bayesian Lee–Carter (PyMC)  ─┐
            → Li–Lee coherent variant      ├→ probabilistic future m(x,t)
                                            │   → probabilistic life tables, e0
            → cohort-component projection ◄─┘   (base = 2024 Census)
            → population by age/sex to 2074 + dependency/age-structure indicators
            → validation, benchmarking, sensitivity, reporting
```

### 5.1 The Lee–Carter model

The Lee–Carter (1992) model decomposes the log central death rate for age `x`
in year `t`:

```
log m(x,t) = a_x + b_x · k_t + ε(x,t)
```

- `a_x` — average age-specific mortality shape (the mean log-rate by age).
- `k_t` — a single time index capturing the overall mortality *level* each year.
- `b_x` — the sensitivity of each age to changes in `k_t`.
- Identifiability constraints: `Σ_x b_x = 1` and `Σ_t k_t = 0`.

Classically `a_x, b_x, k_t` are estimated by **Singular Value Decomposition** of
the centred log-rate matrix (the first left/right singular vectors give `b_x`
and `k_t`), then `k_t` is re-estimated to match observed deaths, and forecast as
a random walk with drift. We implement and report this **SVD baseline** for
transparency and as a reference.

### 5.2 Bayesian Lee–Carter (core model)

The SVD/least-squares route assumes homoskedastic Gaussian errors on log-rates —
a poor assumption where death counts are small and noisy (exactly our case).
Following the Poisson log-bilinear tradition (Brouhns, Denuit & Vermunt 2002;
Czado, Delwarde & Denuit 2005) we instead model **deaths as Poisson** around the
Lee–Carter rate, and estimate all parameters jointly by **MCMC in PyMC**:

```
D(x,t)  ~  Poisson( E(x,t) · exp( a_x + b_x · k_t ) )
k_t     ~  random walk with drift:  k_t = k_{t-1} + d + η_t,  η_t ~ N(0, σ_k²)
priors  on a_x, b_x, d, σ_k ; identifiability imposed by post-hoc rescaling
         or a sum-to-one / sum-to-zero re-parameterisation
```

Where exposures `E(x,t)` for true counts are unavailable (reconstructed-rate
case), we either (a) use WPP exposures, or (b) fit a Gaussian/Student-t
observation model on log-rates with age-varying variance — both variants are run
and compared. **Why Bayesian:** it yields full posterior (and posterior
predictive) distributions, so *parameter* uncertainty and *forecast* uncertainty
are propagated coherently into every downstream quantity (life expectancy,
population, dependency ratios) — the core selling point over point forecasts.

**Convergence & diagnostics:** multiple chains; R-hat < 1.01; sufficient bulk
and tail ESS; divergence checks; posterior predictive checks against held-out
years; trace and rank plots. Reported, not assumed.

### 5.3 Coherent forecasting — the Li–Lee extension

A single small population forecast can drift to implausible long-run mortality.
The **Li–Lee (2005)** coherent method forecasts a group of related populations
with a **common** rate-of-change factor plus population-specific deviations that
are mean-reverting, so members do not diverge without bound. We treat The Gambia
as one member of a **West-African reference group** (e.g. Senegal, Guinea-Bissau,
Mali, Guinea, Sierra Leone — selection criteria documented), borrowing strength
from the group's better-estimated common trend. We **compare** the single-population
and coherent forecasts; the difference is a substantive result (RQ4).

### 5.4 From mortality to life tables

For each posterior draw of future `m(x,t)` we build a full **period life table**
(`m → q → l → L → T → e`), producing a posterior distribution of life
expectancy `e0` and the survivorship needed by the projection. Old-age closure
uses a standard method (e.g. Kannisto/Coale–Kisker extrapolation) documented and
held fixed across models.

### 5.5 Population projection — cohort-component method

Standard demographic accounting (the **Leslie matrix / cohort-component**
method), the same family used by the UN:

- **Base:** 2024 Census population by single age and sex (smoothed/graduated as
  needed; jump-off mortality reconciled to the model).
- **Mortality:** posterior survivorship from §5.2–5.4.
- **Fertility:** age-specific fertility rates with the **TFR projected**
  (declining trajectory consistent with WPP / a Bayesian double-logistic in the
  spirit of Alkema et al. 2011); sex ratio at birth from data.
- **Migration:** net migration scenarios (WPP-based central case + low/high),
  with sensitivity reported given migration's outsized role for The Gambia.

Uncertainty is propagated by running the projection over posterior mortality
(and fertility/migration) draws, yielding **fan charts** for total population,
age structure, and derived indicators to 2074.

### 5.6 Benchmark against the UN's own method

The UN adopted **Bayesian hierarchical projection** for WPP from 2015 (Raftery
et al. 2012; `bayesPop`/`bayesLife`/`bayesTFR`, Ševčíková & Raftery 2016). We do
**not** merely cite this — we (a) reproduce the relevant logic conceptually, (b)
compare our independent projection against the *official* WPP 2024 result for
The Gambia, and (c) explain agreements and divergences. This is how the "same
methodology the UN uses" claim is earned honestly.

## 6. Validation plan

1. **Out-of-sample backtest.** Fit the model to data up to 2010; forecast
   2010–2023; score against the held-out series (interval coverage, MAE/RMSE on
   log-rates and e0, calibration of credible intervals).
2. **Empirical anchor.** Compare modelled age-specific mortality and e0 against
   **Farafenni/Basse HDSS published life tables** for overlapping years/ages.
3. **Independent cross-check.** Compare against **GBD** age-specific mortality
   (a different reconstruction with its own uncertainty).
4. **Benchmark.** Compare population & e0 trajectories against **WPP 2024**
   official projection.
5. **Sanity indicators.** e0, crude death rate, IMR vs. World Bank/WHO.

## 7. Deliverables

- **D1** Research report / paper draft (`reports/`) — full methodology, results,
  uncertainty, limitations, policy implications.
- **D2** Reproducible pipeline (`src/` + `notebooks/`) — raw data → every figure.
- **D3** Figures: mortality surface, `a_x/b_x/k_t` posteriors, e0 fan chart,
  population pyramids 2024 vs 2074, dependency-ratio fan charts, method
  comparison.
- **D4** Plain-language policy brief (2–3 pp) for a general/national-media
  audience.
- **D5** Annotated bibliography (`docs/02-literature.md`).

## 8. Project phases

| Phase | Goal | Key output |
|---|---|---|
| **0. Scoping** ✅ | Frame, decide method, verify data | This protocol |
| **1. Data** | Acquire, harmonise, document every source | Clean `data/processed`, `01-data-sources.md` |
| **2. Lit review** | Read & annotate core + Gambia-specific work | `02-literature.md`, `.bib` |
| **3. EDA** | Mortality surface, trends, data-quality checks | Notebook + figures |
| **4. SVD baseline** | Classical Lee–Carter, reproduced | Notebook + figures |
| **5. Bayesian LC** | PyMC model, diagnostics, posterior life tables | Model code + diagnostics |
| **6. Coherent (Li–Lee)** | West-African reference group variant | Comparison results |
| **7. Projection** | Cohort-component to 2074 + indicators | Projection module + fan charts |
| **8. Validation** | Backtest, HDSS/GBD/WPP comparisons | Validation notebook |
| **9. Reporting** | Report + policy brief + figures | D1, D4 |
| **10. Dissemination** | Repo polish, README, optional preprint/media | Public repo |

## 9. Ethics and data governance

All primary analysis uses **public, aggregate** data. HDSS **microdata**, if
pursued, requires a formal data-sharing request (iSHARE / MRCG@LSHTM) and ethics
clearance; **no restricted microdata is committed to the repository**. Sources
are cited; licences recorded in `01-data-sources.md`.

## 10. Limitations (stated up front)

- Mortality inputs are partly **reconstructions**, not direct observations; the
  forecast inherits their assumptions (mitigated by HDSS anchoring + sensitivity).
- HDSS sites are **not nationally representative** (rural; specific regions) —
  used for validation of *level/shape*, not as the national series.
- **Migration** is volatile and hard to forecast for The Gambia — handled by
  explicit scenarios rather than a single number.
- Single-age detail in the 2024 Census may be **graduated/smoothed**; method
  documented and its effect tested.

## 11. References (working — bibliographic details verified in Phase 2)

- Lee, R. D., & Carter, L. R. (1992). Modeling and Forecasting U.S. Mortality.
  *JASA*, 87(419), 659–671.
- Brouhns, N., Denuit, M., & Vermunt, J. K. (2002). A Poisson log-bilinear
  regression approach to the construction of projected life tables.
  *Insurance: Mathematics and Economics*, 31(3), 373–393.
- Czado, C., Delwarde, A., & Denuit, M. (2005). Bayesian Poisson log-bilinear
  mortality projections. *Insurance: Mathematics and Economics*, 36(3), 260–284.
- Li, N., & Lee, R. D. (2005). Coherent mortality forecasts for a group of
  populations: An extension of the Lee–Carter method. *Demography*, 42(3),
  575–594. https://pmc.ncbi.nlm.nih.gov/articles/PMC1356525/
- Booth, H., & Tickle, L. (2008). Mortality modelling and forecasting: A review
  of methods. *Annals of Actuarial Science*, 3(1–2), 3–43.
- Alkema, L., Raftery, A. E., Gerland, P., et al. (2011). Probabilistic
  projections of the total fertility rate for all countries. *Demography*,
  48(3), 815–839.
- Raftery, A. E., Li, N., Ševčíková, H., Gerland, P., & Heilig, G. K. (2012).
  Bayesian probabilistic population projections for all countries. *PNAS*,
  109(35), 13915–13921. https://www.pnas.org/doi/10.1073/pnas.1211452109
- Ševčíková, H., & Raftery, A. E. (2016). bayesPop: Probabilistic Population
  Projections. *Journal of Statistical Software*, 75(5).
  https://www.jstatsoft.org/article/view/v075i05
- Wilmoth, J., Zureick, S., Canudas-Romo, V., Inoue, M., & Sawyer, C. (2012).
  A flexible two-dimensional mortality model for use in indirect estimation.
  *Population Studies*, 66(1), 1–28.
- Jasseh, M., et al. (2015). Health & Demographic Surveillance System Profile:
  Farafenni HDSS, The Gambia. *Int. J. Epidemiology*, 44(3), 837–847.
  https://academic.oup.com/ije/article/44/3/837/631778
- Basse HDSS Cohort Profile. *Int. J. Epidemiology* (2025).
  https://doi.org/10.1093/ije/dyaf021
- United Nations, DESA, Population Division (2024). *World Population Prospects
  2024.* https://population.un.org/wpp/
- Gambia Bureau of Statistics (2024). *2024 Population and Housing Census,
  Preliminary Report.* https://www.gbosdata.org/

> Methodological touchstones to verify and possibly add in Phase 2: Hyndman &
> Ullah (2007, functional data); Cairns–Blake–Dowd (alternative to Lee–Carter);
> Girosi & King (2008, Bayesian mortality); Sharrow et al. (Bayesian mortality
> for data-sparse settings); UN WPP 2024 Methodology Report.
