# Probabilistic Mortality and Population Projections for The Gambia to 2074

### A Bayesian Lee–Carter Approach Anchored on the 2024 Census

**Abdoulie Balisa** · BSc Statistics, KNUST · Final-year research project
*Draft — reproducible from the public pipeline in this repository.*

---

## Abstract

The Gambia lacks a complete civil-registration and vital-statistics (CRVS)
system, so its mortality and population dynamics must be inferred from censuses,
sample surveys, demographic surveillance and model-based reconstructions. Using
the United Nations *World Population Prospects 2024* (WPP) reconstructed series,
the first fully-digital 2024 Census, and a West-African reference group, we
build an **independent, fully-reproducible, uncertainty-quantified projection of
The Gambia's population to 2074**. We model age-specific mortality with three
methods of increasing sophistication — classical Lee–Carter (SVD), a Bayesian
Poisson Lee–Carter (PyMC/MCMC), and the Li–Lee coherent extension — and feed the
results through a cohort-component projection that we validate against WPP's own
output (reproduced to within 0.9%). Two results stand out. First, **WPP's
population figures for The Gambia are ~13% higher than the 2024 Census**, so
re-basing on the census lowers the projected 2074 population from WPP's 5.35M to
**4.66M (95% CI 4.35–4.98M)**. Second, all three extrapolative mortality models
produce **far narrower forecast intervals than WPP's hierarchical Bayesian
method**, indicating that WPP's uncertainty is dominated by *structural* model
uncertainty that single-method extrapolation does not capture. The projection
shows The Gambia entering its **demographic-dividend window** — the total
dependency ratio falls from 77 (2023) to ~49 (2074) — even as the old-age
dependency ratio triples.

---

## 1. Introduction

Population projections underpin almost every long-range public decision in The
Gambia: classrooms and teachers, the future labour force and pension liability,
clinics and vaccines. Yet the data foundation is thin — there is no national
time series of age-specific death rates, and The Gambia is absent from the Human
Mortality Database. This project confronts that limitation directly: rather than
presenting a single deterministic "answer," it asks **how to build a credible,
uncertainty-honest projection under data scarcity, and how much the conclusions
depend on the method and on the data vintage.**

We make four contributions: (i) an independent projection re-based on the new
2024 Census; (ii) a like-for-like comparison of three mortality-forecasting
methods with honest uncertainty; (iii) a cohort-component engine validated
against the UN's own output; and (iv) a quantification of the gap between the
census and WPP and its consequences for national planning.

## 2. Data

All inputs are public and the acquisition is pinned for exact reproducibility
(see `docs/01-data-sources.md`).

- **UN WPP 2024** — single-age m(x,t) 1950–2100, single-age population 1949–2023,
  life expectancy estimates and the official probabilistic projection (median +
  80/95% bounds), TFR, %ASFR, sex ratio at birth and net migration. Obtained from
  the UN Population Division's official `PPgp/wpp2024` data package, pinned to
  commit `2da7768`, filtered to The Gambia (country code 270) and to a
  West-African reference group (Senegal, Guinea-Bissau, Guinea, Mali, Sierra
  Leone, Mauritania, Burkina Faso).
- **2024 Population and Housing Census (GBoS)** — total population **2,422,712**
  with broad age structure (40.8% under 15, 3.0% aged 65+); used as the
  projection base.
- **Reference data** — DHS (2013, 2019–20), Farafenni/Basse HDSS and IHME/GBD
  are identified for validation; HDSS life-table digitisation remains a data task.

**Data-honesty note.** WPP and GBD series are *model-based reconstructions*, not
raw observations. Where a model is fit to a reconstructed series this is stated,
and the analysis is framed around the *vintage* gap the new census reveals.

## 3. Methods

### 3.1 Life tables
Period life tables are built from m(x) from first principles (Coale–Demeny a(0),
open-interval closeout); the implementation reproduces WPP's published e0(2023) =
**65.86** exactly, validating the engine (`src/lifetable.py`).

### 3.2 Mortality models
**(a) Classical Lee–Carter (SVD).** `log m(x,t) = a_x + b_x k_t`; k_t forecast as
a random walk with drift, uncertainty by simulation (`src/leecarter.py`).

**(b) Bayesian Poisson Lee–Carter.** Deaths modelled as
`D(x,t) ~ Poisson(E(x,t)·exp(a_x + b_x k_t))`, with a_x fixed at the empirical
mean log-rate, `b ~ Dirichlet(1)`, and k_t a random walk with drift; estimated by
NUTS in PyMC (`src/bayes_leecarter.py`). Following Brouhns et al. (2002) and
Czado et al. (2005).

**(c) Li–Lee coherent.** The Gambia is forecast jointly with the reference group:
a common factor (B,K) shared by the group plus a country-specific deviation
forecast as a mean-reverting AR(1), so members cannot diverge (`src/coherent.py`;
Li & Lee 2005).

### 3.3 Population projection
A single-age, two-sex cohort-component (Leslie) projection (Preston, Heuveline &
Guillot 2001): survivorship from the life tables, fertility from TFR × %ASFR,
sex ratio at birth, and net migration. Validated by reproducing WPP's projection
from WPP inputs (`src/projection.py`), then run with our census base, our
mortality uncertainty and WPP fertility uncertainty, vectorised over 1,000
simulations to yield credible intervals (`src/independent_projection.py`).

## 4. Results

### 4.1 The census–WPP gap
WPP estimates The Gambia's 2023 population at **2,728,905**, ~**12.6% above** the
2024 Census total of 2,422,712. WPP also exceeds every historical Gambian census
(1973, 1983, 1993, 2003, 2013). WPP 2024 predates the new census, so its figures
require a downward revision (Figure: `projection_overview.png`).

### 4.2 Mortality
Lee–Carter fits The Gambia exceptionally well — the first singular component
explains **99.3%** of variance (Figure: `lc_parameters.png`). Life expectancy
rose from 31.2 (1950) to 65.9 (2023). **Out-of-sample backtest** (fit ≤2010,
predict 2011–2023): mean absolute e0 error **0.65 years**, with the actual within
the 95% interval in **100%** of years.

Forecast e0 in 2074 by method:

| Method | e0(2074) | 95% interval |
|---|---|---|
| Classical LC | 75.6 | [73.0, 77.8] |
| Bayesian LC | 75.0 | [73.8, 76.2] |
| Li–Lee coherent | 74.8 | [73.0, 76.6] |
| **UN WPP (Bayesian hierarchical)** | **73.2** | **[63.5, 83.4]** |

All three extrapolative models agree closely on the median (~75) but produce
intervals **3–5 years wide**, versus WPP's **~20 years**
(Figures: `bayes_lc_vs_wpp_e0.png`, `coherent_vs_single_vs_wpp.png`). This is a
substantive methodological finding: WPP's width reflects *structural* uncertainty
(about the long-run pace and ceiling of mortality decline) that single-method
extrapolation — even Bayesian or coherent — does not represent.

### 4.3 Population projection
The cohort-component engine reproduces WPP's medium projection to within
**0.3–0.9%** through 2074 when fed WPP inputs, validating it. Re-based on the 2024
Census and driven by our mortality uncertainty:

| Year | Population (median) | 95% credible interval | WPP medium |
|---|---|---|---|
| 2050 | 3.74M | 3.59–3.90M | 4.33M |
| 2074 | **4.66M** | **4.35–4.98M** | 5.35M |

The census re-basing lowers the 2074 figure by ~0.7M relative to WPP
(Figure: `independent_projection.png`).

### 4.4 Age structure and the demographic dividend
The total dependency ratio falls from **77 per 100 working-age (2023) to ~49
(2074)** — the demographic-dividend window. Its composition shifts sharply: the
child dependency ratio falls from 71.7 to **31.8 [27.3–36.5]**, while the old-age
dependency ratio **triples**, from 5.4 to **17.8 [16.5–19.2]**.

## 5. Discussion

The headline policy message is twofold. (1) The Gambia's population will continue
to grow strongly — likely reaching **~4.7 million by 2074** — but the credible
range (4.35–4.98M) and the **~0.7M downward revision implied by the new census**
matter for any plan calibrated to UN figures. (2) The country is entering a
window of falling dependency that can yield a demographic dividend *if* the
growing working-age population is productively employed; simultaneously, the
tripling of old-age dependency signals a future need for pension and health-system
capacity that barely exists today.

Methodologically, the convergence of three independent extrapolative models on a
narrow e0 band — well inside WPP's interval — suggests that conventional mortality
forecasting may understate long-run uncertainty for data-sparse countries, and
that the UN's wider intervals are doing important work by representing structural
uncertainty.

## 6. Limitations
- Mortality inputs are partly **reconstructions**; the forecasts inherit their
  assumptions. HDSS empirical life-table digitisation (Farafenni/Basse) and a GBD
  cross-check are the next validation steps.
- The census base uses the **broad** age structure (the only public detail at
  time of writing) with WPP's within-group shape; detailed single-age census
  tables will refine it.
- **Migration** is treated as a central scenario; given its volatility for The
  Gambia, scenario analysis is a priority extension.

## 7. Conclusion
On open, fully-reproducible data, this project produces the first independent,
uncertainty-quantified population projection for The Gambia off the 2024 Census,
with a validated engine and a transparent comparison of mortality methods. The
Gambia is projected to reach **4.66M (95% CI 4.35–4.98M) by 2074** — materially
below the prevailing UN figure once the new census is taken into account — while
entering a demographic-dividend window that is also the prelude to population
ageing.

## References
See `docs/02-literature.md` and `docs/00-research-protocol.md`. Core: Lee & Carter
(1992); Brouhns, Denuit & Vermunt (2002); Czado, Delwarde & Denuit (2005); Li &
Lee (2005); Raftery et al. (2012); Ševčíková & Raftery (2016); Preston, Heuveline
& Guillot (2001); UN WPP 2024; GBoS 2024 PHC.
