# Research Protocol

**Working title:** *Probabilistic Mortality and Population Projections for The
Gambia to 2074: a Bayesian Lee-Carter approach anchored on the 2024 census and
demographic surveillance data*

**Author:** Abdoulie Balisa, BSc Statistics, KNUST.
**Status:** Draft (Phase 0). Last updated 2026-06-23.

---

## Abstract (working)

The Gambia has no complete civil registration and vital statistics system, so its
mortality and population have to be inferred from a patchwork of censuses, sample
surveys, demographic surveillance, and model-based reconstructions. This project
produces an independent, reproducible projection of the country's population to
2074, with uncertainty stated throughout, taking the 2024 census as the base. I
model age-specific mortality with a Bayesian Lee-Carter specification fitted by
MCMC, carry the parameter and forecast uncertainty into probabilistic life tables,
and extend the single-population model with the Li-Lee coherent method (forecasting
The Gambia together with a West-African reference group so it borrows strength and
cannot drift to implausible long-run values). The population is projected with a
cohort-component (Leslie) method, reporting population by age and sex, the total,
child and old-age dependency ratios, and the school-age and working-age
populations, each with 80% and 95% credible intervals. Results are validated by an
out-of-sample backtest, compared with the Farafenni and Basse surveillance life
tables, cross-checked against the Global Burden of Disease estimates, and
benchmarked against the official UN World Population Prospects 2024 projection,
which itself uses a Bayesian hierarchical method. The central methodological
finding is a measure of how much The Gambia's demographic future depends on the
modelling choices an analyst has to make under data scarcity.

---

## 1. Background and motivation

Population projections are the quiet input behind almost every long-range public
decision: how many classrooms and teachers to fund, how big the future workforce
and pension bill will be, how many clinics and how much vaccine to buy. For The
Gambia, a young and fast-growing country of about 2.4 million (2024 census), these
numbers matter a great deal, yet the data to produce them is thin.

- **No complete civil registration.** Deaths are heavily under-registered, so
  there is no national time series of age-specific death rates. The Gambia is
  therefore absent from the Human Mortality Database, the usual starting point for
  classical mortality forecasting.
- **Sparse, mixed evidence.** Mortality is known mainly through periodic censuses,
  DHS surveys (child mortality directly, adult mortality through sibling-survival
  histories), two long-running surveillance sites, and internationally produced
  reconstructions (UN WPP, IHME/GBD) that themselves ingest all of the above.

This is the normal condition across most of sub-Saharan Africa, and it is exactly
why a careful, uncertainty-honest treatment is worth more than a single
point-estimate answer. The project does not pretend the data is better than it is.
It makes the uncertainty explicit and shows which conclusions hold up.

### Why The Gambia is unusually tractable

- **Farafenni surveillance site.** Continuous demographic surveillance since 1981
  (about 55,000 people by 2015), with cause of death by verbal autopsy. One of the
  longest-running such platforms in Africa.
- **Basse surveillance site.** Since 2007, following about 180,000 people, among
  the largest in the region, reporting falling adult and child mortality.
- **2024 census.** The country's first fully digital census, giving a fresh,
  high-quality base population by age and sex.

These give rare empirical mortality data to anchor and validate a model that would
otherwise lean entirely on reconstructions.

## 2. Problem statement and research questions

**RQ1 (Estimation).** What is the most defensible reconstruction of The Gambia's
age-specific mortality and its trend over time, given only sparse, mixed data, and
with honest uncertainty?

**RQ2 (Forecasting).** Under a Bayesian Lee-Carter model, what is the projected
mortality and life expectancy of The Gambia to 2074, with credible intervals, and
how does a coherent (Li-Lee, West-African group) version change the forecast?

**RQ3 (Projection).** Taking the 2024 census as the base, what is the projected
population by age and sex to 2074, and the resulting total, child and old-age
dependency ratios, school-age (3 to 17) and working-age (15 to 64) populations,
each with credible intervals?

**RQ4 (Robustness).** How sensitive are the headline conclusions to the modelling
choices forced by data scarcity (single versus coherent model, choice of mortality
input, jump-off adjustment), and how do the projections compare with the official
UN WPP 2024 projection?

## 3. Significance and intended contribution

1. **Substantive.** An independent, transparent, reproducible projection for The
   Gambia off the new 2024 census, usable by the Gambia Bureau of Statistics, the
   ministries (Education, Health, Finance), and development partners.
2. **Methodological.** A worked, honest template for probabilistic mortality
   projection in a no-registration setting: how to triangulate sources, where the
   uncertainty really comes from, and how much the method matters.
3. **Communication.** Findings written for both a technical audience (the report)
   and the public (a plain-language brief for national media), with uncertainty
   shown honestly rather than hidden behind a single line.

## 4. Data sources (summary; full detail in `01-data-sources.md`)

| Source | What it provides | Role | Access |
|---|---|---|---|
| UN WPP 2024 | Reconstructed single-age population and life tables 1950 to 2023; official projection to 2100 | Primary mortality input and benchmark | Open |
| 2024 census (GBoS) | Base population by age and sex | Projection base year | Open (preliminary out; detailed tables to come) |
| DHS 2013, 2019-20 | Under-5 mortality (direct); adult mortality (sibling survival); fertility | Empirical anchors and inputs | Open |
| Farafenni and Basse surveillance | Empirical age-specific mortality, multi-decade | Validation | Published life tables open; microdata via request (stretch) |
| IHME / GBD | Independent age-specific mortality with uncertainty | Cross-check | Open |
| World Bank / WHO | Headline indicators | Sanity checks | Open |

**Data honesty note.** The UN WPP and GBD series are model-based reconstructions,
not raw observations. Where the Bayesian Lee-Carter is fit to a reconstructed
series, this is stated, the reconstruction's own assumptions are summarised, and
the surveillance life tables are used as an outside check that the reconstruction
(and my forecast) is not detached from observed Gambian mortality.

## 5. Methodology

### 5.0 Pipeline overview

```
raw sources -> harmonised age-specific death rates m(x,t) and exposures
            -> Bayesian Lee-Carter (PyMC)
            -> Li-Lee coherent variant
            -> probabilistic future m(x,t), then probabilistic life tables and e0
            -> cohort-component projection (base = 2024 census)
            -> population by age and sex to 2074, plus dependency and age-structure
            -> validation, benchmarking, sensitivity, reporting
```

### 5.1 The Lee-Carter model

Lee and Carter (1992) write the log central death rate for age `x` in year `t` as:

```
log m(x,t) = a_x + b_x * k_t + error
```

- `a_x` is the average age pattern of mortality (the mean log rate by age).
- `k_t` is a single time index for the overall level of mortality each year.
- `b_x` is how sensitive each age is to changes in `k_t`.
- Identifiability constraints: the b values sum to 1, and the k values sum to 0.

Classically the parameters come from a singular value decomposition of the centred
log-rate matrix, after which `k_t` is re-estimated to match observed deaths and
then forecast as a random walk with drift. I implement and report this SVD baseline
for transparency and as a reference.

### 5.2 Bayesian Lee-Carter (core model)

The SVD route assumes constant-variance Gaussian errors on log rates, which is a
poor assumption where death counts are small and noisy, exactly the case here.
Following the Poisson log-bilinear tradition (Brouhns, Denuit and Vermunt 2002;
Czado, Delwarde and Denuit 2005), I instead model deaths as Poisson around the
Lee-Carter rate and estimate all parameters jointly by MCMC in PyMC:

```
D(x,t)  ~  Poisson( E(x,t) * exp( a_x + b_x * k_t ) )
k_t     ~  random walk with drift:  k_t = k_{t-1} + d + noise
```

with `a_x` fixed at the empirical mean log rate and `b` constrained to the simplex,
which resolves the identifiability problem cleanly. Where true death counts are not
available, exposures come from the WPP population. The reason for going Bayesian is
that it gives full posterior and posterior-predictive distributions, so parameter
uncertainty and forecast uncertainty both flow through to every downstream quantity
(life expectancy, population, dependency ratios), which is the point that
point-estimate methods miss.

**Convergence and diagnostics:** multiple chains, R-hat below 1.01, healthy bulk
and tail effective sample size, divergence checks, posterior-predictive checks
against held-out years, and trace and rank plots. Reported, not assumed.

### 5.3 Coherent forecasting (the Li-Lee extension)

A single small-population forecast can drift to implausible long-run mortality.
The Li-Lee (2005) coherent method forecasts a group of related populations with a
common rate-of-change factor plus population-specific deviations that are
mean-reverting, so members cannot diverge without bound. I treat The Gambia as one
member of a West-African reference group (Senegal, Guinea-Bissau, Mali, Guinea,
Sierra Leone, and others; selection documented), borrowing strength from the
group's better-estimated common trend. I compare the single-population and coherent
forecasts; the difference is a result in its own right (RQ4).

### 5.4 From mortality to life tables

For each posterior draw of future `m(x,t)` I build a full period life table,
producing a posterior distribution of life expectancy and the survivorship the
projection needs. The open age group is closed with a standard method, documented
and held fixed across models.

### 5.5 Population projection (cohort-component method)

Standard demographic accounting (the Leslie matrix), the same family the UN uses:

- **Base:** 2024 census population by single age and sex (graduated as needed; the
  jump-off mortality reconciled to the model).
- **Mortality:** posterior survivorship from 5.2 to 5.4.
- **Fertility:** age-specific fertility rates with the total fertility rate
  projected on a declining path consistent with WPP (in the spirit of Alkema et al.
  2011); sex ratio at birth from the data.
- **Migration:** net migration scenarios (a WPP-based central case plus low and
  high), with sensitivity reported given how big a role migration plays for The
  Gambia.

Uncertainty is propagated by running the projection over the posterior mortality
(and fertility and migration) draws, giving fan charts for total population, age
structure, and the derived indicators to 2074.

### 5.6 Benchmark against the UN's own method

The UN adopted Bayesian hierarchical projection for WPP from 2015 (Raftery et al.
2012; the bayesPop, bayesLife and bayesTFR packages, Sevcikova and Raftery 2016). I
do not just cite this. I reproduce the relevant logic conceptually, compare my
independent projection against the official WPP 2024 result for The Gambia, and
explain where they agree and differ. That is how the "same methods the UN uses"
claim is earned honestly.

## 6. Validation plan

1. **Out-of-sample backtest.** Fit the model to data up to 2010, forecast 2010 to
   2023, and score against the held-out series (interval coverage, error on log
   rates and life expectancy, calibration of the credible intervals).
2. **Empirical anchor.** Compare modelled death rates and life expectancy against
   the published Farafenni and Basse life tables for overlapping years and ages.
3. **Independent cross-check.** Compare against the GBD age-specific mortality, a
   different reconstruction with its own uncertainty.
4. **Benchmark.** Compare population and life expectancy against the official WPP
   2024 projection.
5. **Sanity checks.** Life expectancy, crude death rate and infant mortality
   against the World Bank and WHO.

## 7. Deliverables

- **D1.** Research report (`reports/`): full methods, results, uncertainty,
  limitations, policy implications.
- **D2.** Reproducible pipeline (`src/`): raw data through to every figure.
- **D3.** Figures: mortality surface, the Lee-Carter parameters, life-expectancy
  fan chart, population pyramids 2024 versus 2074, dependency-ratio fan charts,
  method comparison.
- **D4.** Plain-language policy brief for a general and national-media audience.
- **D5.** Annotated bibliography (`docs/02-literature.md`).

## 8. Project phases

| Phase | Goal | Key output |
|---|---|---|
| 0. Scoping | Frame, decide method, verify data | This protocol |
| 1. Data | Acquire, harmonise, document every source | Clean `data/processed`, `01-data-sources.md` |
| 2. Literature | Read and annotate core and Gambia-specific work | `02-literature.md` |
| 3. EDA | Mortality surface, trends, data-quality checks | Figures |
| 4. SVD baseline | Classical Lee-Carter, reproduced | Figures |
| 5. Bayesian LC | PyMC model, diagnostics, posterior life tables | Model code |
| 6. Coherent (Li-Lee) | West-African reference-group variant | Comparison results |
| 7. Projection | Cohort-component to 2074, plus indicators | Projection module |
| 8. Validation | Backtest and the comparisons above | Validation notes |
| 9. Reporting | Report, policy brief, figures | D1, D4 |
| 10. Dissemination | Repo polish, README, optional preprint/media | Public repo |

## 9. Ethics and data governance

All primary analysis uses public, aggregate data. Surveillance microdata, if
pursued, needs a formal data-sharing request and ethics clearance; no restricted
microdata is committed to the repository. Sources are cited and licences recorded
in `01-data-sources.md`.

## 10. Limitations (stated up front)

- Mortality inputs are partly reconstructions, not direct observations, so the
  forecast inherits their assumptions. This is mitigated by anchoring to the
  surveillance data and by sensitivity analysis.
- The surveillance sites are not nationally representative (rural, specific
  regions), so they are used to validate the level and shape, not as the national
  series.
- Migration is volatile and hard to forecast for The Gambia, so it is handled with
  explicit scenarios rather than a single number.
- The single-age detail in the 2024 census may be graduated; the method is
  documented and its effect tested.

## 11. References (working; bibliographic details verified in Phase 2)

- Lee, R. D., & Carter, L. R. (1992). Modeling and Forecasting U.S. Mortality.
  *JASA*, 87(419), 659 to 671.
- Brouhns, N., Denuit, M., & Vermunt, J. K. (2002). A Poisson log-bilinear
  regression approach to the construction of projected life tables.
  *Insurance: Mathematics and Economics*, 31(3), 373 to 393.
- Czado, C., Delwarde, A., & Denuit, M. (2005). Bayesian Poisson log-bilinear
  mortality projections. *Insurance: Mathematics and Economics*, 36(3), 260 to 284.
- Li, N., & Lee, R. D. (2005). Coherent mortality forecasts for a group of
  populations: an extension of the Lee-Carter method. *Demography*, 42(3), 575 to
  594. https://pmc.ncbi.nlm.nih.gov/articles/PMC1356525/
- Booth, H., & Tickle, L. (2008). Mortality modelling and forecasting: a review of
  methods. *Annals of Actuarial Science*, 3(1-2), 3 to 43.
- Alkema, L., Raftery, A. E., Gerland, P., et al. (2011). Probabilistic projections
  of the total fertility rate for all countries. *Demography*, 48(3), 815 to 839.
- Raftery, A. E., Li, N., Sevcikova, H., Gerland, P., & Heilig, G. K. (2012).
  Bayesian probabilistic population projections for all countries. *PNAS*, 109(35),
  13915 to 13921. https://www.pnas.org/doi/10.1073/pnas.1211452109
- Sevcikova, H., & Raftery, A. E. (2016). bayesPop: Probabilistic Population
  Projections. *Journal of Statistical Software*, 75(5).
  https://www.jstatsoft.org/article/view/v075i05
- Wilmoth, J., Zureick, S., Canudas-Romo, V., Inoue, M., & Sawyer, C. (2012). A
  flexible two-dimensional mortality model for use in indirect estimation.
  *Population Studies*, 66(1), 1 to 28.
- Jasseh, M., et al. (2015). Health and Demographic Surveillance System Profile:
  Farafenni, The Gambia. *Int. J. Epidemiology*, 44(3), 837 to 847.
  https://academic.oup.com/ije/article/44/3/837/631778
- Basse HDSS Cohort Profile. *Int. J. Epidemiology* (2025).
  https://doi.org/10.1093/ije/dyaf021
- United Nations, DESA, Population Division (2024). *World Population Prospects
  2024.* https://population.un.org/wpp/
- Gambia Bureau of Statistics (2024). *2024 Population and Housing Census,
  Preliminary Report.* https://www.gbosdata.org/

> Touchstones to verify and possibly add in Phase 2: Hyndman and Ullah (2007,
> functional data); Cairns-Blake-Dowd (an alternative to Lee-Carter); Girosi and
> King (2008, Bayesian mortality); Sharrow et al. (Bayesian mortality for
> data-sparse settings); the UN WPP 2024 methodology report.
