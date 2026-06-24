# Probabilistic Mortality and Population Projections for The Gambia to 2074

### A Bayesian Lee-Carter approach anchored on the 2024 census

**Abdoulie Balisa. BSc Statistics, KNUST. Final-year research project.**
*Reproducible from the public pipeline in this repository.*

---

## Abstract

The Gambia has no complete system for registering deaths, so its mortality and
population dynamics have to be pieced together from censuses, sample surveys,
demographic surveillance, and model-based reconstructions. This project builds an
independent, reproducible projection of the country's population to 2074, with
uncertainty stated throughout. It uses the United Nations World Population
Prospects 2024 reconstructed series, the first fully digital census from 2024, and
a group of West-African neighbours. Age-specific mortality is modelled three ways:
a classical Lee-Carter fit by singular value decomposition, a Bayesian Poisson
Lee-Carter sampled with MCMC, and the Li-Lee coherent extension. The forecasts run
through a cohort-component projection that reproduces the UN's own output to within
1%. Two results stand out. First, the UN's figures for The Gambia sit about 13%
above the 2024 census, so re-basing on the census lowers the 2074 population from
the UN's 5.35 million to 4.66 million (range 4.35 to 4.98 million). Second, all
three of my mortality models give much tighter forecast ranges than the UN's
hierarchical method, which suggests the UN's width is dominated by structural model
uncertainty that simple extrapolation does not capture. The projection shows the
country entering its demographic-dividend window, with the total dependency ratio
falling from 77 to about 49, even as old-age dependency triples.

---

## 1. Introduction

Population projections sit quietly behind almost every long-range public decision
in The Gambia: how many classrooms and teachers to fund, how big the future
workforce and pension bill will be, how many clinics to build and where. Yet the
data to make them is thin. There is no national record of who dies at what age, so
the country does not appear in the Human Mortality Database, the usual starting
point for this work.

This project takes that limitation as the question rather than working around it.
How do you build a credible projection, with honest uncertainty, for a country
with no death registration, and how much do the conclusions depend on the method
and on which census you trust? It makes four contributions: an independent
projection re-based on the new census; a like-for-like comparison of three
mortality methods with their uncertainty; a projection engine checked against the
UN's own output; and a clear measure of the gap between the census and the UN, with
its consequences for planning.

## 2. Data

All inputs are public, and the way they were collected is pinned so the work
reproduces exactly (see `docs/01-data-sources.md`).

- **UN WPP 2024.** Single-age death rates from 1950 to 2100, single-age population
  from 1949 to 2023, life expectancy estimates, and the UN's official probabilistic
  projection. Pulled from the UN team's own data package, pinned to one commit,
  filtered to The Gambia and to a West-African reference group (Senegal,
  Guinea-Bissau, Guinea, Mali, Sierra Leone, Mauritania, Burkina Faso).
- **2024 census (GBoS).** Total population about 2,422,712, with a broad age
  breakdown (40.8% under 15, 3.0% aged 65 and over). Used as the starting point.
- **For validation.** The DHS surveys (2013, 2019-20) and the Farafenni and Basse
  surveillance sites. Digitising the surveillance life tables is a remaining task.

A note on honesty: the UN and GBD series are reconstructions, not raw
observations. Where a model is fit to a reconstructed series, the report says so,
and the analysis is framed around the gap the new census reveals.

## 3. Methods

### 3.1 Life tables
Period life tables are built from death rates from first principles, using the
Coale-Demeny rule for infant person-years and a standard close-out for the open
age group. The implementation reproduces the UN's published life expectancy for
2023 (65.86) exactly, which is a basic check that the engine is sound.

### 3.2 Mortality models
**Classical Lee-Carter.** The model writes log death rates as `a_x + b_x k_t`,
where `k_t` is a single index of how mortality changes over time. It is fit by SVD
and `k_t` is carried forward as a random walk with drift.

**Bayesian Poisson Lee-Carter.** Deaths are modelled as Poisson around the
Lee-Carter rate, with `a_x` fixed at the average log rate (which is how Lee and
Carter define it), `b_x` constrained to the simplex, and `k_t` a random walk with
drift. It is sampled with MCMC in PyMC. Fixing `a_x` removes an identifiability
problem that otherwise stops the sampler from converging.

**Li-Lee coherent.** The Gambia is forecast together with the reference group. A
common factor shared by the group carries the long-run trend, and a
country-specific deviation is forecast as a mean-reverting process so the country
cannot drift away from its neighbours.

### 3.3 Population projection
A single-age, two-sex cohort-component (Leslie) projection. Survival comes from the
life tables, fertility from the total fertility rate times the age pattern of
childbearing, plus the sex ratio at birth and net migration. The engine is first
checked by reproducing the UN's projection from the UN's inputs, then run on the
census base with my own mortality uncertainty and the UN's fertility uncertainty,
over 1,000 simulations, to give credible ranges.

## 4. Results

### 4.1 The census gap
The UN puts The Gambia at 2,728,905 in 2023, about 13% above the 2024 census total
of 2,422,712. The UN also runs above every historical Gambian census. Its 2024
numbers were finished before the census came out, so they need revising downward.

### 4.2 Mortality
Lee-Carter fits the country's history well: the first component explains 99.3% of
the variation. Life expectancy rose from 31.2 years in 1950 to 65.9 in 2023.
Trained only on data up to 2010 and asked to predict 2011 to 2023, the model was
off by 0.65 years on average, and the truth fell inside its 95% range in every one
of those years.

Forecast life expectancy in 2074, by method:

| Method | e0 in 2074 | 95% range |
|---|---|---|
| Classical Lee-Carter | 75.6 | 73.0 to 77.8 |
| Bayesian Lee-Carter | 75.0 | 73.8 to 76.2 |
| Li-Lee coherent | 74.8 | 73.0 to 76.6 |
| UN WPP (hierarchical Bayesian) | 73.2 | 63.5 to 83.4 |

All three of my models agree closely on the middle value, around 75, but produce
ranges three to five years wide. The UN's range is about 20 years wide. That
difference is a real finding: the UN's width reflects uncertainty about the
long-run pace and ceiling of mortality decline that single-method extrapolation,
even a Bayesian one, does not represent.

### 4.3 Population projection
Fed the UN's own inputs, the engine reproduces the UN's projection to within 0.3 to
0.9% out to 2074, which validates it. Re-based on the 2024 census and driven by my
mortality uncertainty:

| Year | Population (median) | 95% range | UN medium |
|---|---|---|---|
| 2050 | 3.74M | 3.59 to 3.90M | 4.33M |
| 2074 | 4.66M | 4.35 to 4.98M | 5.35M |

Re-basing on the census lowers the 2074 figure by about 0.7 million compared with
the UN.

### 4.4 Age structure and the dividend
The total dependency ratio falls from 77 per 100 working-age adults in 2023 to
about 49 by 2074. Its make-up changes sharply: child dependency drops from 71.7 to
about 31.8, while old-age dependency triples, from 5.4 to about 17.8.

## 5. Discussion

The policy message has two parts. The population will keep growing strongly and
will likely reach about 4.7 million by 2074, but the range (4.35 to 4.98 million)
and the 0.7 million downward revision the census implies both matter for any plan
calibrated to UN figures. And the country is entering a window of falling
dependency that can deliver a real dividend if the growing working-age population
finds productive work. At the same time, old-age dependency tripling points to a
need for pension and health capacity that barely exists today.

On the methods, the fact that three independent models land on a narrow range for
life expectancy, well inside the UN's, suggests that conventional mortality
forecasting may understate long-run uncertainty for data-sparse countries, and
that the UN's wider ranges are doing useful work.

## 6. Limitations

The mortality inputs are partly reconstructions, so the forecasts inherit their
assumptions. The next validation steps are to compare against the Farafenni and
Basse surveillance life tables and against the Global Burden of Disease estimates.
The census base uses the broad age structure that is public so far, with the UN's
within-group shape; detailed single-age census tables will refine it. Migration is
treated as a central scenario, and given how volatile it is for The Gambia, scenario
analysis is a priority.

## 7. Conclusion

On open, reproducible data, and using the methods the UN itself trusts, this
project produces the first independent, uncertainty-aware population projection for
The Gambia built on the 2024 census. It puts the country at 4.66 million by 2074
(range 4.35 to 4.98 million), meaningfully below the figure in circulation once the
census is taken into account, while heading into a demographic-dividend window that
is also the run-up to an ageing society.

## References

See `docs/02-literature.md` and `docs/00-research-protocol.md`. Core references:
Lee and Carter (1992); Brouhns, Denuit and Vermunt (2002); Czado, Delwarde and
Denuit (2005); Li and Lee (2005); Raftery et al. (2012); Sevcikova and Raftery
(2016); Preston, Heuveline and Guillot (2001); UN WPP 2024; GBoS 2024 census.
