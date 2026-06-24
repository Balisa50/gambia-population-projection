# The Gambia 2074: a population projection built on the new census

An independent forecast of The Gambia's population out to 2074. It is built on
the 2024 census, it uses the same mortality methods the UN relies on, and it
comes with honest uncertainty ranges. Everything here runs from public data, so
anyone can check it.

> Final-year research project. BSc Statistics, KNUST. Abdoulie Balisa.
> Read the plain-language write-up: https://balisa50.github.io/research/gambia-2074

---

## What this is

The Gambia has no working death-registration system. Most deaths are never
recorded, so the country does not appear in the
[Human Mortality Database](https://www.mortality.org/), and its future
population is mostly guesswork. The only figures around come from the UN, and
those were finished before the country ran its first digital census in 2024. So
I rebuilt the forecast from scratch on the census, and checked it against the
UN's own numbers.

The interesting part is the constraint. You cannot just fit a mortality model to
a country with no death data and call it a day. So the real question became this:
how do you build a credible, uncertainty-aware projection when the data is thin,
and how much does the answer depend on which method you pick?

## What I found

1. **The UN is about 13% too high.** It puts The Gambia at 2.73M in 2023; the
   2024 census counted about 2.42M. Re-basing on the census brings the 2074
   projection down from 5.35M to **4.66M** (range 4.35 to 4.98M).
2. **Mortality has improved a lot, and the model captures it.** Life expectancy
   rose from 31 years in 1950 to 66 in 2023. Trained only on data up to 2010, the
   model predicted 2011 to 2023 with about 0.65 years of error.
3. **All three methods give tighter ranges than the UN.** That gap is a finding
   in itself. The UN's wide band reflects long-run structural uncertainty that
   simple trend models tend to miss.
4. **The country is entering its demographic dividend.** The dependency ratio
   falls from 77 to 49, even as the number of older people roughly triples.

## The data

- UN WPP 2024 single-age series, plus the UN's official projection to compare against.
- The 2024 census (about 2.42M, the first fully digital one) for the base population.
- DHS surveys (2013 and 2019-20) for child and adult mortality.
- The Farafenni (running since 1981) and Basse (since 2007) surveillance sites,
  rare long-run mortality data, kept for validation.

## How it works

```
Data (census, DHS, surveillance, UN WPP)
   -> age-specific death rates m(x,t)
   -> mortality model, three ways:
        classical Lee-Carter (SVD baseline)
        Bayesian Lee-Carter (PyMC, with uncertainty)
        Li-Lee coherent (forecast alongside West-African neighbours)
   -> probabilistic life tables and life expectancy
   -> cohort-component (Leslie) projection, re-based on the 2024 census
   -> population by age and sex to 2074, with credible ranges
   -> dependency ratios, working-age and school-age populations
Checked against: a held-out backtest, and the UN's own projection (matched to 1%).
```

## Reproduce it

```bash
python src/fetch_wpp.py              # UN data for The Gambia
python src/fetch_refgroup.py         # West-African reference group
python src/eda.py                    # mortality surface
python src/leecarter.py              # classical model + backtest
python src/bayes_leecarter.py        # Bayesian model (PyMC)
python src/coherent.py               # Li-Lee coherent model
python src/projection.py             # engine, validated against the UN
python src/independent_projection.py # the headline: census-based projection
python src/validation.py             # validation summary
```

Figures land in `figures/`, processed tables in `data/processed/`.

## Repository

```
docs/    research protocol, data sources, literature, validation notes
src/     Python modules: data, life tables, the three models, the projection
data/    raw downloads (Gambia extracts committed), processed tables
figures/ generated charts
reports/ the full technical report and a plain-language policy brief
```

## Running it

System Python 3.12 or 3.13 with numpy, pandas, scipy, matplotlib, plus
`pip install pymc arviz` for the Bayesian model. PyMC samples fine without a C
compiler (just slower). The data fetch is pinned to a specific commit of the UN
team's data package, so re-running gives the exact same inputs.

## What is still open

I have not yet compared the model against the published life tables from the
Farafenni and Basse sites, or against the Global Burden of Disease estimates.
Both need data I have to extract by hand, so they are listed as honest next steps
rather than skipped quietly. See [`docs/03-validation.md`](docs/03-validation.md).

## Licence and data

Code is MIT. Each data source keeps its own licence (see
[`docs/01-data-sources.md`](docs/01-data-sources.md)). No restricted microdata is
in this repository.
