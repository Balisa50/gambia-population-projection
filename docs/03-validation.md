# Validation (Phase 8)

How I know the pipeline is trustworthy. Every number here is reproduced by
`python src/validation.py` and `python src/projection.py`.

## 1. Internal consistency (life-table engine)
My from-scratch period life table reproduces WPP's published life expectancy for
2023 (65.86) exactly, to within 0.001. That confirms the life-table code (the
Coale-Demeny rule for a(0), the open-interval close-out) and the path from death
rates to life expectancy that everything downstream relies on.

## 2. Out-of-sample backtest (mortality forecast)
Fit the model only on data up to 2010 and predict 2011 to 2023:
- mean absolute error in life expectancy: 0.65 years;
- the actual value fell inside the 95% range in all 13 held-out years.

See `figures/validation_backtest.png`. The model slightly under-predicts the
post-2010 improvement, because real life expectancy rose a touch faster than the
1950 to 2010 trend. That is the cautious direction to be off in.

## 3. Projection engine (reproduces WPP)
Fed WPP's own inputs (mortality, TFR, age pattern of fertility, sex ratio at
birth, net migration), the cohort-component engine reproduces WPP's published
medium projection to within 0.3 to 0.9% across 2030 to 2074. An engine that can
rebuild the benchmark can be trusted to carry my own census-based inputs.

| Year | my engine | WPP | difference |
|---|---|---|---|
| 2030 | 3,149,576 | 3,160,605 | -0.3% |
| 2040 | 3,738,105 | 3,762,462 | -0.6% |
| 2050 | 4,300,578 | 4,328,827 | -0.7% |
| 2074 | 5,298,317 | 5,346,053 | -0.9% |

## 4. Cross-model agreement
The three mortality methods (classical, Bayesian, and Li-Lee coherent) agree on
life expectancy in 2074 to within about 0.8 years (74.8 to 75.6). Independent
methods landing in the same place is a good sign.

## Remaining checks (data-dependent)
- **HDSS empirical check.** Digitise the published Farafenni and Basse life-table
  values and compare them with the modelled death rates and life expectancy over
  the overlapping years.
- **GBD cross-check.** Compare against the IHME age-specific mortality estimates,
  an independent reconstruction with its own uncertainty.

Both need data I have to extract by hand. They are flagged in
`docs/01-data-sources.md`.
