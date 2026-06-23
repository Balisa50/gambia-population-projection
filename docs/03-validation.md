# Validation (Phase 8)

How we know the pipeline is trustworthy. All numbers are reproduced by
`python src/validation.py` and `python src/projection.py`.

## 1. Internal consistency — life-table engine
Our from-first-principles period life table reproduces **WPP's published
e0(2023) = 65.86 exactly** (|diff| < 0.001). This confirms the life-table code
(Coale–Demeny a(0), open-interval closeout) and the rate→e0 path used everywhere
downstream.

## 2. Out-of-sample backtest — mortality forecast
Fit the model on data **≤ 2010** and predict **2011–2023**:
- mean absolute error in e0 = **0.65 years**;
- the actual e0 falls inside the **95% interval in 100%** of the 13 held-out years.
Figure: `figures/validation_backtest.png`. The model slightly under-predicts the
post-2010 improvement (real e0 rose a touch faster than the 1950–2010 trend),
which is the expected, conservative direction.

## 3. Projection engine — reproduces WPP
Fed WPP's own inputs (mortality, TFR, %ASFR, SRB, net migration), the
cohort-component engine reproduces **WPP's published medium projection to within
0.3–0.9%** across 2030–2074. An engine that reconstructs the benchmark can be
trusted to carry our independent (census-based, Lee–Carter/Bayesian) inputs.

| Year | our engine | WPP | diff |
|---|---|---|---|
| 2030 | 3,149,576 | 3,160,605 | −0.3% |
| 2040 | 3,738,105 | 3,762,462 | −0.6% |
| 2050 | 4,300,578 | 4,328,827 | −0.7% |
| 2074 | 5,298,317 | 5,346,053 | −0.9% |

## 4. Cross-model agreement
Three independent mortality methods (classical LC, Bayesian LC, Li–Lee coherent)
agree on e0(2074) within ~0.8 years (74.8–75.6), a convergent-validity signal.

## Remaining validation tasks (data-dependent)
- **HDSS empirical check** — digitise published Farafenni/Basse life-table values
  and compare modelled age-specific mortality / e0 over overlapping years.
- **GBD cross-check** — compare against IHME age-specific mortality (an
  independent reconstruction with its own uncertainty).
Both require manual data extraction and are flagged in `docs/01-data-sources.md`.
