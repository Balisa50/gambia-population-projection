# Annotated Bibliography

> Each entry has the citation, a one-line summary, and how it is used here.
> Mirrored in `references/references.bib`.

## A. Mortality modelling, foundations

- **Lee & Carter (1992), JASA.** The original single-index log-bilinear model.
  Used for the SVD baseline (section 5.1).
- **Brouhns, Denuit & Vermunt (2002), IME.** Poisson log-bilinear estimation of
  Lee-Carter, which handles small and noisy death counts. Used for the
  observation model (section 5.2).
- **Czado, Delwarde & Denuit (2005), IME.** Bayesian (MCMC) Poisson Lee-Carter.
  Used for the core Bayesian specification.
- **Booth & Tickle (2008), Annals of Actuarial Science.** A review of mortality
  forecasting methods. Used to place the method choice in context.

## B. Coherent / multi-population

- **Li & Lee (2005), Demography.** Coherent forecasts for a group of populations.
  Used for the West-African reference-group variant (section 5.3).

## C. Bayesian projection (the UN approach)

- **Alkema et al. (2011), Demography.** Bayesian probabilistic TFR projection.
  Used for the fertility trajectory logic (section 5.5).
- **Raftery et al. (2012), PNAS.** Bayesian probabilistic population projections
  for all countries, the method the UN adopted. Used to frame the benchmark
  (section 5.6).
- **Sevcikova & Raftery (2016), JStatSoft.** The `bayesPop` package. A reference
  implementation and benchmark.

## D. Data-sparse / indirect estimation

- **Wilmoth et al. (2012), Population Studies.** Log-quadratic model life tables
  for indirect estimation. Used as a structural prior and a plausibility check.
- _To add:_ Sharrow et al.; UN MORTPAK and model life table systems; the Brass
  relational logit; Murray et al. model life tables.

## E. The Gambia and HDSS-specific

- **Jasseh et al. (2015), IJE.** Farafenni HDSS profile.
- **Basse HDSS Cohort Profile (2025), IJE.** dyaf021.
- _To add:_ Gambia DHS 2013 (FR289) and 2019-20; the GBoS 2024 census
  methodology; any published Gambian life tables.

## F. Alternatives considered (for the method-choice discussion)

- The Cairns-Blake-Dowd (CBD) family; the Hyndman-Ullah functional-data approach;
  Girosi & King (2008) on Bayesian mortality. Noted for completeness. Lee-Carter
  was chosen because it is transparent, aligned with what the UN uses, well
  understood, and suited to the data here.
