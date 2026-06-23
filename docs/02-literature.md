# Annotated Bibliography

> Built in Phase 2. Each entry: full citation (verified), one-line summary, and
> how it's used here. Mirrored in `references/references.bib`.

## A. Mortality modelling — foundations

- **Lee & Carter (1992), JASA.** The original single-index log-bilinear model.
  *Use:* SVD baseline (§5.1).
- **Brouhns, Denuit & Vermunt (2002), IME.** Poisson log-bilinear estimation of
  Lee–Carter — handles small/noisy death counts. *Use:* observation model (§5.2).
- **Czado, Delwarde & Denuit (2005), IME.** Bayesian (MCMC) Poisson Lee–Carter.
  *Use:* the core Bayesian specification.
- **Booth & Tickle (2008), Annals of Actuarial Science.** Review of mortality
  forecasting methods. *Use:* situating method choice.

## B. Coherent / multi-population

- **Li & Lee (2005), Demography.** Coherent forecasts for a group of
  populations. *Use:* the West-African reference-group variant (§5.3).

## C. Bayesian projection (the UN approach)

- **Alkema et al. (2011), Demography.** Bayesian probabilistic TFR projection.
  *Use:* fertility trajectory logic (§5.5).
- **Raftery et al. (2012), PNAS.** Bayesian probabilistic population projections
  for all countries — adopted by the UN. *Use:* benchmark framing (§5.6).
- **Ševčíková & Raftery (2016), JStatSoft.** `bayesPop`. *Use:* reference
  implementation / benchmark.

## D. Data-sparse / indirect estimation

- **Wilmoth et al. (2012), Population Studies.** Log-quadratic model life tables
  for indirect estimation. *Use:* structural prior / plausibility checks.
- _To add:_ Sharrow et al.; UN MORTPAK / model life table systems; Brass
  relational logit; Murray et al. model life tables.

## E. The Gambia / HDSS-specific

- **Jasseh et al. (2015), IJE.** Farafenni HDSS profile.
- **Basse HDSS Cohort Profile (2025), IJE.** dyaf021.
- _To add:_ Gambia DHS 2013 (FR289) & 2019–20; GBoS 2024 census methodology;
  any published Gambian life tables.

## F. Alternatives considered (for the "method choice" discussion)

- Cairns–Blake–Dowd (CBD) family; Hyndman–Ullah functional data approach;
  Girosi & King (2008) Bayesian mortality. Noted for completeness; Lee–Carter
  chosen as the transparent, UN-aligned, well-understood baseline appropriate
  to the data.
