# Data Sources — exact locations, access, licences

> Rule: `data/raw/` holds files **exactly as downloaded** (never edited). All
> cleaning happens in code and lands in `data/processed/`. Each download is
> logged below with date and URL so the pipeline is reproducible.

## 1. UN World Population Prospects 2024 (PRIMARY mortality input + benchmark)

- **Portal:** https://population.un.org/wpp/
- **Data Portal (interactive + open API):** https://population.un.org/dataportal/
- **Bulk CSV downloads:** https://population.un.org/wpp/downloads
- **Methodology report:** https://population.un.org/wpp/assets/Files/WPP2024_Methodology-Report_Final.pdf
- **Data sources doc (per-country inputs, incl. The Gambia):**
  https://population.un.org/wpp/assets/Files/WPP2024_Data_Sources.pdf
- **What to pull for The Gambia (ISO code GMB / 270):**
  - Life tables (abridged + single age), both sexes & by sex, 1950–2023 (estimates) and 2024–2100 (medium + probabilistic projection).
  - Age-specific mortality rates `m(x,t)` / `nMx`, and `nqx`.
  - Population by single age & sex (for exposures and as a WPP cross-check on the census base).
  - Total fertility rate, ASFR, sex ratio at birth, net migration.
- **Why both estimates and projection:** estimates 1950–2023 = our fitting
  series; official projection = our benchmark (RQ4).
- **Licence:** UN — free use with attribution.

## 2. 2024 Population and Housing Census — GBoS (PROJECTION BASE)

- **GBoS data portal:** https://www.gbosdata.org/ → Downloads → 2024 PHC reports
  - https://www.gbosdata.org/downloads/142-2024-population-and-housing-census-reports
  - https://www.gbosdata.org/downloads/157-2024-population-and-housing-census
- **Preliminary report (UNFPA mirror):**
  https://gambia.unfpa.org/en/publications/preliminary-report-2024-census-gambia
- **Confirmed headline figures (preliminary):** total ≈ **2.4 million**;
  **40.8%** under 15; **3.0%** aged 65+; female ≈1.24M / male ≈1.19M; +~300k
  since 2013.
- **Needed:** population by **single year of age × sex** (base for the Leslie
  projection). If only broad age groups are public at download time, use WPP
  2024 single-age structure reconciled to the census total, and document the
  reconciliation; upgrade to detailed census tables as GBoS releases them.
- **Licence:** GBoS official statistics — public.

## 3. Demographic and Health Surveys (mortality + fertility inputs/anchors)

- **DHS Program (registration required, free):** https://dhsprogram.com/
- **The Gambia DHS 2013 — final report (FR289):**
  https://dhsprogram.com/pubs/pdf/FR289/FR289.pdf
- **The Gambia DHS 2019–20 — report (SR268 summary; locate full FR):**
  https://dhsprogram.com/pubs/pdf/SR268/SR268.pdf
- **Use:** under-5 mortality (direct, reliable); adult mortality from
  sibling-survival histories (noisier — treat with care); ASFR/TFR.
- **Licence:** DHS data-use agreement; report PDFs are open.

## 4. Health & Demographic Surveillance Systems (EMPIRICAL VALIDATION)

- **Farafenni HDSS profile (IJE 2015):**
  https://academic.oup.com/ije/article/44/3/837/631778
- **Basse HDSS cohort profile (IJE 2025):**
  https://doi.org/10.1093/ije/dyaf021
- **GHDx record (Farafenni):**
  https://ghdx.healthdata.org/record/gambia-farafenni-health-and-demographic-surveillance-system
- **INDEPTH consolidated HDSS data, 29 SSA sites 1990–2018 (GHDx):**
  https://ghdx.healthdata.org/record/indepth-network-consolidated-hdss-data-29-sub-saharan-african-sites-1990-2018
- **iSHARE repository:** https://www.indepth-ishare.org/
- **Open now:** published age-specific mortality / life-table figures inside the
  profile papers → digitise for validation.
- **Stretch (microdata):** formal data-sharing request via iSHARE / MRCG@LSHTM;
  needs a research statement + likely ethics clearance. **Not on the critical
  path** — open data is sufficient for the core project.

## 5. Independent cross-checks

- **IHME / GBD (GHDx):** https://ghdx.healthdata.org/ — age-specific mortality
  with uncertainty intervals (independent reconstruction).
- **World Bank:** https://data.worldbank.org/country/gambia-the — e0, CDR, IMR,
  population (sanity checks).
- **WHO Data:** https://data.who.int/countries/270
- **Our World in Data (Gambia):** https://ourworldindata.org/country/gambia

## Download log

| Date | Source | File(s) | Saved to | Notes |
|---|---|---|---|---|
| _pending_ | WPP 2024 | life tables, pop, fertility (GMB) | `data/raw/wpp2024/` | Phase 1 |
| _pending_ | GBoS 2024 PHC | age×sex tables | `data/raw/census2024/` | Phase 1 |
| _pending_ | DHS | FR289, 2019–20 | `data/raw/dhs/` | Phase 1 |
