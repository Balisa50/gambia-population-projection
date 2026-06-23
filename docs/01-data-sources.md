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

## WPP 2024 acquisition — reproducible route (DONE)

The standard WPP `/data` API endpoint is **token-gated** (HTTP 401 without a
generated bearer token). Instead we use the UN Population Division's official
**`wpp2024` R data package** (`github.com/PPgp/wpp2024`), whose `data-raw/`
directory publishes the WPP 2024 series as plain tab-separated text — open, and
exactly reproducible because we **pin to a commit SHA**.

- Repo: `PPgp/wpp2024`  ·  Pinned SHA: `2da7768ae64fc74105d3f9e98f9a74d37b62f99a` (main @ 2025-06-24)
- Fetcher: [`src/fetch_wpp.py`](../src/fetch_wpp.py) — streams each file, keeps
  only The Gambia (country_code **270**), writes `data/raw/wpp2024/*_GMB.tsv` +
  `_manifest.json` (with per-file SHA-256). Re-run with `python src/fetch_wpp.py`.
- Loader/reshaper: [`src/wpp_data.py`](../src/wpp_data.py) (wide → tidy long).
- Got 31/31 files: m(x,t) single-age 1950–2100 (F/M/Both); population by single
  age 1949–2023 (F/M); e0 estimates + UN probabilistic projection (median +
  80/95% PI); TFR + %ASFR + projections; SRB; migration; UN projected population
  by age (medium); deaths/births/CDR/CBR.

### First integrity + triangulation findings (from `python src/wpp_data.py`)
- m(x,t) fully populated (101 ages × 151 years). Infant m(0): 0.164 (1950) →
  0.030 (2023). e0: 31.2 (1950) → 65.9 (2023).
- UN probabilistic e0: 2050 = 70.2 [62.4, 77.3]; 2074 = 73.2 [63.5, 83.4].
- **Population gap (a key finding):** WPP 2023 ≈ **2,728,905** vs 2024 Census
  ≈ **2,422,712** → WPP is ~**+12.6%** higher. WPP also exceeds every historical
  Gambian census (1973/1983/1993/2003/2013). WPP 2024 predates the new census;
  the projection must reconcile its **jump-off** to the census base (allowing for
  census under-enumeration). This is exactly the kind of result the thesis exists
  to surface.

## Download log

| Date | Source | File(s) | Saved to | Notes |
|---|---|---|---|---|
| 2026-06-23 | WPP 2024 (`PPgp/wpp2024` @ `2da7768`) | 31 GMB series | `data/raw/wpp2024/` | ✅ via `src/fetch_wpp.py` |
| _pending_ | GBoS 2024 PHC | age×sex tables (exact totals) | `data/raw/census2024/` | confirm 2,422,712 + single-age base |
| _pending_ | DHS | FR289, 2019–20 | `data/raw/dhs/` | child + adult mortality anchors |
| _pending_ | HDSS | Farafenni/Basse life tables | `data/raw/hdss/` | digitise from papers (validation) |
