# Data Sources: locations, access, licences

> Rule: `data/raw/` holds files exactly as downloaded (never edited). All
> cleaning happens in code and lands in `data/processed/`. Each download is
> logged below with date and URL so the pipeline reproduces.

## 1. UN World Population Prospects 2024 (primary mortality input and benchmark)

- **Portal:** https://population.un.org/wpp/
- **Data portal (interactive, API):** https://population.un.org/dataportal/
- **Bulk CSV downloads:** https://population.un.org/wpp/downloads
- **Methodology report:** https://population.un.org/wpp/assets/Files/WPP2024_Methodology-Report_Final.pdf
- **Data sources doc (per-country inputs, including The Gambia):**
  https://population.un.org/wpp/assets/Files/WPP2024_Data_Sources.pdf
- **What to pull for The Gambia (ISO code GMB, M49 270):**
  - Life tables (abridged and single age), both sexes and by sex, 1950 to 2023 (estimates) and 2024 to 2100 (medium plus probabilistic projection).
  - Age-specific death rates `m(x,t)` and `nqx`.
  - Population by single age and sex (for exposures and as a cross-check on the census base).
  - Total fertility rate, age pattern of fertility, sex ratio at birth, net migration.
- **Why both:** the 1950 to 2023 estimates are the fitting series; the official projection is the benchmark (RQ4).
- **Licence:** UN, free use with attribution.

## 2. 2024 Population and Housing Census, GBoS (the projection base)

- **GBoS data portal:** https://www.gbosdata.org/ (Downloads, then 2024 PHC reports)
  - https://www.gbosdata.org/downloads/142-2024-population-and-housing-census-reports
  - https://www.gbosdata.org/downloads/157-2024-population-and-housing-census
- **Preliminary report (UNFPA mirror):**
  https://gambia.unfpa.org/en/publications/preliminary-report-2024-census-gambia
- **Headline figures (preliminary):** total about 2.4 million; 40.8% under 15;
  3.0% aged 65 and over; female about 1.24M, male about 1.19M.
- **Needed:** population by single year of age and sex (the base for the Leslie
  projection). If only broad age groups are public, use the WPP 2024 single-age
  shape reconciled to the census total, and document that reconciliation, then
  upgrade to detailed census tables when GBoS releases them.
- **Licence:** GBoS official statistics, public.

## 3. Demographic and Health Surveys (mortality and fertility anchors)

- **DHS Program (free, registration required):** https://dhsprogram.com/
- **The Gambia DHS 2013, final report (FR289):**
  https://dhsprogram.com/pubs/pdf/FR289/FR289.pdf
- **The Gambia DHS 2019-20 (SR268 summary; locate the full FR):**
  https://dhsprogram.com/pubs/pdf/SR268/SR268.pdf
- **Use:** under-5 mortality (direct and reliable); adult mortality from
  sibling-survival histories (noisier, treat with care); fertility.
- **Licence:** DHS data-use agreement; the report PDFs are open.

## 4. Health and Demographic Surveillance Systems (empirical validation)

- **Farafenni HDSS profile (IJE 2015):**
  https://academic.oup.com/ije/article/44/3/837/631778
- **Basse HDSS cohort profile (IJE 2025):**
  https://doi.org/10.1093/ije/dyaf021
- **GHDx record (Farafenni):**
  https://ghdx.healthdata.org/record/gambia-farafenni-health-and-demographic-surveillance-system
- **INDEPTH consolidated HDSS data, 29 SSA sites 1990 to 2018 (GHDx):**
  https://ghdx.healthdata.org/record/indepth-network-consolidated-hdss-data-29-sub-saharan-african-sites-1990-2018
- **iSHARE repository:** https://www.indepth-ishare.org/
- **Open now:** the published death rates and life-table figures inside the
  profile papers, which can be digitised for validation.
- **Stretch goal (microdata):** a formal data-sharing request via iSHARE or
  MRCG@LSHTM, which needs a research statement and probably ethics clearance.
  This is not on the critical path; open data is enough for the core project.

## 5. Independent cross-checks

- **IHME / GBD (GHDx):** https://ghdx.healthdata.org/ for age-specific mortality
  with uncertainty intervals (an independent reconstruction).
- **World Bank:** https://data.worldbank.org/country/gambia-the for life
  expectancy, crude death rate, infant mortality, population.
- **WHO data:** https://data.who.int/countries/270
- **Our World in Data (Gambia):** https://ourworldindata.org/country/gambia

## How the WPP 2024 data was actually pulled (done)

The standard WPP `/data` API endpoint is token-gated (it returns HTTP 401 without
a generated bearer token). So instead I use the UN Population Division's own
`wpp2024` R data package (`github.com/PPgp/wpp2024`), whose `data-raw/` directory
publishes the WPP 2024 series as plain tab-separated text. It is open, and it is
exactly reproducible because the fetch is pinned to a specific commit.

- Repo: `PPgp/wpp2024`. Pinned SHA: `2da7768ae64fc74105d3f9e98f9a74d37b62f99a` (main, 2025-06-24).
- Fetcher: [`src/fetch_wpp.py`](../src/fetch_wpp.py) streams each file, keeps only
  The Gambia (country_code 270), and writes `data/raw/wpp2024/*_GMB.tsv` plus a
  `_manifest.json` with a SHA-256 per file. Re-run with `python src/fetch_wpp.py`.
- Loader: [`src/wpp_data.py`](../src/wpp_data.py) reshapes wide to tidy long.
- Got 31 of 31 files: single-age death rates 1950 to 2100 (female, male, both);
  population by single age 1949 to 2023; life expectancy estimates plus the UN
  probabilistic projection (median and 80/95% intervals); fertility and its
  projection; sex ratio at birth; migration; the UN projected population by age
  (medium); and counts.

### First integrity and triangulation findings (from `python src/wpp_data.py`)
- The death-rate series is fully populated (101 ages by 151 years). Infant
  mortality fell from 0.164 (1950) to 0.030 (2023). Life expectancy rose from
  31.2 (1950) to 65.9 (2023).
- UN probabilistic life expectancy: 70.2 [62.4, 77.3] in 2050; 73.2 [63.5, 83.4]
  in 2074.
- **The population gap (a key finding):** WPP puts 2023 at about 2,728,905 against
  the 2024 census of about 2,422,712, so WPP is roughly 12.6% higher. WPP also
  runs above every historical Gambian census (1973, 1983, 1993, 2003, 2013). Its
  2024 numbers predate the new census, so the projection has to reconcile its
  starting point to the census base (allowing for census under-enumeration). This
  is exactly the kind of result this project exists to surface.

## Download log

| Date | Source | File(s) | Saved to | Notes |
|---|---|---|---|---|
| 2026-06-23 | WPP 2024 (`PPgp/wpp2024` @ `2da7768`) | 31 GMB series | `data/raw/wpp2024/` | done, via `src/fetch_wpp.py` |
| pending | GBoS 2024 PHC | age by sex tables (exact totals) | `data/raw/census2024/` | confirm 2,422,712 and the single-age base |
| pending | DHS | FR289, 2019-20 | `data/raw/dhs/` | child and adult mortality anchors |
| pending | HDSS | Farafenni/Basse life tables | `data/raw/hdss/` | digitise from the papers (validation) |
