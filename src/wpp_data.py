"""
Load and tidy the WPP 2024 Gambia extracts produced by ``fetch_wpp.py``.

The raw extracts are wide, tab-separated tables (one row per age, one column per
year). This module reshapes them into tidy long-format DataFrames that the
mortality model and projection consume, and exposes small convenience loaders.

All functions return pandas DataFrames with explicit, documented columns.
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "wpp2024"

_YEAR_RE = re.compile(r"^\d{4}$")          # single-year column, e.g. "2024"
_PERIOD_RE = re.compile(r"^\d{4}-\d{4}$")  # abridged period, e.g. "2020-2025"


def _year_columns(df: pd.DataFrame) -> list[str]:
    return [c for c in df.columns if _YEAR_RE.match(str(c)) or _PERIOD_RE.match(str(c))]


def _read_wide(fname: str) -> pd.DataFrame:
    path = RAW / fname
    if not path.exists():
        raise FileNotFoundError(f"{path} - run `python src/fetch_wpp.py` first.")
    return pd.read_csv(path, sep="\t")


def _melt(fname: str, id_vars: list[str], value_name: str) -> pd.DataFrame:
    df = _read_wide(fname)
    years = _year_columns(df)
    long = df.melt(id_vars=id_vars, value_vars=years, var_name="year", value_name=value_name)
    # single-year columns -> int year; keep period strings as-is
    if all(_YEAR_RE.match(str(y)) for y in years):
        long["year"] = long["year"].astype(int)
    long[value_name] = pd.to_numeric(long[value_name], errors="coerce")
    return long


# --- mortality ---------------------------------------------------------------
def load_mx(sex: str = "both") -> pd.DataFrame:
    """Age-specific central death rate m(x, t). Columns: age, year, mx."""
    fname = {"both": "mxB_GMB.tsv", "female": "mxF_GMB.tsv", "male": "mxM_GMB.tsv"}[sex]
    long = _melt(fname, id_vars=["country_code", "country", "age"], value_name="mx")
    return long[["age", "year", "mx"]].sort_values(["year", "age"]).reset_index(drop=True)


def mx_matrix(sex: str = "both", year_max: int | None = 2023) -> pd.DataFrame:
    """Wide age x year matrix of m(x,t) (rows=age, cols=year), for model fitting.

    Defaults to estimates only (<= 2023); pass year_max=None to include the WPP
    medium projection years too.
    """
    long = load_mx(sex)
    if year_max is not None:
        long = long[long["year"] <= year_max]
    return long.pivot(index="age", columns="year", values="mx")


# --- population --------------------------------------------------------------
def load_pop(sex: str = "female") -> pd.DataFrame:
    """Population by single age (mid-year estimates). Columns: age, year, pop."""
    fname = {"female": "popF_GMB.tsv", "male": "popM_GMB.tsv"}[sex]
    long = _melt(fname, id_vars=["country_code", "country", "age"], value_name="pop")
    return long[["age", "year", "pop"]].sort_values(["year", "age"]).reset_index(drop=True)


# --- single-series indicators (e0, tfr, migration, counts) -------------------
def load_series(fname: str, value_name: str) -> pd.DataFrame:
    """A one-row-per-country indicator (e.g. e0, tfr). Columns: year, <value>."""
    long = _melt(fname, id_vars=["country_code", "country"], value_name=value_name)
    return long[["year", value_name]].sort_values("year").reset_index(drop=True)


def _verify() -> None:
    """Print an integrity + triangulation summary on the fetched data."""
    print("=== WPP 2024 - The Gambia: data integrity check ===\n")

    mxB = load_mx("both")
    print(f"m(x,t): ages {mxB.age.min()}-{mxB.age.max()}, "
          f"years {mxB.year.min()}-{mxB.year.max()}, "
          f"{mxB.mx.notna().sum()} non-null cells")
    # infant mortality decline
    m0 = mxB[mxB.age == 0].set_index("year")["mx"]
    print(f"  infant m(0): 1950={m0.loc[1950]:.4f}  2000={m0.loc[2000]:.4f}  2023={m0.loc[2023]:.4f}")

    e0 = load_series("e0B_GMB.tsv", "e0")
    print("\nLife expectancy at birth (both sexes), WPP estimates:")
    for y in [1950, 1980, 2000, 2023]:
        if y in set(e0.year):
            print(f"  e0({y}) = {e0.set_index('year').loc[y, 'e0']:.2f}")

    e0p = load_series("e0Bproj_GMB.tsv", "e0")
    e0lo = load_series("e0Bproj95l_GMB.tsv", "e0")
    e0hi = load_series("e0Bproj95u_GMB.tsv", "e0")
    proj = e0p.merge(e0lo, on="year", suffixes=("", "_lo")).merge(e0hi, on="year", suffixes=("", "_hi"))
    print("\nUN probabilistic e0 projection (median [95% PI]):")
    for y in [2050, 2074, 2100]:
        row = proj[proj.year == y]
        if len(row):
            r = row.iloc[0]
            print(f"  {y}: {r['e0']:.2f}  [{r['e0_lo']:.2f}, {r['e0_hi']:.2f}]")

    # --- TRIANGULATION: WPP 2024 population vs the 2024 Census -----------------
    # WPP 2024 was finalised before the 2024 Gambian census results, so this gap
    # is a substantive finding, not a bug. The exact census total is still to be
    # confirmed from the GBoS primary report (~2.4M per preliminary reporting).
    # GBoS 2024 PHC preliminary total, as reported in Gambian media (Alkamba
    # Times / Kerr Fatou): 2,422,712 (F 51% / M 49%). Confirm against the GBoS
    # primary report before final use - a slightly different figure also circulates.
    CENSUS_2024_APPROX = 2_422_712
    popF = load_pop("female"); popM = load_pop("male")
    tot = (popF.groupby("year")["pop"].sum() + popM.groupby("year")["pop"].sum())  # thousands
    latest = int(tot.index.max())
    print(f"\nTotal population (WPP estimates, persons):")
    for y in [1973, 1983, 1993, 2003, 2013, latest]:
        if y in tot.index:
            print(f"  {y}: {tot.loc[y]*1000:,.0f}")
    print(f"  2024 Census (GBoS): ~{CENSUS_2024_APPROX:,.0f} (preliminary; exact TBC)")
    diff = tot.loc[latest] * 1000 - CENSUS_2024_APPROX
    print(f"  WPP {latest} ({tot.loc[latest]*1000:,.0f}) vs Census ~2.4M: "
          f"gap ~= {diff:,.0f} ({diff/CENSUS_2024_APPROX*100:+.1f}%)")
    print("  -> WPP appears to overstate the population; the new census implies a")
    print("     downward revision. A jump-off reconciliation is needed (see protocol).")


if __name__ == "__main__":
    _verify()
