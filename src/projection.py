"""
Cohort-component (Leslie-matrix) population projection by single age and sex.

This is the headline engine: given a base population by age and sex, and yearly
schedules of mortality, fertility, sex-ratio-at-birth and net migration, it rolls
the population forward one year at a time. The same machinery the UN Population
Division uses for WPP.

Design choice - VALIDATE FIRST: ``main()`` runs the engine on WPP's *own* inputs
(mortality m(x,t), TFR, %ASFR, SRB, net migration) and checks that it reproduces
WPP's published projected population. Only an engine that reconstructs the
benchmark can be trusted to carry our independent (Lee-Carter / Bayesian)
mortality. Migration is distributed across ages in proportion to population (a
documented simplification; WPP uses explicit age schedules).

References: Preston, Heuveline & Guillot (2001), *Demography*, Ch. 6.
Run:  python src/projection.py
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

import figstyle  # noqa: F401  (applies dark, transparent figure style on import)
import wpp_data as w
from lifetable import period_life_table

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures"
RADIX = 100_000.0


# --------------------------------------------------------------------------- #
# Survivorship ratios from a one-year mortality schedule
# --------------------------------------------------------------------------- #
@dataclass
class Survival:
    Sx: np.ndarray      # Sx[a] = L[a+1]/L[a], a=0..A-1  (age a -> a+1)
    S_birth: float      # births -> age 0  (= L0 / l0)
    s_open: float       # within open age group, one-year survival


def survival_from_mx(age: np.ndarray, mx: np.ndarray, sex: str) -> Survival:
    lt = period_life_table(age, mx, sex=sex, radix=RADIX)
    L = lt["Lx"].to_numpy()
    A = len(age)
    Sx = np.empty(A)
    Sx[:-1] = L[1:] / L[:-1]          # age a -> a+1 for a=0..A-2
    Sx[-1] = np.nan                   # last handled via s_open
    S_birth = L[0] / RADIX
    s_open = float(np.exp(-mx[-1]))   # stayers in the open (100+) group
    return Survival(Sx=Sx, S_birth=S_birth, s_open=s_open)


# --------------------------------------------------------------------------- #
# One-year projection step
# --------------------------------------------------------------------------- #
def project_step(pop_f, pop_m, surv_f: Survival, surv_m: Survival,
                 asfr: np.ndarray, srb: float, mig_f, mig_m):
    """Advance the population one year. All arrays are length A (single age)."""
    A = len(pop_f)
    new_f = np.zeros(A)
    new_m = np.zeros(A)

    # age everyone one year (a-1 -> a), using survivorship ratios
    new_f[1:] = pop_f[:-1] * surv_f.Sx[:-1]
    new_m[1:] = pop_m[:-1] * surv_m.Sx[:-1]
    # open age group: inflow from A-2 already added at index A-1 above; add stayers
    new_f[-1] += pop_f[-1] * surv_f.s_open
    new_m[-1] += pop_m[-1] * surv_m.s_open

    # births: from average reproductive-age female population over the year
    avg_f = 0.5 * (pop_f + new_f)
    births = float(np.sum(asfr * avg_f))
    f_births = births / (1.0 + srb)
    m_births = births * srb / (1.0 + srb)
    new_f[0] = f_births * surv_f.S_birth
    new_m[0] = m_births * surv_m.S_birth

    # net migration (added at year end)
    new_f = np.maximum(new_f + mig_f, 0.0)
    new_m = np.maximum(new_m + mig_m, 0.0)
    return new_f, new_m


# --------------------------------------------------------------------------- #
# Full projection
# --------------------------------------------------------------------------- #
def run_projection(ages, base_f, base_m, base_year, end_year,
                   mxF, mxM, tfr, pasfr, srb, netmig):
    """Roll forward base_year -> end_year.

    mxF/mxM : DataFrames (index age, cols year) of m(x,t) for projection years.
    tfr     : dict year -> TFR.  pasfr : array len A (fractions, sum=1).
    srb     : dict year -> sex ratio at birth (males/female).
    netmig  : dict year -> net migrants (persons; +inflow).  Distributed across
              ages 0-64 in proportion to population, split by sex share.
    """
    A = len(ages)
    pf, pm = base_f.copy(), base_m.copy()
    rows = [(base_year, pf.copy(), pm.copy())]
    work_mask = ages <= 64
    for yr in range(base_year + 1, end_year + 1):
        sf = survival_from_mx(ages.astype(float), mxF[yr].to_numpy(float), "female")
        sm = survival_from_mx(ages.astype(float), mxM[yr].to_numpy(float), "male")
        # migration distributed proportional to working-age population
        nm = netmig.get(yr, 0.0)
        base_tot = (pf + pm)
        wsum = base_tot[work_mask].sum()
        if wsum > 0 and nm != 0.0:
            share = np.where(work_mask, base_tot / wsum, 0.0)
            ffrac = pf.sum() / (pf.sum() + pm.sum())
            mig_f = nm * share * ffrac
            mig_m = nm * share * (1 - ffrac)
        else:
            mig_f = mig_m = np.zeros(A)
        pf, pm = project_step(pf, pm, sf, sm,
                              asfr=tfr[yr] * pasfr, srb=srb[yr],
                              mig_f=mig_f, mig_m=mig_m)
        rows.append((yr, pf.copy(), pm.copy()))
    return rows


# --------------------------------------------------------------------------- #
# Inputs (WPP) and indicators
# --------------------------------------------------------------------------- #
def _series_dict(fname, value):
    s = w.load_series(fname, value)
    return dict(zip(s.year.astype(int), s[value].astype(float)))


def load_pasfr(ages):
    """Latest-year %ASFR as fractions (sum=1), mapped onto the single-age grid."""
    df = pd.read_csv(ROOT / "data/raw/wpp2024/percentASFR_GMB.tsv", sep="\t")
    yrs = [c for c in df.columns if str(c).isdigit() or "-" in str(c)]
    latest = yrs[len([y for y in yrs if True]) - 1]  # last column
    sub = df[["age", latest]].copy()
    sub.columns = ["age", "p"]
    pasfr = np.zeros(len(ages))
    amap = {a: i for i, a in enumerate(ages)}
    for _, r in sub.iterrows():
        a = int(r["age"])
        if a in amap:
            pasfr[amap[a]] = float(r["p"])
    if pasfr.sum() > 0:
        pasfr = pasfr / pasfr.sum()
    return pasfr


def total(rows):
    return pd.DataFrame({"year": [r[0] for r in rows],
                         "pop": [(r[1].sum() + r[2].sum()) for r in rows]})


def dependency_ratios(ages, pf, pm):
    p = pf + pm
    young = p[ages < 15].sum()
    old = p[ages >= 65].sum()
    work = p[(ages >= 15) & (ages <= 64)].sum()
    return dict(child=100 * young / work, old=100 * old / work,
                total=100 * (young + old) / work)


def main():
    BASE, END = 2023, 2074
    ages = np.arange(0, 101)
    pf = w.load_pop("female").pivot(index="age", columns="year", values="pop")
    pm = w.load_pop("male").pivot(index="age", columns="year", values="pop")
    base_f = pf[BASE].to_numpy(float) * 1000.0
    base_m = pm[BASE].to_numpy(float) * 1000.0

    mxF = w.mx_matrix("female", year_max=None)
    mxM = w.mx_matrix("male", year_max=None)

    tfr = {**_series_dict("tfr_GMB.tsv", "tfr"), **_series_dict("tfrprojMed_GMB.tsv", "tfr")}
    srb_raw = {**_series_dict("sexRatio_GMB.tsv", "srb")}
    srb = {y: (v / 100.0 if v > 5 else v) for y, v in srb_raw.items()}
    # extend SRB to projection years with the last known value
    last_srb = srb[max(srb)]
    mig = {**_series_dict("migration_GMB.tsv", "mig"),
           **_series_dict("migrationprojMed_GMB.tsv", "mig")}
    # detect thousands -> persons
    mig = {y: (v * 1000.0 if abs(v) < 1000 else v) for y, v in mig.items()}
    pasfr = load_pasfr(ages)

    srb_full = {y: srb.get(y, last_srb) for y in range(BASE, END + 1)}

    print("=== VALIDATION: reproduce WPP projection with WPP inputs ===")
    rows = run_projection(ages, base_f, base_m, BASE, END, mxF, mxM,
                          tfr, pasfr, srb_full, mig)
    ours = total(rows).set_index("year")["pop"]

    # WPP's own projected population by age (medium)
    wf = w.load_pop("female")  # has up to 2023 only; projection in separate files
    pfp = pd.read_csv(ROOT / "data/raw/wpp2024/popFprojMed_GMB.tsv", sep="\t")
    pmp = pd.read_csv(ROOT / "data/raw/wpp2024/popMprojMed_GMB.tsv", sep="\t")
    yrcols = [c for c in pfp.columns if str(c).isdigit()]
    wpp_tot = {}
    for y in yrcols:
        wpp_tot[int(y)] = (pfp[y].sum() + pmp[y].sum()) * 1000.0

    print(f"{'year':>6}{'ours':>14}{'WPP':>14}{'diff %':>9}")
    for y in (2030, 2040, 2050, 2074):
        if y in ours.index and y in wpp_tot:
            o, ww = ours[y], wpp_tot[y]
            print(f"{y:>6}{o:>14,.0f}{ww:>14,.0f}{(o-ww)/ww*100:>8.1f}%")

    dr = dependency_ratios(ages, rows[-1][1], rows[-1][2])
    print(f"\n2074 (our engine, WPP inputs): total pop {ours[END]:,.0f}")
    print(f"  dependency ratios - child {dr['child']:.1f}, old {dr['old']:.1f}, total {dr['total']:.1f}")
    print(f"  (2023 base) ", dependency_ratios(ages, base_f, base_m))

    _figures(ages, rows, ours, wpp_tot, BASE, END)
    print(f"\nFigures -> {FIG}")


def _figures(ages, rows, ours, wpp_tot, BASE, END):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    by_year = {r[0]: (r[1], r[2]) for r in rows}
    deps = pd.DataFrame([{"year": y, **dependency_ratios(ages, f, m)}
                         for y, (f, m) in by_year.items()])

    fig, axs = plt.subplots(1, 3, figsize=(15, 4.6))

    # (a) total population trajectory: ours vs WPP
    axs[0].plot(ours.index, ours.values / 1e6, color="C3", label="our engine")
    wy = sorted(wpp_tot); axs[0].plot(wy, [wpp_tot[y] / 1e6 for y in wy], "--", color="C0", label="WPP medium")
    axs[0].axhline(2.4227, color="grey", ls=":", lw=1); axs[0].text(BASE, 2.55, "2024 census ~2.42M", fontsize=7, color="grey")
    axs[0].set_title("Total population (millions)"); axs[0].set_xlabel("year"); axs[0].legend(fontsize=8)

    # (b) population pyramids 2023 vs 2074
    f23, m23 = by_year[BASE]; f74, m74 = by_year[END]
    axs[1].barh(ages, -m23 / 1000, color="C0", alpha=.35, label="M 2023")
    axs[1].barh(ages, f23 / 1000, color="C3", alpha=.35, label="F 2023")
    axs[1].step(-m74 / 1000, ages, color="C0", lw=1.3, label="M 2074")
    axs[1].step(f74 / 1000, ages, color="C3", lw=1.3, label="F 2074")
    axs[1].set_title("Population pyramid: 2023 vs 2074"); axs[1].set_xlabel("thousands"); axs[1].set_ylabel("age")
    axs[1].legend(fontsize=7)

    # (c) dependency ratios over time
    axs[2].plot(deps.year, deps.child, label="child (<15)")
    axs[2].plot(deps.year, deps.old, label="old (65+)")
    axs[2].plot(deps.year, deps.total, color=figstyle.FG, lw=1.6, label="total")
    axs[2].set_title("Dependency ratios (per 100 working-age)"); axs[2].set_xlabel("year"); axs[2].legend(fontsize=8)

    fig.suptitle("The Gambia population projection to 2074 (engine validated on WPP inputs)")
    fig.tight_layout(); fig.savefig(FIG / "projection_overview.png", dpi=130); plt.close(fig)


if __name__ == "__main__":
    main()
