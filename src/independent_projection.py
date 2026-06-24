"""
THE HEADLINE: an independent probabilistic population projection for The Gambia
to 2074, based on the 2024 Census.

Pulls together everything:
  * BASE  : 2024 Census total (2,422,712) + its broad age structure (40.8% <15,
            3.0% 65+), distributed to single ages/sexes using the WPP 2023 shape.
            This re-bases off the new census (WPP 2023 was ~13% too high).
  * MORTALITY : Lee-Carter, sex-specific, with forecast uncertainty propagated by
                simulating k_t paths (random walk + drift).
  * FERTILITY : WPP TFR medium with uncertainty from its 95% projection bounds.
  * MIGRATION : WPP medium net migration (central scenario).
  * ENGINE    : the cohort-component engine validated in projection.py, here
                vectorised over simulations so we get full credible intervals.

Outputs population by age/sex to 2074 with 80%/95% credible intervals, plus
dependency ratios and a fan-chart figure.

Run:  python src/independent_projection.py
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd

import figstyle  # noqa: F401  (applies dark, transparent figure style on import)
import wpp_data as w
from leecarter import fit_lee_carter, forecast_kt, _a0_coale_demeny

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures"
CENSUS_2024 = 2_422_712
PROPS = {"u15": 0.408, "mid": 0.562, "o65": 0.030}   # GBoS 2024 PHC broad ages


# --- vectorised life-table person-years L(x) for survival ratios -------------
def Lx_from_mx(age, mx, sex):
    """Return L(x) with radix 1 for each column of mx (ages x sims)."""
    age = np.asarray(age, float); mx = np.atleast_2d(mx)
    if mx.shape[0] != age.shape[0]:
        mx = mx.T
    A, S = mx.shape
    n = np.empty(A); n[:-1] = np.diff(age); n[-1] = np.inf
    nc = n[:, None]
    ax = np.where(np.isinf(nc), 0.0, nc / 2.0) * np.ones((A, S))
    if age[0] == 0:
        ax[0, :] = _a0_coale_demeny(mx[0, :], sex)
    closed = ~np.isinf(n)
    qx = np.empty((A, S))
    qx[closed] = (nc[closed] * mx[closed]) / (1 + (nc[closed] - ax[closed]) * mx[closed])
    qx = np.clip(qx, 0, 1); qx[-1, :] = 1.0
    lx = np.empty((A, S)); lx[0, :] = 1.0
    lx[1:, :] = np.cumprod(1 - qx[:-1, :], axis=0)
    dx = lx * qx
    Lx = np.empty((A, S))
    Lx[:-1, :] = nc[:-1] * lx[1:, :] + ax[:-1, :] * dx[:-1, :]
    Lx[-1, :] = lx[-1, :] / mx[-1, :]
    return Lx


def census_base(ages):
    """2024 base population by single age & sex, reconciled to the census."""
    pf = w.load_pop("female").pivot(index="age", columns="year", values="pop")[2023].to_numpy(float) * 1000
    pm = w.load_pop("male").pivot(index="age", columns="year", values="pop")[2023].to_numpy(float) * 1000
    grp = np.where(ages < 15, "u15", np.where(ages <= 64, "mid", "o65"))
    base_f, base_m = pf.copy(), pm.copy()
    for g, prop in PROPS.items():
        mask = grp == g
        wpp_g = pf[mask].sum() + pm[mask].sum()
        scale = (CENSUS_2024 * prop) / wpp_g
        base_f[mask] *= scale; base_m[mask] *= scale
    return base_f, base_m


def mortality_paths(sex, horizon, n_sims, seed):
    """Simulated m(x,t) futures: returns array (horizon, A, n_sims)."""
    mx = w.mx_matrix(sex, year_max=2023)
    fit = fit_lee_carter(mx)
    _, kt = forecast_kt(fit, horizon, n_sims=n_sims, seed=seed)   # (n_sims, horizon)
    out = np.empty((horizon, len(fit.ages), n_sims))
    for h in range(horizon):
        out[h] = np.exp(fit.ax[:, None] + np.outer(fit.bx, kt[:, h]))
    return fit.ages.astype(float), out


def main():
    BASE, END = 2024, 2074
    H = END - BASE
    S = 1000
    ages = np.arange(0, 101)
    rng = np.random.default_rng(0)

    base_f, base_m = census_base(ages)
    print(f"2024 census-reconciled base: {base_f.sum()+base_m.sum():,.0f} "
          f"(target {CENSUS_2024:,.0f})")

    agesF, mxF = mortality_paths("female", H, S, seed=1)
    _, mxM = mortality_paths("male", H, S, seed=2)

    # fertility: TFR medium + uncertainty from 95% bounds
    tfr_med = dict(zip(*[w.load_series("tfrprojMed_GMB.tsv", "tfr")[c] for c in ["year", "tfr"]]))
    tfr_lo = dict(zip(*[w.load_series("tfrproj95l_GMB.tsv", "tfr")[c] for c in ["year", "tfr"]]))
    tfr_hi = dict(zip(*[w.load_series("tfrproj95u_GMB.tsv", "tfr")[c] for c in ["year", "tfr"]]))
    # %ASFR (constant), SRB, migration (medium)
    from projection import load_pasfr, _series_dict
    pasfr = load_pasfr(ages)
    srb_raw = _series_dict("sexRatio_GMB.tsv", "srb"); last_srb = srb_raw[max(srb_raw)]
    srb = (last_srb / 100 if last_srb > 5 else last_srb)
    mig = _series_dict("migrationprojMed_GMB.tsv", "mig")
    mig = {y: (v * 1000 if abs(v) < 1000 else v) for y, v in mig.items()}

    # sample TFR paths (S,) per year
    work = ages <= 64
    pf = np.repeat(base_f[:, None], S, axis=1)
    pm = np.repeat(base_m[:, None], S, axis=1)
    tot_paths = np.empty((H + 1, S)); tot_paths[0] = pf.sum(0) + pm.sum(0)
    store = {BASE: (pf.copy(), pm.copy())}

    for h in range(H):
        yr = BASE + 1 + h
        Lf = Lx_from_mx(ages, mxF[h], "female"); Lm = Lx_from_mx(ages, mxM[h], "male")
        Sf = Lf[1:] / Lf[:-1]; Sm = Lm[1:] / Lm[:-1]
        sopen_f = np.exp(-mxF[h][-1]); sopen_m = np.exp(-mxM[h][-1])
        sbirth_f = Lf[0]; sbirth_m = Lm[0]

        new_f = np.zeros((101, S)); new_m = np.zeros((101, S))
        new_f[1:] = pf[:-1] * Sf; new_m[1:] = pm[:-1] * Sm
        new_f[-1] += pf[-1] * sopen_f; new_m[-1] += pm[-1] * sopen_m

        # fertility with uncertainty
        med, lo, hi = tfr_med.get(yr), tfr_lo.get(yr), tfr_hi.get(yr)
        sd = max((hi - lo) / (2 * 1.96), 1e-3)
        tfr = np.clip(rng.normal(med, sd, size=S), 0.5, None)
        avg_f = 0.5 * (pf + new_f)
        births = (pasfr[:, None] * avg_f).sum(0) * tfr            # (S,)
        f_b = births / (1 + srb); m_b = births * srb / (1 + srb)
        new_f[0] = f_b * sbirth_f; new_m[0] = m_b * sbirth_m

        # migration (central), distributed over working ages by population share
        nm = mig.get(yr, 0.0)
        bt = pf + pm; wsum = bt[work].sum(0)
        share = np.where(work[:, None], bt / np.where(wsum > 0, wsum, 1), 0.0)
        ffrac = pf.sum(0) / (pf.sum(0) + pm.sum(0))
        new_f = np.maximum(new_f + nm * share * ffrac, 0)
        new_m = np.maximum(new_m + nm * share * (1 - ffrac), 0)

        pf, pm = new_f, new_m
        tot_paths[h + 1] = pf.sum(0) + pm.sum(0)
        if yr in (2050, END):
            store[yr] = (pf.copy(), pm.copy())

    # --- headline numbers ---
    def ci(a):
        return np.median(a), np.percentile(a, 2.5), np.percentile(a, 97.5), np.percentile(a, 10), np.percentile(a, 90)

    print("\n=== INDEPENDENT PROJECTION (2024 census base, LC mortality + fertility uncertainty) ===")
    for yr in (2050, END):
        t = tot_paths[yr - BASE]
        m, lo, hi, l80, h80 = ci(t)
        print(f"  {yr} population: {m:,.0f}  [80% {l80:,.0f}-{h80:,.0f}] [95% {lo:,.0f}-{hi:,.0f}]")

    pf_end, pm_end = store[END]
    p = pf_end + pm_end
    odr = 100 * p[ages >= 65].sum(0) / p[(ages >= 15) & (ages <= 64)].sum(0)
    cdr = 100 * p[ages < 15].sum(0) / p[(ages >= 15) & (ages <= 64)].sum(0)
    print(f"  2074 old-age dependency: {np.median(odr):.1f} [{np.percentile(odr,2.5):.1f}, {np.percentile(odr,97.5):.1f}]")
    print(f"  2074 child   dependency: {np.median(cdr):.1f} [{np.percentile(cdr,2.5):.1f}, {np.percentile(cdr,97.5):.1f}]")

    pd.DataFrame({
        "year": np.arange(BASE, END + 1),
        "pop_med": np.median(tot_paths, axis=1),
        "pop_lo95": np.percentile(tot_paths, 2.5, axis=1),
        "pop_hi95": np.percentile(tot_paths, 97.5, axis=1),
        "pop_lo80": np.percentile(tot_paths, 10, axis=1),
        "pop_hi80": np.percentile(tot_paths, 90, axis=1),
    }).to_csv(ROOT / "data" / "processed" / "independent_projection.csv", index=False)

    _figure(ages, np.arange(BASE, END + 1), tot_paths, store, BASE, END)
    print("\nSaved -> data/processed/independent_projection.csv + figures/independent_projection.png")


def _figure(ages, years, tot_paths, store, BASE, END):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    med = np.median(tot_paths, axis=1) / 1e6
    f, axs = plt.subplots(1, 2, figsize=(13, 4.8))
    axs[0].fill_between(years, np.percentile(tot_paths, 2.5, axis=1) / 1e6,
                        np.percentile(tot_paths, 97.5, axis=1) / 1e6, alpha=.2, color="C3", label="95% CI")
    axs[0].fill_between(years, np.percentile(tot_paths, 10, axis=1) / 1e6,
                        np.percentile(tot_paths, 90, axis=1) / 1e6, alpha=.3, color="C3", label="80% CI")
    axs[0].plot(years, med, color="C3", label="median (independent)")
    # WPP medium total (for comparison)
    pfp = pd.read_csv(ROOT / "data/raw/wpp2024/popFprojMed_GMB.tsv", sep="\t")
    pmp = pd.read_csv(ROOT / "data/raw/wpp2024/popMprojMed_GMB.tsv", sep="\t")
    yc = [c for c in pfp.columns if str(c).isdigit() and BASE <= int(c) <= END]
    wy = [int(c) for c in yc]; wv = [(pfp[c].sum() + pmp[c].sum()) * 1000 / 1e6 for c in yc]
    axs[0].plot(wy, wv, "--", color="C0", label="WPP medium")
    axs[0].set_title("The Gambia: independent population projection to 2074")
    axs[0].set_xlabel("year"); axs[0].set_ylabel("population (millions)"); axs[0].legend(fontsize=8)

    f24, m24 = store[BASE]; f74, m74 = store[END]
    axs[1].barh(ages, -m24.mean(1) / 1000, color="C0", alpha=.35, label="M 2024")
    axs[1].barh(ages, f24.mean(1) / 1000, color="C3", alpha=.35, label="F 2024")
    axs[1].step(-np.median(m74, 1) / 1000, ages, color="C0", lw=1.3, label="M 2074")
    axs[1].step(np.median(f74, 1) / 1000, ages, color="C3", lw=1.3, label="F 2074")
    axs[1].set_title("Population pyramid: 2024 vs 2074"); axs[1].set_xlabel("thousands"); axs[1].set_ylabel("age"); axs[1].legend(fontsize=7)
    f.tight_layout(); f.savefig(FIG / "independent_projection.png", dpi=130); plt.close(f)


if __name__ == "__main__":
    main()
