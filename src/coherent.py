"""
Li-Lee (2005) coherent mortality forecast for The Gambia within a West-African
reference group.

A single-population Lee-Carter can drift to implausible long-run mortality and,
as we showed, understates uncertainty relative to the UN. The coherent method
decomposes each country's log-rates into a COMMON factor shared by the group
plus a country-specific deviation that is forecast as a mean-reverting AR(1) —
so members cannot diverge without bound and The Gambia borrows strength from its
neighbours' better-estimated common trend.

    log m_i(x,t) = a_i(x) + B(x) K(t) + b_i(x) k_i(t)
      common  (B, K): SVD of the group-average centred log-rates; K ~ RW+drift
      specific(b_i,k_i): SVD of the residual; k_i ~ AR(1) (mean-reverting)

Reference: Li & Lee (2005), Demography 42(3):575-594.
Run:  python src/coherent.py
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd

from leecarter import e0_from_mx, fit_lee_carter, forecast_e0
import figstyle  # noqa: F401  (applies dark, transparent figure style on import)
import wpp_data as w

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "data" / "raw" / "wpp2024_refgroup"
FIG = ROOT / "figures"
GAMBIA = 270


def load_group(sex="B", year_max=2023):
    df = pd.read_csv(REF / f"mx{sex}_ref.tsv", sep="\t")
    yrs = [c for c in df.columns if str(c).isdigit() and int(c) <= year_max]
    ages = np.sort(df["age"].unique())
    mats = {}
    for code, g in df.groupby("country_code"):
        g = g.set_index("age").loc[ages, yrs]
        mats[int(code)] = np.log(g.to_numpy(float))   # log m(x,t)
    years = np.array([int(y) for y in yrs])
    return ages.astype(float), years, mats


def _svd1(M):
    """Rank-1 factor of matrix M (ages x years): returns (beta sum=1, kappa)."""
    U, S, Vt = np.linalg.svd(M, full_matrices=False)
    b = U[:, 0]; k = S[0] * Vt[0, :]
    s = b.sum()
    b, k = b / s, k * s
    if b.sum() < 0:
        b, k = -b, -k
    return b, k


def fit_lilee(ages, years, mats):
    codes = list(mats)
    A, T = len(ages), len(years)
    a = {i: mats[i].mean(axis=1) for i in codes}               # a_i(x)
    centred = {i: mats[i] - a[i][:, None] for i in codes}
    Ybar = np.mean([centred[i] for i in codes], axis=0)        # group average
    B, K = _svd1(Ybar)                                         # common factor
    b, k = {}, {}
    for i in codes:
        R = centred[i] - np.outer(B, K)
        b[i], k[i] = _svd1(R)                                  # country residual
    return dict(ages=ages, years=years, a=a, B=B, K=K, b=b, k=k)


def _ar1(series):
    """Fit k_t = c + phi k_{t-1} + eps. Returns (c, phi, sigma, last)."""
    y, x = series[1:], series[:-1]
    X = np.vstack([np.ones_like(x), x]).T
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    c, phi = coef
    resid = y - X @ coef
    sigma = resid.std(ddof=2) if len(resid) > 2 else resid.std()
    return float(c), float(np.clip(phi, -0.999, 0.999)), float(sigma), float(series[-1])


def forecast_gambia_e0(fit, horizon, sex="both", n_sims=4000, seed=3):
    ages, years = fit["ages"], fit["years"]
    a, B, K, b, k = fit["a"][GAMBIA], fit["B"], fit["K"], fit["b"][GAMBIA], fit["k"][GAMBIA]
    rng = np.random.default_rng(seed)
    T = len(K)
    # common K ~ random walk with drift (+ drift uncertainty)
    dK = (K[-1] - K[0]) / (T - 1)
    sK = np.std(np.diff(K) - dK, ddof=1)
    drift = rng.normal(dK, sK / np.sqrt(T - 1), size=n_sims)
    Kf = K[-1] + np.cumsum(drift[:, None] + rng.normal(0, sK, size=(n_sims, horizon)), axis=1)
    # country-specific k ~ AR(1) mean-reverting
    c, phi, sig, klast = _ar1(k)
    kf = np.empty((n_sims, horizon)); prev = np.full(n_sims, klast)
    for h in range(horizon):
        prev = c + phi * prev + rng.normal(0, sig, size=n_sims)
        kf[:, h] = prev

    yrs = np.arange(int(years[-1]) + 1, int(years[-1]) + 1 + horizon)
    recs = []
    for h in range(horizon):
        logm = a[None, :] + np.outer(Kf[:, h], B) + np.outer(kf[:, h], b)   # (S, A)
        e0 = e0_from_mx(ages, np.exp(logm).T, sex)
        recs.append((yrs[h], np.median(e0), np.percentile(e0, 2.5), np.percentile(e0, 97.5),
                     np.percentile(e0, 10), np.percentile(e0, 90)))
    return pd.DataFrame(recs, columns=["year", "e0_med", "e0_lo95", "e0_hi95", "e0_lo80", "e0_hi80"])


def main():
    FIG.mkdir(exist_ok=True)
    ages, years, mats = load_group("B")
    print(f"Li-Lee coherent group ({len(mats)} countries), ages {ages.min():.0f}-{ages.max():.0f}, "
          f"years {years.min()}-{years.max()}")
    fit = fit_lilee(ages, years, mats)
    c, phi, sig, _ = _ar1(fit["k"][GAMBIA])
    print(f"  common K drift = {(fit['K'][-1]-fit['K'][0])/(len(fit['K'])-1):.3f}/yr | "
          f"Gambia-specific k AR(1) phi = {phi:.2f} (mean-reverting)")

    fc = forecast_gambia_e0(fit, horizon=2074 - int(years[-1]))
    for yr in (2050, 2074):
        r = fc[fc.year == yr].iloc[0]
        print(f"  coherent e0({yr}) = {r.e0_med:.2f}  [95% PI {r.e0_lo95:.2f}, {r.e0_hi95:.2f}]")
    print("  Compare 2074 — single-pop LC 75.6 [73.0,77.8]; Bayesian 75.0 [73.8,76.2]; "
          "WPP 73.2 [63.5,83.4]")
    fc.to_csv(ROOT / "data" / "processed" / "coherent_e0_forecast.csv", index=False)

    # figure: coherent vs single-pop vs WPP
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    mxG = w.mx_matrix("both", year_max=2023)
    sp = forecast_e0(fit_lee_carter(mxG), horizon=2074 - 2023)
    e0h = w.load_series("e0B_GMB.tsv", "e0")
    p = w.load_series("e0Bproj_GMB.tsv", "e0")
    lo = w.load_series("e0Bproj95l_GMB.tsv", "e0"); hi = w.load_series("e0Bproj95u_GMB.tsv", "e0")
    m = p.merge(lo, on="year", suffixes=("", "_lo")).merge(hi, on="year", suffixes=("", "_hi"))
    m = m[m.year <= 2074]
    f, ax = plt.subplots(figsize=(8.5, 4.8))
    ax.plot(e0h.year, e0h.e0, color=figstyle.FG, lw=1.3, label="WPP estimates")
    ax.fill_between(m.year, m.e0_lo, m.e0_hi, alpha=0.12, color="C0", label="WPP 95% PI")
    ax.plot(m.year, m.e0, color="C0", ls="--", label="WPP median")
    ax.fill_between(sp.year, sp.e0_lo95, sp.e0_hi95, alpha=0.15, color="C2", label="single-pop LC 95%")
    ax.plot(sp.year, sp.e0_med, color="C2", label="single-pop LC")
    ax.fill_between(fc.year, fc.e0_lo95, fc.e0_hi95, alpha=0.22, color="C3", label="coherent 95%")
    ax.plot(fc.year, fc.e0_med, color="C3", lw=1.6, label="Li-Lee coherent")
    ax.set_title("e0 forecast — Li-Lee coherent vs single-pop vs WPP (The Gambia)")
    ax.set_xlabel("year"); ax.set_ylabel("e0 (years)"); ax.legend(fontsize=8)
    f.tight_layout(); f.savefig(FIG / "coherent_vs_single_vs_wpp.png", dpi=130); plt.close(f)
    print(f"\nFigure -> figures/coherent_vs_single_vs_wpp.png")


if __name__ == "__main__":
    main()
