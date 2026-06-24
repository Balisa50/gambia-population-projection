"""
Bayesian Poisson Lee-Carter for The Gambia (PyMC / MCMC).

Why Bayesian rather than the SVD baseline (src/leecarter.py): for a small,
data-sparse population the homoskedastic-Gaussian SVD assumption understates
uncertainty (we showed its e0 interval is far narrower than the UN's). Modelling
deaths as Poisson around the Lee-Carter rate and sampling the full posterior
propagates parameter AND forecast uncertainty into every downstream quantity
(rates, life expectancy, eventually the population projection).

Model
-----
    D(x,t) ~ Poisson( E(x,t) * exp(a_x + b_x * k_t) )
    b ~ Dirichlet(1)                      # sum_x b_x = 1, b_x > 0  (resolves the
                                          #   b/k scale identifiability cleanly)
    k_t = k_{t-1} + drift + sigma_k * z_t,   k_0 = 0   (random walk w/ drift)
    a_x ~ Normal(emp_a_x, 1)              # emp_a_x = mean_t log m(x,t)

Deaths/exposure are reconstructed from WPP rates x population (documented as a
reconstruction in the protocol; the structure is what matters here).

Run:  python src/bayes_leecarter.py [draws] [tune]
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pymc as pm
import pytensor.tensor as pt
import arviz as az

import figstyle  # noqa: F401  (applies dark, transparent figure style on import)
import wpp_data as w
from leecarter import e0_from_mx

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures"


def load_fit_data(sex: str = "both", year_max: int = 2023):
    """Aligned arrays for fitting. Returns ages, years, D (int), E, emp_a."""
    mx = w.mx_matrix(sex=sex, year_max=year_max)            # index age, cols year
    popsex = "female" if sex in ("both", "female") else "male"
    # population by single age; for 'both' use F+M
    pf = w.load_pop("female").pivot(index="age", columns="year", values="pop")
    pmm = w.load_pop("male").pivot(index="age", columns="year", values="pop")
    pop = (pf + pmm) if sex == "both" else (pf if sex == "female" else pmm)

    years = [y for y in mx.columns if y in pop.columns and y <= year_max]
    ages = mx.index.to_numpy(int)
    mx = mx.loc[ages, years].to_numpy(float)
    E = pop.loc[ages, years].to_numpy(float) * 1000.0       # thousands -> persons
    D = np.rint(mx * E).astype(int)
    D = np.clip(D, 0, None)
    emp_a = np.log(mx).mean(axis=1)
    return ages, np.array(years, int), D, E, emp_a


def build_model(D: np.ndarray, E: np.ndarray, emp_a: np.ndarray) -> pm.Model:
    """a_x is fixed at the empirical mean log-rate (the classical Lee-Carter
    definition); this removes the a/k identifiability ridge and lets NUTS mix.
    k_t has a free level (k0) plus a random walk with drift, so its level is
    pinned by the (a_x-fixed) likelihood rather than an arbitrary anchor."""
    A, T = D.shape
    with pm.Model() as model:
        b = pm.Dirichlet("b", a=np.ones(A))                 # sum_x b_x = 1, b_x>0
        k0 = pm.Normal("k0", 0.0, 10.0)
        drift = pm.Normal("drift", 0.0, 1.0)
        sigma_k = pm.HalfNormal("sigma_k", 1.0)
        z = pm.Normal("z", 0.0, 1.0, shape=T - 1)           # non-centred innovations
        steps = drift + sigma_k * z
        k = pm.Deterministic("k", k0 + pt.concatenate([[0.0], pt.cumsum(steps)]))
        log_rate = emp_a[:, None] + pt.outer(b, k)          # (A, T)
        pm.Poisson("D", mu=pt.exp(log_rate) * E, observed=D)
    return model


def _post(idata, name):
    arr = idata.posterior[name].values
    return arr.reshape(-1, *arr.shape[2:])                  # (chains*draws, ...)


def forecast_e0_bayes(idata, ages, emp_a, last_year, horizon, sex="both", seed=7):
    b = _post(idata, "b")
    drift = _post(idata, "drift"); sigma_k = _post(idata, "sigma_k")
    k = _post(idata, "k")
    S = b.shape[0]
    k_last = k[:, -1]
    rng = np.random.default_rng(seed)
    innov = rng.normal(0, 1, size=(S, horizon)) * sigma_k[:, None]
    steps = drift[:, None] + innov
    k_future = k_last[:, None] + np.cumsum(steps, axis=1)    # (S, horizon)

    years = np.arange(last_year + 1, last_year + 1 + horizon)
    recs = []
    for h in range(horizon):
        log_rate = emp_a[None, :] + b * k_future[:, h:h + 1]   # (S, A)
        e0 = e0_from_mx(ages.astype(float), np.exp(log_rate).T, sex)  # (S,)
        recs.append((years[h], np.median(e0),
                     np.percentile(e0, 2.5), np.percentile(e0, 97.5),
                     np.percentile(e0, 10), np.percentile(e0, 90)))
    return pd.DataFrame(recs, columns=["year", "e0_med", "e0_lo95", "e0_hi95",
                                       "e0_lo80", "e0_hi80"])


def main():
    draws = int(sys.argv[1]) if len(sys.argv) > 1 else 500
    tune = int(sys.argv[2]) if len(sys.argv) > 2 else 500
    sex = "both"
    ages, years, D, E, emp_a = load_fit_data(sex)
    print(f"Fitting Bayesian Poisson Lee-Carter: ages {ages.min()}-{ages.max()}, "
          f"years {years.min()}-{years.max()}, {D.size} cells, total deaths/yr ~ {D.sum()/len(years):,.0f}")

    model = build_model(D, E, emp_a)
    with model:
        idata = pm.sample(draws=draws, tune=tune, chains=4, cores=1, init="adapt_diag",
                          target_accept=0.95, progressbar=False, random_seed=11)

    s = az.summary(idata, var_names=["drift", "sigma_k", "k0"])
    print(s[["mean", "sd", "r_hat", "ess_bulk"]].to_string())
    full = az.summary(idata, var_names=["b", "k", "drift", "sigma_k", "k0"])
    rhat_max = float(pd.to_numeric(full["r_hat"], errors="coerce").max())
    ess_min = float(pd.to_numeric(full["ess_bulk"], errors="coerce").min())
    print(f"max r-hat = {rhat_max:.3f} (want <1.01) | min ESS = {ess_min:.0f}")

    fc = forecast_e0_bayes(idata, ages, emp_a, int(years.max()),
                           horizon=2074 - int(years.max()), sex=sex)
    for yr in (2050, 2074):
        r = fc[fc.year == yr].iloc[0]
        print(f"  Bayesian e0({yr}) = {r.e0_med:.2f}  [95% PI {r.e0_lo95:.2f}, {r.e0_hi95:.2f}]")
    print("  Compare — classical LC 2074: 75.6 [73.0,77.8]; WPP Bayesian: 73.2 [63.5,83.4]")

    fc.to_csv(ROOT / "data" / "processed" / "bayes_lc_e0_forecast.csv", index=False)
    print(f"\nSaved e0 forecast -> data/processed/bayes_lc_e0_forecast.csv")

    # persist trace (gitignored) + a comparison figure
    traces = ROOT / "data" / "processed" / "traces"; traces.mkdir(exist_ok=True)
    try:
        idata.to_netcdf(str(traces / f"bayes_lc_{sex}.nc"))
    except Exception as e:  # noqa: BLE001
        print("  (trace not saved:", e, ")")
    _make_figure(fc, sex, rhat_max)
    return idata, fc


def _make_figure(fc: pd.DataFrame, sex: str, rhat_max: float) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    e0_hist = w.load_series("e0B_GMB.tsv", "e0")
    f, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(e0_hist.year, e0_hist.e0, color=figstyle.FG, lw=1.4, label="WPP estimates (1950-2023)")
    ax.plot(fc.year, fc.e0_med, color="C3", label="Bayesian LC (median)")
    ax.fill_between(fc.year, fc.e0_lo95, fc.e0_hi95, alpha=0.20, color="C3", label="95% credible")
    ax.fill_between(fc.year, fc.e0_lo80, fc.e0_hi80, alpha=0.30, color="C3")
    # WPP probabilistic benchmark
    p = w.load_series("e0Bproj_GMB.tsv", "e0")
    lo = w.load_series("e0Bproj95l_GMB.tsv", "e0"); hi = w.load_series("e0Bproj95u_GMB.tsv", "e0")
    m = p.merge(lo, on="year", suffixes=("", "_lo")).merge(hi, on="year", suffixes=("", "_hi"))
    m = m[m.year <= fc.year.max()]
    ax.plot(m.year, m.e0, color="C0", ls="--", label="WPP Bayesian (median)")
    ax.fill_between(m.year, m.e0_lo, m.e0_hi, alpha=0.12, color="C0", label="WPP 95% PI")
    ax.set_title(f"e0 forecast — Bayesian LC vs UN WPP ({sex}); max r-hat={rhat_max:.3f}")
    ax.set_xlabel("year"); ax.set_ylabel("e0 (years)"); ax.legend(fontsize=8)
    f.tight_layout(); f.savefig(FIG / "bayes_lc_vs_wpp_e0.png", dpi=130); plt.close(f)


if __name__ == "__main__":
    main()
