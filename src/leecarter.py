"""
Classical Lee-Carter mortality model via SVD, with random-walk-with-drift
forecasting and Monte-Carlo uncertainty.

This is the transparent, well-understood baseline (Lee & Carter 1992) that the
Bayesian Poisson version (Phase 5) is benchmarked against. The model is

    log m(x,t) = a_x + b_x * k_t + eps(x,t),    sum_x b_x = 1,  sum_t k_t = 0

estimated by SVD of the centred log-rate matrix. The time index k_t is forecast
as a random walk with drift; uncertainty (innovation + drift) is propagated by
simulation into the death-rate surface and into life expectancy e0.

Run:  python src/leecarter.py   (fits Gambia both-sex rates, writes figures/)
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

import figstyle  # noqa: F401  (applies dark, transparent figure style on import)
import wpp_data as w

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures"


# --------------------------------------------------------------------------- #
# Vectorised life expectancy at birth from a death-rate matrix.
# --------------------------------------------------------------------------- #
def _a0_coale_demeny(m0: np.ndarray, sex: str) -> np.ndarray:
    """Vectorised Coale-Demeny (West) a(0). m0 is a 1-D array over scenarios."""
    m0 = np.asarray(m0, dtype=float)
    if sex.startswith("m"):
        return np.where(m0 < 0.107, 0.045 + 2.684 * m0, 0.330)
    if sex.startswith("f"):
        return np.where(m0 < 0.107, 0.053 + 2.800 * m0, 0.350)
    return 0.5 * (_a0_coale_demeny(m0, "m") + _a0_coale_demeny(m0, "f"))


def e0_from_mx(age: np.ndarray, mx: np.ndarray, sex: str = "both") -> np.ndarray:
    """Life expectancy at birth for each column of ``mx`` (ages x scenarios).

    Mirrors src.lifetable.period_life_table but fully vectorised over scenarios,
    so thousands of simulated mortality schedules cost one array computation.
    """
    age = np.asarray(age, dtype=float)
    mx = np.atleast_2d(np.asarray(mx, dtype=float))
    if mx.shape[0] != age.shape[0]:
        mx = mx.T  # accept (scenarios, ages) too
    A, S = mx.shape

    n = np.empty(A)
    n[:-1] = np.diff(age)
    n[-1] = np.inf
    n_col = n[:, None]

    ax = np.where(np.isinf(n_col), 0.0, n_col / 2.0) * np.ones((A, S))
    if age[0] == 0:
        ax[0, :] = _a0_coale_demeny(mx[0, :], sex)

    closed = ~np.isinf(n)
    qx = np.empty((A, S))
    qx[closed] = (n_col[closed] * mx[closed]) / (1.0 + (n_col[closed] - ax[closed]) * mx[closed])
    qx = np.clip(qx, 0.0, 1.0)
    qx[-1, :] = 1.0

    px = 1.0 - qx
    lx = np.empty((A, S))
    lx[0, :] = 1.0
    lx[1:, :] = np.cumprod(px[:-1, :], axis=0)
    dx = lx * qx

    Lx = np.empty((A, S))
    Lx[:-1, :] = n_col[:-1] * lx[1:, :] + ax[:-1, :] * dx[:-1, :]
    Lx[-1, :] = lx[-1, :] / mx[-1, :]
    e0 = Lx.sum(axis=0) / lx[0, :]
    return e0


# --------------------------------------------------------------------------- #
# Lee-Carter fit + forecast
# --------------------------------------------------------------------------- #
@dataclass
class LeeCarterFit:
    ages: np.ndarray
    years: np.ndarray
    ax: np.ndarray          # (A,)
    bx: np.ndarray          # (A,)  sum = 1
    kt: np.ndarray          # (T,)  sum = 0
    logm: np.ndarray        # (A, T) observed
    svd_explained: float    # variance share of 1st component

    def drift(self) -> float:
        return (self.kt[-1] - self.kt[0]) / (len(self.kt) - 1)

    def innovation_sd(self) -> float:
        d = self.drift()
        eps = np.diff(self.kt) - d
        return float(np.std(eps, ddof=1))


def fit_lee_carter(mx_wide: pd.DataFrame) -> LeeCarterFit:
    """Fit Lee-Carter by SVD. ``mx_wide`` is a DataFrame (index=age, cols=year)."""
    ages = mx_wide.index.to_numpy(dtype=float)
    years = mx_wide.columns.to_numpy(dtype=float)
    mx = mx_wide.to_numpy(dtype=float)
    if np.any(mx <= 0) or not np.all(np.isfinite(mx)):
        raise ValueError("m(x,t) must be finite and strictly positive for log.")

    logm = np.log(mx)
    ax = logm.mean(axis=1)
    centred = logm - ax[:, None]

    U, S, Vt = np.linalg.svd(centred, full_matrices=False)
    explained = float(S[0] ** 2 / np.sum(S ** 2))

    b_raw = U[:, 0]
    k_raw = S[0] * Vt[0, :]
    s = b_raw.sum()
    bx = b_raw / s
    kt = k_raw * s

    # Sign convention: with sum(bx)=1>0, mortality decline => kt declines.
    # If the first component came out inverted, flip both so bx sums to +1.
    if bx.sum() < 0:
        bx, kt = -bx, -kt

    return LeeCarterFit(ages=ages, years=years, ax=ax, bx=bx, kt=kt,
                        logm=logm, svd_explained=explained)


def forecast_kt(fit: LeeCarterFit, horizon: int, n_sims: int = 4000,
                seed: int = 42) -> tuple[np.ndarray, np.ndarray]:
    """Simulate future k_t as a random walk with drift.

    Returns (future_years, kt_paths) with kt_paths shape (n_sims, horizon).
    Uncertainty includes both the innovation variance and drift estimation
    uncertainty (drift se = sigma / sqrt(T-1)).
    """
    rng = np.random.default_rng(seed)
    T = len(fit.kt)
    d_hat = fit.drift()
    sigma = fit.innovation_sd()
    drift_se = sigma / np.sqrt(T - 1)

    last_year = int(fit.years[-1])
    future_years = np.arange(last_year + 1, last_year + 1 + horizon)

    drifts = rng.normal(d_hat, drift_se, size=n_sims)
    innov = rng.normal(0.0, sigma, size=(n_sims, horizon))
    steps = drifts[:, None] + innov
    kt_paths = fit.kt[-1] + np.cumsum(steps, axis=1)
    return future_years, kt_paths


def forecast_e0(fit: LeeCarterFit, horizon: int, sex: str = "both",
                n_sims: int = 4000, seed: int = 42) -> pd.DataFrame:
    """Forecast life expectancy at birth with prediction intervals."""
    future_years, kt_paths = forecast_kt(fit, horizon, n_sims, seed)
    # log m(x,t) = a_x + b_x k_t  ->  rates (ages x [sims*years]) then e0
    recs = []
    for j, yr in enumerate(future_years):
        kt = kt_paths[:, j]                       # (n_sims,)
        logm = fit.ax[:, None] + np.outer(fit.bx, kt)   # (A, n_sims)
        e0 = e0_from_mx(fit.ages, np.exp(logm), sex)    # (n_sims,)
        recs.append((yr, np.median(e0),
                     np.percentile(e0, 2.5), np.percentile(e0, 97.5),
                     np.percentile(e0, 10), np.percentile(e0, 90)))
    return pd.DataFrame(recs, columns=["year", "e0_med", "e0_lo95", "e0_hi95",
                                       "e0_lo80", "e0_hi80"])


def backtest(mx_wide: pd.DataFrame, cutoff: int, sex: str = "both") -> pd.DataFrame:
    """Fit on years <= cutoff, forecast the rest, compare e0 to actual."""
    train = mx_wide.loc[:, mx_wide.columns <= cutoff]
    test_years = [y for y in mx_wide.columns if y > cutoff]
    fit = fit_lee_carter(train)
    fc = forecast_e0(fit, horizon=len(test_years), sex=sex)
    actual = {int(y): e0_from_mx(mx_wide.index.to_numpy(float),
                                 mx_wide[y].to_numpy(float), sex)[0]
              for y in test_years}
    fc["e0_actual"] = fc["year"].map(actual)
    fc["in_95PI"] = (fc["e0_actual"] >= fc["e0_lo95"]) & (fc["e0_actual"] <= fc["e0_hi95"])
    fc["abs_err"] = (fc["e0_med"] - fc["e0_actual"]).abs()
    return fc


def _demo() -> None:
    FIG.mkdir(exist_ok=True)
    sex = "both"
    mx = w.mx_matrix(sex=sex, year_max=2023)
    fit = fit_lee_carter(mx)

    print("=== Lee-Carter (SVD) - The Gambia, both sexes, 1950-2023 ===")
    print(f"  1st singular component explains {fit.svd_explained*100:.1f}% of variance")
    print(f"  k_t drift = {fit.drift():.4f}/yr, innovation sd = {fit.innovation_sd():.4f}")
    print(f"  k_t: {fit.kt[0]:.2f} (1950) -> {fit.kt[-1]:.2f} (2023)")

    e0_hist = e0_from_mx(fit.ages, mx[2023].to_numpy(float), sex)[0]
    print(f"  fitted-data e0(2023) = {e0_hist:.2f}  (WPP estimate 65.86)")

    fc = forecast_e0(fit, horizon=2074 - 2023, sex=sex)
    for yr in (2050, 2074):
        r = fc[fc.year == yr].iloc[0]
        print(f"  LC e0({yr}) = {r.e0_med:.2f}  [95% PI {r.e0_lo95:.2f}, {r.e0_hi95:.2f}]")
    print("  (WPP probabilistic benchmark: 2050=70.2 [62.4,77.3], 2074=73.2 [63.5,83.4])")

    print("\n=== Backtest: fit<=2010, predict 2011-2023 ===")
    bt = backtest(mx, cutoff=2010, sex=sex)
    print(f"  mean |e0 error| = {bt.abs_err.mean():.2f} yrs | "
          f"95% PI coverage = {bt.in_95PI.mean()*100:.0f}% of {len(bt)} years")
    last = bt.iloc[-1]
    print(f"  2023: predicted {last.e0_med:.2f} [{last.e0_lo95:.2f}, {last.e0_hi95:.2f}]"
          f" vs actual {last.e0_actual:.2f}")

    _make_figures(fit, fc)
    print(f"\nFigures written to {FIG}")


def _make_figures(fit: LeeCarterFit, fc: pd.DataFrame) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # a_x, b_x, k_t
    f, axs = plt.subplots(1, 3, figsize=(13, 3.6))
    axs[0].plot(fit.ages, fit.ax); axs[0].set_title("$a_x$ (avg log-mortality)"); axs[0].set_xlabel("age")
    axs[1].plot(fit.ages, fit.bx); axs[1].set_title("$b_x$ (sensitivity)"); axs[1].set_xlabel("age")
    axs[2].plot(fit.years, fit.kt); axs[2].set_title("$k_t$ (mortality index)"); axs[2].set_xlabel("year")
    f.suptitle("Lee-Carter parameters - The Gambia (both sexes, 1950-2023)")
    f.tight_layout(); f.savefig(FIG / "lc_parameters.png", dpi=130); plt.close(f)

    # e0 fan chart
    e0_hist = np.array([e0_from_mx(fit.ages, np.exp(fit.ax + fit.bx * k), "both")[0]
                        for k in fit.kt])
    f, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(fit.years, e0_hist, color=figstyle.FG, label="fitted (1950-2023)")
    ax.plot(fc.year, fc.e0_med, color="C0", label="LC forecast (median)")
    ax.fill_between(fc.year, fc.e0_lo95, fc.e0_hi95, alpha=0.2, color="C0", label="95% PI")
    ax.fill_between(fc.year, fc.e0_lo80, fc.e0_hi80, alpha=0.3, color="C0")
    ax.set_title("Life expectancy at birth - The Gambia (Lee-Carter)")
    ax.set_xlabel("year"); ax.set_ylabel("e0 (years)"); ax.legend()
    f.tight_layout(); f.savefig(FIG / "lc_e0_forecast.png", dpi=130); plt.close(f)


if __name__ == "__main__":
    _demo()
