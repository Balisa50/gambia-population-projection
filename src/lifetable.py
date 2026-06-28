"""
Period life-table construction from age-specific central death rates.

This is the deterministic demographic core of the project. Every mortality
forecast - classical SVD Lee-Carter or Bayesian - ultimately produces a matrix
of central death rates m(x, t); this module turns each such schedule into a full
period life table, from which life expectancy e(x) and the survivorship needed
by the cohort-component projection are read off.

Implemented from first principles (no external demography library) so the
methodology is fully transparent and auditable, per the research protocol.

Conventions
-----------
- ``age``  : left edge of each age interval (e.g. 0, 1, 2, ... for single year).
- ``mx``   : central death rate  m(x) = deaths / person-years in the interval.
- The final interval is treated as open-ended (x and above).
- ``ax``   : average person-years lived in [x, x+n) by those who die there.
             Defaults to n/2 for closed intervals; age 0 uses the Coale-Demeny
             (West) rule, which depends on the infant mortality level and sex.

References
----------
Preston, Heuveline & Guillot (2001), *Demography: Measuring and Modeling
Population Processes*, Ch. 3. Coale & Demeny (1983) regression rules for a(0).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

__all__ = ["a0_coale_demeny", "period_life_table"]


def a0_coale_demeny(m0: float, sex: str = "both") -> float:
    """Average years lived by infants who die before age 1 (Coale-Demeny West).

    Infant deaths cluster in the first days/weeks of life, so a(0) is well below
    0.5. These are the standard Coale-Demeny regression rules used by the UN and
    in Preston et al. (2001), Table 3.3.
    """
    sex = sex.lower()
    if sex.startswith("m"):
        return 0.045 + 2.684 * m0 if m0 < 0.107 else 0.330
    if sex.startswith("f"):
        return 0.053 + 2.800 * m0 if m0 < 0.107 else 0.350
    # both sexes: average of the male and female rules
    return 0.5 * (a0_coale_demeny(m0, "m") + a0_coale_demeny(m0, "f"))


def period_life_table(
    age,
    mx,
    sex: str = "both",
    radix: float = 100_000.0,
    ax=None,
) -> pd.DataFrame:
    """Build a period life table from central death rates.

    Parameters
    ----------
    age : array-like
        Left edges of the age intervals, strictly increasing. The last entry is
        the open age interval (x+). Single-year (0,1,2,...) or abridged
        (0,1,5,10,...) spacing both work; interval widths are inferred.
    mx : array-like
        Central death rate for each interval. ``mx`` for the open interval must
        be > 0 (it sets the closeout: L = l / m).
    sex : {"both","male","female"}
        Controls the a(0) rule only.
    radix : float
        l(0), the starting cohort size. Cosmetic; e(x) is invariant to it.
    ax : array-like, optional
        Override the average years lived in each interval by decedents. If None,
        defaults are used (n/2 for closed intervals, Coale-Demeny for age 0).

    Returns
    -------
    pandas.DataFrame
        Columns: age, n, mx, ax, qx, lx, dx, Lx, Tx, ex.
    """
    age = np.asarray(age, dtype=float)
    mx = np.asarray(mx, dtype=float)
    if age.ndim != 1 or age.shape != mx.shape:
        raise ValueError("age and mx must be 1-D arrays of the same length.")
    if np.any(np.diff(age) <= 0):
        raise ValueError("age must be strictly increasing.")
    if np.any(mx < 0) or not np.all(np.isfinite(mx)):
        raise ValueError("mx must be finite and non-negative.")
    if mx[-1] <= 0:
        raise ValueError("mx for the open age interval must be > 0 (closeout).")

    n_ages = len(age)
    # Interval widths: difference between successive left edges; open interval
    # width is conventionally infinite (handled separately below).
    n = np.empty(n_ages)
    n[:-1] = np.diff(age)
    n[-1] = np.inf

    # a(x): average years lived in interval by those who die in it.
    if ax is None:
        ax = np.where(np.isinf(n), 0.0, n / 2.0)  # n/2 for closed intervals
        if age[0] == 0:
            ax[0] = a0_coale_demeny(mx[0], sex)
        # Abridged tables: the 1-4 interval is also front-loaded. Apply a simple
        # Coale-Demeny-style adjustment if such an interval is present.
        if n_ages > 1 and age[0] == 0 and np.isclose(age[1], 1) and not np.isinf(n[1]) and np.isclose(n[1], 4):
            if sex.startswith("m"):
                ax[1] = 1.651 - 2.816 * mx[0] if mx[0] < 0.107 else 1.352
            elif sex.startswith("f"):
                ax[1] = 1.522 - 1.518 * mx[0] if mx[0] < 0.107 else 1.361
            else:
                ax[1] = 0.5 * (
                    (1.651 - 2.816 * mx[0] if mx[0] < 0.107 else 1.352)
                    + (1.522 - 1.518 * mx[0] if mx[0] < 0.107 else 1.361)
                )
    else:
        ax = np.asarray(ax, dtype=float).copy()

    # qx: probability of dying in the interval.  Closed: q = n*m/(1+(n-a)m).
    qx = np.empty(n_ages)
    closed = ~np.isinf(n)
    qx[closed] = (n[closed] * mx[closed]) / (1.0 + (n[closed] - ax[closed]) * mx[closed])
    qx = np.clip(qx, 0.0, 1.0)
    qx[-1] = 1.0  # everyone in the open interval eventually dies

    # lx, dx
    lx = np.empty(n_ages)
    lx[0] = radix
    for i in range(1, n_ages):
        lx[i] = lx[i - 1] * (1.0 - qx[i - 1])
    dx = lx * qx
    # ensure d sums exactly to radix (closes the table)
    dx[-1] = lx[-1]

    # Lx: person-years lived in each interval.
    Lx = np.empty(n_ages)
    Lx[:-1] = n[:-1] * lx[1:] + ax[:-1] * dx[:-1]
    Lx[-1] = lx[-1] / mx[-1]  # open interval closeout

    # Tx, ex
    Tx = np.cumsum(Lx[::-1])[::-1]
    ex = Tx / lx

    return pd.DataFrame(
        {
            "age": age.astype(int) if np.all(age == age.astype(int)) else age,
            "n": n,
            "mx": mx,
            "ax": ax,
            "qx": qx,
            "lx": lx,
            "dx": dx,
            "Lx": Lx,
            "Tx": Tx,
            "ex": ex,
        }
    )


def _self_test() -> None:
    """Sanity checks against known demographic regularities.

    We don't have a Gambian gold-standard table yet, so we verify internal
    consistency and plausibility on a synthetic Gompertz-like schedule plus
    elevated infant mortality (typical of a high-mortality setting):
      * e(x) decreases then the table closes correctly (T = sum of L),
      * d sums to the radix,
      * a flat hazard m gives e0 ~= 1/m (exponential-survival check).
    """
    # --- Check 1: constant hazard => exponential survival, e0 ~ 1/m ---
    ages = np.arange(0, 111)
    m_const = np.full_like(ages, 0.02, dtype=float)
    lt = period_life_table(ages, m_const, sex="both")
    # For a constant hazard with a(x)=0.5, e0 is close to 1/m - 0.5; check order.
    assert abs(lt["ex"].iloc[0] - (1 / 0.02 - 0.5)) < 1.0, lt["ex"].iloc[0]

    # --- Check 2: d sums to radix, T0 == sum(L) ---
    assert abs(lt["dx"].sum() - 100_000) < 1e-6
    assert abs(lt["Tx"].iloc[0] - lt["Lx"].sum()) < 1e-6

    # --- Check 3: a high-mortality schedule yields a plausible e0 ---
    # Brass-style: high infant mortality, low adult hazard rising with age.
    m = 0.00003 * np.exp(0.095 * ages)      # Gompertz adult component
    m[0] = 0.06                              # ~60/1000 infant mortality
    m[1:5] += 0.004                          # elevated early childhood
    lt2 = period_life_table(ages, m, sex="both")
    e0 = lt2["ex"].iloc[0]
    assert 45 < e0 < 75, f"implausible e0={e0:.2f}"
    assert (lt2["qx"] >= 0).all() and (lt2["qx"] <= 1).all()
    assert (lt2["lx"].diff().dropna() <= 1e-9).all(), "lx must be non-increasing"

    print("lifetable self-test passed.")
    print(f"  constant-hazard e0 = {lt['ex'].iloc[0]:.2f} (expect ~{1/0.02 - 0.5:.1f})")
    print(f"  high-mortality   e0 = {e0:.2f}, infant q0 = {lt2['qx'].iloc[0]:.4f}")


if __name__ == "__main__":
    _self_test()
