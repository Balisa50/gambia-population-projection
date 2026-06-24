"""
Validation summary (Phase 8): consolidates the quantitative checks and draws the
out-of-sample backtest calibration figure.

Checks reported:
  1. Internal consistency  — our life table reproduces WPP e0(2023).
  2. Out-of-sample backtest — fit <=2010, predict 2011-2023; error + PI coverage.
  3. Projection-engine validation — reproduced in src/projection.py (<1% vs WPP).
HDSS (Farafenni/Basse) and GBD cross-checks require digitising published life
tables and are listed as the next data tasks in docs/03-validation.md.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import figstyle  # noqa: F401  (applies dark, transparent figure style on import)
import wpp_data as w
from leecarter import e0_from_mx, backtest

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures"


def main():
    FIG.mkdir(exist_ok=True)
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

    mx = w.mx_matrix("both", year_max=2023)
    ages = mx.index.to_numpy(float)

    # 1. internal consistency
    e0_2023 = e0_from_mx(ages, mx[2023].to_numpy(float), "both")[0]
    wpp_e0 = w.load_series("e0B_GMB.tsv", "e0").set_index("year").loc[2023, "e0"]
    print(f"[1] internal: our e0(2023)={e0_2023:.2f} vs WPP {wpp_e0:.2f} "
          f"(|diff|={abs(e0_2023-wpp_e0):.3f})")

    # 2. backtest
    bt = backtest(mx, cutoff=2010, sex="both")
    print(f"[2] backtest fit<=2010: mean|e0 err|={bt.abs_err.mean():.2f}y, "
          f"95% PI coverage={bt.in_95PI.mean()*100:.0f}% over {len(bt)} years")

    f, ax = plt.subplots(figsize=(7.5, 4.5))
    ax.fill_between(bt.year, bt.e0_lo95, bt.e0_hi95, alpha=.25, color="C0", label="forecast 95% PI")
    ax.plot(bt.year, bt.e0_med, color="C0", label="forecast median")
    ax.plot(bt.year, bt.e0_actual, "o-", color=figstyle.FG, ms=4, label="actual (WPP)")
    ax.set_title("Backtest — fit ≤2010, predict 2011–2023 (The Gambia, e0)")
    ax.set_xlabel("year"); ax.set_ylabel("e0 (years)"); ax.legend()
    f.tight_layout(); f.savefig(FIG / "validation_backtest.png", dpi=130); plt.close(f)
    print("[fig] figures/validation_backtest.png")
    print("[3] engine validation: <1% vs WPP (see src/projection.py)")


if __name__ == "__main__":
    main()
