"""
Exploratory data analysis figures for the mortality data (Phase 3).

Produces the two views a mortality analyst always wants first:
  1. the log-mortality SURFACE (age x year heatmap), and
  2. age schedules log m(x) at selected years, showing the decline and the
     characteristic accident-hump / infant-mortality shape.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import figstyle  # noqa: F401  (applies dark, transparent figure style on import)
import wpp_data as w

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures"


def main():
    FIG.mkdir(exist_ok=True)
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    mx = w.mx_matrix("both", year_max=2023)
    ages = mx.index.to_numpy(); years = mx.columns.to_numpy()
    logm = np.log(mx.to_numpy(float))

    f, axs = plt.subplots(1, 2, figsize=(13, 4.8))

    im = axs[0].imshow(logm, aspect="auto", origin="lower", cmap="magma",
                       extent=[years.min(), years.max(), ages.min(), ages.max()])
    axs[0].set_title("log m(x,t) surface — The Gambia, 1950-2023")
    axs[0].set_xlabel("year"); axs[0].set_ylabel("age")
    f.colorbar(im, ax=axs[0], label="log mortality rate")

    for yr in (1950, 1980, 2000, 2023):
        axs[1].plot(ages, logm[:, list(years).index(yr)], label=str(yr))
    axs[1].set_title("Age schedule of log-mortality, selected years")
    axs[1].set_xlabel("age"); axs[1].set_ylabel("log m(x)"); axs[1].legend(title="year")

    f.tight_layout(); f.savefig(FIG / "eda_mortality_surface.png", dpi=130); plt.close(f)
    print("Saved -> figures/eda_mortality_surface.png")
    print(f"  log m(x,t): ages {ages.min()}-{ages.max()}, years {years.min()}-{years.max()}")
    print(f"  overall mortality decline: mean log m fell from "
          f"{logm[:, 0].mean():.2f} (1950) to {logm[:, -1].mean():.2f} (2023)")


if __name__ == "__main__":
    main()
