"""
Shared dark, transparent Matplotlib style so figures blend into a dark site
(the portfolio) instead of showing as white blocks. Importing this module
applies the style globally; expose FG for any explicit foreground colours.
"""
import matplotlib as mpl
from cycler import cycler

mpl.use("Agg")

FG = "#e6edf3"          # near-white foreground (text, key lines)
MUTED = "#9aa4b2"       # secondary
GRID = "#262b36"
EDGE = "#3a4150"
# accent cycle matching the site: cyan, pink, violet, amber, green
CYCLE = ["#22d3ee", "#f472b6", "#a78bfa", "#fbbf24", "#34d399"]

mpl.rcParams.update({
    "savefig.transparent": True,
    "figure.facecolor": "none",
    "axes.facecolor": "none",
    "savefig.facecolor": "none",
    "savefig.edgecolor": "none",
    "text.color": FG,
    "axes.labelcolor": FG,
    "axes.titlecolor": FG,
    "axes.edgecolor": EDGE,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "axes.grid": True,
    "grid.color": GRID,
    "grid.alpha": 0.5,
    "grid.linewidth": 0.6,
    "legend.frameon": False,
    "legend.labelcolor": FG,
    "axes.prop_cycle": cycler(color=CYCLE),
    "font.size": 11,
    "axes.titlesize": 12,
    "figure.dpi": 130,
})
