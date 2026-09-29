"""Bad-cell rate versus refinement depth (reads badrate.csv written by
badrate_sweep.sh).  One colour + marker per curve (fixed categorical order,
palette validated for CVD separation); solid = tangency threshold + repair
rounds, dashed = baseline.  Output: badrate.pdf"""
import csv
from collections import defaultdict
import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm",
                     "font.size": 9, "axes.linewidth": 0.6})

curves = [("conic", "conic ($g=0$)", "#2a78d6", "o"),
          ("elliptic", "elliptic ($g=1$)", "#eb6834", "s"),
          ("fermat_cubic", "Fermat cubic ($g=1$)", "#1baf7a", "^"),
          ("fermat_quartic", "Fermat quartic ($g=3$)", "#eda100", "D")]

data = defaultdict(list)
for r in csv.DictReader(open("badrate.csv")):
    ok, bad = int(r["ok"]), int(r["bad"])
    data[(r["curve"], r["mode"])].append((int(r["depth"]), 100.0 * bad / (ok + bad)))

fig, ax = plt.subplots(figsize=(5.2, 3.3))
for name, label, color, marker in curves:
    for mode, ls, lw, fill in (("baseline", (0, (4, 2)), 1.0, "white"), ("combined", "-", 2.0, color)):
        pts = sorted(data[(name, mode)])
        if not pts:
            continue
        x, y = zip(*pts)
        ax.plot(x, y, ls=ls, lw=lw, color=color, marker=marker, ms=5,
                mfc=fill, mec=color, mew=1.2, zorder=3)
    pts = sorted(data[(name, "combined")])
    if pts:  # direct label at the last point of the solid line
        dy = {"conic": 4, "elliptic": -4}.get(name, 0)
        ax.annotate(label, pts[-1], xytext=(7, dy), textcoords="offset points",
                    va="center", fontsize=8, color="0.15")

ax.set_yscale("log")
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter
ax.yaxis.set_major_locator(FixedLocator([0.2, 0.5, 1, 2, 5, 10, 20, 50]))
ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
ax.yaxis.set_minor_formatter(NullFormatter())
ax.set_ylim(0.3, 45)
ax.set_xlabel("refinement depth")
ax.set_ylabel("bad cells (% of touched cells)")
depths = sorted({d for v in data.values() for d, _ in v})
ax.set_xticks(depths)
ax.set_xlim(min(depths) - 0.5, max(depths) + 4.5)
ax.grid(True, which="major", color="0.88", lw=0.5)
ax.grid(True, which="minor", axis="y", color="0.94", lw=0.4)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
handles = [Line2D([], [], color="0.3", ls=(0, (4, 2)), lw=1.0, marker="o", mfc="white", mec="0.3", label="baseline"),
           Line2D([], [], color="0.3", ls="-", lw=2.0, marker="o", mfc="0.3", mec="0.3",
                  label="tangency threshold + repair")]
ax.legend(handles=handles, loc="upper right", frameon=False, fontsize=8)
fig.tight_layout()
fig.savefig("badrate.pdf")
print("wrote badrate.pdf")
