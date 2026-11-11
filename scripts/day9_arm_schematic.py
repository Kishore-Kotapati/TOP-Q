"""Schematic (not data): the Day 9-10 bridge architecture -- one MZI of the mesh with an
array element in each internal arm. Writes figures/fig_bridge_schematic.{png,pdf}."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

ROOT = Path(__file__).resolve().parents[1]
fig, ax = plt.subplots(figsize=(7.4, 2.9))
ax.set_xlim(0, 16); ax.set_ylim(-0.6, 5.4); ax.axis("off")
for y in (1.0, 4.0):
    ax.plot([0.3, 2.2], [y, y], color="0.4", lw=1.4)
    ax.plot([13.8, 15.7], [y, y], color="0.4", lw=1.4)
for x0 in (2.2, 12.2):
    ax.add_patch(FancyBboxPatch((x0, 0.7), 1.6, 3.6, boxstyle="round,pad=0.05", fc="#eee", ec="0.4"))
    ax.text(x0 + 0.8, 2.5, "coupler\n(σ_η)", ha="center", va="center", fontsize=7)
for y, lab in ((4.0, "top arm"), (1.0, "bottom arm")):
    ax.plot([3.8, 4.6], [y, y], color="0.4", lw=1.4); ax.plot([11.4, 12.2], [y, y], color="0.4", lw=1.4)
    ax.add_patch(FancyBboxPatch((4.6, y - 0.75), 6.8, 1.5, boxstyle="round,pad=0.04", fc="#dce8f7", ec="#2a78d6"))
    for k in range(-3, 4):
        ax.plot([4.9, 11.1], [y + 0.18 * k, y + 0.18 * k], color="#2a78d6" if k == 0 else "0.6",
                lw=1.6 if k == 0 else 0.6)
    ax.text(8.0, y + 0.95 if y > 2 else y - 1.15, f"{lab}: arm element T, D or W (own random errors) → t = |t| e^(iδφ)",
            ha="center", fontsize=7)
ax.text(3.75, 3.35, "θ heater", fontsize=6.5, color="#a33")
ax.text(0.3, 4.3, "in", fontsize=7); ax.text(15.2, 4.3, "out", fontsize=7)
ax.set_title("SCHEMATIC of one MZI in the Day 9-10 mesh (drawn by us; not built). "
             "Light enters and leaves each array through its centre guide.", loc="left", fontsize=7.8, color="#a33")
for ext in ("png", "pdf"):
    fig.savefig(ROOT / "figures" / f"fig_bridge_schematic.{ext}", dpi=300, bbox_inches="tight")
