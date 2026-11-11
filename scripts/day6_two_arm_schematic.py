"""Schematic (not data): the 'two identical arrays on one chip' thought experiment of Day 6.
Writes figures/fig_two_arm_schematic.{png,pdf}."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ROOT = Path(__file__).resolve().parents[1]
fig, ax = plt.subplots(figsize=(7.4, 2.6))
ax.set_xlim(0, 15); ax.set_ylim(0, 5); ax.axis("off")


def box(x, y, w, h, text, fc):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.05", fc=fc, ec="0.3", lw=0.8))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=7.5)


def arrow(x0, y0, x1, y1, c="#c22"):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=10, color=c, lw=1.2))


ax.text(0.1, 2.5, "1550 nm\nlaser", fontsize=7.5, va="center", color="#c22")
arrow(1.2, 2.5, 2.0, 2.5)
box(2.0, 2.0, 1.6, 1.0, "50:50\nsplitter", "#eee")
arrow(3.6, 2.8, 4.6, 3.8); arrow(3.6, 2.2, 4.6, 1.2)
for y, lab in ((3.4, "arm A: array (T or D)\nlocal errors set A"), (0.8, "arm B: same design\nlocal errors set B")):
    box(4.6, y, 4.6, 0.9, lab, "#dce8f7")
ax.text(6.9, 4.55, "same chip → common (wafer-scale) error cancels in φA − φB", ha="center", fontsize=7, color="0.3")
arrow(9.2, 3.85, 10.3, 2.8); arrow(9.2, 1.25, 10.3, 2.2)
box(10.3, 2.0, 1.6, 1.0, "50:50\ncombiner", "#eee")
arrow(11.9, 2.8, 12.8, 3.4, "#555"); arrow(11.9, 2.2, 12.8, 1.6, "#555")
ax.text(12.9, 3.5, "detector 1", fontsize=7, va="center"); ax.text(12.9, 1.5, "detector 2", fontsize=7, va="center")
ax.text(7.5, 0.15, "output power depends on Δφ = φA − φB  →  Day 6 computes the statistics of Δφ",
        ha="center", fontsize=7.5)
ax.set_title("SCHEMATIC of the Day-6 numerical experiment (not built, not measured)", loc="left", fontsize=8.5, color="#a33")
for ext in ("png", "pdf"):
    fig.savefig(ROOT / "figures" / f"fig_two_arm_schematic.{ext}", dpi=300, bbox_inches="tight")
