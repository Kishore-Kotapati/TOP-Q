"""Schematic (not data): how real fabrication errors map onto the Day-4 disorder models.
Writes figures/fig_error_schematic.{png,pdf}. Errors are drawn EXAGGERATED for visibility."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

ROOT = Path(__file__).resolve().parents[1]
W, H = 450, 220                      # nominal Si strip (nm), Blanco-Redondo et al. 2016
GAPS = [166, 294, 294, 166]          # long-long wall around the centre

def design_x():
    xs, x = [], 0.0
    for i in range(5):
        xs.append(x)
        x += W + (GAPS[i] if i < 4 else 0)
    return xs


def row(ax, actual, title, model, highlight, gap_mark=None):
    """actual: list of (x_left, width, height); dashed outlines show the design."""
    xd = design_x()
    right = xd[-1] + W
    ax.add_patch(Rectangle((-300, -150), right + 600, 150, color="#d9e6f2", lw=0))
    for i, (x0, w, h) in enumerate(actual):
        ax.add_patch(Rectangle((x0, 0), w, h, color="#c0392b" if i in highlight else "#555", lw=0))
    for x0 in xd:
        ax.add_patch(Rectangle((x0, 0), W, H, fill=False, ls="--", lw=0.7, ec="0.75"))
    if gap_mark is not None:
        x_a, x_b = gap_mark
        ax.annotate("", (x_a, H + 50), (x_b, H + 50), arrowprops=dict(arrowstyle="<->", lw=0.7, color="#c0392b"))
        ax.text((x_a + x_b) / 2, H + 80, "narrower gap", ha="center", fontsize=6.5, color="#c0392b")
    ax.set_xlim(-300, right + 300); ax.set_ylim(-170, 380); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(title, loc="left", fontsize=8.5)
    ax.text(right + 340, 110, model, fontsize=7.5, va="center")


xd = design_x()
fig, axs = plt.subplots(3, 1, figsize=(7.4, 4.6), constrained_layout=True)
g1 = [(x, W, H) for x in xd]
g1[0] = (xd[0] + 90, W, H)                     # outer guide placed too close: one gap shrinks by 90 nm
row(axs[0], g1, "(1) gap error: one gap narrower than designed",
    "→ that coupling κ grows\n   (log-normal coupling disorder;\n   keeps chiral symmetry)", [0],
    gap_mark=(xd[0] + 90 + W, xd[1]))
g2 = [(x, W, H) for x in xd]
g2[2] = (xd[2] - 45, W + 90, H)                # centre guide 90 nm wider, centre position fixed
row(axs[1], g2, "(2) width error: centre guide wider (centre position fixed by mask)",
    "→ its β grows AND both\n   neighbouring gaps shrink\n   (correlated width model;\n   breaks chiral symmetry)", [2])
row(axs[2], [(x, W, H + 40) for x in xd], "(3) wafer-scale thickness error: every guide thicker",
    "→ every β shifts together\n   (global shift; cancels in a\n   relative-phase measurement)", [0, 1, 2, 3, 4])
fig.suptitle("SCHEMATIC — errors exaggerated; dashed = design. Not measured data.", fontsize=8, x=0.02, ha="left", color="#a33")
for ext in ("png", "pdf"):
    fig.savefig(ROOT / "figures" / f"fig_error_schematic.{ext}", dpi=300)
