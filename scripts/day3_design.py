"""Day 3: comparison design (T vs D), first-order disorder predictions, and a
dimensioned hardware SCHEMATIC of the Blanco-Redondo et al. (2016) silicon arrays.

Writes data/day3_design.npz (+ .json) and, FROM THAT FILE,
figures/fig_design.{png,pdf}  and  figures/fig_hardware_schematic.{png,pdf}.
Deterministic. Run from repo root:  python scripts/day3_design.py
"""
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from topq.arrays import (  # noqa: E402
    centre, defect_array, matched_defect_detuning, short_short_interface,
    topological_interface, wall_couplings,
)
from topq.disorder import first_order_std, first_order_std_width  # noqa: E402
from topq.topology import edge_mode, ipr  # noqa: E402

DATA, META = ROOT / "data" / "day3_design.npz", ROOT / "data" / "day3_design.json"
K1, K2, M = 1 / 3, 1.0, 41            # units of k2; r = 1/3 (decision D2)
KBAR = 0.5 * (K1 + K2)
# Blanco-Redondo et al., PRL 116, 163901 (2016) -- values as reported (nm)
HW = {"width": 450, "height": 220, "gap_narrow_LL": 166, "gap_wide_LL": 294,
      "gap_narrow_SS": 182, "gap_wide_SS": 324, "gap_uniform": 230,
      "length_um": 500, "n_guides": 203, "wavelength_nm": 1550}


def compute():
    delta = matched_defect_detuning(K1, K2)
    HT = topological_interface(M, K1, K2)
    HD = defect_array(M, KBAR, delta)
    HS = short_short_interface(M, K1, K2)
    eT, vT = np.linalg.eigh(HT); iT = int(np.argmin(np.abs(eT)))
    eD, vD = np.linalg.eigh(HD); iD = int(np.argmax(eD))
    eS = np.linalg.eigvalsh(HS)
    psiT, psiD = vT[:, iT], vD[:, iD]
    cT, cD = wall_couplings(M, K1, K2), np.full(M - 1, KBAR)
    pred = {}
    for name, psi, c in (("T", psiT, cT), ("D", psiD, cD)):
        pred[name] = {"coupling": first_order_std(psi, c, sigma_kappa=1.0),
                      "onsite": first_order_std(psi, c, sigma_beta=1.0),
                      "width": first_order_std_width(psi, c, 1.0, 1.0)}
    np.savez(DATA, eT=eT, eD=eD, eS=eS, psiT=psiT, psiD=psiD, delta=delta,
             ipr_T=ipr(psiT), ipr_D=ipr(psiD), bT=eT[iT], bD=eD[iD])
    meta = {"script": "scripts/day3_design.py",
            "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "python": platform.python_version(), "numpy": np.__version__,
            "params": {"k1": K1, "k2": K2, "M": M, "kbar": KBAR, "Delta_matched": delta},
            "mode": {"b_T": float(eT[iT]), "b_D": float(eD[iD]),
                     "IPR_T": ipr(psiT), "IPR_D": ipr(psiD)},
            "first_order_std_per_unit_sigma": pred,
            "note_width_model": "sigma_beta = sigma_kappa = 1 per unit width error (ratio unknown; Day 4 scans it)",
            "hardware_parameters_source": "Blanco-Redondo et al., PRL 116, 163901 (2016)",
            "hardware_parameters": HW, "randomness": "none"}
    META.write_text(json.dumps(meta, indent=2))
    return meta


def plot_design(meta):
    D = np.load(DATA)
    C = {"T": "#2a78d6", "D": "#eb6834"}
    plt.rcParams.update({"font.size": 8.5, "axes.spines.top": False, "axes.spines.right": False})
    fig, ax = plt.subplots(1, 3, figsize=(7.4, 2.9), constrained_layout=True)
    n = np.arange(M) - centre(M)

    a = ax[0]
    a.semilogy(n, D["psiT"] ** 2 + 1e-300, "o", ms=3.2, mfc="white", mec=C["T"], mew=1, label="T: topological interface")
    a.semilogy(n, D["psiD"] ** 2, "s", ms=2.8, color=C["D"], label=r"D: detuned defect ($\Delta$ matched)")
    a.set_ylim(1e-12, 2); a.set_xlim(-12, 12)
    a.set_xlabel("guide position from centre"); a.set_ylabel(r"mode power $|\psi_n|^2$")
    a.set_title(f"(a) same IPR = {float(D['ipr_T']):.3f}", loc="left", fontsize=9)
    a.legend(fontsize=6.5, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.2))

    b = ax[1]
    idx = np.arange(1, M + 1)
    b.plot(idx, D["eT"], "o", ms=2.8, mfc="white", mec=C["T"], mew=0.9, label="T spectrum")
    b.plot(idx, D["eD"], "s", ms=2.4, color=C["D"], label="D spectrum")
    b.annotate("T mode, b = 0", (M // 2 + 1, 0), xytext=(2, 1.3), fontsize=7,
               arrowprops=dict(arrowstyle="->", lw=0.6))
    b.annotate(f"D mode, b = {float(D['bD']):.2f}", (M, float(D["bD"])), xytext=(12, 2.05), fontsize=7,
               arrowprops=dict(arrowstyle="->", lw=0.6))
    b.set_xlabel("eigenvalue index"); b.set_ylabel(r"$b/\kappa_2$")
    b.set_title("(b) where the two modes sit", loc="left", fontsize=9)
    b.legend(fontsize=6.5, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.2))

    c = ax[2]
    p = meta["first_order_std_per_unit_sigma"]
    kinds = ["coupling", "onsite", "width"]
    x = np.arange(3)
    c.bar(x - 0.19, [p["T"][k] for k in kinds], 0.36, color=C["T"], label="T")
    c.bar(x + 0.19, [p["D"][k] for k in kinds], 0.36, color=C["D"], label="D")
    for xi, k in zip(x, kinds):
        c.text(xi - 0.19, p["T"][k] + 0.02, f"{p['T'][k]:.2f}", ha="center", fontsize=6.5)
        c.text(xi + 0.19, p["D"][k] + 0.02, f"{p['D'][k]:.2f}", ha="center", fontsize=6.5)
    c.set_xticks(x, ["coupling\n(per $\\sigma_\\kappa$)", "on-site\n(per $\\sigma_\\beta/\\kappa_2$)", "width model\n($\\sigma_\\beta=\\sigma_\\kappa$)"])
    c.set_ylabel(r"predicted std$(\delta b)/\kappa_2$ per unit $\sigma$")
    c.set_title("(c) first-order PREDICTION", loc="left", fontsize=9)
    c.legend(fontsize=6.5, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.32), ncol=2)
    for ext in ("png", "pdf"):
        fig.savefig(ROOT / "figures" / f"fig_design.{ext}", dpi=300)
    plt.close(fig)


def plot_hardware():
    """Dimensioned schematic. Cross-sections to scale (nm); top view NOT to scale."""
    w, h = HW["width"], HW["height"]
    fig = plt.figure(figsize=(7.4, 5.6), constrained_layout=True)
    gs = fig.add_gridspec(3, 1, height_ratios=[1, 1, 1.25])

    def cross_section(axx, gaps_from_centre, title, centre_width=None, note=None):
        n_side = len(gaps_from_centre)
        widths = [w] * (2 * n_side + 1)
        if centre_width:
            widths[n_side] = centre_width
        xs = [0.0]
        order = list(reversed(gaps_from_centre)) + list(gaps_from_centre)
        for j, g in enumerate(order):
            xs.append(xs[-1] + widths[j] + g)
        xs = np.array(xs); xs -= xs[n_side] + widths[n_side] / 2
        axx.add_patch(Rectangle((xs[0] - 600, -400), xs[-1] - xs[0] + w + 1200, 400, color="#d9e6f2", lw=0))
        axx.text(xs[0] - 560, -330, "buried oxide (SiO$_2$) / Si substrate (SOI)", fontsize=6.5, color="#345")
        axx.text(xs[0] - 560, h + 60, "air cladding", fontsize=6.5, color="#555")
        for j, x0 in enumerate(xs):
            is_c = j == n_side
            axx.add_patch(Rectangle((x0, 0), widths[j], h, color="#555" if not is_c else "#2a78d6", lw=0))
        for j, g in enumerate(order[n_side - 2:n_side + 2]):
            jj = n_side - 2 + j
            x_left = xs[jj] + widths[jj]
            axx.annotate("", (x_left, h + 25), (x_left + g, h + 25),
                         arrowprops=dict(arrowstyle="<->", lw=0.6))
            axx.text(x_left + g / 2, h + 45, f"{g}", ha="center", fontsize=6.5)
        axx.annotate("", (xs[n_side] - 10, -40), (xs[n_side] + widths[n_side] + 10, -40),
                     arrowprops=dict(arrowstyle="<->", lw=0.6))
        axx.text(xs[n_side] + widths[n_side] / 2, -95, f"{widths[n_side]}" if not centre_width else "wider (Δ)",
                 ha="center", fontsize=6.5)
        axx.text(xs[-1] + w + 60, h / 2, f"{h} nm", fontsize=6.5, va="center")
        axx.set_xlim(xs[0] - 600, xs[-1] + w + 600); axx.set_ylim(-560 if note else -420, h + 160)
        axx.set_aspect("equal"); axx.axis("off")
        axx.set_title(title, loc="left", fontsize=8.5)
        if note:
            axx.text(xs[0] - 560, -470, note, fontsize=6.3, color="#a33", va="top")

    g_LL = [HW["gap_wide_LL"], HW["gap_narrow_LL"], HW["gap_wide_LL"], HW["gap_narrow_LL"]]
    cross_section(fig.add_subplot(gs[0]), g_LL,
                  "(a) T: 'long-long' topological interface, cross-section to scale (dimensions in nm, from Blanco-Redondo et al. 2016)")
    cross_section(fig.add_subplot(gs[1]), [HW["gap_uniform"]] * 4,
                  "(b) D: conventional defect — uniform 230 nm gaps, wider centre guide (our trivial reference)",
                  centre_width=560, note="centre width for the matched Δ NOT computed yet (needs a mode solver) — drawn illustratively")

    t = fig.add_subplot(gs[2])
    ys = np.arange(-5, 6) * 1.0
    for y in ys:
        t.plot([0, 10], [y, y], color="#2a78d6" if y == 0 else "#666", lw=2.2 if y == 0 else 1.2)
    t.add_patch(FancyArrowPatch((-1.8, 0), (-0.1, 0), arrowstyle="-|>", mutation_scale=12, color="#c22"))
    t.text(-1.9, 0.45, "1550 nm laser\ninto centre guide", fontsize=6.5, ha="left", color="#c22")
    t.add_patch(FancyArrowPatch((10.1, 0), (11.6, 0), arrowstyle="-|>", mutation_scale=12, color="#c22"))
    t.text(10.2, 0.45, "output facet → NIR camera /\nfibre + power meter", fontsize=6.5, color="#c22")
    t.annotate("", (0, -6.2), (10, -6.2), arrowprops=dict(arrowstyle="<->", lw=0.7))
    t.text(5, -7.0, "propagation length z = 500 μm (not to scale)", ha="center", fontsize=6.8)
    t.text(5, 5.6, "array of 203 guides; domain wall at the centre (only 11 drawn)", ha="center", fontsize=6.8)
    t.set_xlim(-2.2, 14); t.set_ylim(-7.6, 6.4); t.axis("off")
    t.set_title("(c) top view of the measurement idea in [6] — SCHEMATIC, not our hardware", loc="left", fontsize=8.5)
    for ext in ("png", "pdf"):
        fig.savefig(ROOT / "figures" / f"fig_hardware_schematic.{ext}", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    meta = compute()
    plot_design(meta)
    plot_hardware()
    print(json.dumps({k: meta[k] for k in ("params", "mode", "first_order_std_per_unit_sigma")}, indent=2))
