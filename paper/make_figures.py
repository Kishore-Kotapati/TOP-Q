"""Figures for the TOP-Q term paper, drawn from the stored simulation data (data/*.npz).
Plain scientific style: black/grey lines, serif font, white background.
Usage: python3 make_figures.py <path-to-topq-data-dir> <output-dir>
"""
import sys, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch

DATA, OUT = sys.argv[1], sys.argv[2]
plt.rcParams.update({
    "font.family": "serif", "mathtext.fontset": "cm", "font.size": 8,
    "axes.labelsize": 8, "axes.titlesize": 8, "legend.fontsize": 6.5,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.linewidth": 0.6,
    "lines.linewidth": 1.0, "lines.markersize": 3.5, "legend.frameon": False,
    "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
    "axes.grid": True, "grid.color": "0.88", "grid.linewidth": 0.4, "axes.axisbelow": True,
})
W = 7.1
K = "black"; G = "0.55"


def load(n):
    return np.load(f"{DATA}/{n}.npz")


def tag(ax, s, *args, **kw):
    """Panel letter placed outside the plot frame, above its top-left corner."""
    if ax.axison:
        ax.set_title(s, loc="left", fontsize=8, fontweight="bold", pad=3)
    else:
        ax.text(0.0, 1.0, s, transform=ax.transAxes, fontsize=8, fontweight="bold", va="bottom", ha="left")


# ----------------------------------------------------------------------------- Fig 1
def fig_schematic():
    fig, axs = plt.subplots(1, 3, figsize=(W, 1.66), gridspec_kw={"width_ratios": [1.1, 1.15, 1.3]})
    for a in axs:
        a.set_axis_off()

    # (a) SSH array cross-section
    a = axs[0]; a.set_xlim(0, 10); a.set_ylim(-0.3, 5.2)
    w, x, xs = 0.62, 0.5, []
    gaps = [1.0, 0.45, 1.0, 0.45, 1.0]
    for i in range(6):
        a.add_patch(Rectangle((x, 1.5), w, 0.9, fc="0.75", ec=K, lw=0.6))
        a.text(x + w / 2, 1.05, "A" if i % 2 == 0 else "B", ha="center", fontsize=7)
        xs.append(x)
        if i < 5:
            lab = r"$\kappa_1$" if i % 2 == 0 else r"$\kappa_2$"
            a.annotate("", xy=(x + w + gaps[i], 2.75), xytext=(x + w, 2.75),
                       arrowprops=dict(arrowstyle="<->", lw=0.5))
            a.text(x + w + gaps[i] / 2, 2.95, lab, ha="center", fontsize=7)
            x += w + gaps[i]
    a.plot([xs[0] - 0.1, xs[1] + w + 0.1], [0.75, 0.75], color=K, lw=0.6)
    a.text((xs[0] + xs[1] + w) / 2, 0.2, "one unit cell", ha="center", fontsize=6.5)
    a.text(5.3, 4.1, "wide gap: weak $\\kappa_1$\nnarrow gap: strong $\\kappa_2$", ha="center", fontsize=6.5)
    tag(a, "(a)", 0.0, 0.98)

    # (b) arrays T and D around the centre
    a = axs[1]; a.set_xlim(0, 10); a.set_ylim(-0.3, 5.2)

    def row(y, gaps, widths, label):
        x = 0.9
        a.text(0.35, y + 0.38, label, ha="center", va="center", fontsize=7.5)
        for i, wd in enumerate(widths):
            centre = (i == len(widths) // 2)
            a.add_patch(Rectangle((x, y), wd, 0.75, fc="0.35" if centre else "0.75", ec=K, lw=0.6))
            x += wd + (gaps[i] if i < len(gaps) else 0)
        pass
    row(2.9, [0.35, 0.85, 0.35, 0.85, 0.85, 0.35, 0.85, 0.35], [0.45] * 9, "T")
    row(1.0, [0.6] * 8, [0.45, 0.45, 0.45, 0.45, 0.62, 0.45, 0.45, 0.45, 0.45], "D")
    a.text(5.2, 4.25, "T: weak coupling on both\nsides of the centre guide", ha="center", fontsize=6.5)
    a.text(5.0, 0.1, "D: equal gaps; centre guide detuned by $\\Delta$", ha="center", fontsize=6.5)
    tag(a, "(b)", 0.0, 0.98)

    # (c) one MZI with an arm element in each arm
    a = axs[2]; a.set_xlim(0, 10); a.set_ylim(-0.3, 5.2)
    yu, yl = 3.0, 1.2
    for y in (yu, yl):
        a.plot([0.2, 1.3], [y, y], color=K, lw=0.9)
        a.plot([2.3, 7.7], [y, y], color=K, lw=0.9)
        a.plot([8.7, 9.8], [y, y], color=K, lw=0.9)
    for x0 in (1.3, 7.7):
        a.add_patch(Rectangle((x0, yl - 0.25), 1.0, yu - yl + 0.5, fc="white", ec=K, lw=0.6, zorder=3))
        a.text(x0 + 0.5, (yu + yl) / 2, "50:50", ha="center", va="center", fontsize=6, rotation=90, zorder=4)
    a.add_patch(Rectangle((2.7, yu - 0.28), 0.9, 0.56, fc="white", ec=K, lw=0.6, zorder=3)); a.text(3.15, yu, r"$\theta$", ha="center", va="center", fontsize=7, zorder=4)
    a.add_patch(Rectangle((0.35, yu - 0.28), 0.7, 0.56, fc="white", ec=K, lw=0.6, zorder=3)); a.text(0.7, yu, r"$\phi$", ha="center", va="center", fontsize=7, zorder=4)
    for y in (yu, yl):
        a.add_patch(Rectangle((4.0, y - 0.45), 3.3, 0.9, fc="0.8", ec=K, lw=0.6, zorder=3))
        a.text(5.65, y, "arm element", ha="center", va="center", fontsize=5.6, zorder=4)
    a.text(5.65, 4.25, "arm element = T, D or W,\neach with its own random errors", ha="center", fontsize=6.5)
    a.text(1.8, 0.25, "coupler", ha="center", fontsize=6.3); a.text(8.2, 0.25, "coupler", ha="center", fontsize=6.3); a.text(3.15, 2.2, "phase\nshifter", ha="center", va="center", fontsize=6)
    tag(a, "(c)", 0.0, 0.98)
    fig.subplots_adjust(wspace=0.05)
    fig.savefig(f"{OUT}/fig1_schematic.pdf"); plt.close(fig)


# ----------------------------------------------------------------------------- Fig 2
def fig_theory():
    D1, D2 = load("day1_bands"), load("day2_topology")
    fig, axs = plt.subplots(1, 4, figsize=(W, 1.68))
    a = axs[0]
    k = D1["k"]
    a.plot(k, D1["band_plus_2"], color=K, label=r"$\kappa_1/\kappa_2=1/2$ or $2$")
    a.plot(k, D1["band_minus_2"], color=K)
    a.plot(k, D1["band_plus_1"], color=G, ls="--", label=r"$\kappa_1=\kappa_2$")
    a.plot(k, D1["band_minus_1"], color=G, ls="--")
    kr = D1["ring_kabs_2"]; e = D1["ring_eigs_2"]
    kk = np.concatenate([kr, -kr]); ee = np.concatenate([e, e])
    for sgn in (1, -1):
        sel = np.sign(ee) == sgn; o = np.argsort(kk[sel])
        a.plot(kk[sel][o], ee[sel][o], color="0.35", ls="-.", lw=0.8, label="ring, 48 guides" if sgn > 0 else None)
    a.set_xticks([-np.pi, 0, np.pi]); a.set_xticklabels([r"$-\pi$", "0", r"$\pi$"])
    a.set_xlabel(r"Bloch phase $k$ (rad)"); a.set_ylabel(r"$b/\bar\kappa$")
    a.set_ylim(-2.3, 4.6); a.legend(loc="upper center", handlelength=1.6, borderpad=0.1, labelspacing=0.2)
    tag(a, "(a)", 0.03, 0.12)

    a = axs[1]
    n = np.arange(1, 41)
    a.plot(n, D1["open_eigs_2"], color=K, label=r"$\kappa_1<\kappa_2$")
    a.plot(n, D1["open_eigs_0"], color=G, ls="--", label=r"$\kappa_1>\kappa_2$")
    a.axhspan(-2 / 3, 2 / 3, color="0.92", lw=0)
    a.set_xlabel("mode number"); a.set_ylabel(r"$b/\bar\kappa$")
    a.legend(loc="upper left", handletextpad=0.1, borderpad=0.1); tag(a, "(b)", 0.80, 0.14)

    a = axs[2]
    a.plot(D2["loop_dx_0"], D2["loop_dy_0"], color=K, label=r"$\kappa_1<\kappa_2$: $W=1$")
    a.plot(D2["loop_dx_1"], D2["loop_dy_1"], color=G, ls="--", label=r"$\kappa_1>\kappa_2$: $W=0$")
    a.axhline(0, color="0.8", lw=0.4, zorder=0); a.axvline(0, color="0.8", lw=0.4, zorder=0)
    a.set_xlim(-1.5, 3.5)
    a.set_xlabel(r"$d_x/\kappa_2$"); a.set_ylabel(r"$d_y/\kappa_2$")
    a.set_ylim(-1.5, 3.5); a.set_aspect("equal"); a.legend(loc="upper center", handlelength=1.6, borderpad=0.1, labelspacing=0.2)
    tag(a, "(c)", 0.03, 0.12)

    a = axs[3]
    rs = [1 / 3, 0.5, 0.8]; mk = ["o", "s", "^"]; lab = ["1/3", "1/2", "0.8"]
    for i, r in enumerate(rs):
        p = D2[f"prof_power_{i}"][::2]
        cell = np.arange(p.size)
        m = p > 1e-15
        a.semilogy(cell, (1 - r * r) * r ** (2 * cell), color=G, lw=2.2, alpha=0.6)
        a.semilogy(cell[m], p[m], color=K, ls=["-", "--", "-."][i], lw=0.9, label=rf"$r={lab[i]}$")
    a.set_ylim(1e-15, 5); a.set_xlim(-0.5, 20.5)
    a.set_xlabel("cell number $n-1$"); a.set_ylabel("power on A guide")
    a.legend(loc="lower left", handletextpad=0.1, borderpad=0.1, labelspacing=0.2); tag(a, "(d)", 0.80, 0.93)
    fig.subplots_adjust(wspace=0.62)
    fig.savefig(f"{OUT}/fig2_theory.pdf"); plt.close(fig)


# ----------------------------------------------------------------------------- Fig 3
def fig_propagation():
    D3, D5 = load("day3_design"), load("day5_propagation")
    fig, axs = plt.subplots(1, 4, figsize=(W, 1.68), gridspec_kw={"width_ratios": [1, 1, 1, 1.05]})
    a = axs[0]
    M = D3["psiT"].size; n = np.arange(M) - M // 2
    pT, pD = D3["psiT"] ** 2, D3["psiD"] ** 2
    m = pT > 1e-14
    a.semilogy(n[m], pT[m], color=K, label="T (interface)")
    a.semilogy(n, pD, color=G, ls="--", label="D (defect)")
    a.set_xlim(-10.5, 10.5); a.set_ylim(1e-11, 3e4)
    a.set_xlabel("guide number from centre"); a.set_ylabel("mode power")
    a.legend(loc="upper center", ncol=1, handletextpad=0.1, borderpad=0.1, labelspacing=0.15); tag(a, "(a)", 0.03, 0.12)

    zs = D5["zs"]; Lc = np.pi / 2
    for a, key, t, lb in ((axs[1], "map_T_centre", "array T", "(b)"), (axs[2], "map_uniform_centre", "uniform array", "(c)")):
        P = D5[key]; Mg = P.shape[1]
        im = a.imshow(P.T ** 0.5, origin="lower", aspect="auto", cmap="gray_r", vmin=0, vmax=1,
                      extent=[0, zs[-1] / Lc, -(Mg // 2) - 0.5, Mg // 2 + 0.5], interpolation="nearest")
        a.grid(False); a.set_ylim(-22, 22); a.set_xlabel("$z$ (coupling lengths)"); a.set_ylabel("guide number")
        a.text(0.97, 0.93, t, transform=a.transAxes, ha="right", va="top", fontsize=7,
               bbox=dict(fc="white", ec="none", pad=1)); tag(a, lb, 0.03, 0.14, bbox=dict(fc="white", ec="none", pad=1))
    a = axs[3]
    a.plot(zs / Lc, D5["self_re_T_clean"], color=K, label="T, clean")
    a.plot(zs / Lc, D5["self_re_T_coupling06"], color=K, ls=":", label=r"T, $\sigma_\kappa=0.6$")
    a.plot(zs / Lc, D5["self_re_uniform"], color=G, ls="--", label="uniform")
    a.axhline(0, color="0.6", lw=0.4)
    a.axhline(0.6, color="0.6", lw=0.4, ls="-."); a.text(9.9, 0.42, r"bound $2w_0-1$", ha="right", fontsize=6, color="0.3")
    a.set_ylim(-0.6, 2.3); a.set_xlabel("$z$ (coupling lengths)"); a.set_ylabel(r"centre field $a_c(z)$")
    a.legend(loc="upper right", ncol=1, handlelength=1.8, borderpad=0.1, labelspacing=0.1)
    tag(a, "(d)", 0.03, 0.12)
    fig.subplots_adjust(wspace=0.6)
    fig.savefig(f"{OUT}/fig3_propagation.pdf", dpi=300); plt.close(fig)


# ----------------------------------------------------------------------------- Fig 4
def fig_montecarlo():
    D4, D6 = load("day4_mc_eigen"), load("day6_mc_phase")
    s4 = D4["sigmas"]; s6 = D6["sigmas"]; iL = 9; L = float(D6["Ls"][iL])
    fig, axs = plt.subplots(1, 4, figsize=(W, 1.68))
    a = axs[0]
    a.loglog(s4, D4["stats_D_coupling"][:, 0], color=K, ls="--", label="D")
    a.loglog(s4, 0.506 * s4, color=G, lw=0.7, label=r"$0.51\,\sigma_\kappa$")
    a.loglog(s4, np.maximum(D4["stats_T_coupling"][:, 0], 1e-17), color=K, label="T (round-off)")
    a.set_ylim(1e-17, 30); a.set_xlabel(r"gap-error strength $\sigma_\kappa$"); a.set_ylabel(r"std of $\delta b/\kappa_2$")
    a.legend(loc="center right", handletextpad=0.2, borderpad=0.1, labelspacing=0.2); tag(a, "(a)")

    a = axs[1]
    a.loglog(s4, D4["stats_D_onsite"][:, 0], color=K, ls="--", lw=1.6, label="D")
    a.loglog(s4, D4["stats_T_onsite"][:, 0], color=K, lw=0.8, label="T")
    a.loglog(s4, 0.81 * s4, color=G, lw=0.7, label=r"$0.81\,\sigma_\beta$")
    a.set_xlabel(r"on-site error $\sigma_\beta/\kappa_2$"); a.set_ylabel(r"std of $\delta b/\kappa_2$")
    a.legend(loc="lower right", handletextpad=0.2, borderpad=0.1, labelspacing=0.2); tag(a, "(b)")

    def col(arr, kind, q):
        return np.array([D6[f"{arr}_{kind}_{i}_{q}"][iL] for i in range(s6.size)])
    a = axs[2]
    a.loglog(s6, col("D", "coupling", "modal_cstd"), color=K, label="D, gap errors")
    a.loglog(s6, col("D", "onsite", "modal_cstd"), color=K, ls="--", lw=1.6, label="D, on-site")
    a.loglog(s6, col("T", "onsite", "modal_cstd"), color=K, ls=":", lw=1.0, label="T, on-site")
    a.loglog(s6, np.sqrt(2) * 0.506 * s6 * L, color=G, lw=0.7)
    a.loglog(s6, np.sqrt(2) * 0.81 * s6 * L, color=G, lw=0.7, label="first order")
    a.set_ylim(2e-2, 400); a.set_xlabel(r"error strength $\sigma$"); a.set_ylabel(r"arm-to-arm phase std (rad)")
    a.legend(loc="upper left", handletextpad=0.2, borderpad=0.1, labelspacing=0.15); tag(a, "(c)", 0.82, 0.14)

    a = axs[3]
    a.semilogx(s6, col("T", "coupling", "guide_R"), color=K, label="T, gap")
    a.semilogx(s6, col("D", "coupling", "guide_R"), color=K, ls="--", label="D, gap")
    a.semilogx(s6, col("T", "onsite", "guide_R"), color=G, lw=1.6, label="T, on-site")
    a.semilogx(s6, col("D", "onsite", "guide_R"), color=G, ls=":", lw=1.3, label="D, on-site")
    a.set_ylim(-0.05, 1.08); a.set_xlabel(r"error strength $\sigma$"); a.set_ylabel(r"fringe contrast $R$")
    a.legend(loc="lower left", handletextpad=0.2, borderpad=0.1, labelspacing=0.15); tag(a, "(d)", 0.03, 0.55)
    fig.subplots_adjust(wspace=0.66)
    fig.savefig(f"{OUT}/fig4_montecarlo.pdf"); plt.close(fig)


# ----------------------------------------------------------------------------- Fig 5
def fig_circuits():
    D8, D9 = load("day8_mzi_mesh"), load("day9_bridge")
    coef = json.load(open(f"{DATA}/day9_bridge.json"))["arm_coefficients_per_unit_sigma"]
    fig, axs = plt.subplots(1, 4, figsize=(W, 1.68))
    a = axs[0]
    th = D8["theta"]
    a.semilogy(th, np.maximum(D8["bar_ideal"], 1e-6), color=K, label=r"$\eta=0.50$")
    a.semilogy(th, D8["bar_eta045"], color=K, ls="--", label=r"$\eta=0.45$")
    a.semilogy(th, D8["bar_eta040"], color=G, ls=":", lw=1.3, label=r"$\eta=0.40$")
    a.set_xticks([0, np.pi, 2 * np.pi]); a.set_xticklabels(["0", r"$\pi$", r"$2\pi$"])
    a.set_ylim(1e-4, 60); a.set_xlabel(r"internal phase $\theta$ (rad)"); a.set_ylabel("bar-port power")
    a.legend(loc="upper center", ncol=1, borderpad=0.1, labelspacing=0.15, handlelength=1.6); tag(a, "(a)", 0.03, 0.12)

    a = axs[1]
    s = D8["sigmas"]; mk = {4: "v", 6: "o", 8: "s"}
    for N in (4, 6, 8):
        a.loglog(s, D8[f"inf_phase_{N}"], color=K, ls={4: "-", 6: "--", 8: ":"}[N], label=f"$N={N}$")
        a.loglog(s, (N - 1) * s ** 2, color=G, lw=0.7)
    a.set_ylim(1e-6, 3); a.set_xlabel(r"phase error $\sigma_\phi$ (rad)"); a.set_ylabel(r"mesh infidelity $1-F$")
    a.legend(loc="upper left", handletextpad=0.2, borderpad=0.1, labelspacing=0.15); tag(a, "(b)", 0.82, 0.14)

    a = axs[2]
    s9 = D9["sigmas"]

    def m(key):
        return np.array([D9[f"{key}_{i}_inf"].mean() for i in range(s9.size)])
    a.loglog(s9, m("T_coupling_cpl"), color=K, label="T, gap")
    a.loglog(s9, m("D_coupling_cpl"), color=K, ls="--", label="D, gap")
    a.loglog(s9, m("W_coupling_cpl"), color=K, ls=":", lw=1.3, label="W, gap")
    a.loglog(s9, m("T_onsite_arms"), color=G, lw=1.6, label="T, on-site")
    a.loglog(s9, m("W_onsite_arms"), color=G, ls=":", lw=1.3, label="W, on-site")
    a.axhline(0.01, color="0.6", lw=0.4, ls="-.")
    a.set_ylim(3e-6, 3e8); a.set_xlabel(r"error strength $\sigma$"); a.set_ylabel(r"mesh infidelity $1-F$")
    a.legend(loc="upper left", ncol=1, handletextpad=0.3, borderpad=0.1, labelspacing=0.05, handlelength=1.6, fontsize=5.6)
    tag(a, "(c)", 0.82, 0.14)

    a = axs[3]
    N = 6
    sb = np.logspace(-4, -1.5, 240); sk = np.logspace(-3.3, -0.5, 240)
    SB, SK = np.meshgrid(sb, sk)
    for el, ls, c in (("T", "-", K), ("D", "--", K), ("W", ":", K)):
        cc = coef[el]
        var = (cc["onsite"]["phase_per_sigma"] * SB) ** 2 + (cc["onsite"]["logamp_per_sigma"] * SB) ** 2 \
            + (cc["coupling"]["phase_per_sigma"] * SK) ** 2 + (cc["coupling"]["logamp_per_sigma"] * SK) ** 2
        Z = (N - 1) ** 2 / N * var + 2 * (N - 1) * (np.pi / 4 * SK) ** 2
        a.contour(SB, SK, Z, levels=[0.01], colors=[c], linestyles=[ls], linewidths=1.0)
        a.plot([], [], color=c, ls=ls, label=el)
    a.set_xscale("log"); a.set_yscale("log")
    a.set_xlabel(r"on-site error $\sigma_\beta/\kappa_2$"); a.set_ylabel(r"gap error $\sigma_\kappa$")
    a.set_ylim(3e-4, 3); a.text(1.3e-4, 0.2, r"$F\geq0.99$" + "\nbelow each line", fontsize=6)
    a.legend(loc="lower left", borderpad=0.1, labelspacing=0.1, handlelength=1.8); tag(a, "(d)", 0.82, 0.93)
    fig.subplots_adjust(wspace=0.66)
    fig.savefig(f"{OUT}/fig5_circuits.pdf"); plt.close(fig)


# ----------------------------------------------------------------------------- Fig 6
def fig_hom():
    D = load("day11_hom")
    fig, axs = plt.subplots(1, 3, figsize=(W, 1.66), gridspec_kw={"width_ratios": [1, 1, 1.35]})
    a = axs[0]
    a.plot(D["eta"], D["V_eta"], color=K, label=r"$2RT/(R^2+T^2)$")
    a.plot(D["eta_fock"], D["V_fock"], color=G, ls="--", lw=1.4, label="number state check")
    a.set_xlabel(r"coupler bar fraction $\eta$"); a.set_ylabel("HOM visibility $V$")
    a.set_ylim(0, 1.3); a.legend(loc="upper center", borderpad=0.1, labelspacing=0.15); tag(a, "(a)", 0.03, 0.12)

    a = axs[1]
    t = D["tau"]
    a.plot(t, D["Pc_pure"], color=K, label="identical, pure")
    a.plot(t, D["Pc_detuned"], color=K, ls="--", label="colour mismatch")
    a.plot(t, D["Pc_jitter"], color=G, ls=":", lw=1.3, label="timing jitter")
    a.set_ylim(0, 0.72); a.set_xlabel(r"delay $\tau$ (units of $1/\sigma_\omega$)"); a.set_ylabel(r"coincidence probability $P_c$")
    a.legend(loc="upper center", ncol=1, borderpad=0.1, labelspacing=0.15, handlelength=1.8); tag(a, "(b)", 0.03, 0.12)

    a = axs[2]
    data, labs = [], []
    for kind, i, name in (("coupling", 1, "gap\n0.03"), ("onsite", 0, "on-site\n0.01"), ("onsite", 1, "on-site\n0.03")):
        for el in "TDW":
            data.append(np.clip(1 - D[f"V_{el}_{kind}_{i}"], 1e-6, None)); labs.append(el)
    pos = [1, 2, 3, 4.6, 5.6, 6.6, 8.2, 9.2, 10.2]
    a.boxplot(data, positions=pos, widths=0.6, whis=(5, 95), showfliers=False,
              medianprops=dict(color=K, lw=0.9), boxprops=dict(lw=0.6), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6))
    a.set_yscale("log"); a.set_ylim(1e-5, 100)
    a.set_xticks(pos); a.set_xticklabels(labs)
    for x, s in ((2, "gap\n0.03"), (5.6, "on-site\n0.01"), (9.2, "on-site\n0.03")):
        a.text(x, 4, s, ha="center", fontsize=6.3)
    a.text(3, 3e-5, "W: 0", ha="center", fontsize=6.3); a.set_ylabel(r"$1-V$ of MZI splitter"); tag(a, "(c)", 0.03, 0.12)
    fig.subplots_adjust(wspace=0.42)
    fig.savefig(f"{OUT}/fig6_hom.pdf"); plt.close(fig)


if __name__ == "__main__":
    fig_schematic(); fig_theory(); fig_propagation(); fig_montecarlo(); fig_circuits(); fig_hom()
    print("ok")
