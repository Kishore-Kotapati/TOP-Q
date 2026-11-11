"""Day 2: winding number, bulk-boundary correspondence, edge-mode localisation,
edge-mode hybridisation. Computes -> data/day2_topology.npz (+ .json), then plots
figures/fig_topology.{png,pdf} FROM THE SAVED FILE.  Deterministic (no randomness).

Run from the repo root:  python scripts/day2_topology.py
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

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from topq.ssh import ssh_couplings, ssh_hamiltonian  # noqa: E402
from topq.topology import (  # noqa: E402
    even_chain_splitting, fit_localization_length, winding_number,
    xi_analytic, zero_mode_from_couplings,
)

DATA = ROOT / "data" / "day2_topology.npz"
META = ROOT / "data" / "day2_topology.json"
FIG = ROOT / "figures" / "fig_topology"

K2 = 1.0                                   # everything in units of k2
R_LOOPS = (0.5, 2.0)                       # k1/k2 for the d-vector loops
R_PROFILE = (1 / 3, 0.5, 0.8)              # 1/3 ~ Blanco-Redondo 2016 hopping ratio (0.31-0.34)
M_PROFILE = 41                             # odd chain: exactly one zero mode
N_BB = 100                                 # cells for bulk-boundary count
R_SPLIT = (0.5, 0.7)
N_SPLIT = np.arange(2, 31)


def compute():
    out = {}
    # (a) d-vector loops
    k = np.linspace(0, 2 * np.pi, 400)
    for i, r in enumerate(R_LOOPS):
        out[f"loop_dx_{i}"] = r * K2 + K2 * np.cos(k)
        out[f"loop_dy_{i}"] = K2 * np.sin(k)
    # (b) bulk-boundary: W vs in-gap modes of an even open chain (per edge)
    r_grid = np.r_[np.linspace(0.05, 0.95, 37), np.linspace(1.05, 3.0, 40)]
    W = np.array([winding_number(r * K2, K2) for r in r_grid])
    n_edge = []
    for r in r_grid:
        e = np.linalg.eigvalsh(ssh_hamiltonian(2 * N_BB, r * K2, K2))
        n_edge.append(np.sum(np.abs(e) < 0.5 * abs(r * K2 - K2)) / 2)
    out.update(r_grid=r_grid, W=W, n_edge=np.array(n_edge))
    # (c) edge profiles: numerical eigenvector (eigh) and analytic exp decay
    xi_table = {}
    for i, r in enumerate(R_PROFILE):
        H = ssh_hamiltonian(M_PROFILE, r * K2, K2)
        e, v = np.linalg.eigh(H)
        j = np.argmin(np.abs(e))
        psi = v[:, j] * np.sign(v[0, j])
        out[f"prof_b_{i}"] = e[j]
        out[f"prof_power_{i}"] = np.abs(psi) ** 2
        xi_a = xi_analytic(r * K2, K2)
        # fit on the NUMERICAL eigenvector (eigh) of a long chain: independent of the formula
        e201, v201 = np.linalg.eigh(ssh_hamiltonian(201, r * K2, K2))
        xi_f = fit_localization_length(v201[:, np.argmin(np.abs(e201))], floor=1e-10)
        xi_table[f"{r:.4f}"] = {"xi_analytic_cells": xi_a, "xi_fit_cells": xi_f,
                                "rel_err": xi_f / xi_a - 1, "b_numeric": float(e[j])}
    # (d) even-chain splitting vs N
    for i, r in enumerate(R_SPLIT):
        num = [np.sort(np.abs(np.linalg.eigvalsh(ssh_hamiltonian(2 * N, r * K2, K2))))[0] for N in N_SPLIT]
        out[f"split_num_{i}"] = np.array(num)
        out[f"split_th_{i}"] = np.array([even_chain_splitting(r * K2, K2, N) for N in N_SPLIT])
    out["N_split"] = N_SPLIT
    np.savez(DATA, **out)
    meta = {
        "script": "scripts/day2_topology.py",
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "python": platform.python_version(), "numpy": np.__version__,
        "matplotlib": matplotlib.__version__, "units": "all b, kappa in units of kappa2",
        "bulk_boundary_mismatches": int(np.sum(np.abs(np.array(n_edge) - np.round(W)) > 0)),
        "xi": xi_table, "M_profile": M_PROFILE, "N_bulk_boundary": N_BB,
        "randomness": "none (deterministic)",
    }
    META.write_text(json.dumps(meta, indent=2))
    return meta


def plot():
    D = np.load(DATA)
    C = ["#2a78d6", "#eb6834", "#1baf7a"]           # blue, orange, aqua (fixed order)
    plt.rcParams.update({"font.size": 8.5, "axes.spines.top": False, "axes.spines.right": False})
    fig, ax = plt.subplots(2, 2, figsize=(7.2, 6.2), constrained_layout=True)

    a = ax[0, 0]
    for i, (r, ls) in enumerate(zip(R_LOOPS, ["-", "--"])):
        W = winding_number(r, 1.0)
        a.plot(D[f"loop_dx_{i}"], D[f"loop_dy_{i}"], ls, color=C[i], lw=1.6,
               label=fr"$\kappa_1/\kappa_2={r:g}$:  $W={int(round(W)) + 0}$")
        a.annotate("", xy=(D[f"loop_dx_{i}"][60], D[f"loop_dy_{i}"][60]),
                   xytext=(D[f"loop_dx_{i}"][50], D[f"loop_dy_{i}"][50]),
                   arrowprops=dict(arrowstyle="-|>", color=C[i], lw=1.2))
    a.plot(0, 0, "k+", ms=10, mew=1.5)
    a.annotate("origin", (0, 0), xytext=(0.15, -0.35), fontsize=7)
    a.set_aspect("equal")
    a.set_xlabel(r"$d_x/\kappa_2$"); a.set_ylabel(r"$d_y/\kappa_2$")
    a.set_title(r"(a) $\mathbf{d}(k)$ as $k$ goes $0\to2\pi$", loc="left", fontsize=9)
    a.legend(fontsize=7, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.26), ncol=1)

    b = ax[0, 1]
    b.plot(D["r_grid"], np.round(D["W"]), "-", color="0.35", lw=1.2, drawstyle="steps-mid",
           label="winding number $W$ (bulk)")
    b.plot(D["r_grid"], D["n_edge"], "o", ms=3.5, mfc="white", mec=C[0], mew=1.0,
           label=f"in-gap modes per edge (open chain, {2*N_BB} guides)")
    b.axvline(1, color="0.8", lw=0.8, ls=":")
    b.set_xlabel(r"$\kappa_1/\kappa_2$"); b.set_ylabel("count")
    b.set_yticks([0, 1])
    b.set_title("(b) bulk-boundary correspondence", loc="left", fontsize=9)
    b.legend(fontsize=7, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.18))

    c = ax[1, 0]
    for i, r in enumerate(R_PROFILE):
        p = D[f"prof_power_{i}"]
        n_wg = np.arange(1, p.size + 1)
        A = slice(0, None, 2)
        xi = xi_analytic(r, 1.0)
        cells = np.arange(p[A].size)
        c.semilogy(n_wg[A], p[A], "o", ms=3.2, mfc="white", mec=C[i], mew=1.0)
        c.semilogy(n_wg[A], p[0] * np.exp(-2 * cells / xi), "-", color=C[i], lw=1.0,
                   label=fr"$\kappa_1/\kappa_2={r:.2f}$, $\xi={xi:.2f}$ cells")
    c.set_ylim(1e-16, 2)
    c.set_xlabel("waveguide number (A guides shown; B guides carry 0 power)")
    c.set_ylabel(r"power $|a_n|^2$")
    c.set_title(f"(c) zero-mode profile, odd chain M={M_PROFILE}", loc="left", fontsize=9)
    c.legend(fontsize=7, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.18),
             title="circles: numerical eigenvector; lines: exp(-2n/ξ)", title_fontsize=6.5)

    d = ax[1, 1]
    N = D["N_split"]
    for i, r in enumerate(R_SPLIT):
        d.semilogy(N, D[f"split_num_{i}"], "o", ms=3.2, mfc="white", mec=C[i], mew=1.0)
        d.semilogy(N, D[f"split_th_{i}"], "-", color=C[i], lw=1.0,
                   label=fr"$\kappa_1/\kappa_2={r}$: $\kappa_2(1-r^2)r^N$")
    d.set_xlabel("number of unit cells N (even chain, 2N guides)")
    d.set_ylabel(r"$|b|/\kappa_2$ of the two edge modes")
    d.set_title("(d) even chain: edge modes are NOT at exactly 0", loc="left", fontsize=9)
    d.legend(fontsize=7, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.18),
             title="circles: numerical; lines: first-order formula", title_fontsize=6.5)

    for ext in ("png", "pdf"):
        fig.savefig(f"{FIG}.{ext}", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    meta = compute()
    plot()
    print(json.dumps({k: meta[k] for k in ("bulk_boundary_mismatches", "xi")}, indent=2))
    print("saved data/day2_topology.npz/.json, figures/fig_topology.png/.pdf")
