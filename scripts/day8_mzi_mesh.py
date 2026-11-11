"""Day 8: MZI behaviour and error propagation through rectangular MZI meshes.

Writes data/day8_mzi_mesh.npz + .json, then FROM THE FILE figures/fig_mzi_mesh.{png,pdf};
also draws figures/fig_mesh_schematic.{png,pdf} (a schematic of the N = 6 layout).
Run from repo root:  python scripts/day8_mzi_mesh.py   (~1 min)
"""
import json
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from topq.disorder import rng_for  # noqa: E402
from topq.mzi import (clements_layout, fidelity, intensity_error, mesh_unitary, mzi,  # noqa: E402
                      random_mesh_params)

DATA, META = ROOT / "data" / "day8_mzi_mesh.npz", ROOT / "data" / "day8_mzi_mesh.json"
SEED = 2026_11_03
NS = (4, 6, 8)
SIGMAS = np.geomspace(0.001, 0.3, 10)
N_SAMPLES = 1000
N_TARGETS = 5                     # different random target meshes per N (results averaged)


def compute():
    t0 = time.time()
    out = {"sigmas": SIGMAS, "Ns": np.array(NS)}
    th = np.linspace(0, 2 * np.pi, 361)
    out["theta"] = th
    for tag, eta in (("ideal", 0.5), ("eta045", 0.45), ("eta040", 0.40)):
        T = np.array([mzi(t, 0.0, eta, eta) for t in th])
        out[f"bar_{tag}"] = np.abs(T[:, 0, 0]) ** 2
        out[f"cross_{tag}"] = np.abs(T[:, 1, 0]) ** 2
    for N in NS:
        n = len(clements_layout(N))
        for err in ("phase", "split"):
            inf = np.zeros((N_TARGETS, SIGMAS.size, N_SAMPLES))
            ie = np.zeros_like(inf)
            for t in range(N_TARGETS):
                p = random_mesh_params(N, rng_for(SEED, 1000 * N + t))
                U0 = mesh_unitary(N, **p)
                for si, s in enumerate(SIGMAS):
                    for k in range(N_SAMPLES):
                        rng = rng_for(SEED + (1 if err == "phase" else 2), 10**6 * N + 10**4 * t + 10 * si + k * 997)
                        if err == "phase":
                            U = mesh_unitary(N, p["theta"] + s * rng.standard_normal(n),
                                             p["phi"] + s * rng.standard_normal(n),
                                             p["out"] + s * rng.standard_normal(N))
                        else:
                            e1 = np.clip(0.5 + s * rng.standard_normal(n), 0, 1)
                            e2 = np.clip(0.5 + s * rng.standard_normal(n), 0, 1)
                            U = mesh_unitary(N, **p, eta1=e1, eta2=e2)
                        inf[t, si, k] = 1 - fidelity(U0, U)
                        ie[t, si, k] = intensity_error(U0, U)
            out[f"inf_{err}_{N}"] = inf.mean(axis=(0, 2))
            out[f"inf_sem_{err}_{N}"] = inf.mean(axis=2).std(axis=0) / np.sqrt(N_TARGETS)   # spread over targets
            out[f"ie_{err}_{N}"] = ie.mean(axis=(0, 2))
    np.savez(DATA, **out)
    meta = {"script": "scripts/day8_mzi_mesh.py",
            "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "python": platform.python_version(), "numpy": np.__version__, "seed": SEED,
            "Ns": list(NS), "sigmas": SIGMAS.tolist(), "samples_per_point": N_SAMPLES,
            "targets_per_N": N_TARGETS,
            "first_order": {"phase": "1-F = (N-1) sigma_phi^2  (N^2 phase shifters)",
                            "split": "1-F = 2 (N-1) sigma_eta^2  (N(N-1) couplers)"},
            "ratio_MC_to_first_order_at_smallest_sigma": {
                f"{err}_N{N}": float(out[f"inf_{err}_{N}"][0] / ((1 if err == "phase" else 2) * (N - 1) * SIGMAS[0] ** 2))
                for err in ("phase", "split") for N in NS},
            "infidelity": {f"{err}_N{N}": out[f"inf_{err}_{N}"].tolist() for err in ("phase", "split") for N in NS},
            "intensity_error": {f"{err}_N{N}": out[f"ie_{err}_{N}"].tolist() for err in ("phase", "split") for N in NS},
            "max_cross_eta045": float(out["cross_eta045"].max()), "max_cross_eta040": float(out["cross_eta040"].max()),
            "runtime_s": round(time.time() - t0, 1)}
    META.write_text(json.dumps(meta, indent=1))
    return meta


def plot(meta):
    D = np.load(DATA)
    S = D["sigmas"]
    C3 = ["#2a78d6", "#eb6834", "#1baf7a"]
    plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False})
    fig, ax = plt.subplots(1, 3, figsize=(7.6, 2.9), constrained_layout=True)
    a = ax[0]
    th = D["theta"] / np.pi
    for tag, c, lab in (("ideal", "0.2", "ideal 50:50"), ("eta045", C3[1], "η = 0.45"), ("eta040", C3[2], "η = 0.40")):
        a.plot(th, D[f"cross_{tag}"], "-", color=c, lw=1.2, label=f"cross, {lab}")
        a.plot(th, D[f"bar_{tag}"], ":", color=c, lw=1.0)
    a.set_xlabel(r"internal phase $\theta/\pi$"); a.set_ylabel("output power")
    a.set_title("(a) one MZI (dotted: bar)", loc="left", fontsize=8.5)
    a.legend(fontsize=6, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=1)
    for j, (err, lab, k) in enumerate((("phase", r"phase error $\sigma_\phi$ (rad)", 1),
                                       ("split", r"coupler error $\sigma_\eta$", 2))):
        a = ax[j + 1]
        for i, N in enumerate(meta["Ns"]):
            a.loglog(S, D[f"inf_{err}_{N}"], "o", ms=3.3, color=C3[i], mfc="white", label=f"N={N}: simulated")
            a.loglog(S, k * (N - 1) * S**2, "-", color=C3[i], lw=0.9, label=f"N={N}: {k}(N−1)σ²")
        a.set_ylim(1e-7, 1.5)
        a.set_xlabel(lab); a.set_ylabel(r"infidelity $1-F$")
        a.set_title(f"({'bc'[j]}) mesh: {err} errors", loc="left", fontsize=8.5)
        a.legend(fontsize=5.6, frameon=False, loc="lower right", ncol=1)
    for ext in ("png", "pdf"):
        fig.savefig(ROOT / "figures" / f"fig_mzi_mesh.{ext}", dpi=300)
    plt.close(fig)


def plot_mesh_schematic(N: int = 6):
    lay = clements_layout(N)
    fig, ax = plt.subplots(figsize=(7.4, 2.8))
    for m in range(N):
        ax.plot([-0.6, N + 0.6], [m, m], color="0.55", lw=1.0)
        ax.text(-0.75, m, f"in {m + 1}", ha="right", va="center", fontsize=7)
        ax.text(N + 1.15, m, f"out {m + 1}", ha="left", va="center", fontsize=7)
        ax.add_patch(FancyBboxPatch((N + 0.15, m - 0.12), 0.3, 0.24, boxstyle="round,pad=0.02", fc="#f3d9b1", ec="0.4", lw=0.6))
    for c, m in lay:
        x = c + 0.1
        ax.add_patch(FancyBboxPatch((x, m - 0.18), 0.8, 1.36, boxstyle="round,pad=0.03", fc="#dce8f7", ec="#2a78d6", lw=0.8))
        ax.text(x + 0.4, m + 0.5, "MZI", ha="center", va="center", fontsize=6.5)
        ax.text(x + 0.4, m + 0.18, "θ, φ", ha="center", va="center", fontsize=5.5, color="0.3")
    ax.set_xlim(-1.6, N + 2.1); ax.set_ylim(-0.6, N - 0.3); ax.invert_yaxis(); ax.axis("off")
    ax.set_title(f"SCHEMATIC: rectangular (Clements-type) mesh, N = {N}: {len(lay)} MZIs + {N} output phase shifters "
                 "(drawn by us)", loc="left", fontsize=8)
    for ext in ("png", "pdf"):
        fig.savefig(ROOT / "figures" / f"fig_mesh_schematic.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    meta = compute()
    plot(meta)
    plot_mesh_schematic()
    print(json.dumps({k: meta[k] for k in ("ratio_MC_to_first_order_at_smallest_sigma", "max_cross_eta045",
                                          "max_cross_eta040", "runtime_s")}, indent=1))
