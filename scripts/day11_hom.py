"""Day 11: Hong-Ou-Mandel extension (bounded scope: two photons, one 2x2 device).

Writes data/day11_hom.npz + .json and FROM THEM figures/fig_hom.{png,pdf}.
Run from repo root:  python scripts/day11_hom.py   (~1 min)
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
from topq.bridge import arm_pool  # noqa: E402
from topq.disorder import rng_for  # noqa: E402
from topq.hom import (bs_visibility, coincidence, fock_coincidence_bs, gaussian_overlap_sq,  # noqa: E402
                      jitter_averaged_overlap_sq, visibility)
from topq.mzi import coupler, phase  # noqa: E402

DATA, META = ROOT / "data" / "day11_hom.npz", ROOT / "data" / "day11_hom.json"
SEED = 2026_11_05
L = 5 * np.pi
N_ARMS = 4000
SIG_ARMS = (0.01, 0.03)


def mzi_splitter(t_top, t_bot, eta=0.5):
    """MZI set to 50:50 (theta = pi/2) whose arms carry transmissions t_top, t_bot."""
    return coupler(eta) @ np.diag([t_top, t_bot]) @ phase(np.pi / 2) @ coupler(eta)


def compute():
    out, meta = {}, {}
    # (a) visibility vs splitting ratio, analytic + Fock-space points
    eta = np.linspace(0.05, 0.95, 181)
    out["eta"], out["V_eta"] = eta, np.array([bs_visibility(e) for e in eta])
    kl = np.linspace(0.1, 1.4, 9)
    eta_f = np.cos(kl) ** 2
    pc_f = np.array([fock_coincidence_bs(k) for k in kl])
    out["eta_fock"], out["V_fock"] = eta_f, 1 - pc_f / (eta_f**2 + (1 - eta_f) ** 2)
    # coupler gap error -> visibility (Monte Carlo over random couplers)
    sk = np.geomspace(0.003, 0.3, 9)
    Vmc = []
    for i, s in enumerate(sk):
        e = 0.5 - (np.pi / 4) * s * rng_for(SEED, i).standard_normal(20000)
        Vmc.append(np.mean([bs_visibility(np.clip(x, 0, 1)) for x in e]))
    out["sk"], out["V_gap_mc"] = sk, np.array(Vmc)
    out["V_gap_pred"] = 1 - (np.pi**2 / 2) * sk**2
    # (b) dip vs delay: pure, detuned (pure but distinguishable), jittered (mixed)
    tau = np.linspace(-4, 4, 401)                      # in units of 1/sigma_w
    out["tau"] = tau
    out["Pc_pure"] = np.array([0.5 - 0.5 * g for g in gaussian_overlap_sq(tau, 1.0)])
    det = 2 * np.sqrt(np.log(2))                        # chosen so that |J|^2 at tau = 0 is 1/2
    out["Pc_detuned"] = 0.5 - 0.5 * gaussian_overlap_sq(tau, 1.0, det)
    jit = np.sqrt((4 - 1) / 2)                          # chosen so that 1/sqrt(1 + 2 j^2) = 1/2
    out["Pc_jitter"] = 0.5 - 0.5 * jitter_averaged_overlap_sq(tau, 1.0, jit)
    meta["dip_models"] = {"detuning_delta_over_sigma": det, "jitter_over_inv_sigma": jit,
                          "V_detuned": float(1 - out["Pc_detuned"][200] / 0.5),
                          "V_jitter": float(1 - out["Pc_jitter"][200] / 0.5)}
    def fwhm(pc):
        depth = 0.5 - pc
        half = depth.max() / 2
        idx = np.where(depth >= half)[0]
        return float(tau[idx[-1]] - tau[idx[0]])
    meta["dip_fwhm"] = {k: fwhm(out[f"Pc_{k}"]) for k in ("pure", "detuned", "jitter")}
    # (c) MZI-based splitter whose arms are T, D or W elements
    res = {}
    for el in "TDW":
        for kind in ("coupling", "onsite"):
            for si, s in enumerate(SIG_ARMS):
                pool = arm_pool(el, kind, s, 2 * N_ARMS, int(np.random.SeedSequence([SEED, ord(el), len(kind), si]).generate_state(1)[0]), L)
                V = np.array([visibility(mzi_splitter(pool[2 * k], pool[2 * k + 1]), 1.0) for k in range(N_ARMS)])
                out[f"V_{el}_{kind}_{si}"] = V
                res[f"{el}_{kind}_sigma{s}"] = {"mean": float(V.mean()), "p05": float(np.quantile(V, 0.05))}
    meta["mzi_splitter_visibility"] = res
    meta["phase_insensitivity_check"] = float(max(abs(visibility(coupler(0.5) @ np.diag([np.exp(1j * p), 1]), 1.0) - 1)
                                                  for p in np.linspace(0, 2 * np.pi, 13)))
    np.savez(DATA, **out)
    meta.update({"script": "scripts/day11_hom.py",
                 "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                 "python": platform.python_version(), "numpy": np.__version__, "seed": SEED,
                 "max_abs_V_fock_minus_formula": float(np.max(np.abs(out["V_fock"] - np.array([bs_visibility(e) for e in eta_f])))),
                 "V_gap": {"sigma_kappa": sk.tolist(), "mc": out["V_gap_mc"].tolist(), "pred": out["V_gap_pred"].tolist()}})
    META.write_text(json.dumps(meta, indent=1))
    return meta


def plot(meta):
    D = np.load(DATA)
    C = {"T": "#2a78d6", "D": "#eb6834", "W": "#555555"}
    plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False})
    fig, ax = plt.subplots(1, 3, figsize=(7.6, 2.9), constrained_layout=True)
    a = ax[0]
    a.plot(D["eta"], D["V_eta"], "-", color="0.2", lw=1.2, label=r"$2RT/(R^2+T^2)$")
    a.plot(D["eta_fock"], D["V_fock"], "o", ms=3.5, mfc="white", mec="#2a78d6", label="Fock-space check")
    a.set_xlabel(r"coupler bar fraction $\eta$"); a.set_ylabel("HOM visibility V")
    a.set_title("(a) splitting ratio", loc="left", fontsize=8.5)
    a.legend(fontsize=6.3, frameon=False, loc="lower center")
    a = ax[1]
    t = D["tau"]
    a.plot(t, D["Pc_pure"], "-", color="0.2", lw=1.2, label="identical, pure")
    a.plot(t, D["Pc_detuned"], "--", color="#eb6834", lw=1.2, label="pure, detuned (distinguishable)")
    a.plot(t, D["Pc_jitter"], ":", color="#2a78d6", lw=1.6, label="identical, timing jitter (mixed)")
    a.set_xlabel(r"delay $\tau\,\sigma_\omega$"); a.set_ylabel("coincidence probability")
    a.set_ylim(0, 0.55)
    a.set_title("(b) same depth, different width", loc="left", fontsize=8.5)
    a.legend(fontsize=6, frameon=False, loc="lower left")
    a = ax[2]
    pos = 0
    ticks, labels = [], []
    for kind in ("coupling", "onsite"):
        for si, s in enumerate((0.01, 0.03)):
            for el in "TDW":
                V = D[f"V_{el}_{kind}_{si}"]
                a.boxplot(1 - V, positions=[pos], widths=0.6, showfliers=False,
                          patch_artist=True, boxprops=dict(facecolor=C[el], alpha=0.5, lw=0.6),
                          medianprops=dict(color="k", lw=0.8), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6))
                pos += 1
            ticks.append(pos - 2); labels.append(f"{kind[:4]}.\nσ={s}")
            pos += 0.8
    a.set_yscale("symlog", linthresh=1e-6); a.set_ylim(-1e-7, 1)
    a.set_xticks(ticks, labels, fontsize=6.5)
    a.set_ylabel(r"$1-V$ (MZI as 50:50 splitter)")
    for el in "TDW":
        a.bar(0, 0, color=C[el], alpha=0.5, label=el)
    a.legend(fontsize=6.3, frameon=False, loc="upper left", ncol=3)
    a.set_title("(c) arms made of T, D or W", loc="left", fontsize=8.5)
    for ext in ("png", "pdf"):
        fig.savefig(ROOT / "figures" / f"fig_hom.{ext}", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    meta = compute()
    plot(meta)
    print(json.dumps({k: meta[k] for k in ("dip_models", "dip_fwhm", "mzi_splitter_visibility",
                                          "phase_insensitivity_check", "max_abs_V_fock_minus_formula", "V_gap")}, indent=1))
