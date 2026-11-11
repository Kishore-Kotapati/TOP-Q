"""Day 6: Monte Carlo study 2 -- PHASE and TRANSMISSION statistics.

Two nominally identical arrays on one chip (independent local errors; any chip-wide
common error cancels exactly in their difference) -> differential phase, the quantity an
interferometer built from them would see. Both phase definitions (modal, output-guide).

Writes data/day6_mc_phase.npz + .json and figures/fig_mc_phase.{png,pdf} FROM THE FILE.
Run from repo root:  python scripts/day6_mc_phase.py     (~2 min)
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

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from topq.disorder import first_order_std, first_order_std_width  # noqa: E402
from topq.montecarlo import (clean_amplitudes, clean_mode, make_spec, phase_stats,  # noqa: E402
                             run_propagation_ensemble)

DATA, META = ROOT / "data" / "day6_mc_phase.npz", ROOT / "data" / "day6_mc_phase.json"
SEED = 2026_11_02
M = 81
N_CHIPS = 2000                         # -> 1000 independent arm pairs
LC = np.pi / 2                         # coupling length pi/(2 k2)
LS = LC * np.arange(1, 21)             # 1 ... 20 coupling lengths
I_MAIN = 9                             # 10 coupling lengths
SIGMAS = np.geomspace(0.003, 0.3, 9)
KINDS = ("coupling", "onsite", "width")
SIGMA_LSCAN = 0.03


def seed_for(*k):
    return int(np.random.SeedSequence([SEED, *k]).generate_state(1)[0])


def analyse(r, ref):
    """Per-length statistics for one ensemble."""
    out = {}
    for q in ("modal", "guide"):
        err = np.angle(r[q] / ref[q])                       # phase error of each chip, (n, nL)
        diff = err[0::2] - err[1::2]                         # arm A minus arm B, same chip
        st = [phase_stats(diff[:, j]) for j in range(diff.shape[1])]
        out[f"{q}_R"] = np.array([s["R"] for s in st])
        out[f"{q}_cstd"] = np.array([s["circ_std"] for s in st])
        out[f"{q}_pflip"] = np.array([s["p_flip"] for s in st])
    p = np.abs(r["guide"]) ** 2                              # power left in the centre guide
    out["T_mean"], out["T_std"] = p.mean(0), p.std(0)
    pa, pb = p[0::2], p[1::2]
    out["imbalance_std"] = np.std(pa - pb, axis=0)           # arm-to-arm power difference
    return out


def compute():
    t0 = time.time()
    out = {"Ls": LS, "sigmas": SIGMAS}
    specs = {a: make_spec(a, M=M) for a in "TD"}
    pred = {}
    for a in "TD":
        _, psi = clean_mode(specs[a])
        c = specs[a].couplings
        pred[a] = {"coupling": first_order_std(psi, c, sigma_kappa=1),
                   "onsite": first_order_std(psi, c, sigma_beta=1),
                   "width": first_order_std_width(psi, c, 1, 1)}
    for a in "TD":
        ref = clean_amplitudes(specs[a], LS)
        out[f"clean_guide_power_{a}"] = np.abs(ref["guide"]) ** 2
        for ki, kind in enumerate(KINDS):
            for si, s in enumerate(SIGMAS):
                r = run_propagation_ensemble(specs[a], kind, s, N_CHIPS, seed_for(ord(a), ki, si), LS)
                for k, v in analyse(r, ref).items():
                    out[f"{a}_{kind}_{si}_{k}"] = v
    np.savez(DATA, **out)

    summary = {}
    for a in "TD":
        for kind in KINDS:
            row = {}
            for q in ("modal", "guide"):
                cs = np.array([out[f"{a}_{kind}_{si}_{q}_cstd"][I_MAIN] for si in range(SIGMAS.size)])
                R = np.array([out[f"{a}_{kind}_{si}_{q}_R"][I_MAIN] for si in range(SIGMAS.size)])
                pf = np.array([out[f"{a}_{kind}_{si}_{q}_pflip"][I_MAIN] for si in range(SIGMAS.size)])
                row[q] = {"circ_std_rad": cs.tolist(), "R": R.tolist(), "p_flip": pf.tolist()}
            row["predicted_modal_diff_std_rad"] = (np.sqrt(2) * pred[a][kind] * SIGMAS * LS[I_MAIN]).tolist()
            row["centre_power_mean"] = [float(out[f"{a}_{kind}_{si}_T_mean"][I_MAIN]) for si in range(SIGMAS.size)]
            row["arm_power_imbalance_std"] = [float(out[f"{a}_{kind}_{si}_imbalance_std"][I_MAIN]) for si in range(SIGMAS.size)]
            summary[f"{a}_{kind}"] = row
    meta = {"script": "scripts/day6_mc_phase.py",
            "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "python": platform.python_version(), "numpy": np.__version__,
            "seed": SEED, "chips_per_point": N_CHIPS, "arm_pairs_per_point": N_CHIPS // 2, "M": M,
            "main_length_coupling_lengths": int(I_MAIN + 1), "sigmas": SIGMAS.tolist(),
            "first_order_db_per_sigma": pred, "summary_at_main_length": summary,
            "note": "chip-wide common errors (global shift) cancel exactly in the differential phase",
            "runtime_s": round(time.time() - t0, 1)}
    META.write_text(json.dumps(meta, indent=1))
    return meta


def plot(meta):
    D = np.load(DATA)
    S, Ls = D["sigmas"], D["Ls"] / LC
    C = {"T": "#2a78d6", "D": "#eb6834"}
    plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False})
    fig, ax = plt.subplots(2, 3, figsize=(7.6, 5.4), constrained_layout=True)
    names = {"coupling": "gap errors (coupling)", "onsite": "width errors (on-site)", "width": "realistic mix (width model)"}
    for j, kind in enumerate(KINDS):
        a = ax[0, j]
        for arr in "TD":
            row = meta["summary_at_main_length"][f"{arr}_{kind}"]
            for q, mk, ls in (("modal", "o", "-"), ("guide", "^", ":")):
                y = np.array(row[q]["circ_std_rad"])
                a.loglog(S, np.maximum(y, 1e-16), mk + ls, ms=3.2, lw=0.9, color=C[arr],
                         mfc="white" if arr == "T" else C[arr], label=f"{arr}, {q} phase")
            if max(row["predicted_modal_diff_std_rad"]) > 1e-10:
                a.loglog(S, row["predicted_modal_diff_std_rad"], "-", color=C[arr], lw=2.2, alpha=0.25)
        a.axhline(np.pi / 4, color="0.6", lw=0.6, ls="--")
        a.text(S[0], np.pi / 4 * 1.25, "π/4", fontsize=6.5, color="0.4")
        a.set_ylim(1e-16, 10)
        a.set_xlabel(r"disorder strength $\sigma$")
        if j == 0:
            a.set_ylabel("arm-to-arm phase error\n(circular std, rad)")
        a.set_title(f"({'abc'[j]}) {names[kind]}", loc="left", fontsize=8)
        if j == 0:
            a.legend(fontsize=5.8, frameon=False, loc="center right")

    a = ax[1, 0]
    for kind, ls in (("coupling", "-"), ("onsite", "--"), ("width", ":")):
        si = int(np.argmin(np.abs(S - SIGMA_LSCAN)))
        for arr in "TD":
            y = D[f"{arr}_{kind}_{si}_modal_cstd"]
            a.plot(Ls, y, ls, color=C[arr], lw=1.1, label=f"{arr}, {kind}")
    a.set_xlabel("device length (coupling lengths)")
    a.set_ylabel("modal phase error (rad)")
    a.set_title(f"(d) grows with length (σ≈{S[int(np.argmin(np.abs(S - SIGMA_LSCAN)))]:.2g})", loc="left", fontsize=8)
    a.legend(fontsize=5.6, frameon=False, ncol=2, loc="upper left")
    a.set_ylim(-0.05, 2.2)

    a = ax[1, 1]
    for arr in "TD":
        for kind, ls in (("coupling", "-"), ("onsite", "--")):
            R = meta["summary_at_main_length"][f"{arr}_{kind}"]["guide"]["R"]
            a.semilogx(S, R, ls, marker="o" if arr == "T" else "s", ms=3, color=C[arr],
                       mfc="white" if arr == "T" else C[arr], lw=1, label=f"{arr}, {kind}")
    a.set_ylim(0, 1.03)
    a.set_xlabel(r"$\sigma$"); a.set_ylabel(r"$R=|\langle e^{i\Delta\varphi}\rangle|$ (output-guide)")
    a.set_title("(e) fringe contrast an MZI would keep", loc="left", fontsize=8)
    a.legend(fontsize=5.8, frameon=False, loc="lower left")

    a = ax[1, 2]
    for arr in "TD":
        for kind, ls in (("coupling", "-"), ("onsite", "--")):
            sd = meta["summary_at_main_length"][f"{arr}_{kind}"]["arm_power_imbalance_std"]
            a.loglog(S, sd, ls, color=C[arr], lw=1.1, marker="o" if arr == "T" else "s", ms=3,
                     mfc="white" if arr == "T" else C[arr], label=f"{arr}, {kind}")
    a.set_xlabel(r"$\sigma$"); a.set_ylabel("std of arm-to-arm\noutput-power difference")
    a.set_title("(f) power imbalance: T is WORSE", loc="left", fontsize=8)
    a.legend(fontsize=5.8, frameon=False, loc="upper left")
    for ext in ("png", "pdf"):
        fig.savefig(ROOT / "figures" / f"fig_mc_phase.{ext}", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    meta = compute()
    plot(meta)
    s = meta["summary_at_main_length"]
    for k in s:
        print(k, "modal cstd", np.round(s[k]["modal"]["circ_std_rad"], 4).tolist())
        print(k, "guide cstd", np.round(s[k]["guide"]["circ_std_rad"], 4).tolist(), "pflip", s[k]["guide"]["p_flip"])
        print(k, "pred", np.round(s[k]["predicted_modal_diff_std_rad"], 4).tolist())
        print(k, "power", np.round(s[k]["centre_power_mean"], 3).tolist(), "imb", np.round(s[k]["arm_power_imbalance_std"], 3).tolist())
    print("runtime", meta["runtime_s"])
