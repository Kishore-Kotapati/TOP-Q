"""Days 9-10: does the array-level protection survive inside an MZI mesh?

Every internal MZI arm of an N = 6 rectangular mesh contains one arm element
(T: topological interface array, D: matched defect array, W: plain single waveguide),
each with independent fabrication errors, length L = 10 coupling lengths.
Optionally the couplers carry the same gap errors (sigma_eta = pi/4 sigma_kappa).

Writes data/day9_bridge.npz + .json and, FROM THEM, figures/fig_bridge.{png,pdf}
and figures/fig_design_map.{png,pdf}.   Run from repo root (~6 min).
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
from topq.bridge import arm_pool, clean_insertion_loss, run_bridge  # noqa: E402
from topq.mzi import clements_layout  # noqa: E402

DATA, META = ROOT / "data" / "day9_bridge.npz", ROOT / "data" / "day9_bridge.json"
SEED = 2026_11_04
N = 6
L = 5 * np.pi                        # 10 coupling lengths
N_MESH = 400
SIGMAS = np.array([0.001, 0.003, 0.01, 0.03, 0.1])
ELEMENTS = ("T", "D", "W")
KINDS = ("coupling", "onsite", "width")
SIGMA_CAL = 0.003                    # small sigma used to calibrate first-order arm coefficients
N_CAL = 4000


def seed_for(*k):
    return int(np.random.SeedSequence([SEED, *k]).generate_state(1)[0])


def compute():
    t0 = time.time()
    n_arms = N_MESH * len(clements_layout(N)) * 2
    out = {"sigmas": SIGMAS}
    summ = {}
    for ei, el in enumerate(ELEMENTS):
        for ki, kind in enumerate(KINDS):
            for si, s in enumerate(SIGMAS):
                sd = seed_for(ei, ki, si)
                pool = arm_pool(el, kind, s, n_arms, sd, L)
                for cp in (False, True):
                    if cp and kind == "onsite":
                        continue                      # on-site errors do not touch couplers
                    r = run_bridge(el, kind, s, N, N_MESH, sd, L, with_couplers=cp, pool=pool)
                    tag = f"{el}_{kind}_{'cpl' if cp else 'arms'}_{si}"
                    for q in ("inf", "ie", "power"):
                        out[f"{tag}_{q}"] = r[q]
                    summ[tag] = {"sigma": float(s), "inf_mean": float(r["inf"].mean()),
                                 "inf_sem": float(r["inf"].std() / np.sqrt(N_MESH)),
                                 "pred_first_order": r["pred"], "ie_mean": float(r["ie"].mean()),
                                 "arm_phase_std": r["arm_phase_std"], "arm_logamp_std": r["arm_logamp_std"]}
    # --- first-order arm coefficients (per unit sigma) for the design map -----------------
    coef = {}
    for el in ELEMENTS:
        coef[el] = {}
        for kind in ("coupling", "onsite"):
            t = arm_pool(el, kind, SIGMA_CAL, N_CAL, seed_for(99, ord(el), len(kind)), L)
            coef[el][kind] = {"phase_per_sigma": float(np.std(np.angle(t)) / SIGMA_CAL),
                              "logamp_per_sigma": float(np.std(np.log(np.abs(t))) / SIGMA_CAL)}
    np.savez(DATA, **out)
    meta = {"script": "scripts/day9_bridge.py",
            "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "python": platform.python_version(), "numpy": np.__version__, "seed": SEED,
            "N": N, "L_coupling_lengths": 10, "meshes_per_point": N_MESH, "sigmas": SIGMAS.tolist(),
            "architecture": "each internal MZI arm = one arm element with independent errors; centre-guide in/out",
            "clean_insertion_loss": {el: clean_insertion_loss(el, L) for el in ELEMENTS},
            "arm_coefficients_per_unit_sigma": coef, "summary": summ,
            "runtime_s": round(time.time() - t0, 1)}
    META.write_text(json.dumps(meta, indent=1))
    return meta


def infidelity_model(coef, el, sb, sk, couplers=True):
    """First-order 1-F for independent on-site (sb, units k2) and gap (sk) errors."""
    c = coef[el]
    var = (c["onsite"]["phase_per_sigma"] * sb) ** 2 + (c["onsite"]["logamp_per_sigma"] * sb) ** 2 \
        + (c["coupling"]["phase_per_sigma"] * sk) ** 2 + (c["coupling"]["logamp_per_sigma"] * sk) ** 2
    return (N - 1) ** 2 / N * var + (2 * (N - 1) * (np.pi / 4 * sk) ** 2 if couplers else 0)


def plot(meta):
    D = np.load(DATA)
    S = D["sigmas"]
    C = {"T": "#2a78d6", "D": "#eb6834", "W": "#555555"}
    MK = {"T": "o", "D": "s", "W": "^"}
    plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False})
    fig, ax = plt.subplots(1, 3, figsize=(7.6, 2.9), constrained_layout=True)
    titles = {"coupling": "(a) gap errors only", "onsite": "(b) width (on-site) errors only",
              "width": r"(c) realistic mix ($\sigma_\beta=\sigma_\kappa$)"}
    sm = meta["summary"]
    for j, kind in enumerate(KINDS):
        a = ax[j]
        for el in ELEMENTS:
            for cp, ls in ((False, "-"), (True, "--")):
                if cp and kind == "onsite":
                    continue
                tags = [f"{el}_{kind}_{'cpl' if cp else 'arms'}_{si}" for si in range(S.size)]
                y = np.array([sm[t]["inf_mean"] for t in tags])
                a.loglog(S, np.maximum(y, 1e-9), MK[el] + ls, color=C[el], ms=3.3, lw=1,
                         mfc="white" if not cp else C[el],
                         label=f"{el}{' + coupler errors' if cp else ''}")
        a.axhline(0.01, color="0.7", lw=0.7, ls=":")
        a.text(S[0], 0.013, "F = 0.99", fontsize=6.3, color="0.4")
        a.set_ylim(1e-9, 2)
        a.set_xlabel(r"disorder strength $\sigma$")
        if j == 0:
            a.set_ylabel(r"mesh infidelity $1-F$ (N = 6)")
        a.set_title(titles[kind], loc="left", fontsize=8.3)
        a.legend(fontsize=5.6, frameon=False, loc="lower right")
    for ext in ("png", "pdf"):
        fig.savefig(ROOT / "figures" / f"fig_bridge.{ext}", dpi=300)
    plt.close(fig)

    # --- design map ------------------------------------------------------------------------
    coef = meta["arm_coefficients_per_unit_sigma"]
    sb = np.geomspace(1e-4, 0.1, 300)
    sk = np.geomspace(1e-4, 0.3, 300)
    SB, SK = np.meshgrid(sb, sk)
    fig, ax = plt.subplots(1, 2, figsize=(7.4, 3.1), constrained_layout=True)
    for j, cp in enumerate((False, True)):
        a = ax[j]
        for el in ELEMENTS:
            Z = infidelity_model(coef, el, SB, SK, couplers=cp)
            a.contour(SB, SK, Z, levels=[0.01], colors=[C[el]], linewidths=1.6,
                      linestyles="-" if el != "W" else "--")
            a.plot([], [], "-" if el != "W" else "--", color=C[el], lw=1.6, label=f"{el}: F = 0.99 boundary")
        a.set_xscale("log"); a.set_yscale("log")
        a.set_xlim(1e-4, 2e-2)
        a.set_xlabel(r"on-site (width) error $\sigma_\beta/\kappa_2$")
        a.set_ylabel(r"gap (coupling) error $\sigma_\kappa$")
        a.set_title(f"({'ab'[j]}) {'arms only' if not cp else 'arms + coupler gap errors'}: "
                    "allowed region is below/left", loc="left", fontsize=8)
        a.legend(fontsize=6.3, frameon=False, loc="upper right")
    for ext in ("png", "pdf"):
        fig.savefig(ROOT / "figures" / f"fig_design_map.{ext}", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    meta = compute()
    plot(meta)
    print(json.dumps({"insertion": meta["clean_insertion_loss"], "coef": meta["arm_coefficients_per_unit_sigma"],
                      "runtime": meta["runtime_s"]}, indent=1))
    for k, v in meta["summary"].items():
        print(k, "1-F %.3g (pred %.3g)  ie %.3g" % (v["inf_mean"], v["pred_first_order"], v["ie_mean"]))
