"""Day 4: Monte Carlo study 1 -- eigenvalue (propagation-constant) statistics of the
T and D modes under coupling, on-site and width-model disorder. Tests H2-H5.

Writes data/day4_mc_eigen.npz + .json, then figures/fig_mc_eigen.{png,pdf} FROM THE FILE.
Run from repo root:  python scripts/day4_mc_eigen.py      (~1-2 min on a laptop)
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
from topq.arrays import centre  # noqa: E402
from topq.disorder import coupling_disorder, first_order_std, first_order_std_width, rng_for  # noqa: E402
from topq.montecarlo import bootstrap_std, clean_mode, make_spec, run_ensemble  # noqa: E402
from topq.topology import zero_mode_from_couplings  # noqa: E402

DATA, META = ROOT / "data" / "day4_mc_eigen.npz", ROOT / "data" / "day4_mc_eigen.json"
SEED = 2026_10_31
N = 2000                                  # samples per (array, kind, sigma)
SIGMAS = np.geomspace(0.003, 0.3, 11)
KINDS = ("coupling", "onsite", "width")   # width uses ratio = sigma_beta/sigma_kappa = 1
RATIOS = np.array([0.1, 0.3, 1.0, 3.0, 10.0])
SIGMA_RATIO_SCAN = 0.05
N_CONV = np.array([125, 250, 500, 1000, 2000, 4000])
SIGMA_CONV = 0.05
H2_SIGMAS = (0.1, 0.3, 0.6)
H2_M = 201
H2_N = 2000


def seed_for(*keys):
    """Distinct, documented seed per experiment block (avoids reusing random numbers)."""
    return int(np.random.SeedSequence([SEED, *keys]).generate_state(1)[0])


def compute():
    t0 = time.time()
    out, meta_stats = {"sigmas": SIGMAS, "ratios": RATIOS, "n_conv": N_CONV}, {}
    specs = {a: make_spec(a) for a in "TD"}
    psi = {a: clean_mode(specs[a])[1] for a in "TD"}
    pred = {}
    for a in "TD":
        c = specs[a].couplings
        pred[a] = {"coupling": first_order_std(psi[a], c, sigma_kappa=1),
                   "onsite": first_order_std(psi[a], c, sigma_beta=1),
                   "width": first_order_std_width(psi[a], c, 1, 1)}
    # --- main grid -------------------------------------------------------------
    for a in "TD":
        for ki, kind in enumerate(KINDS):
            stats = np.zeros((SIGMAS.size, 8))   # std, lo, hi, mean, sem, ipr_mean, overlap_mean, overlap_min
            for si, s in enumerate(SIGMAS):
                r = run_ensemble(specs[a], kind, s, N, seed_for(ord(a), ki, si))
                sd, lo, hi = bootstrap_std(r["b"], seed=si)
                stats[si] = [sd, lo, hi, r["b"].mean(), r["b"].std(ddof=1) / np.sqrt(N),
                             r["ipr"].mean(), r["overlap"].mean(), r["overlap"].min()]
            out[f"stats_{a}_{kind}"] = stats
    # --- width-model ratio scan (H5) -------------------------------------------------
    for a in "TD":
        rs = []
        for ri, rho in enumerate(RATIOS):
            r = run_ensemble(specs[a], "width", SIGMA_RATIO_SCAN, N, seed_for(ord(a), 9, ri), ratio=rho)
            rs.append(bootstrap_std(r["b"], seed=ri))
        out[f"ratio_scan_{a}"] = np.array(rs)
    # --- convergence in sample count (T, on-site) ---------------------------------------
    conv = []
    for ni, n in enumerate(N_CONV):
        r = run_ensemble(specs["T"], "onsite", SIGMA_CONV, int(n), seed_for(7, ni))
        conv.append(bootstrap_std(r["b"], seed=ni))
    out["conv"] = np.array(conv)
    # --- H2: localisation of the exact T zero mode under coupling disorder --------------
    from topq.arrays import wall_couplings
    c_clean = wall_couplings(H2_M, 1 / 3, 1.0)
    c0 = centre(H2_M)
    dist = np.arange(0, 30)                       # cells from centre (right side)
    out["h2_dist"] = dist
    for hi, s in enumerate((0.0,) + H2_SIGMAS):
        logs = np.empty((H2_N if s > 0 else 1, dist.size))
        for k in range(logs.shape[0]):
            c = coupling_disorder(c_clean, s, rng_for(seed_for(11, hi), k)) if s > 0 else c_clean
            z = zero_mode_from_couplings(c)
            # exact product formula; on the centre sublattice, right of the centre
            amp = np.abs(z[c0::2][:dist.size]) / np.abs(z[c0])
            logs[k] = np.log(amp)
        out[f"h2_meanlog_{hi}"] = logs.mean(0)
        out[f"h2_stdlog_{hi}"] = logs.std(0)
    # --- strong coupling disorder: does the zero mode stay at the launch site? -----------
    strong = {}
    for si, s in enumerate((0.3, 0.45, 0.6, 0.8)):
        r = run_ensemble(specs["T"], "coupling", s, N, seed_for(13, si))
        strong[str(s)] = {"fraction_tracked_mode_not_zero_mode": float(np.mean(np.abs(r["b"]) > 1e-10)),
                          "mean_overlap": float(r["overlap"].mean()),
                          "mean_centre_power": float(r["p_centre"].mean())}
    np.savez(DATA, **out)

    for a in "TD":
        for kind in KINDS:
            st = out[f"stats_{a}_{kind}"]
            small = SIGMAS <= 0.03
            meta_stats[f"{a}_{kind}"] = {
                "predicted_slope_first_order": pred[a][kind],
                "measured_std_over_sigma_small_sigma_mean": float(np.mean(st[small, 0] / SIGMAS[small])),
                "std_at_sigma_0.3": float(st[-1, 0]), "mean_shift_at_sigma_0.3": float(st[-1, 3]),
                "mean_shift_sem_at_sigma_0.3": float(st[-1, 4]),
                "mean_overlap_at_sigma_0.3": float(st[-1, 6]), "min_overlap_at_sigma_0.3": float(st[-1, 7])}
    rT, rD = out["ratio_scan_T"][:, 0], out["ratio_scan_D"][:, 0]
    meta = {"script": "scripts/day4_mc_eigen.py",
            "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "python": platform.python_version(), "numpy": np.__version__,
            "base_seed": SEED, "samples_per_point": N, "sigmas": SIGMAS.tolist(),
            "arrays": "T: long-long interface, M=41, k1/k2=1/3; D: uniform kbar=2/3, matched Delta=16/9",
            "units": "b and sigma_beta in units of k2; sigma_kappa dimensionless (log-normal)",
            "first_order_predictions_per_unit_sigma": pred,
            "summary": meta_stats,
            "width_ratio_scan": {"sigma_kappa": SIGMA_RATIO_SCAN, "ratios": RATIOS.tolist(),
                                 "std_T": rT.tolist(), "std_D": rD.tolist(),
                                 "T_over_D": (rT / rD).tolist()},
            "convergence_T_onsite_sigma_0.05": {"n": N_CONV.tolist(), "std": out["conv"][:, 0].tolist(),
                                                "ci_low": out["conv"][:, 1].tolist(),
                                                "ci_high": out["conv"][:, 2].tolist()},
            "h2_slope_per_cell": {str(s): float(np.polyfit(dist[:15], out[f"h2_meanlog_{i}"][:15], 1)[0])
                                  for i, s in enumerate((0.0,) + H2_SIGMAS)},
            "h2_expected_slope_clean": float(np.log(1 / 3)),
            "T_strong_coupling_disorder": strong,
            "runtime_s": round(time.time() - t0, 1)}
    META.write_text(json.dumps(meta, indent=2))
    return meta


def plot(meta):
    D = np.load(DATA)
    S = D["sigmas"]
    C = {"T": "#2a78d6", "D": "#eb6834"}
    pred = meta["first_order_predictions_per_unit_sigma"]
    plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False})
    fig, ax = plt.subplots(2, 3, figsize=(7.4, 5.4), constrained_layout=True)
    titles = {"coupling": "(a) coupling disorder (chiral)", "onsite": "(b) on-site disorder",
              "width": r"(c) width model ($\sigma_\beta=\sigma_\kappa$)"}
    for j, kind in enumerate(("coupling", "onsite", "width")):
        a = ax[0, j]
        for arr, mk in (("T", "o"), ("D", "s")):
            st = D[f"stats_{arr}_{kind}"]
            y = np.maximum(st[:, 0], 1e-17)
            a.errorbar(S, y, yerr=[y - np.maximum(st[:, 1], 1e-17), st[:, 2] - y], fmt=mk, ms=3.2,
                       mfc="white" if arr == "T" else C[arr], mec=C[arr], ecolor=C[arr], lw=0.8,
                       label=f"{arr}: Monte Carlo")
            if pred[arr][kind] > 1e-10:
                a.plot(S, pred[arr][kind] * S, "-", color=C[arr], lw=0.9, alpha=0.8,
                       label=f"{arr}: 1st order, {pred[arr][kind]:.2f}σ")
        a.set_xscale("log"); a.set_yscale("log")
        a.set_ylim(1e-17, 1)
        a.set_xlabel(r"disorder strength $\sigma$")
        if j == 0:
            a.set_ylabel(r"std of $\delta b$  (units of $\kappa_2$)")
            a.annotate("T: round-off level\n(exactly pinned)", (0.01, 1e-15), xytext=(0.004, 1e-11), fontsize=6.5,
                       arrowprops=dict(arrowstyle="->", lw=0.6))
        a.set_title(titles[kind], loc="left", fontsize=8.5)
        a.legend(fontsize=6, frameon=False, loc="center right" if kind == "coupling" else "lower right")

    a = ax[1, 0]
    for arr, mk in (("T", "o"), ("D", "s")):
        for kind, ls in (("coupling", "-"), ("onsite", "--"), ("width", ":")):
            st = D[f"stats_{arr}_{kind}"]
            a.plot(S, st[:, 6], ls, marker=mk, ms=2.6, color=C[arr], lw=0.9,
                   mfc="white" if arr == "T" else C[arr], label=f"{arr} {kind}")
    a.set_xscale("log"); a.set_xlabel(r"$\sigma$"); a.set_ylabel(r"mean overlap $|\langle\psi_0|\psi\rangle|^2$")
    a.set_title("(d) mode SHAPE change", loc="left", fontsize=8.5)
    a.legend(fontsize=5.6, frameon=False, ncol=2, loc="lower left")

    a = ax[1, 1]
    dist = D["h2_dist"]
    cols = ["0.2", "#2a78d6", "#1baf7a", "#eb6834"]
    for i, s in enumerate((0.0, 0.1, 0.3, 0.6)):
        m, sd = D[f"h2_meanlog_{i}"], D[f"h2_stdlog_{i}"]
        a.plot(dist, m, "-", color=cols[i], lw=1.1, label=fr"$\sigma_\kappa={s}$")
        if s > 0:
            a.fill_between(dist, m - sd, m + sd, color=cols[i], alpha=0.15, lw=0)
    a.set_xlabel("cells from the centre"); a.set_ylabel(r"$\ln|\psi_n/\psi_0|$ (T zero mode)")
    a.set_title("(e) H2: T localisation", loc="left", fontsize=8.5)
    a.legend(fontsize=6, frameon=False, loc="lower left", title="mean ± 1 std (2000 samples)", title_fontsize=6)

    a = ax[1, 2]
    rs = meta["width_ratio_scan"]
    R = np.array(rs["ratios"])
    a.plot(R, rs["std_T"], "o-", color=C["T"], mfc="white", ms=3.5, lw=1, label="T")
    a.plot(R, rs["std_D"], "s-", color=C["D"], ms=3.2, lw=1, label="D")
    a2 = a.twinx()
    a2.plot(R, rs["T_over_D"], "k:", lw=1)
    a2.set_ylabel("T / D (dotted)", fontsize=7); a2.set_ylim(0, 1.05)
    a2.spines["top"].set_visible(False)
    a.set_xscale("log"); a.set_yscale("log")
    a.set_xlabel(r"$\sigma_\beta/\sigma_\kappa$  (units of $\kappa_2$)")
    a.set_ylabel(r"std $\delta b$ at $\sigma_\kappa=0.05$")
    a.set_title("(f) H5: width model, unknown ratio", loc="left", fontsize=8.5)
    a.legend(fontsize=6.5, frameon=False, loc="upper left")
    for ext in ("png", "pdf"):
        fig.savefig(ROOT / "figures" / f"fig_mc_eigen.{ext}", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    meta = compute()
    plot(meta)
    print(json.dumps({k: meta[k] for k in ("summary", "width_ratio_scan", "convergence_T_onsite_sigma_0.05",
                                          "h2_slope_per_cell", "h2_expected_slope_clean",
                                          "T_strong_coupling_disorder", "runtime_s")}, indent=1))
