"""Day 1: compute SSH bands (analytic + ring numerics) and open-chain spectra,
save them to data/, then draw figures/fig_bands.{png,pdf} FROM THE SAVED FILE.

Parametrisation: mean coupling kbar and dimerisation delta,
    k1 = kbar (1 - delta),  k2 = kbar (1 + delta),
so delta > 0  <=>  k2 > k1 (topological for an A-terminated chain, Day 2) and
delta -> -delta swaps k1 <-> k2 while keeping the band width 2(k1+k2) = 4 kbar fixed.

Run from the repo root:  python scripts/day1_bands.py
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
from topq.ssh import bloch_bands, ssh_hamiltonian  # noqa: E402

DATA = ROOT / "data" / "day1_bands.npz"
META = ROOT / "data" / "day1_bands.json"
FIG = ROOT / "figures" / "fig_bands"

KBAR = 1.0                      # all propagation constants in units of kbar
DELTAS = (-1 / 3, 0.0, +1 / 3)  # k1/k2 = 2, 1, 0.5
N_RING = 24                     # unit cells in the periodic ring (even -> contains k = pi)
M_OPEN = 40                     # waveguides in the open chain (even: 20 complete cells)


def compute():
    k = np.linspace(-np.pi, np.pi, 801)
    km = 2 * np.pi * np.arange(N_RING) / N_RING
    km = np.where(km > np.pi, km - 2 * np.pi, km)          # fold into (-pi, pi]
    out = {"k": k, "k_ring": km, "deltas": np.array(DELTAS)}
    for i, d in enumerate(DELTAS):
        k1, k2 = KBAR * (1 - d), KBAR * (1 + d)
        bm, bp = bloch_bands(k, k1, k2)
        out[f"band_minus_{i}"], out[f"band_plus_{i}"] = bm, bp
        # Ring eigenvalues with momenta assigned INDEPENDENTLY of the analytic
        # bands: for each real eigenvector psi, the one-cell translation T
        # (site j -> j+2) gives <psi|T|psi> = cos k, so |k| = arccos(<T>).
        ev, vec = np.linalg.eigh(ssh_hamiltonian(2 * N_RING, k1, k2, periodic=True))
        cos_k = np.einsum("jm,jm->m", vec, np.roll(vec, -2, axis=0))
        kabs = np.arccos(np.clip(cos_k, -1.0, 1.0))
        out[f"ring_eigs_{i}"], out[f"ring_kabs_{i}"] = ev, kabs
        an = np.sort(np.concatenate(bloch_bands(km, k1, k2)))
        out[f"ring_maxerr_{i}"] = np.max(np.abs(np.sort(ev) - an))
        # point-wise check: each numerical (|k|, b) pair lies on an analytic band
        out[f"ring_pointerr_{i}"] = np.max(np.abs(np.abs(ev) - bloch_bands(kabs, k1, k2)[1]))
        out[f"open_eigs_{i}"] = np.sort(np.linalg.eigvalsh(ssh_hamiltonian(M_OPEN, k1, k2)))
    np.savez(DATA, **out)
    meta = {
        "script": "scripts/day1_bands.py",
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "matplotlib": matplotlib.__version__,
        "kbar": KBAR, "deltas": list(DELTAS), "N_ring": N_RING, "M_open": M_OPEN,
        "ring_max_abs_error_sorted": {str(d): float(out[f"ring_maxerr_{i}"]) for i, d in enumerate(DELTAS)},
        "ring_max_abs_error_pointwise": {str(d): float(out[f"ring_pointerr_{i}"]) for i, d in enumerate(DELTAS)},
        "units": "propagation-constant shift b in units of kbar (mean coupling)",
        "randomness": "none (deterministic)",
    }
    META.write_text(json.dumps(meta, indent=2))
    return meta


def plot():
    D = np.load(DATA)
    colors = ["#eb6834", "#2a78d6", "#1baf7a"]          # orange, blue, aqua (fixed order)
    styles = ["-", ":", "--"]
    labels = [r"$\delta=-1/3$  ($\kappa_1/\kappa_2=2$)",
              r"$\delta=0$  (uniform)",
              r"$\delta=+1/3$  ($\kappa_1/\kappa_2=1/2$)"]

    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(7.2, 3.7), constrained_layout=True)

    for i in range(3):
        lw = 2.4 if i == 0 else 1.4
        ax.plot(D["k"], D[f"band_plus_{i}"], styles[i], color=colors[i], lw=lw, label=labels[i])
        ax.plot(D["k"], D[f"band_minus_{i}"], styles[i], color=colors[i], lw=lw)
    # ring eigenvalues of the delta=+1/3 case: they fall on the delta=-1/3 curve too
    kk = np.r_[D["ring_kabs_2"], -D["ring_kabs_2"]]
    bb = np.r_[D["ring_eigs_2"], D["ring_eigs_2"]]
    ax.plot(kk, bb, "o", ms=3.5, mfc="white", mec="#1baf7a", mew=1.0,
            label=f"ring eigenvalues, $\\delta=+1/3$ (N={len(D['k_ring'])} cells)")
    ax.set_xlabel(r"Bloch momentum $k=q\Lambda$ (rad)")
    ax.set_ylabel(r"propagation-constant shift $b/\bar\kappa$")
    ax.set_xticks([-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi],
                  [r"$-\pi$", r"$-\pi/2$", "0", r"$\pi/2$", r"$\pi$"])
    ax.axhline(0, color="0.8", lw=0.6, zorder=0)
    ax.set_title("(a) bulk bands: identical for $\\pm\\delta$", loc="left", fontsize=9)
    ax.legend(fontsize=6.5, frameon=False, loc="upper center", ncol=2, bbox_to_anchor=(0.5, -0.2))
    ax.set_ylim(-2.2, 2.2)

    n = np.arange(1, len(D["open_eigs_0"]) + 1)
    bx.plot(n, D["open_eigs_0"], "s", ms=3.5, color=colors[0], label=labels[0])
    bx.plot(n, D["open_eigs_2"], "o", ms=3.5, mfc="white", mec=colors[2], mew=1.0, label=labels[2])
    bx.axhline(0, color="0.8", lw=0.6, zorder=0)
    bx.set_xlabel("eigenvalue index (sorted)")
    bx.set_ylabel(r"$b/\bar\kappa$")
    bx.set_title(f"(b) open chain, M={len(n)} waveguides", loc="left", fontsize=9)
    bx.legend(fontsize=6.5, frameon=False, loc="upper center", ncol=1, bbox_to_anchor=(0.5, -0.2))
    mid = len(n) // 2
    bx.annotate("two in-gap modes\nnear $b=0$ (Day 2)", xy=(mid + 0.5, 0.0), xytext=(mid + 3, -1.4),
                fontsize=7, arrowprops=dict(arrowstyle="->", lw=0.7))

    for ext in ("png", "pdf"):
        fig.savefig(f"{FIG}.{ext}", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    meta = compute()
    plot()
    print(json.dumps({k: meta[k] for k in meta if k.startswith("ring_")}, indent=2))
    print(f"saved {DATA.relative_to(ROOT)}, {META.relative_to(ROOT)}, figures/fig_bands.png/.pdf")
