"""Day 5: light propagation in T and D arrays, two phase definitions, and loss.

Writes data/day5_propagation.npz (+ .json), then FROM THAT FILE
figures/fig_propagation.{png,pdf} and figures/fig_loss.{png,pdf}.
Run from repo root:  python scripts/day5_propagation.py    (~20 s)
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
from topq.arrays import centre  # noqa: E402
from topq.disorder import rng_for  # noqa: E402
from topq.montecarlo import clean_mode, disordered_H, make_spec  # noqa: E402
from topq.propagation import launch, propagate, with_loss  # noqa: E402
from topq.ssh import hamiltonian_from_couplings  # noqa: E402

DATA, META = ROOT / "data" / "day5_propagation.npz", ROOT / "data" / "day5_propagation.json"
M = 81                          # wide enough that light does not reach the ends before z = L
L = 5 * np.pi                   # = 10 coupling lengths pi/(2 k2); units 1/k2
ZS = np.linspace(0, L, 301)
SEED = 2026_11_01
LOSS_UNIFORM = 0.05             # power loss per unit length (units k2)
LOSS_SIGMAS = np.array([0.01, 0.03, 0.1, 0.3])   # random-loss spread (uniform on [0, 2*mean])
N_LOSS = 500


def H_of(spec):
    return hamiltonian_from_couplings(spec.couplings, onsite=spec.onsite)


def compute():
    out = {"zs": ZS}
    specs = {a: make_spec(a, M=M) for a in "TD"}
    c0 = centre(M)
    meta = {}
    # --- clean propagation maps -------------------------------------------------------
    uni = hamiltonian_from_couplings(np.full(M - 1, 2 / 3))
    cases = {"T_centre": (H_of(specs["T"]), c0), "D_centre": (H_of(specs["D"]), c0),
             "T_neighbour": (H_of(specs["T"]), c0 + 1), "uniform_centre": (uni, c0)}
    for key, (H, g) in cases.items():
        a = propagate(H, launch(M, g), ZS)
        out[f"map_{key}"] = np.abs(a) ** 2
        out[f"edge_{key}"] = np.sum(np.abs(a[:, :5]) ** 2 + np.abs(a[:, -5:]) ** 2, axis=1)
    # --- phases of the clean modes ------------------------------------------------------
    for name in "TD":
        b0, psi0 = clean_mode(specs[name])
        a = propagate(H_of(specs[name]), launch(M, c0), ZS)
        modal = a @ psi0                                   # <psi0|a(z)>, psi0 real
        out[f"modal_power_{name}"] = np.abs(modal) ** 2
        out[f"modal_resid_{name}"] = np.unwrap(np.angle(modal * np.exp(-1j * b0 * ZS)))
        out[f"guide_power_{name}"] = np.abs(a[:, c0]) ** 2
        out[f"guide_resid_{name}"] = np.unwrap(np.angle(a[:, c0] * np.exp(-1j * b0 * ZS)))
        meta[f"{name}_b0"] = b0
        meta[f"{name}_modal_power_at_L"] = float(out[f"modal_power_{name}"][-1])
        meta[f"{name}_centre_guide_power_at_L"] = float(out[f"guide_power_{name}"][-1])
        meta[f"{name}_guide_phase_resid_at_L_rad"] = float(out[f"guide_resid_{name}"][-1])
        # --- one disordered chip of each kind: phase drift relative to the clean chip --
        for kind in ("coupling", "onsite"):
            H = disordered_H(specs[name], kind, 0.1, rng_for(SEED, 1 if kind == "coupling" else 2))
            ad = propagate(H, launch(M, c0), ZS)
            out[f"drift_modal_{name}_{kind}"] = np.unwrap(np.angle((ad @ psi0) / modal))
            out[f"drift_guide_{name}_{kind}"] = np.unwrap(np.angle(ad[:, c0] / a[:, c0]))
    # --- chiral symmetry makes the self-amplitude real; topology keeps it positive ---------
    Ht6 = disordered_H(specs["T"], "coupling", 0.6, rng_for(SEED, 3))
    for key, H in (("uniform", uni), ("T_clean", H_of(specs["T"])), ("T_coupling06", Ht6)):
        a = propagate(H, launch(M, c0), ZS)[:, c0]
        out[f"self_re_{key}"] = a.real
        meta[f"self_amp_{key}_max_abs_imag"] = float(np.max(np.abs(a.imag)))
        meta[f"self_amp_{key}_min_real"] = float(np.min(a.real))
    e6, v6 = np.linalg.eigh(Ht6)
    meta["T_coupling06_zero_mode_weight_on_centre"] = float(v6[c0, np.argmin(np.abs(e6))] ** 2)
    # --- loss ----------------------------------------------------------------------------
    for name in "TD":
        H = H_of(specs[name])
        a0 = propagate(H, launch(M, c0), ZS)
        al = propagate(with_loss(H, LOSS_UNIFORM), launch(M, c0), ZS)
        out[f"uniform_loss_power_{name}"] = np.sum(np.abs(al) ** 2, axis=1)
        out[f"uniform_loss_phase_change_{name}"] = np.angle(al[:, c0] / a0[:, c0])
        b0, _ = clean_mode(specs[name])
        re_std, im_mean, im_std = [], [], []
        for si, s in enumerate(LOSS_SIGMAS):
            re, im = [], []
            for k in range(N_LOSS):
                rng = rng_for(SEED + 100 + si, k)
                alpha = 2 * s * rng.random(M)          # mean s, spread uniform in [0, 2s]
                e = np.linalg.eigvals(with_loss(H, alpha))
                i = int(np.argmin(np.abs(e - b0)))
                re.append(e[i].real - b0)
                im.append(e[i].imag)
            re_std.append(np.std(re)); im_mean.append(np.mean(im)); im_std.append(np.std(im))
        out[f"loss_re_std_{name}"] = np.array(re_std)
        out[f"loss_im_mean_{name}"] = np.array(im_mean)
        out[f"loss_im_std_{name}"] = np.array(im_std)
    out["loss_sigmas"] = LOSS_SIGMAS
    np.savez(DATA, **out)
    meta.update({
        "script": "scripts/day5_propagation.py",
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "python": platform.python_version(), "numpy": np.__version__,
        "M": M, "L_over_inv_k2": L, "L_in_coupling_lengths": L / (np.pi / 2),
        "max_edge_power_before_L": {k[5:]: float(np.max(out[k])) for k in out if k.startswith("edge_")},
        "random_loss": {"mean_alpha": LOSS_SIGMAS.tolist(), "samples": N_LOSS,
                        "T_std_Re_db": out["loss_re_std_T"].tolist(),
                        "D_std_Re_db": out["loss_re_std_D"].tolist(),
                        "T_mean_Im_b": out["loss_im_mean_T"].tolist(),
                        "D_mean_Im_b": out["loss_im_mean_D"].tolist()},
        "uniform_loss": {"alpha": LOSS_UNIFORM,
                         "max_abs_phase_change_T": float(np.max(np.abs(out["uniform_loss_phase_change_T"]))),
                         "max_abs_phase_change_D": float(np.max(np.abs(out["uniform_loss_phase_change_D"])))},
        "seed": SEED})
    META.write_text(json.dumps(meta, indent=2))
    return meta


def plot(meta):
    D = np.load(DATA)
    z = D["zs"]
    c0 = M // 2
    C = {"T": "#2a78d6", "D": "#eb6834"}
    plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False})
    fig = plt.figure(figsize=(8.2, 7.0), constrained_layout=True)
    gs = fig.add_gridspec(2, 12, height_ratios=[1.15, 1])
    titles = {"T_centre": "(a) T, light into centre", "D_centre": "(b) D, light into centre",
              "T_neighbour": "(c) T, light into neighbour", "uniform_centre": "(d) uniform array"}
    for j, key in enumerate(titles):
        a = fig.add_subplot(gs[0, 3 * j:3 * j + 3])
        m = D[f"map_{key}"][:, c0 - 20:c0 + 21]
        im = a.imshow(m.T, origin="lower", aspect="auto", extent=[0, z[-1], -20.5, 20.5],
                      cmap="magma", vmin=0, vmax=0.25)
        a.set_title(titles[key], loc="left", fontsize=8)
        a.set_xlabel(r"$z\ (1/\kappa_2)$")
        if j == 0:
            a.set_ylabel("guide (0 = centre)")
    fig.colorbar(im, ax=fig.axes[:4], shrink=0.8, pad=0.01, label=r"power (clip 0.25)")

    a = fig.add_subplot(gs[1, 0:4])
    for name, mk in (("T", "-"), ("D", "--")):
        a.plot(z, D[f"guide_power_{name}"], mk, color=C[name], lw=1.2, label=f"{name}: centre-guide power")
        a.plot(z, D[f"modal_power_{name}"], mk, color=C[name], lw=0.8, alpha=0.6, label=f"{name}: power in the mode")
    a.set_ylim(0, 1.02); a.set_xlabel(r"$z\ (1/\kappa_2)$"); a.set_ylabel("power")
    a.set_title("(e) clean arrays: what stays at the centre", loc="left", fontsize=8.5)
    a.legend(fontsize=5.8, frameon=False, ncol=1, loc="lower left")

    a = fig.add_subplot(gs[1, 4:8])
    for name in "TD":
        for kind, ls in (("coupling", "-"), ("onsite", ":")):
            a.plot(z, D[f"drift_modal_{name}_{kind}"], ls, color=C[name], lw=1.2,
                   label=f"{name}, one chip, {kind} σ=0.1")
    a.axhline(0, color="0.8", lw=0.6)
    a.set_xlabel(r"$z\ (1/\kappa_2)$"); a.set_ylabel("modal phase error (rad)")
    a.set_title("(f) one chip, σ=0.1: modal phase drift", loc="left", fontsize=8.5)
    a.legend(fontsize=5.8, frameon=False, loc="upper left")

    a = fig.add_subplot(gs[1, 8:12])
    a.plot(z, D["self_re_uniform"], color="0.5", lw=1.1, label="uniform array (chiral, no zero mode)")
    a.plot(z, D["self_re_T_clean"], color=C["T"], lw=1.3, label="T, clean")
    a.plot(z, D["self_re_T_coupling06"], "--", color=C["T"], lw=1.1, label=r"T, coupling $\sigma_\kappa$=0.6")
    a.axhline(0, color="0.8", lw=0.6)
    a.set_xlabel(r"$z\ (1/\kappa_2)$"); a.set_ylabel(r"centre-guide field $a_c(z)$ (it is real)")
    a.set_title("(g) phase is 0 or π; T never flips", loc="left", fontsize=8.5)
    a.legend(fontsize=5.8, frameon=False, loc="lower right")
    for ext in ("png", "pdf"):
        fig.savefig(ROOT / "figures" / f"fig_propagation.{ext}", dpi=300)
    plt.close(fig)

    fig, ax = plt.subplots(1, 2, figsize=(7.4, 2.8), constrained_layout=True)
    a = ax[0]
    for name, ls in (("T", "-"), ("D", "--")):
        a.semilogy(z, D[f"uniform_loss_power_{name}"], ls, color=C[name], lw=1.3, label=f"{name}: total power")
    a.semilogy(z, np.exp(-0.05 * z), "k:", lw=1, label=r"$e^{-\alpha z}$")
    a.set_xlabel(r"$z\ (1/\kappa_2)$"); a.set_ylabel("power")
    a.set_title(f"(a) uniform loss α={meta['uniform_loss']['alpha']}: phase change ≤ "
                f"{max(meta['uniform_loss']['max_abs_phase_change_T'], meta['uniform_loss']['max_abs_phase_change_D']):.0e} rad",
                loc="left", fontsize=8)
    a.legend(fontsize=6.5, frameon=False)
    a = ax[1]
    s = D["loss_sigmas"]
    for name, mk in (("T", "o"), ("D", "s")):
        a.loglog(s, np.maximum(D[f"loss_re_std_{name}"], 1e-17), mk + "-", color=C[name], ms=4,
                 mfc="white" if name == "T" else C[name], label=f"{name}: std of Re δb (phase)")
        a.loglog(s, D[f"loss_im_std_{name}"], mk + ":", color=C[name], ms=3, alpha=0.6,
                 mfc="white" if name == "T" else C[name], label=f"{name}: std of Im b (amplitude)")
    a.set_xlabel(r"mean random loss $\bar\alpha$ (spread $0$–$2\bar\alpha$)")
    a.set_ylabel(r"std (units of $\kappa_2$)")
    a.set_title("(b) random loss: amplitude noise, not phase", loc="left", fontsize=8)
    a.legend(fontsize=6, frameon=False, loc="center right")
    for ext in ("png", "pdf"):
        fig.savefig(ROOT / "figures" / f"fig_loss.{ext}", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    meta = compute()
    plot(meta)
    print(json.dumps(meta, indent=1))
