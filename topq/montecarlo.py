"""Monte Carlo ensembles for the T (topological interface) and D (matched defect) arrays.

For every random sample we diagonalise the full H and TRACK the mode of interest as the
eigenvector with the largest overlap |<psi_clean|psi>|^2 with the clean mode. We record:
  b        : its propagation-constant shift (units of k2)
  overlap  : |<psi_clean|psi>|^2  (1 = shape unchanged)
  ipr      : sum |psi|^4  (confinement)
  p_centre : power on the centre guide
Disorder kinds: "coupling" (sigma = sigma_kappa), "onsite" (sigma = sigma_beta / k2),
"width" (sigma = sigma_kappa, sigma_beta = ratio * sigma), "global" (sigma = sigma_beta).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .arrays import centre, matched_defect_detuning, wall_couplings
from .disorder import (coupling_disorder, global_shift, onsite_disorder, rng_for,
                       width_disorder)
from .ssh import hamiltonian_from_couplings

KINDS = ("coupling", "onsite", "width", "global")


@dataclass(frozen=True)
class ArraySpec:
    name: str              # "T" or "D"
    M: int
    couplings: np.ndarray
    onsite: np.ndarray     # clean on-site terms (Delta on the centre for D)


def make_spec(name: str, M: int = 41, k1: float = 1 / 3, k2: float = 1.0) -> ArraySpec:
    if name == "T":
        return ArraySpec("T", M, wall_couplings(M, k1, k2), np.zeros(M))
    if name == "D":
        kbar = 0.5 * (k1 + k2)
        on = np.zeros(M)
        on[centre(M)] = matched_defect_detuning(k1, k2, kbar)
        return ArraySpec("D", M, np.full(M - 1, kbar), on)
    raise ValueError("array name must be 'T' or 'D'")


def clean_mode(spec: ArraySpec):
    e, v = np.linalg.eigh(hamiltonian_from_couplings(spec.couplings, onsite=spec.onsite))
    i = int(np.argmin(np.abs(e))) if spec.name == "T" else int(np.argmax(e))
    psi = v[:, i] * np.sign(v[centre(spec.M), i])
    return float(e[i]), psi


def disordered_H(spec: ArraySpec, kind: str, sigma: float, rng, ratio: float = 1.0):
    c, on = spec.couplings, spec.onsite.copy()
    if kind == "coupling":
        c = coupling_disorder(c, sigma, rng)
    elif kind == "onsite":
        on = on + onsite_disorder(spec.M, sigma, rng)
    elif kind == "width":
        c, d = width_disorder(c, ratio * sigma, sigma, rng)
        on = on + d
    elif kind == "global":
        on = on + global_shift(spec.M, sigma, rng)
    else:
        raise ValueError(f"unknown disorder kind {kind!r}")
    return hamiltonian_from_couplings(c, onsite=on)


def run_ensemble(spec: ArraySpec, kind: str, sigma: float, n: int, seed: int,
                 ratio: float = 1.0) -> dict:
    """n independent samples; sample s uses rng_for(seed, s) -> fully reproducible."""
    b0, psi0 = clean_mode(spec)
    c0 = centre(spec.M)
    out = {k: np.empty(n) for k in ("b", "overlap", "ipr", "p_centre")}
    for s in range(n):
        e, v = np.linalg.eigh(disordered_H(spec, kind, sigma, rng_for(seed, s), ratio))
        ov = (psi0 @ v) ** 2
        i = int(np.argmax(ov))
        p = v[:, i] ** 2
        out["b"][s] = e[i] - b0                 # shift from the clean value
        out["overlap"][s] = ov[i]
        out["ipr"][s] = np.sum(p**2)
        out["p_centre"][s] = p[c0]
    return out


def bootstrap_std(x, n_boot: int = 1000, seed: int = 0, ci: float = 0.95):
    """(std, low, high): percentile bootstrap confidence interval of the std."""
    x = np.asarray(x)
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, x.size, size=(n_boot, x.size))
    boots = x[idx].std(axis=1, ddof=1)
    a = (1 - ci) / 2
    return float(x.std(ddof=1)), float(np.quantile(boots, a)), float(np.quantile(boots, 1 - a))


def run_propagation_ensemble(spec: ArraySpec, kind: str, sigma: float, n: int, seed: int,
                             lengths, ratio: float = 1.0) -> dict:
    """Centre-guide launch; for every sample and every length L returns
       modal : <psi_clean | U(L) | c>   (what an ideal mode filter sees)
       guide : <c | U(L) | c>           (field in the centre output guide)
    as complex arrays of shape (n, len(lengths)). Phases are relative to beta0.
    Sample s uses rng_for(seed, s): the same seed reproduces the same chips."""
    Ls = np.asarray(lengths, dtype=float)
    c0 = centre(spec.M)
    _, psi0 = clean_mode(spec)
    modal = np.empty((n, Ls.size), dtype=complex)
    guide = np.empty((n, Ls.size), dtype=complex)
    for s in range(n):
        b, V = np.linalg.eigh(disordered_H(spec, kind, sigma, rng_for(seed, s), ratio))
        vc = V[c0, :]                         # <m|c>
        ph = np.exp(1j * np.outer(Ls, b))     # (nL, M)
        modal[s] = ph @ ((psi0 @ V) * vc)
        guide[s] = ph @ (vc * vc)
    return {"modal": modal, "guide": guide}


def clean_amplitudes(spec: ArraySpec, lengths) -> dict:
    """Same quantities for the clean array (reference for phase errors)."""
    r = run_propagation_ensemble(spec, "onsite", 0.0, 1, 0, lengths)
    return {"modal": r["modal"][0], "guide": r["guide"][0]}


def phase_stats(dphi) -> dict:
    """Circular statistics of phase errors (radians, any range):
       R    = |<exp(i dphi)>|  (1 = perfectly stable; it is the average fringe visibility
              an interferometer would keep if this were its arm-phase error)
       circ_std = sqrt(-2 ln R)  (equals the ordinary std for small Gaussian errors)
       p_flip   = fraction with |dphi| > pi/2 (wrapped)."""
    d = np.angle(np.exp(1j * np.asarray(dphi)))
    R = float(np.abs(np.mean(np.exp(1j * d))))
    return {"R": R, "circ_std": float(np.sqrt(-2 * np.log(max(R, 1e-300)))),
            "std_wrapped": float(np.std(d)), "p_flip": float(np.mean(np.abs(d) > np.pi / 2))}
