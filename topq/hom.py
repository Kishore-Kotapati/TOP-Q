"""Hong-Ou-Mandel (HOM) two-photon interference.

One photon enters each input of a lossless or lossy 2x2 device with matrix T
(same coupler convention as topq.mzi). Each photon also has an internal state
(spectrum, arrival time, spatial profile, polarisation) described by a density matrix.
For inputs in internal states rho1, rho2 the coincidence probability is
    P_c = |T11 T22|^2 + |T12 T21|^2 + 2 Re(T11 T22 conj(T12 T21)) * Tr(rho1 rho2)
(derived in docs/derivations/day11_theory.pdf). For a beam splitter B(eta):
    P_c = R^2 + T^2 - 2 R T Tr(rho1 rho2),   T = eta, R = 1 - eta.
Visibility: V = 1 - P_c(actual) / P_c(fully distinguishable, Tr = 0).
"""
from __future__ import annotations

import numpy as np
from scipy.linalg import expm

__all__ = ["coincidence", "visibility", "bs_visibility", "gaussian_overlap_sq",
           "jitter_averaged_overlap_sq", "fock_coincidence_bs", "purity_from_spectrum"]


def coincidence(T: np.ndarray, overlap_sq: float) -> float:
    """Coincidence probability for one photon in each input; overlap_sq = Tr(rho1 rho2)
    (= |<psi1|psi2>|^2 for pure states)."""
    a, b, c, d = T[0, 0], T[1, 1], T[0, 1], T[1, 0]
    return float(abs(a * b) ** 2 + abs(c * d) ** 2 + 2 * np.real(a * b * np.conj(c * d)) * overlap_sq)


def visibility(T: np.ndarray, overlap_sq: float = 1.0) -> float:
    dist = coincidence(T, 0.0)
    return float(1 - coincidence(T, overlap_sq) / dist) if dist > 0 else 0.0


def bs_visibility(eta: float, overlap_sq: float = 1.0) -> float:
    """Closed form for a beam splitter: V = 2RT |J|^2 / (R^2 + T^2)."""
    R, Tt = 1 - eta, eta
    return float(2 * R * Tt * overlap_sq / (R**2 + Tt**2))


def gaussian_overlap_sq(tau, sigma_w: float, delta_w: float = 0.0):
    """|<psi1|psi2>|^2 for two photons with Gaussian spectral intensity of rms width sigma_w,
    relative delay tau and centre-frequency difference delta_w:
        exp(-sigma_w^2 tau^2) * exp(-delta_w^2 / (4 sigma_w^2))."""
    tau = np.asarray(tau, dtype=float)
    return np.exp(-(sigma_w * tau) ** 2) * np.exp(-delta_w**2 / (4 * sigma_w**2))


def jitter_averaged_overlap_sq(tau, sigma_w: float, jitter: float):
    """Average of |J|^2 when the delay fluctuates from shot to shot by a Gaussian of rms
    `jitter` (a MIXED state: decoherence-like averaging):
        <exp(-s^2 (tau + x)^2)> = exp(-s^2 tau^2 / (1 + 2 s^2 j^2)) / sqrt(1 + 2 s^2 j^2)."""
    tau = np.asarray(tau, dtype=float)
    q = 1 + 2 * (sigma_w * jitter) ** 2
    return np.exp(-(sigma_w * tau) ** 2 / q) / np.sqrt(q)


def fock_coincidence_bs(kappa_l: float) -> float:
    """INDEPENDENT check: build the two-mode Fock space (0..2 photons per mode), evolve
    |1,1> with H = kappa_l (a^dag b + b^dag a), return |<1,1|U|1,1>|^2."""
    n = 3
    a = np.diag(np.sqrt(np.arange(1, n)), 1)
    I = np.eye(n)
    A, B = np.kron(a, I), np.kron(I, a)
    H = kappa_l * (A.conj().T @ B + B.conj().T @ A)
    U = expm(1j * H)
    s11 = np.zeros(n * n); s11[1 * n + 1] = 1.0
    return float(abs(s11 @ U @ s11) ** 2)


def purity_from_spectrum(weights) -> float:
    """Purity Tr(rho^2) of a photon whose state is an incoherent mixture of orthogonal
    spectral components with probabilities `weights`."""
    w = np.asarray(weights, dtype=float)
    w = w / w.sum()
    return float(np.sum(w**2))
