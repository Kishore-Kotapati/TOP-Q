"""Directional couplers, phase shifters, Mach-Zehnder interferometers (MZIs) and
rectangular (Clements-type) MZI meshes, with error models and figures of merit.

Conventions (see docs/derivations/day8_theory.pdf)
* Two-mode amplitudes (top, bottom). Same sign convention as the arrays: a lossless
  directional coupler is the solution of da/dz = i*kappa*[[0,1],[1,0]] a over its length:
      B(eta) = [[sqrt(eta), i sqrt(1-eta)], [i sqrt(1-eta), sqrt(eta)]]
  eta = bar (straight-through) power fraction; eta = 1/2 is an ideal 50:50 coupler.
* Phase shifter on the TOP arm: P(phi) = diag(exp(i phi), 1).
* MZI(theta, phi) = B(eta2) P(theta) B(eta1) P(phi): external phase phi first, then
  coupler, internal phase theta, coupler. With ideal couplers:
      |T11|^2 = sin^2(theta/2) (bar),   |T21|^2 = cos^2(theta/2) (cross).
"""
from __future__ import annotations

import numpy as np

__all__ = ["coupler", "phase", "mzi", "embed", "clements_layout", "mesh_unitary",
           "random_mesh_params", "fidelity", "intensity_error", "unitarity_error"]


def coupler(eta: float = 0.5) -> np.ndarray:
    if not 0.0 <= eta <= 1.0:
        raise ValueError("coupler bar power fraction eta must be in [0, 1]")
    t, k = np.sqrt(eta), np.sqrt(1.0 - eta)
    return np.array([[t, 1j * k], [1j * k, t]])


def phase(phi: float) -> np.ndarray:
    return np.diag([np.exp(1j * phi), 1.0 + 0j])


def mzi(theta: float, phi: float, eta1: float = 0.5, eta2: float = 0.5) -> np.ndarray:
    return coupler(eta2) @ phase(theta) @ coupler(eta1) @ phase(phi)


def embed(T2: np.ndarray, N: int, m: int) -> np.ndarray:
    """N x N identity with the 2x2 block T2 acting on modes (m, m+1)."""
    U = np.eye(N, dtype=complex)
    U[m:m + 2, m:m + 2] = T2
    return U


def clements_layout(N: int):
    """Rectangular mesh: N columns; column c holds MZIs on modes (m, m+1) for
    m = c % 2, c % 2 + 2, ... -> N(N-1)/2 MZIs in total (Clements et al. 2016 layout)."""
    return [(c, m) for c in range(N) for m in range(c % 2, N - 1, 2)]


def random_mesh_params(N: int, rng: np.random.Generator) -> dict:
    """Random settings for the clean (target) mesh: theta, phi uniform in [0, 2pi)
    for every MZI, plus N output phases. (Not Haar-distributed: stated as a choice.)"""
    n = len(clements_layout(N))
    return {"theta": rng.uniform(0, 2 * np.pi, n), "phi": rng.uniform(0, 2 * np.pi, n),
            "out": rng.uniform(0, 2 * np.pi, N)}


def mesh_unitary(N: int, theta, phi, out, eta1=None, eta2=None, arm_amp=None) -> np.ndarray:
    """Transfer matrix of the mesh. Light passes the columns in order, so later MZIs
    multiply from the LEFT. Optional per-MZI coupler ratios eta1/eta2 (default 1/2) and
    per-MZI complex amplitude factors arm_amp = (top_arm, bottom_arm) placed between the
    couplers (used for arm loss / arm amplitude errors; default 1, 1)."""
    lay = clements_layout(N)
    n = len(lay)
    e1 = np.full(n, 0.5) if eta1 is None else np.asarray(eta1)
    e2 = np.full(n, 0.5) if eta2 is None else np.asarray(eta2)
    U = np.eye(N, dtype=complex)
    for k, (_, m) in enumerate(lay):
        inner = phase(theta[k])
        if arm_amp is not None:
            inner = np.diag(arm_amp[k]) @ inner
        T2 = coupler(e2[k]) @ inner @ coupler(e1[k]) @ phase(phi[k])
        U = embed(T2, N, m) @ U
    return np.diag(np.exp(1j * np.asarray(out))) @ U


def fidelity(U0: np.ndarray, U: np.ndarray) -> float:
    """Normalised unitary fidelity F = |Tr(U0^dag U)|^2 / (N Tr(U^dag U)).
    F = 1 iff U = (complex number) * U0 (a global phase and a uniform loss do not count);
    F <= 1 by Cauchy-Schwarz."""
    N = U0.shape[0]
    return float(np.abs(np.trace(U0.conj().T @ U)) ** 2 / (N * np.real(np.trace(U.conj().T @ U))))


def intensity_error(U0: np.ndarray, U: np.ndarray) -> float:
    """Mean over the N single-port inputs of the total-variation distance between the
    normalised output power distributions: 0 = identical, 1 = completely different.
    Read it as 'fraction of the light that ends up in the wrong output ports'."""
    P0 = np.abs(U0) ** 2
    P = np.abs(U) ** 2
    P0 = P0 / P0.sum(axis=0, keepdims=True)
    P = P / P.sum(axis=0, keepdims=True)
    return float(0.5 * np.mean(np.sum(np.abs(P - P0), axis=0)))


def unitarity_error(U: np.ndarray) -> float:
    return float(np.linalg.norm(U.conj().T @ U - np.eye(U.shape[0])))
