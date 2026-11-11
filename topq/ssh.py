"""SSH coupled-waveguide array: real-space Hamiltonians and Bloch bands.

Conventions (fixed for the whole project; see docs/derivations/day1_theory.pdf)
------------------------------------------------------------------------------
* Field  E(x,y,z,t) = sum_n a_n(z) e_n(x,y) exp(i(beta0 z - omega t)).
* Coupled-mode equations:  da/dz = i H a   =>   a(z) = expm(i H z) a(0).
* An eigenvalue b of H is a propagation-constant SHIFT: the supermode
  propagates as exp(i (beta0 + b) z).
* Couplings k1 (intra-cell, A_n -> B_n) and k2 (inter-cell, B_n -> A_{n+1})
  are real, positive, nearest-neighbour. Units: any inverse length
  (e.g. mm^-1); every eigenvalue comes out in the same unit.
* The open chain always starts on an A site with k1.
* Lattice momentum k is dimensionless (k = q * Lambda), k in (-pi, pi].
"""
from __future__ import annotations

import numpy as np

__all__ = [
    "ssh_couplings",
    "hamiltonian_from_couplings",
    "ssh_hamiltonian",
    "bloch_hamiltonian",
    "bloch_bands",
    "band_gap",
]


def _check_coupling(name: str, value: float) -> float:
    value = float(value)
    if not np.isfinite(value) or value < 0.0:
        raise ValueError(f"{name} must be a finite, non-negative coupling; got {value!r}")
    return value


def ssh_couplings(M: int, k1: float, k2: float) -> np.ndarray:
    """Nearest-neighbour couplings of an open chain of M waveguides.

    Returns the length-(M-1) array [k1, k2, k1, k2, ...]. Element j couples
    waveguide j to waveguide j+1 (0-based).
    """
    if int(M) != M or M < 2:
        raise ValueError(f"M must be an integer >= 2; got {M!r}")
    k1 = _check_coupling("k1", k1)
    k2 = _check_coupling("k2", k2)
    c = np.empty(int(M) - 1)
    c[0::2] = k1
    c[1::2] = k2
    return c


def hamiltonian_from_couplings(couplings, onsite=None, closing: float | None = None) -> np.ndarray:
    """Build a real symmetric tridiagonal coupled-mode matrix.

    Parameters
    ----------
    couplings : array, shape (M-1,)
        kappa_{j, j+1}. Must be non-negative.
    onsite : array, shape (M,), optional
        On-site propagation-constant detuning delta_beta_j (same units).
        None means zero detuning (chiral-symmetric for bipartite couplings).
    closing : float, optional
        Coupling between the last and first waveguides (periodic ring).

    This general builder will be reused for disordered arrays (Day 3).
    """
    c = np.asarray(couplings, dtype=float)
    if c.ndim != 1 or c.size < 1:
        raise ValueError("couplings must be a 1-D array with at least one entry")
    if np.any(~np.isfinite(c)) or np.any(c < 0):
        raise ValueError("couplings must be finite and non-negative")
    M = c.size + 1
    H = np.diag(c, 1) + np.diag(c, -1)
    if onsite is not None:
        d = np.asarray(onsite, dtype=float)
        if d.shape != (M,):
            raise ValueError(f"onsite must have shape ({M},); got {d.shape}")
        H = H + np.diag(d)
    if closing is not None:
        kc = _check_coupling("closing", closing)
        if M < 3:
            raise ValueError("a ring needs at least 3 waveguides")
        H[0, -1] += kc
        H[-1, 0] += kc
    return H


def ssh_hamiltonian(M: int, k1: float, k2: float, onsite=None, periodic: bool = False) -> np.ndarray:
    """Real-space SSH coupled-mode matrix for M waveguides.

    Open chain (periodic=False): A1 -k1- B1 -k2- A2 -k1- B2 ...
    Ring (periodic=True): requires even M = 2N; the last B is joined to the
    first A with k2, so the ring has exactly N complete unit cells.
    """
    if periodic and (int(M) != M or M % 2 != 0 or M < 4):
        raise ValueError("a periodic SSH ring needs an even number M >= 4 of waveguides")
    c = ssh_couplings(M, k1, k2)
    return hamiltonian_from_couplings(c, onsite=onsite, closing=float(k2) if periodic else None)


def bloch_hamiltonian(k, k1: float, k2: float) -> np.ndarray:
    """2x2 Bloch Hamiltonian H(k) = [[0, h], [h*, 0]], h = k1 + k2 exp(-i k).

    k may be a scalar or an array; the result has shape k.shape + (2, 2).
    """
    k1 = _check_coupling("k1", k1)
    k2 = _check_coupling("k2", k2)
    k = np.asarray(k, dtype=float)
    h = k1 + k2 * np.exp(-1j * k)
    H = np.zeros(k.shape + (2, 2), dtype=complex)
    H[..., 0, 1] = h
    H[..., 1, 0] = np.conj(h)
    return H


def bloch_bands(k, k1: float, k2: float):
    """Analytic bands b_-(k), b_+(k) = -/+ sqrt(k1^2 + k2^2 + 2 k1 k2 cos k)."""
    k1 = _check_coupling("k1", k1)
    k2 = _check_coupling("k2", k2)
    k = np.asarray(k, dtype=float)
    # clip guards against tiny negative round-off when k1 == k2 and k == pi
    b = np.sqrt(np.clip(k1**2 + k2**2 + 2.0 * k1 * k2 * np.cos(k), 0.0, None))
    return -b, b


def band_gap(k1: float, k2: float) -> float:
    """Full bulk gap 2|k1 - k2| (located at k = pi)."""
    return 2.0 * abs(_check_coupling("k1", k1) - _check_coupling("k2", k2))
