"""Chiral symmetry, winding number and edge (zero) modes of SSH arrays.

Same conventions as topq.ssh: da/dz = i H a, b = propagation-constant shift,
chain starts on an A site, A sites are indices 0, 2, 4, ... (0-based).
"""
from __future__ import annotations

import numpy as np

__all__ = [
    "chiral_operator",
    "chiral_violation",
    "winding_number",
    "zero_mode_from_couplings",
    "edge_mode",
    "sublattice_weights",
    "ipr",
    "fit_localization_length",
    "xi_analytic",
    "even_chain_splitting",
]


# --- chiral symmetry --------------------------------------------------------
def chiral_operator(M: int) -> np.ndarray:
    """Gamma = diag(+1, -1, +1, -1, ...): +1 on A guides, -1 on B guides."""
    g = np.ones(int(M))
    g[1::2] = -1.0
    return np.diag(g)


def chiral_violation(H: np.ndarray) -> float:
    """|| Gamma H Gamma + H ||_F / || H ||_F.  Zero  <=>  H is chiral-symmetric.

    Relative (divided by ||H||) so the number does not depend on the units of kappa.
    """
    H = np.asarray(H)
    G = chiral_operator(H.shape[0])
    nH = np.linalg.norm(H)
    if nH == 0:
        return 0.0
    return float(np.linalg.norm(G @ H @ G + H) / nH)


# --- winding number ---------------------------------------------------------
def winding_number(k1: float, k2: float, nk: int = 2001, gap_tol: float = 1e-9) -> float:
    """Winding of d(k) = (dx, dy) = (k1 + k2 cos k, k2 sin k) around the origin.

    Convention: W = (1/2pi) * sum of wrapped phase steps of (dx + i dy) as k goes
    once around the Brillouin zone in the + direction. For the A-first unit cell,
    W = 1 if k2 > k1 and W = 0 if k2 < k1. Returns nan when the gap is closed
    (|d| ~ 0 somewhere), because W is then undefined.
    """
    k = np.linspace(0.0, 2 * np.pi, int(nk))           # closed loop: k[-1] == k[0] + 2 pi
    z = (k1 + k2 * np.cos(k)) + 1j * (k2 * np.sin(k))
    if np.min(np.abs(z)) < gap_tol * max(abs(k1), abs(k2), 1e-300):
        return float("nan")
    dphi = np.angle(z[1:] / z[:-1])                    # each step wrapped into (-pi, pi]
    return float(np.sum(dphi) / (2 * np.pi))


# --- zero modes -------------------------------------------------------------
def zero_mode_from_couplings(couplings) -> np.ndarray:
    """EXACT b = 0 eigenvector of an ODD chain (A at both ends), any positive couplings.

    Zero-mode condition on every B site j:  c[2j] a_j + c[2j+1] a_{j+1} = 0
    =>  a_{j+1} = -(c[2j] / c[2j+1]) a_j , and all B amplitudes are 0.
    Works for disordered couplings too (each cell may have its own k1, k2).
    """
    c = np.asarray(couplings, dtype=float)
    M = c.size + 1
    if M % 2 == 0:
        raise ValueError("the exact zero mode formula needs an odd number of waveguides")
    if np.any(c <= 0):
        raise ValueError("couplings must be strictly positive")
    nA = (M + 1) // 2
    a = np.empty(nA)
    a[0] = 1.0
    ratios = -c[0::2] / c[1::2]                        # length nA - 1
    a[1:] = np.cumprod(ratios)
    psi = np.zeros(M)
    psi[0::2] = a
    return psi / np.linalg.norm(psi)


def edge_mode(H: np.ndarray):
    """Eigenpair with the smallest |b| (the in-gap / zero mode candidate)."""
    e, v = np.linalg.eigh(H)
    i = int(np.argmin(np.abs(e)))
    return float(e[i]), v[:, i]


def sublattice_weights(psi):
    """(power on A guides, power on B guides) for a normalised mode."""
    p = np.abs(np.asarray(psi)) ** 2
    return float(p[0::2].sum()), float(p[1::2].sum())


def ipr(psi) -> float:
    """Inverse participation ratio sum|psi|^4 / (sum|psi|^2)^2.
    1 = all power in one guide; ~1/M = spread over the whole array."""
    p = np.abs(np.asarray(psi)) ** 2
    return float(np.sum(p ** 2) / np.sum(p) ** 2)


def xi_analytic(k1: float, k2: float) -> float:
    """Amplitude localisation length in UNIT CELLS: |a_n| ~ exp(-n/xi), xi = 1/ln(k2/k1)."""
    if not (0 < k1 < k2):
        raise ValueError("a left-edge zero mode exists only for 0 < k1 < k2")
    return 1.0 / np.log(k2 / k1)


def fit_localization_length(psi, n_fit: int | None = None, floor: float = 1e-12) -> float:
    """Fit ln|a_n| = c - n/xi on the A-sublattice amplitudes (n = cell index).

    Points below `floor` (relative to the largest) are dropped: they are round-off,
    not physics. Returns xi in unit cells.
    """
    a = np.abs(np.asarray(psi)[0::2])
    a = a / a.max()
    n = np.arange(a.size)
    keep = a > floor
    if n_fit is not None:
        keep &= n < n_fit
    if keep.sum() < 3:
        raise ValueError("too few points above the floor to fit a decay")
    slope = np.polyfit(n[keep], np.log(a[keep]), 1)[0]
    return float(-1.0 / slope)


def even_chain_splitting(k1: float, k2: float, N: int) -> float:
    """First-order |b| of the two hybridised edge modes of an even chain (N cells):
    |t| = k2 (1 - r^2) r^N,  r = k1/k2 < 1   (derived in day2_theory.pdf)."""
    r = k1 / k2
    return float(k2 * (1 - r**2) * r**N)
