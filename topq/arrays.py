"""The two arrays compared in TOP-Q (Day 3 design, see docs/derivations/day3_theory.pdf).

T  : topological interface ("long-long" domain wall, after Blanco-Redondo et al. 2016).
     Odd number M of identical guides. The centre guide is joined to BOTH neighbours by
     the weak coupling k1; moving outwards the couplings alternate k2, k1, k2, ...
     Chiral-symmetric; exact zero mode localised on the centre guide's sublattice.
D  : conventional (trivial) defect. Same M, uniform coupling kbar = (k1 + k2)/2
     everywhere, centre guide detuned by Delta (e.g. a wider guide). Not chiral-symmetric.
     Delta is chosen so that the defect mode has the SAME inverse participation ratio
     (IPR) as the T mode  ->  same first-order sensitivity to on-site disorder.
SS : "short-short" wall (strong bonds next to the centre) built with couplings only.
     Kept for discussion: in a nearest-neighbour-only model it is ALSO a domain wall.

All propagation constants in the same (arbitrary) inverse-length unit as the couplings.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

from .ssh import hamiltonian_from_couplings

__all__ = [
    "centre", "wall_couplings", "topological_interface", "short_short_interface",
    "defect_array", "defect_q", "ipr_geometric", "matched_defect_detuning",
    "defect_mode_analytic", "interface_mode_analytic",
]


def centre(M: int) -> int:
    if int(M) != M or M < 5 or M % 2 == 0:
        raise ValueError("M must be an odd integer >= 5 so that there is a centre guide")
    return int(M) // 2


def wall_couplings(M: int, k_next_to_centre: float, k_other: float) -> np.ndarray:
    """Couplings of a symmetric domain wall: bonds at distance 0, 2, 4, ... from the
    centre get k_next_to_centre, bonds at distance 1, 3, ... get k_other."""
    c0 = centre(M)
    j = np.arange(M - 1)                       # bond j joins guides j and j+1
    dist = np.where(j < c0, c0 - 1 - j, j - c0)
    return np.where(dist % 2 == 0, float(k_next_to_centre), float(k_other))


def topological_interface(M: int, k1: float, k2: float, onsite=None) -> np.ndarray:
    """T array: weak k1 on both sides of the centre ("long-long" gaps). Requires k1 < k2."""
    if not 0 < k1 < k2:
        raise ValueError("the long-long interface needs 0 < k1 < k2")
    return hamiltonian_from_couplings(wall_couplings(M, k1, k2), onsite=onsite)


def short_short_interface(M: int, k1: float, k2: float) -> np.ndarray:
    """SS wall: strong k2 on both sides of the centre ("short-short" gaps)."""
    return hamiltonian_from_couplings(wall_couplings(M, k2, k1))


def defect_array(M: int, kbar: float, delta: float, onsite=None) -> np.ndarray:
    """D array: uniform coupling kbar, centre on-site detuning delta (+ optional extra onsite)."""
    c0 = centre(M)
    d = np.zeros(M)
    d[c0] = delta
    if onsite is not None:
        d = d + np.asarray(onsite, dtype=float)
    return hamiltonian_from_couplings(np.full(M - 1, float(kbar)), onsite=d)


# --- analytic results for the D array (infinite chain) ------------------------------
def defect_q(kbar: float, delta: float) -> float:
    """Amplitude ratio between neighbouring guides of the bound defect mode:
    a_n ~ q^|n|,  q = (sqrt(delta^2 + 4 kbar^2) - |delta|) / (2 kbar),  0 < q < 1."""
    return float((np.sqrt(delta**2 + 4 * kbar**2) - abs(delta)) / (2 * kbar))


def ipr_geometric(q: float) -> float:
    """IPR of a symmetric profile with POWER p_n = q^(2|n|), n = ..., -1, 0, 1, ...:
    IPR = (1 + q^4)(1 - q^2) / (1 + q^2)^3.  The T mode has the same form with q -> r,
    because its zeros on the other sublattice do not change the sums."""
    q2 = q * q
    return float((1 + q2 * q2) * (1 - q2) / (1 + q2) ** 3)


def matched_defect_detuning(k1: float, k2: float, kbar: float | None = None) -> float:
    """Delta such that IPR(D mode) = IPR(T mode) in the infinite-chain limit.
    Since the IPR depends only on q (resp. r), matching IPR <=> q = r = k1/k2, giving
    Delta = kbar (1/r - r).  Solved numerically as a cross-check of that closed form."""
    kbar = 0.5 * (k1 + k2) if kbar is None else kbar
    r = k1 / k2
    target = ipr_geometric(r)
    f = lambda dlt: ipr_geometric(defect_q(kbar, dlt)) - target   # noqa: E731
    return float(brentq(f, 1e-9 * kbar, 1e3 * kbar, xtol=1e-14, rtol=1e-14))


def defect_mode_analytic(kbar: float, delta: float):
    """(b_D, q) of the bound defect mode: b_D = sign(delta) sqrt(delta^2 + 4 kbar^2)."""
    return float(np.sign(delta) * np.sqrt(delta**2 + 4 * kbar**2)), defect_q(kbar, delta)


def interface_mode_analytic(M: int, k1: float, k2: float) -> np.ndarray:
    """Exact normalised zero mode of the T array: centre amplitude 1, then -(k1/k2) per
    cell on the centre's sublattice in both directions, zero on the other sublattice."""
    c0 = centre(M)
    n = np.arange(M) - c0
    psi = np.where(n % 2 == 0, (-(k1 / k2)) ** (np.abs(n) // 2), 0.0)
    return psi / np.linalg.norm(psi)
