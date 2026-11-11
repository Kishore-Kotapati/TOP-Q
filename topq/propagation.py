"""Coupled-mode propagation, phase extraction and loss.

Convention (unchanged): da/dz = i H a, so a(z) = expm(i H z) a(0). The fast common factor
exp(i beta0 z) is NOT included: all phases here are measured relative to the reference
propagation constant beta0 (a common shift of beta0 adds the same phase to every guide).

Loss: amplitude loss alpha_n/2 per unit length in guide n enters as H -> H + i diag(alpha_n/2),
because i * (i alpha/2) = -alpha/2. Power then decays as exp(-alpha z) for uniform alpha.
"""
from __future__ import annotations

import numpy as np
from scipy.linalg import expm

__all__ = ["with_loss", "propagator", "propagate", "launch", "modal_amplitude",
           "guide_field", "power"]


def with_loss(H: np.ndarray, alpha) -> np.ndarray:
    """Return H + i*diag(alpha/2). `alpha` is the POWER loss per unit length (scalar or per guide)."""
    H = np.asarray(H, dtype=complex)
    a = np.broadcast_to(np.asarray(alpha, dtype=float), (H.shape[0],))
    if np.any(a < 0):
        raise ValueError("loss must be non-negative (no gain)")
    return H + 1j * np.diag(a / 2)


def propagator(H: np.ndarray, z: float) -> np.ndarray:
    """U(z) = expm(i H z). Hermitian H: exact eigen-decomposition (fast, unitary to
    round-off). Non-Hermitian H (loss): scipy.linalg.expm (Pade approximation)."""
    H = np.asarray(H)
    if np.allclose(H, H.conj().T, atol=0, rtol=0):
        b, V = np.linalg.eigh(H)
        return (V * np.exp(1j * b * z)) @ V.conj().T
    return expm(1j * H * z)


def propagate(H: np.ndarray, a0, zs) -> np.ndarray:
    """Field at each z in zs; returns array of shape (len(zs), M)."""
    a0 = np.asarray(a0, dtype=complex)
    H = np.asarray(H)
    if np.allclose(H, H.conj().T, atol=0, rtol=0):
        b, V = np.linalg.eigh(H)
        c = V.conj().T @ a0                                  # modal coefficients
        return (np.exp(1j * np.outer(zs, b)) * c) @ V.T
    return np.array([propagator(H, z) @ a0 for z in zs])


def launch(M: int, guide: int) -> np.ndarray:
    """Unit power in one guide (an ideal butt-coupled input)."""
    a = np.zeros(M, dtype=complex)
    a[guide] = 1.0
    return a


def modal_amplitude(psi_ref, a) -> complex:
    """<psi_ref|a>: what an ideal filter for the reference mode would see.
    Its argument is the MODAL PHASE, its |.|^2 the power in that mode."""
    return complex(np.vdot(psi_ref, a))


def guide_field(a, guide: int) -> complex:
    """Complex field in one output guide: its argument is the OUTPUT-GUIDE PHASE."""
    return complex(np.asarray(a)[guide])


def power(a) -> float:
    return float(np.sum(np.abs(a) ** 2))
