"""Days 9-10: the BRIDGE between array physics and circuit performance.

Architecture (stated assumption A-bridge): in every MZI of a rectangular mesh, each of the
two internal arms contains one copy of the same 'arm element', and each copy has its own
independent fabrication errors. The arm element is one of
    "T" : topological interface array (Day 3), light in/out through the centre guide
    "D" : matched conventional-defect array (Day 3), same ports
    "W" : a plain single waveguide of the same length (no neighbours at all)
The arm acts on the light as a complex transmission t = a_out / a_out(clean):
its phase is the arm phase error, its modulus the arm amplitude error.
The clean insertion loss |a_out(clean)|^2 (light left behind in other guides of the array)
is common to every arm, so it does not change the normalised fidelity; it is reported
separately because on a real chip it is lost light.
Couplers can additionally carry gap errors: sigma_eta = (pi/4) sigma_kappa (Day 8).
"""
from __future__ import annotations

import numpy as np

from .disorder import rng_for
from .montecarlo import clean_amplitudes, make_spec, run_propagation_ensemble
from .mzi import clements_layout, fidelity, intensity_error, mesh_unitary, random_mesh_params

__all__ = ["arm_pool", "clean_insertion_loss", "mesh_with_arm_errors", "predicted_infidelity", "run_bridge"]


def arm_pool(element: str, kind: str, sigma: float, n: int, seed: int, L: float,
             M: int = 81, ratio: float = 1.0) -> np.ndarray:
    """n independent arm transmissions t (complex, clean value = 1)."""
    if element in ("T", "D"):
        spec = make_spec(element, M=M)
        r = run_propagation_ensemble(spec, kind, sigma, n, seed, [L], ratio)["guide"][:, 0]
        return r / clean_amplitudes(spec, [L])["guide"][0]
    if element == "W":
        if kind == "coupling":
            return np.ones(n, dtype=complex)            # a lone guide has no couplings
        sb = sigma * (ratio if kind == "width" else 1.0)
        g = np.array([rng_for(seed, s).standard_normal() for s in range(n)])
        return np.exp(1j * sb * g * L)
    raise ValueError("element must be 'T', 'D' or 'W'")


def clean_insertion_loss(element: str, L: float, M: int = 81) -> float:
    """Fraction of the input power that reaches the arm's output guide in a clean arm."""
    if element == "W":
        return 1.0
    return float(np.abs(clean_amplitudes(make_spec(element, M=M), [L])["guide"][0]) ** 2)


def mesh_with_arm_errors(N: int, params: dict, t_arms: np.ndarray, eta1=None, eta2=None):
    """Mesh transfer matrix with arm transmissions t_arms of shape (n_MZI, 2)."""
    return mesh_unitary(N, params["theta"], params["phi"], params["out"],
                        eta1=eta1, eta2=eta2, arm_amp=t_arms)


def predicted_infidelity(N: int, t_pool: np.ndarray, sigma_eta: float = 0.0) -> float:
    """First-order prediction (Day 8 rules): every arm error acts on one mode, and a small
    log-amplitude error counts exactly like a phase error of the same size:
        1 - F ~ (N-1)^2/N * (var(arg t) + var(ln|t|))   [N(N-1) arms]
              + 2 (N-1) sigma_eta^2                     [N(N-1) couplers]."""
    var_arm = np.var(np.angle(t_pool)) + np.var(np.log(np.abs(t_pool)))
    return float((N - 1) ** 2 / N * var_arm + 2 * (N - 1) * sigma_eta**2)


def run_bridge(element: str, kind: str, sigma: float, N: int, n_mesh: int, seed: int, L: float,
               with_couplers: bool = False, n_targets: int = 4, M: int = 81, ratio: float = 1.0,
               pool: np.ndarray | None = None) -> dict:
    """Monte Carlo over meshes. Returns per-mesh 1-F, intensity error, mean transmitted power,
    plus the first-order prediction from the same arm pool. Pass `pool` to reuse arm samples
    (must hold at least n_mesh * n_MZI * 2 values)."""
    n_mzi = len(clements_layout(N))
    if pool is None:
        pool = arm_pool(element, kind, sigma, n_mesh * n_mzi * 2, seed, L, M, ratio)
    pool = pool[: n_mesh * n_mzi * 2]
    arms = pool.reshape(n_mesh, n_mzi, 2)
    s_eta = (np.pi / 4) * sigma if (with_couplers and kind in ("coupling", "width")) else 0.0
    inf, ie, pw = np.empty(n_mesh), np.empty(n_mesh), np.empty(n_mesh)
    for k in range(n_mesh):
        params = random_mesh_params(N, rng_for(seed + 1, k % n_targets))
        U0 = mesh_unitary(N, **params)
        e1 = e2 = None
        if s_eta > 0:
            rng = rng_for(seed + 2, k)
            e1 = np.clip(0.5 + s_eta * rng.standard_normal(n_mzi), 0, 1)
            e2 = np.clip(0.5 + s_eta * rng.standard_normal(n_mzi), 0, 1)
        U = mesh_with_arm_errors(N, params, arms[k], e1, e2)
        inf[k] = 1 - fidelity(U0, U)
        ie[k] = intensity_error(U0, U)
        pw[k] = np.real(np.trace(U.conj().T @ U)) / N
    return {"inf": inf, "ie": ie, "power": pw,
            "pred": predicted_infidelity(N, pool, s_eta),
            "arm_phase_std": float(np.std(np.angle(pool))),
            "arm_logamp_std": float(np.std(np.log(np.abs(pool))))}
