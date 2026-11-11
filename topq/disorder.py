"""Reproducible disorder generators and first-order (perturbation-theory) predictions.

Disorder models (all strengths are dimensionless or in units of the couplings):
  coupling : kappa_j -> kappa_j * exp(sigma_kappa * g_j)               (keeps chiral symmetry)
             Justification: kappa ~ exp(-gamma * gap); a Gaussian gap error dg gives
             d(ln kappa) = -gamma dg, i.e. EXACTLY log-normal kappa, always positive.
  onsite   : delta_beta_n = sigma_beta * g_n                            (breaks chiral symmetry)
  width    : one width error dw_n per guide (unit-variance Gaussian), used twice:
             delta_beta_n = sigma_beta * dw_n                (wider guide -> larger beta)
             kappa_{n,n+1} -> kappa * exp(sigma_kappa * (dw_n + dw_{n+1}) / 2)
             (centres fixed: wider guides -> narrower gaps -> larger kappa).
             Correlated, mixed: this is the realistic fabrication case (risk R3).
             Neglected: the direct change of kappa with mode confinement (Limitation L5).
  global   : one common offset sigma_global * g added to every beta (wafer-scale error).

Random numbers: rng_for(seed, sample) builds an independent, reproducible generator for
each Monte Carlo sample from numpy's SeedSequence, so results do not depend on the order
in which samples are run.
"""
from __future__ import annotations

import numpy as np

__all__ = ["rng_for", "coupling_disorder", "onsite_disorder", "width_disorder",
           "global_shift", "first_order_shift", "first_order_std", "first_order_std_width"]


def rng_for(seed: int, sample: int = 0) -> np.random.Generator:
    return np.random.default_rng(np.random.SeedSequence([int(seed), int(sample)]))


def coupling_disorder(couplings, sigma_kappa: float, rng: np.random.Generator) -> np.ndarray:
    c = np.asarray(couplings, dtype=float)
    return c * np.exp(float(sigma_kappa) * rng.standard_normal(c.shape))


def onsite_disorder(M: int, sigma_beta: float, rng: np.random.Generator) -> np.ndarray:
    return float(sigma_beta) * rng.standard_normal(int(M))


def width_disorder(couplings, sigma_beta: float, sigma_kappa: float, rng: np.random.Generator):
    """Returns (new_couplings, onsite) from ONE set of width errors."""
    c = np.asarray(couplings, dtype=float)
    dw = rng.standard_normal(c.size + 1)
    onsite = float(sigma_beta) * dw
    c_new = c * np.exp(float(sigma_kappa) * 0.5 * (dw[:-1] + dw[1:]))
    return c_new, onsite


def global_shift(M: int, sigma_global: float, rng: np.random.Generator) -> np.ndarray:
    return np.full(int(M), float(sigma_global) * rng.standard_normal())


# --- first-order perturbation theory (predictions to be tested on Day 4) ---------------
def first_order_shift(psi, d_onsite=None, d_couplings=None) -> float:
    """delta_b = <psi| dH |psi> for real psi:
    sum_n psi_n^2 d_beta_n + 2 sum_j psi_j psi_{j+1} d_kappa_j."""
    psi = np.asarray(psi, dtype=float)
    db = 0.0
    if d_onsite is not None:
        db += float(np.sum(psi**2 * np.asarray(d_onsite)))
    if d_couplings is not None:
        db += float(2 * np.sum(psi[:-1] * psi[1:] * np.asarray(d_couplings)))
    return db


def first_order_std(psi, couplings, sigma_beta: float = 0.0, sigma_kappa: float = 0.0) -> float:
    """Predicted std of delta_b for INDEPENDENT on-site (std sigma_beta) and multiplicative
    coupling (d_kappa_j ~ kappa_j sigma_kappa g_j, small sigma) disorder:
        var = sigma_beta^2 * sum psi^4  +  4 sigma_kappa^2 * sum (kappa_j psi_j psi_{j+1})^2.
    The first sum is the IPR: on-site spread = sigma_beta * sqrt(IPR)."""
    psi = np.asarray(psi, dtype=float)
    c = np.asarray(couplings, dtype=float)
    var = sigma_beta**2 * np.sum(psi**4) + 4 * sigma_kappa**2 * np.sum((c * psi[:-1] * psi[1:]) ** 2)
    return float(np.sqrt(var))


def first_order_std_width(psi, couplings, sigma_beta: float, sigma_kappa: float) -> float:
    """Predicted std of delta_b for the correlated WIDTH model (one dw_n per guide).
    delta_b = sum_n w_n dw_n with
        w_n = sigma_beta psi_n^2 + sigma_kappa (c_{n-1} psi_{n-1} psi_n + c_n psi_n psi_{n+1})."""
    psi = np.asarray(psi, dtype=float)
    c = np.asarray(couplings, dtype=float)
    bond = c * psi[:-1] * psi[1:]                  # c_j psi_j psi_{j+1}
    w = sigma_beta * psi**2
    w[:-1] += sigma_kappa * bond
    w[1:] += sigma_kappa * bond
    return float(np.sqrt(np.sum(w**2)))
