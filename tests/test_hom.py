"""Day-11 tests for topq.hom."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from topq.hom import (bs_visibility, coincidence, fock_coincidence_bs, gaussian_overlap_sq,  # noqa: E402
                      jitter_averaged_overlap_sq, purity_from_spectrum, visibility)
from topq.mzi import coupler, phase  # noqa: E402


@pytest.mark.parametrize("kl", [0.2, np.pi / 4, 0.7, 1.2])
def test_fock_space_agrees_with_formula(kl):
    eta = np.cos(kl) ** 2
    assert abs(fock_coincidence_bs(kl) - coincidence(coupler(eta), 1.0)) < 1e-12
    assert abs(coincidence(coupler(eta), 1.0) - (2 * eta - 1) ** 2) < 1e-12


def _first_quantised_pc(T, psi1, psi2):
    """Independent check with explicit symmetrised two-photon states (path x internal)."""
    d = psi1.size
    def ket(path, psi):
        v = np.zeros(2 * d, dtype=complex)
        v[path * d:(path + 1) * d] = psi
        return v
    s = (np.kron(ket(0, psi1), ket(1, psi2)) + np.kron(ket(1, psi2), ket(0, psi1))) / np.sqrt(2)
    U1 = np.kron(T, np.eye(d))
    s = np.kron(U1, U1) @ s
    P0 = np.kron(np.diag([1, 0]), np.eye(d))
    P1 = np.kron(np.diag([0, 1]), np.eye(d))
    Pi = np.kron(P0, P1) + np.kron(P1, P0)
    return float(np.real(np.vdot(s, Pi @ s)))


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_general_formula_against_first_quantised_calculation(seed):
    rng = np.random.default_rng(seed)
    d = 4
    psi1 = rng.standard_normal(d) + 1j * rng.standard_normal(d); psi1 /= np.linalg.norm(psi1)
    psi2 = rng.standard_normal(d) + 1j * rng.standard_normal(d); psi2 /= np.linalg.norm(psi2)
    T = coupler(0.37) @ phase(0.8) @ coupler(0.55)
    J2 = abs(np.vdot(psi1, psi2)) ** 2
    assert abs(_first_quantised_pc(T, psi1, psi2) - coincidence(T, J2)) < 1e-12


def test_bs_visibility_closed_form_and_limits():
    for eta in (0.5, 0.4, 0.2):
        assert abs(visibility(coupler(eta), 0.7) - bs_visibility(eta, 0.7)) < 1e-12
    assert abs(bs_visibility(0.5, 1.0) - 1) < 1e-15
    assert bs_visibility(0.5, 0.0) == 0.0
    assert abs(bs_visibility(0.3) - bs_visibility(0.7)) < 1e-15        # R <-> T symmetric


def test_gaussian_overlap_limits():
    assert gaussian_overlap_sq(0.0, 2.0) == 1.0
    assert abs(gaussian_overlap_sq(1.0, 1.0) - np.exp(-1)) < 1e-15
    # jitter -> 0 recovers the pure result; jitter averaging checked by Monte Carlo
    assert abs(jitter_averaged_overlap_sq(0.3, 1.5, 0.0) - gaussian_overlap_sq(0.3, 1.5)) < 1e-15
    x = np.random.default_rng(3).normal(0, 0.4, 400000)
    mc = np.mean(gaussian_overlap_sq(0.2 + x, 1.5))
    assert abs(mc - jitter_averaged_overlap_sq(0.2, 1.5, 0.4)) < 2e-3


def test_purity():
    assert purity_from_spectrum([1]) == 1.0
    assert abs(purity_from_spectrum([1, 1, 1, 1]) - 0.25) < 1e-15


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q", "-p", "no:faulthandler", "--capture=sys", "--color=no"]))
