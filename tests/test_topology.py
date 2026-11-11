"""Day-2 validation tests for topq.topology."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from topq.ssh import hamiltonian_from_couplings, ssh_couplings, ssh_hamiltonian  # noqa: E402
from topq.topology import (  # noqa: E402
    chiral_operator, chiral_violation, edge_mode, even_chain_splitting,
    fit_localization_length, ipr, sublattice_weights, winding_number,
    xi_analytic, zero_mode_from_couplings,
)

RNG_SEED = 20261023


# --- chiral symmetry ----------------------------------------------------------
def test_clean_and_coupling_disordered_chains_are_chiral():
    rng = np.random.default_rng(RNG_SEED)
    assert chiral_violation(ssh_hamiltonian(41, 0.5, 1.0)) < 1e-15
    c = ssh_couplings(41, 0.5, 1.0) * (1 + 0.3 * rng.standard_normal(40)).clip(0.05)
    assert chiral_violation(hamiltonian_from_couplings(c)) < 1e-15


@pytest.mark.parametrize("kind", ["random_onsite", "uniform_shift", "nnn"])
def test_symmetry_breaking_perturbations_are_detected(kind):
    rng = np.random.default_rng(RNG_SEED)
    H = ssh_hamiltonian(41, 0.5, 1.0)
    if kind == "random_onsite":
        H = H + np.diag(0.05 * rng.standard_normal(41))
    elif kind == "uniform_shift":
        H = H + 0.05 * np.eye(41)
    else:  # next-nearest-neighbour coupling joins A-A and B-B
        H = H + 0.05 * (np.eye(41, k=2) + np.eye(41, k=-2))
    assert chiral_violation(H) > 1e-3


def test_chiral_partner_has_opposite_eigenvalue():
    H = ssh_hamiltonian(30, 0.4, 1.0)
    G = chiral_operator(30)
    e, v = np.linalg.eigh(H)
    for j in (0, 5, 17):
        w = G @ v[:, j]
        assert np.allclose(H @ w, -e[j] * w, atol=1e-12)


# --- winding number -------------------------------------------------------------
@pytest.mark.parametrize("k1,k2,W", [(0.5, 1.0, 1), (0.2, 1.3, 1), (1.0, 0.5, 0), (2.0, 0.1, 0)])
def test_winding_values(k1, k2, W):
    assert abs(winding_number(k1, k2) - W) < 1e-9
    assert abs(winding_number(k1, k2, nk=301) - W) < 1e-9    # grid-independent


def test_winding_undefined_at_gap_closing():
    assert np.isnan(winding_number(1.0, 1.0))


def test_winding_depends_on_unit_cell_choice():
    """Choosing the cell (B_n, A_{n+1}) swaps k1 <-> k2 and flips W: W is only
    meaningful together with a stated cell / termination."""
    assert abs(winding_number(0.5, 1.0) - 1) < 1e-9
    assert abs(winding_number(1.0, 0.5) - 0) < 1e-9


# --- exact zero mode of odd chains ----------------------------------------------------
@pytest.mark.parametrize("k1,k2", [(0.5, 1.0), (1.0, 0.5), (0.9, 1.0)])
def test_zero_mode_formula_is_exact(k1, k2):
    c = ssh_couplings(41, k1, k2)
    psi = zero_mode_from_couplings(c)
    H = hamiltonian_from_couplings(c)
    assert np.linalg.norm(H @ psi) < 1e-13
    assert sublattice_weights(psi)[1] == 0.0


def test_zero_mode_survives_coupling_disorder_exactly():
    rng = np.random.default_rng(RNG_SEED)
    for _ in range(20):
        c = ssh_couplings(41, 0.5, 1.0) * np.exp(0.3 * rng.standard_normal(40))
        H = hamiltonian_from_couplings(c)
        b, _ = edge_mode(H)
        assert abs(b) < 1e-12 * np.abs(c).max()
        assert np.linalg.norm(H @ zero_mode_from_couplings(c)) < 1e-12


def test_zero_mode_side_depends_on_dimerisation():
    left = zero_mode_from_couplings(ssh_couplings(41, 0.5, 1.0))
    right = zero_mode_from_couplings(ssh_couplings(41, 1.0, 0.5))
    assert np.abs(left[0]) > 0.8 and np.abs(right[-1]) > 0.8


def test_numeric_edge_mode_matches_formula():
    c = ssh_couplings(41, 0.5, 1.0)
    b, v = edge_mode(hamiltonian_from_couplings(c))
    assert abs(abs(v @ zero_mode_from_couplings(c)) - 1) < 1e-12


# --- localisation length ----------------------------------------------------------------
@pytest.mark.parametrize("r", [0.3, 0.5, 0.8])
def test_xi_fit_matches_analytic(r):
    c = ssh_couplings(201, r, 1.0)
    xi_fit = fit_localization_length(zero_mode_from_couplings(c))
    assert abs(xi_fit / xi_analytic(r, 1.0) - 1) < 0.01


def test_xi_analytic_rejects_trivial_side():
    with pytest.raises(ValueError):
        xi_analytic(1.0, 0.5)


# --- even chains: hybridised edge modes --------------------------------------------------
@pytest.mark.parametrize("r,N", [(0.3, 10), (0.5, 20), (0.7, 30)])
def test_even_chain_splitting_formula(r, N):
    e = np.linalg.eigvalsh(ssh_hamiltonian(2 * N, r, 1.0))
    bmin = np.sort(np.abs(e))[0]
    assert abs(bmin / even_chain_splitting(r, 1.0, N) - 1) < 1e-4


# --- bulk-boundary correspondence ---------------------------------------------------------
@pytest.mark.parametrize("r", [0.3, 0.6, 0.9, 1.1, 1.6, 3.0])
def test_bulk_boundary_even_chain(r):
    """Number of in-gap modes of an even open chain = 2 W (one per edge)."""
    k1, k2 = r, 1.0
    N = 100
    e = np.linalg.eigvalsh(ssh_hamiltonian(2 * N, k1, k2))
    n_in_gap = int(np.sum(np.abs(e) < 0.5 * abs(k1 - k2)))
    assert n_in_gap == 2 * round(winding_number(k1, k2))


# --- IPR ----------------------------------------------------------------------------------------
def test_ipr_limits():
    assert ipr(np.eye(10)[3]) == 1.0
    assert abs(ipr(np.ones(10)) - 0.1) < 1e-15


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q", "-p", "no:faulthandler", "--capture=sys", "--color=no"]))
