"""Day-3 validation tests: comparison arrays (topq.arrays) and disorder (topq.disorder)."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from topq.arrays import (  # noqa: E402
    centre, defect_array, defect_mode_analytic, interface_mode_analytic, ipr_geometric,
    matched_defect_detuning, short_short_interface, topological_interface, wall_couplings,
)
from topq.disorder import (  # noqa: E402
    coupling_disorder, first_order_shift, first_order_std, global_shift,
    onsite_disorder, rng_for, width_disorder,
)
from topq.ssh import hamiltonian_from_couplings  # noqa: E402
from topq.topology import chiral_violation, edge_mode, ipr  # noqa: E402

K1, K2, M = 1 / 3, 1.0, 41
KBAR = 0.5 * (K1 + K2)
SEED = 20261025


def top_mode(H):
    e, v = np.linalg.eigh(H)
    i = int(np.argmax(e))
    return e[i], v[:, i]


# --- array geometry -------------------------------------------------------------------------
def test_wall_pattern_long_long():
    c = wall_couplings(M, K1, K2)
    c0 = centre(M)
    assert c[c0 - 1] == K1 and c[c0] == K1            # weak bonds on both sides of the centre
    assert c[c0 - 2] == K2 and c[c0 + 1] == K2
    assert np.array_equal(c, c[::-1])                 # mirror symmetric


def test_centre_requires_odd_M():
    with pytest.raises(ValueError):
        centre(40)


def test_T_has_exact_centre_zero_mode():
    H = topological_interface(M, K1, K2)
    assert chiral_violation(H) == 0.0
    b, psi = edge_mode(H)
    assert abs(b) < 1e-14
    assert abs(abs(psi @ interface_mode_analytic(M, K1, K2)) - 1) < 1e-12
    assert np.argmax(psi**2) == centre(M)


def test_SS_wall_is_also_a_domain_wall_in_NN_model():
    """Couplings-only SS wall: still chiral, still has near-zero modes; none of them sits
    ON the centre guide (they live on its neighbours / the ends)."""
    H = short_short_interface(M, K1, K2)
    assert chiral_violation(H) == 0.0
    e, v = np.linalg.eigh(H)
    near_zero = np.abs(e) < 1e-3                       # wall + two end modes, weakly hybridised
    assert near_zero.sum() == 3
    assert np.all(v[centre(M), near_zero] ** 2 < 1e-9)   # the centre guide stays dark


# --- matched trivial defect -----------------------------------------------------------------
def test_matched_detuning_closed_form():
    assert abs(matched_defect_detuning(K1, K2) - KBAR * (K2 / K1 - K1 / K2)) < 1e-12


def test_D_mode_matches_analytic_and_IPR():
    delta = matched_defect_detuning(K1, K2)
    b, psi = top_mode(defect_array(M, KBAR, delta))
    b_an, q = defect_mode_analytic(KBAR, delta)
    assert abs(b - b_an) < 1e-10
    assert abs(q - K1 / K2) < 1e-12
    bT, psiT = edge_mode(topological_interface(M, K1, K2))
    assert abs(ipr(psi) - ipr(psiT)) < 1e-8
    assert abs(ipr(psi) - ipr_geometric(q)) < 1e-8


def test_D_is_not_chiral():
    assert chiral_violation(defect_array(M, KBAR, 1.0)) > 0.1


# --- disorder generators --------------------------------------------------------------------
def test_rng_reproducible_and_independent():
    a = rng_for(SEED, 7).standard_normal(5)
    b = rng_for(SEED, 7).standard_normal(5)
    c = rng_for(SEED, 8).standard_normal(5)
    assert np.array_equal(a, b) and not np.allclose(a, c)


def test_coupling_disorder_positive_lognormal_and_chiral():
    c = wall_couplings(M, K1, K2)
    rng = rng_for(SEED)
    samples = np.array([coupling_disorder(c, 0.2, rng) for _ in range(4000)])
    assert np.all(samples > 0)
    logs = np.log(samples / c)
    assert abs(logs.mean()) < 0.005 and abs(logs.std() - 0.2) < 0.005
    assert chiral_violation(hamiltonian_from_couplings(samples[0])) == 0.0


def test_onsite_and_global_break_chirality():
    rng = rng_for(SEED)
    c = wall_couplings(M, K1, K2)
    assert chiral_violation(hamiltonian_from_couplings(c, onsite=onsite_disorder(M, 0.05, rng))) > 1e-3
    g = global_shift(M, 0.05, rng)
    assert np.all(g == g[0]) and chiral_violation(hamiltonian_from_couplings(c, onsite=g)) > 1e-3


def test_width_model_is_correlated():
    """A wide guide raises its own beta AND both neighbouring couplings."""
    c = np.ones(M - 1)
    rng = rng_for(SEED)
    corr = []
    for _ in range(2000):
        cn, d = width_disorder(c, 0.1, 0.1, rng)
        corr.append(np.corrcoef(d[1:-1], np.log(cn[1:]) + np.log(cn[:-1]))[0, 1])
    assert np.mean(corr) > 0.6


def test_zero_sigma_is_identity():
    c = wall_couplings(M, K1, K2)
    rng = rng_for(SEED)
    assert np.array_equal(coupling_disorder(c, 0.0, rng), c)
    cn, d = width_disorder(c, 0.0, 0.0, rng)
    assert np.array_equal(cn, c) and np.all(d == 0)


# --- first-order perturbation theory --------------------------------------------------------
def test_first_order_coupling_shift_of_T_mode_is_zero():
    """psi_T lives on one sublattice, so every psi_j psi_{j+1} = 0."""
    psi = interface_mode_analytic(M, K1, K2)
    dk = rng_for(SEED).standard_normal(M - 1)
    assert first_order_shift(psi, d_couplings=dk) == 0.0


def test_first_order_matches_exact_for_small_perturbation():
    delta = matched_defect_detuning(K1, K2)
    H0 = defect_array(M, KBAR, delta)
    b0, psi = top_mode(H0)
    d = 1e-6 * rng_for(SEED).standard_normal(M)
    b1, _ = top_mode(H0 + np.diag(d))
    assert abs((b1 - b0) - first_order_shift(psi, d_onsite=d)) < 1e-10


def test_onsite_first_order_std_equals_sigma_sqrt_ipr():
    psi = interface_mode_analytic(M, K1, K2)
    s = first_order_std(psi, wall_couplings(M, K1, K2), sigma_beta=0.01)
    assert abs(s - 0.01 * np.sqrt(ipr(psi))) < 1e-15


def test_width_first_order_std_against_monte_carlo():
    """Small-sigma Monte Carlo of the exact top-mode shift agrees with the linear prediction."""
    from topq.disorder import first_order_std_width
    delta = matched_defect_detuning(K1, K2)
    c = np.full(M - 1, KBAR)
    H0 = defect_array(M, KBAR, delta)
    b0, psi = top_mode(H0)
    sb = sk = 1e-4
    pred = first_order_std_width(psi, c, sb, sk)
    db = []
    for s in range(3000):
        cn, d = width_disorder(c, sb, sk, rng_for(SEED, s))
        d[centre(M)] += delta
        db.append(top_mode(hamiltonian_from_couplings(cn, onsite=d))[0] - b0)
    assert abs(np.std(db) / pred - 1) < 0.05          # 3000 samples: ~1.3 % statistical error


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q", "-p", "no:faulthandler", "--capture=sys", "--color=no"]))
