"""Day-5 tests for topq.propagation."""
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import solve_ivp

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from topq.arrays import centre, interface_mode_analytic, topological_interface  # noqa: E402
from topq.disorder import onsite_disorder, rng_for  # noqa: E402
from topq.montecarlo import clean_mode, make_spec  # noqa: E402
from topq.propagation import (guide_field, launch, modal_amplitude, power,  # noqa: E402
                              propagate, propagator, with_loss)

M, K1, K2 = 41, 1 / 3, 1.0
H_T = topological_interface(M, K1, K2)


def test_propagator_is_unitary_and_composes():
    U1, U2 = propagator(H_T, 1.3), propagator(H_T, 2.1)
    assert np.allclose(U1.conj().T @ U1, np.eye(M), atol=1e-12)
    assert np.allclose(U2 @ U1, propagator(H_T, 3.4), atol=1e-12)


def test_power_conserved_without_loss():
    a = propagate(H_T, launch(M, centre(M)), np.linspace(0, 30, 61))
    assert np.allclose(np.sum(np.abs(a) ** 2, axis=1), 1.0, atol=1e-12)


@pytest.mark.parametrize("lossy", [False, True])
def test_matches_independent_ode_solver(lossy):
    rng = rng_for(5)
    H = H_T + np.diag(onsite_disorder(M, 0.1, rng))
    if lossy:
        H = with_loss(H, 0.05 * rng.random(M))
    a0 = launch(M, centre(M) + 1)
    z = 7.0
    sol = solve_ivp(lambda _z, a: 1j * (H @ a), (0, z), a0, rtol=1e-11, atol=1e-13, method="DOP853")
    assert np.allclose(sol.y[:, -1], propagate(H, a0, [z])[0], atol=1e-8)


def test_eigenmode_only_picks_up_phase_bz():
    b, psi = clean_mode(make_spec("D"))
    z = 4.2
    a = propagate(make_spec_H("D"), psi, [z])[0]
    assert np.allclose(a, np.exp(1j * b * z) * psi, atol=1e-12)


def make_spec_H(name):
    from topq.ssh import hamiltonian_from_couplings
    s = make_spec(name)
    return hamiltonian_from_couplings(s.couplings, onsite=s.onsite)


def test_zero_mode_has_constant_phase_and_power():
    psi = interface_mode_analytic(M, K1, K2)
    for z in (1.0, 10.0, 50.0):
        amp = modal_amplitude(psi, propagate(H_T, psi, [z])[0])
        assert abs(amp - 1) < 1e-12


def test_uniform_loss_is_pure_exponential_and_phase_free():
    a0 = launch(M, centre(M))
    z, alpha = 9.0, 0.08
    a_l = propagate(with_loss(H_T, alpha), a0, [z])[0]
    a_0 = propagate(H_T, a0, [z])[0]
    assert abs(power(a_l) - np.exp(-alpha * z)) < 1e-10
    assert np.allclose(a_l, np.exp(-alpha * z / 2) * a_0, atol=1e-10)   # same phases


def test_launch_efficiency_into_matched_modes_is_equal():
    """Centre-guide launch puts |psi(centre)|^2 = (1-r^2)/(1+r^2) = 0.8 into either mode."""
    for name in "TD":
        _, psi = clean_mode(make_spec(name))
        assert abs(psi[centre(M)] ** 2 - 0.8) < 1e-8


def test_neighbour_launch_does_not_excite_T_zero_mode():
    psi = interface_mode_analytic(M, K1, K2)
    a = propagate(H_T, launch(M, centre(M) + 1), [5.0])[0]
    assert abs(modal_amplitude(psi, a)) < 1e-12


@pytest.mark.parametrize("sigma", [0.1, 0.6])
def test_self_amplitude_is_real_in_any_chiral_array(sigma):
    """<c|exp(iHz)|c> = sum_m |psi_m(c)|^2 exp(i b_m z); chiral pairs (b, -b) carry equal
    weight on c, so the sum is REAL: the phase can only be 0 or pi."""
    from topq.montecarlo import disordered_H
    spec = make_spec("T")
    H = disordered_H(spec, "coupling", sigma, rng_for(11, 0))
    a = propagate(H, launch(M, centre(M)), np.linspace(0, 20, 41))[:, centre(M)]
    assert np.max(np.abs(a.imag)) < 1e-12


def test_zero_mode_weight_bounds_the_self_amplitude():
    """Re a_c(z) >= 2 w0 - 1, with w0 = zero-mode weight on the launch guide."""
    from topq.montecarlo import disordered_H
    spec = make_spec("T")
    for k in range(20):
        H = disordered_H(spec, "coupling", 0.3, rng_for(12, k))
        e, v = np.linalg.eigh(H)
        w0 = v[centre(M), np.argmin(np.abs(e))] ** 2
        a = propagate(H, launch(M, centre(M)), np.linspace(0, 30, 61))[:, centre(M)]
        assert np.min(a.real) >= 2 * w0 - 1 - 1e-12


def test_random_loss_keeps_T_zero_mode_real_part_pinned():
    """Gamma (H0 + iD) Gamma = -(H0 + iD)^dagger: eigenvalues pair as (b, -b*), so an isolated
    mode has Re b = 0 exactly, for ANY diagonal loss on a chiral H0."""
    rng = rng_for(13)
    for _ in range(50):
        e = np.linalg.eigvals(with_loss(H_T, 0.3 * rng.random(M)))
        b = e[np.argmin(np.abs(e.real))]
        assert abs(b.real) < 1e-12 and b.imag > 0


def test_negative_loss_rejected():
    with pytest.raises(ValueError):
        with_loss(H_T, -0.1)


def test_guide_field_helper():
    a = launch(M, 3)
    assert guide_field(a, 3) == 1 and power(a) == 1


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q", "-p", "no:faulthandler", "--capture=sys", "--color=no"]))
