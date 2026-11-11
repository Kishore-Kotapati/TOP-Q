"""Day-8 tests for topq.mzi."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from topq.mzi import (clements_layout, coupler, embed, fidelity, intensity_error,  # noqa: E402
                      mesh_unitary, mzi, phase, random_mesh_params, unitarity_error)
from topq.propagation import propagator  # noqa: E402


def test_coupler_is_the_coupled_mode_solution():
    """B(eta) = expm(i kappa l sigma_x) with eta = cos^2(kappa l): same physics as Day 1."""
    kl = 0.37
    U = propagator(np.array([[0.0, 1.0], [1.0, 0.0]]), kl)
    assert np.allclose(U, coupler(np.cos(kl) ** 2), atol=1e-12)


@pytest.mark.parametrize("eta", [0.0, 0.3, 0.5, 1.0])
def test_coupler_and_mzi_unitary(eta):
    assert unitarity_error(coupler(eta)) < 1e-12
    assert unitarity_error(mzi(1.1, -0.4, eta, 0.5)) < 1e-12


@pytest.mark.parametrize("theta", np.linspace(0, 2 * np.pi, 9))
def test_ideal_mzi_closed_form(theta):
    T = mzi(theta, 0.7)
    assert abs(abs(T[0, 0]) ** 2 - np.sin(theta / 2) ** 2) < 1e-12
    assert abs(abs(T[1, 0]) ** 2 - np.cos(theta / 2) ** 2) < 1e-12
    # closed form of B P(theta) B (without the external phase)
    inner = coupler() @ phase(theta) @ coupler()
    closed = 1j * np.exp(1j * theta / 2) * np.array([[np.sin(theta / 2), np.cos(theta / 2)],
                                                    [np.cos(theta / 2), -np.sin(theta / 2)]])
    assert np.allclose(inner, closed, atol=1e-12)


def test_imperfect_couplers_limit_extinction():
    """Equal couplers eta1 = eta2 = eta: full bar still reachable, max cross = 4 eta (1 - eta)."""
    eta = 0.45
    cross_max = max(abs(mzi(t, 0, eta, eta)[1, 0]) ** 2 for t in np.linspace(0, 2 * np.pi, 2001))
    bar_min_cross = min(abs(mzi(t, 0, eta, eta)[1, 0]) ** 2 for t in np.linspace(0, 2 * np.pi, 2001))
    assert abs(cross_max - 4 * eta * (1 - eta)) < 1e-6
    assert bar_min_cross < 1e-6


def test_layout_counts():
    for N in (2, 3, 4, 7, 8):
        lay = clements_layout(N)
        assert len(lay) == N * (N - 1) // 2
        assert all(0 <= m < N - 1 for _, m in lay)


def test_mesh_is_unitary_and_matches_manual_product():
    rng = np.random.default_rng(1)
    N = 4
    p = random_mesh_params(N, rng)
    U = mesh_unitary(N, **p)
    assert unitarity_error(U) < 1e-12
    V = np.eye(N, dtype=complex)
    for k, (_, m) in enumerate(clements_layout(N)):
        V = embed(mzi(p["theta"][k], p["phi"][k]), N, m) @ V
    assert np.allclose(np.diag(np.exp(1j * p["out"])) @ V, U, atol=1e-12)


def test_fidelity_properties():
    rng = np.random.default_rng(2)
    U0 = mesh_unitary(5, **random_mesh_params(5, rng))
    assert abs(fidelity(U0, U0) - 1) < 1e-12
    assert abs(fidelity(U0, 0.7 * np.exp(1j * 1.3) * U0) - 1) < 1e-12   # global phase, uniform loss
    U1 = mesh_unitary(5, **random_mesh_params(5, rng))
    assert 0 <= fidelity(U0, U1) < 1
    assert intensity_error(U0, U0) < 1e-12 and 0 < intensity_error(U0, U1) <= 1


@pytest.mark.parametrize("N", [4, 8])
def test_phase_error_first_order_fidelity(N):
    """Independent phase errors sigma on all N^2 phase shifters: 1 - F ~ (N - 1) sigma^2."""
    rng = np.random.default_rng(3)
    p = random_mesh_params(N, rng)
    U0 = mesh_unitary(N, **p)
    s, n = 0.01, len(clements_layout(N))
    inf = [1 - fidelity(U0, mesh_unitary(N, p["theta"] + s * rng.standard_normal(n),
                                         p["phi"] + s * rng.standard_normal(n),
                                         p["out"] + s * rng.standard_normal(N))) for _ in range(2000)]
    assert abs(np.mean(inf) / ((N - 1) * s**2) - 1) < 0.1


@pytest.mark.parametrize("N", [4, 8])
def test_splitting_error_first_order_fidelity(N):
    """Independent coupler errors eta = 1/2 + eps (std s) on all N(N-1) couplers:
    1 - F ~ 2 (N - 1) s^2."""
    rng = np.random.default_rng(4)
    p = random_mesh_params(N, rng)
    U0 = mesh_unitary(N, **p)
    s, n = 0.01, len(clements_layout(N))
    inf = [1 - fidelity(U0, mesh_unitary(N, **p, eta1=0.5 + s * rng.standard_normal(n),
                                         eta2=0.5 + s * rng.standard_normal(n))) for _ in range(2000)]
    assert abs(np.mean(inf) / (2 * (N - 1) * s**2) - 1) < 0.1


def test_bad_eta_rejected():
    with pytest.raises(ValueError):
        coupler(1.2)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q", "-p", "no:faulthandler", "--capture=sys", "--color=no"]))
