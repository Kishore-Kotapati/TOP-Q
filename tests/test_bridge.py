"""Days 9-10 tests for topq.bridge."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from topq.bridge import arm_pool, mesh_with_arm_errors, predicted_infidelity, run_bridge  # noqa: E402
from topq.disorder import rng_for  # noqa: E402
from topq.mzi import clements_layout, fidelity, mesh_unitary, random_mesh_params  # noqa: E402

L = 5 * np.pi


@pytest.mark.parametrize("el", ["T", "D", "W"])
def test_clean_arm_is_identity(el):
    t = arm_pool(el, "onsite", 0.0, 3, 1, L, M=41)
    assert np.allclose(t, 1, atol=1e-12)


def test_plain_waveguide_ignores_coupling_errors():
    assert np.all(arm_pool("W", "coupling", 0.2, 10, 1, L) == 1)


def test_T_arm_phase_immune_to_coupling_errors():
    t = arm_pool("T", "coupling", 0.05, 50, 2, L, M=41)
    assert np.max(np.abs(np.angle(t))) < 1e-12


def test_unit_arms_reproduce_target():
    N = 5
    p = random_mesh_params(N, rng_for(3))
    U = mesh_with_arm_errors(N, p, np.ones((len(clements_layout(N)), 2)))
    assert np.allclose(U, mesh_unitary(N, **p), atol=1e-12)


def test_arm_amplitude_error_counts_like_phase_error():
    """A log-amplitude error of size s costs the same fidelity as a phase error of size s."""
    N, s = 6, 0.01
    n = len(clements_layout(N))
    rng = np.random.default_rng(4)
    p = random_mesh_params(N, rng)
    U0 = mesh_unitary(N, **p)
    f_amp = np.mean([1 - fidelity(U0, mesh_with_arm_errors(N, p, np.exp(s * rng.standard_normal((n, 2)))))
                     for _ in range(1500)])
    f_ph = np.mean([1 - fidelity(U0, mesh_with_arm_errors(N, p, np.exp(1j * s * rng.standard_normal((n, 2)))))
                    for _ in range(1500)])
    assert abs(f_amp / f_ph - 1) < 0.1
    assert abs(f_ph / ((N - 1) ** 2 / N * s**2) - 1) < 0.1


def test_run_bridge_matches_first_order_at_small_sigma():
    r = run_bridge("D", "coupling", 0.003, 4, 150, 5, L, M=41)
    assert abs(r["inf"].mean() / r["pred"] - 1) < 0.15


def test_predicted_infidelity_zero_for_perfect_arms():
    assert predicted_infidelity(6, np.ones(100)) == 0.0


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q", "-p", "no:faulthandler", "--capture=sys", "--color=no"]))
