"""Day-4 tests for topq.montecarlo."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from topq.montecarlo import bootstrap_std, clean_mode, make_spec, run_ensemble  # noqa: E402


@pytest.mark.parametrize("name", ["T", "D"])
def test_zero_disorder_gives_clean_mode(name):
    r = run_ensemble(make_spec(name), "onsite", 0.0, 5, seed=1)
    assert np.allclose(r["b"], 0, atol=1e-13) and np.allclose(r["overlap"], 1, atol=1e-12)


def test_clean_mode_values():
    assert abs(clean_mode(make_spec("T"))[0]) < 1e-14
    assert abs(clean_mode(make_spec("D"))[0] - 20 / 9) < 1e-10


def test_ensemble_is_reproducible():
    a = run_ensemble(make_spec("D"), "width", 0.1, 50, seed=42)
    b = run_ensemble(make_spec("D"), "width", 0.1, 50, seed=42)
    assert all(np.array_equal(a[k], b[k]) for k in a)


def test_T_pinned_under_coupling_disorder():
    r = run_ensemble(make_spec("T"), "coupling", 0.3, 300, seed=3)
    assert np.max(np.abs(r["b"])) < 1e-12


def test_exact_zero_mode_exists_even_for_strong_disorder():
    """At sigma = 0.6 the TRACKED mode is sometimes a different one, but a b = 0 mode
    always exists (counting argument). Pinning != staying at the launch site."""
    from topq.disorder import rng_for
    from topq.montecarlo import disordered_H
    spec = make_spec("T")
    for s in range(200):
        e = np.linalg.eigvalsh(disordered_H(spec, "coupling", 0.6, rng_for(3, s)))
        assert np.min(np.abs(e)) < 1e-12


def test_global_shift_moves_every_mode_equally():
    r = run_ensemble(make_spec("T"), "global", 0.1, 200, seed=5)
    assert abs(r["b"].std() - 0.1) < 0.015        # b follows the common offset 1:1
    assert np.allclose(r["overlap"], 1, atol=1e-12)  # shape untouched


def test_bootstrap_ci_brackets_std():
    x = np.random.default_rng(0).standard_normal(2000)
    sd, lo, hi = bootstrap_std(x)
    assert lo < sd < hi and lo < 1.0 < hi



def test_propagation_ensemble_matches_direct_propagation():
    from topq.arrays import centre
    from topq.disorder import rng_for
    from topq.montecarlo import disordered_H, run_propagation_ensemble
    from topq.propagation import launch, propagate
    spec = make_spec("D", M=41)
    Ls = [2.0, 7.5]
    r = run_propagation_ensemble(spec, "width", 0.1, 3, seed=9, lengths=Ls)
    _, psi0 = clean_mode(spec)
    for s in range(3):
        a = propagate(disordered_H(spec, "width", 0.1, rng_for(9, s)), launch(41, centre(41)), Ls)
        assert np.allclose(r["guide"][s], a[:, centre(41)], atol=1e-12)
        assert np.allclose(r["modal"][s], a @ psi0, atol=1e-12)


def test_T_differential_phase_zero_under_coupling_disorder():
    from topq.montecarlo import clean_amplitudes, run_propagation_ensemble
    spec = make_spec("T", M=41)
    Ls = [5 * np.pi]
    r = run_propagation_ensemble(spec, "coupling", 0.1, 40, seed=4, lengths=Ls)
    ref = clean_amplitudes(spec, Ls)
    for q in ("modal", "guide"):
        assert np.max(np.abs(np.angle(r[q][:, 0] / ref[q][0]))) < 1e-12


def test_phase_stats_limits():
    from topq.montecarlo import phase_stats
    assert phase_stats(np.zeros(10))["R"] == 1.0
    x = 0.2 * np.random.default_rng(1).standard_normal(200000)
    st = phase_stats(x)
    assert abs(st["R"] - np.exp(-0.02)) < 2e-3 and abs(st["circ_std"] - 0.2) < 5e-3
    assert phase_stats(np.full(5, np.pi))["p_flip"] == 1.0


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q", "-p", "no:faulthandler", "--capture=sys", "--color=no"]))
