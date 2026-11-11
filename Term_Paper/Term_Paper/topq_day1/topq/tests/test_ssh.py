"""Day-1 validation tests for topq.ssh.

Each test checks the code against something it was NOT built from:
an analytic formula, a symmetry, or an exactly solvable limit.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

# Make the repo root importable even when this file is run directly (e.g. F5 in IDLE).
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from topq.ssh import (
    band_gap,
    bloch_bands,
    bloch_hamiltonian,
    hamiltonian_from_couplings,
    ssh_couplings,
    ssh_hamiltonian,
)

RATIOS = [(0.5, 1.0), (1.0, 0.5), (0.75, 0.75), (0.2, 1.3)]
TOL = 1e-12


# 1. Ring spectrum == analytic Bloch bands at the allowed momenta ------------
@pytest.mark.parametrize("k1,k2", RATIOS)
def test_ring_matches_analytic_bands(k1, k2):
    N = 50                                   # unit cells
    H = ssh_hamiltonian(2 * N, k1, k2, periodic=True)
    numeric = np.sort(np.linalg.eigvalsh(H))
    km = 2 * np.pi * np.arange(N) / N        # periodic boundary: exp(i k N) = 1
    bm, bp = bloch_bands(km, k1, k2)
    analytic = np.sort(np.concatenate([bm, bp]))
    assert np.max(np.abs(numeric - analytic)) < TOL


# 2. Hermiticity (lossless => real symmetric) ---------------------------------
@pytest.mark.parametrize("periodic", [False, True])
def test_hamiltonian_is_real_symmetric(periodic):
    H = ssh_hamiltonian(40, 0.4, 1.1, periodic=periodic)
    assert np.isrealobj(H)
    assert np.array_equal(H, H.T)


# 3. Spectrum symmetric about b = 0 for an open chain (chiral symmetry) -------
@pytest.mark.parametrize("M", [40, 41])
@pytest.mark.parametrize("k1,k2", RATIOS)
def test_open_chain_spectrum_is_symmetric(M, k1, k2):
    e = np.sort(np.linalg.eigvalsh(ssh_hamiltonian(M, k1, k2)))
    assert np.max(np.abs(e + e[::-1])) < TOL


# 4. Bulk gap = 2|k1 - k2| (ring with even N contains k = pi exactly) --------
@pytest.mark.parametrize("k1,k2", RATIOS)
def test_ring_gap(k1, k2):
    N = 60                                    # even => k = pi is an allowed momentum
    e = np.linalg.eigvalsh(ssh_hamiltonian(2 * N, k1, k2, periodic=True))
    gap = e[e > 1e-14].min() - e[e < -1e-14].max() if k1 != k2 else 0.0
    assert abs(gap - band_gap(k1, k2)) < 1e-10


# 5. Input validation ----------------------------------------------------------
def test_bad_inputs_raise():
    with pytest.raises(ValueError):
        ssh_hamiltonian(10, -0.1, 1.0)            # negative coupling
    with pytest.raises(ValueError):
        ssh_hamiltonian(11, 0.5, 1.0, periodic=True)  # odd ring breaks cells
    with pytest.raises(ValueError):
        ssh_hamiltonian(1, 0.5, 1.0)              # no couplings
    with pytest.raises(ValueError):
        ssh_hamiltonian(4, 0.5, 1.0, onsite=[0, 0, 0])  # wrong onsite length
    with pytest.raises(ValueError):
        bloch_bands(0.0, np.nan, 1.0)


# 6. Bloch matrix: Hermitian, 2pi-periodic, eigenvalues == bands ---------------
def test_bloch_hamiltonian_consistency():
    k1, k2 = 0.6, 1.0
    k = np.linspace(-np.pi, np.pi, 201)
    Hk = bloch_hamiltonian(k, k1, k2)
    assert np.allclose(Hk, np.conj(np.swapaxes(Hk, -1, -2)), atol=TOL)
    assert np.allclose(Hk, bloch_hamiltonian(k + 2 * np.pi, k1, k2), atol=TOL)
    ev = np.linalg.eigvalsh(Hk)
    bm, bp = bloch_bands(k, k1, k2)
    assert np.allclose(ev[:, 0], bm, atol=TOL) and np.allclose(ev[:, 1], bp, atol=TOL)


# 7. Exactly solvable limits (independent of Bloch theory) ---------------------
def test_uniform_open_chain_limit():
    """k1 = k2 = kappa: open uniform chain, b_j = 2 kappa cos(pi j / (M+1))."""
    M, kap = 37, 0.8
    e = np.sort(np.linalg.eigvalsh(ssh_hamiltonian(M, kap, kap)))
    j = np.arange(1, M + 1)
    assert np.allclose(e, np.sort(2 * kap * np.cos(np.pi * j / (M + 1))), atol=TOL)


def test_fully_dimerised_limits():
    """k2 = 0: isolated dimers (+-k1). k1 = 0: dimers (+-k2) plus 2 free end sites at b = 0."""
    N = 10
    e = np.sort(np.linalg.eigvalsh(ssh_hamiltonian(2 * N, 0.7, 0.0)))
    assert np.allclose(e, np.r_[-0.7 * np.ones(N), 0.7 * np.ones(N)], atol=TOL)
    e = np.sort(np.linalg.eigvalsh(ssh_hamiltonian(2 * N, 0.0, 0.9)))
    expected = np.sort(np.r_[-0.9 * np.ones(N - 1), 0.0, 0.0, 0.9 * np.ones(N - 1)])
    assert np.allclose(e, expected, atol=TOL)


def test_couplings_pattern_and_general_builder():
    c = ssh_couplings(6, 0.3, 0.8)
    assert np.array_equal(c, [0.3, 0.8, 0.3, 0.8, 0.3])
    H = hamiltonian_from_couplings(c, onsite=np.arange(6.0))
    assert np.array_equal(np.diag(H), np.arange(6.0))
    assert np.array_equal(np.diag(H, 1), c)


if __name__ == "__main__":  # F5 in IDLE / "python tests/test_ssh.py"
    # IDLE has no real stderr/stdout file descriptors, so disable pytest's
    # faulthandler plugin and use Python-level (not OS-level) output capture.
    raise SystemExit(pytest.main([__file__, "-q", "-p", "no:faulthandler", "--capture=sys"]))
