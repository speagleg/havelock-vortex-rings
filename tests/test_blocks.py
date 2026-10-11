from fractions import Fraction

import flint
import mpmath as mp
import pytest
import sympy as sp

import havelock as hv
from havelock.blocks import X, block_exprs, blocks_exact, blocks_mp, numerators
from havelock.models import get_model
from havelock.models.mt import bn_exact

from fd_hessian import block_spectrum, fd_spectrum


def _sym(p: flint.fmpz_poly):
    return sum(sp.Integer(int(c)) * X**i for i, c in enumerate(p.coeffs()))


@pytest.mark.parametrize("N", [2, 3, 4, 7, 8, 13])
def test_mt_numerators_equal_paper_formula(N):
    for m in range(N):
        nm = numerators("mt", N, m)
        a, b, g = block_exprs("mt", N, m)
        assert sp.cancel(a * _sym(nm.den_a) - _sym(nm.A)) == 0
        assert sp.cancel(b * _sym(nm.den_b) - _sym(nm.B)) == 0
        assert g == 0 and nm.G == 0
        assert nm.D == nm.A * nm.B


def test_mode_range_and_family():
    with pytest.raises(hv.InvalidInput):
        numerators("mt", 5, 5)
    with pytest.raises(hv.InvalidInput):
        numerators("kirchhoff", 5, 1)


def test_exact_values_and_threshold_zero():
    # a_m(b) = 2(n-1)(1+b^2)/(1+b)^2 - m(n-m): at b = 0, n = 7: C1 = 12 = floor(49/4), so a_3 = 0
    assert blocks_exact("mt", 7, 3, Fraction(0)) == (0, 12, 0)
    assert blocks_exact("mt", 10, 5, Fraction(-1, 7)) == (0, 25, 0)  # b_10 = -1/7
    assert blocks_exact("mt", 4, 1, Fraction(1, 2)) == (Fraction(10, 3) - 3, 3, 0)
    with mp.workdps(50):
        for n in (4, 5, 8, 9, 16, 25):
            b = mp.mpf(str(sp.N(bn_exact(n), 60)))
            a, bb, g = blocks_mp("mt", n, n // 2, b)
            assert abs(a) < mp.mpf(10) ** -45 and bb == (n // 2) * (n - n // 2)
            assert all(blocks_mp("mt", n, m, b)[0] > 0 for m in range(1, n // 2))


# Independent cross-check (spec 9.2): brute-force FD Hessian of the full 2N system vs the block eigenvalues.
# Verified while planning (dps 40, h = 1e-10): max |error| / max(1, max |eigenvalue|) <= 1.3e-19 for every case
# below (absolute errors 1.8e-20 .. 1.2e-16); the tolerance 1e-15 relative leaves a margin of > 1000.
@pytest.mark.parametrize("model,N,v", [
    ("sphere", 5, "0.3"), ("sphere", 4, "0.9"), ("plane", 7, "0"), ("sinh", 9, "-0.3"), ("mt", 8, "-0.0627"),
    ("sinh", 12, "-0.8"),
])
def test_fd_hessian_matches_blocks(model, N, v):
    mdl = get_model(model)
    with mp.workdps(40):
        E, omega, grad = fd_spectrum(mdl, N, mp.mpf(v), mp.mpf("1e-10"))
        P = block_spectrum(mdl.family, N, mp.mpf(v))
        assert grad < mp.mpf("1e-25")            # F = H - Omega J3 is critical at the ring
        scale = max(1, max(abs(p) for p in P))
        assert max(abs(a - b) for a, b in zip(E, P)) < mp.mpf("1e-15") * scale


def test_fd_negative_control_wrong_momentum_multiplier():
    mdl = get_model("sinh")
    with mp.workdps(40):
        E, _, grad = fd_spectrum(mdl, 9, mp.mpf("-0.3"), mp.mpf("1e-10"), omega_shift=mp.mpf("1e-3"))
        P = block_spectrum("mt", 9, mp.mpf("-0.3"))
        assert grad > mp.mpf("1e-6")
        assert max(abs(a - b) for a, b in zip(E, P)) > mp.mpf("1e-6")
