from fractions import Fraction

import flint
import mpmath as mp
import pytest
import sympy as sp

import havelock as hv
from havelock.blocks import X, block_exprs, blocks_exact, numerators
from havelock.models import get_model
from havelock.stability import stability

from fd_hessian import block_spectrum, fd_spectrum


def _sym(p: flint.fmpz_poly):
    return sum(sp.Integer(int(c)) * X**i for i, c in enumerate(p.coeffs()))


@pytest.mark.parametrize("N", [3, 5, 8, 9])
def test_numerators_equal_note_theorem_2(N):
    for m in range(N):
        nm = numerators("tanh", N, m)
        a, b, g = block_exprs("tanh", N, m)
        assert sp.cancel(a * _sym(nm.den_a) - _sym(nm.A)) == 0, m
        assert sp.cancel(b * _sym(nm.den_b) - _sym(nm.B)) == 0, m
        assert sp.cancel(g * _sym(nm.den_b) - _sym(nm.G)) == 0, m


def test_mode_zero_is_rotation_and_symmetries():
    N = 9
    assert numerators("tanh", N, 0).B.is_zero() and numerators("tanh", N, 0).G.is_zero()
    for m in range(1, N):
        assert numerators("tanh", N, m).A == numerators("tanh", N, N - m).A
        assert numerators("tanh", N, m).G == -numerators("tanh", N, N - m).G
    assert numerators("tanh", 8, 4).G.is_zero()


@pytest.mark.parametrize("N,m", [(5, 2), (8, 3), (9, 4)])
def test_reciprocity(N, m):
    a, b, g = block_exprs("tanh", N, m)
    inv = {X: 1 / X}
    assert sp.cancel(a.subs(inv) - a) == 0 and sp.cancel(b.subs(inv) - b) == 0 and sp.cancel(g.subs(inv) + g) == 0


def test_model_conversions():
    m = get_model("tanh")
    with mp.workdps(30):
        assert mp.almosteq(m.radius(mp.mpf("0.0673562230040812622684374")), mp.mpf("0.5312104560"), 1e-9)
        assert mp.almosteq(m.cosh_r(mp.mpf("0.0673562230040812622684374")), mp.mpf("1.1444414784"), 1e-9)
    with pytest.raises(hv.InvalidInput):
        m.check_value(Fraction(1))


# Spec 9.2 cross-check for the tanh model. Verified while planning (dps 40, h = 1e-10): absolute errors
# 5.8e-20 (N=5, 0.5), 7.2e-19 (N=8, 0.0673), 3.9e-19 (N=9, 0.2), 2.1e-18 (N=12, 0.9); relative <= 9e-20.
@pytest.mark.parametrize("N,v", [(5, "0.5"), (8, "0.0673"), (9, "0.2"), (12, "0.9")])
def test_fd_hessian_matches_tanh_blocks(N, v):
    mdl = get_model("tanh")
    with mp.workdps(40):
        E, omega, grad = fd_spectrum(mdl, N, mp.mpf(v), mp.mpf("1e-10"))
        P = block_spectrum("tanh", N, mp.mpf(v))
        assert grad < mp.mpf("1e-25")
        scale = max(1, max(abs(p) for p in P))
        assert max(abs(a - b) for a, b in zip(E, P)) < mp.mpf("1e-15") * scale


def test_fd_negative_control_sinh_blocks_do_not_fit_tanh():
    # a Boatto-type modelling error (wrong pair law) must be caught: the sinh blocks at the same xi do not match
    mdl = get_model("tanh")
    with mp.workdps(40):
        E, _, _ = fd_spectrum(mdl, 8, mp.mpf("0.3"), mp.mpf("1e-10"))
        wrong = block_spectrum("mt", 8, mp.mpf("-0.3"))
        # sinh blocks use log|.|^2, i.e. twice the tanh normalisation; recorded mismatch 3.9
        assert max(abs(a - b / 2) for a, b in zip(E, wrong)) > mp.mpf("1e-3")


@pytest.mark.parametrize("N,xi,verdict", [
    (5, "0.5", "EM-stable"), (7, "0.01", "EM-stable"), (8, "0.05", "unstable"), (8, "0.1", "EM-stable"),
    (9, "0.1128710", "unstable"), (9, "0.11287108", "spectrally-stable-only"), (9, "0.1128711", "EM-stable"),
])
def test_tanh_stability(N, xi, verdict):
    s = stability("tanh", N, xi)
    assert s.verdict == verdict and s.status == hv.Status.CERTIFIED


def test_krein_window_is_explained():
    s = stability("tanh", 9, "0.11287108")      # xi*(9) = 0.1128710632... < 0.11287108 < xi_EM(9) = 0.1128710954...
    assert "Krein window" in s.basis
    k = s.blocks[3]
    assert k.a > 0 and k.b > 0 and k.a * k.b - k.g ** 2 < 0
    assert all(abs(e.real) < 1e-25 for e in k.eigenvalues)
