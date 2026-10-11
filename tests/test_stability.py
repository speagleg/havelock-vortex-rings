from fractions import Fraction

import pytest

import havelock as hv
from havelock.stability import stability


@pytest.mark.parametrize("model,N,value,verdict", [
    ("sinh", 9, "-0.3", "EM-stable"), ("sinh", 9, "-0.05", "unstable"),
    ("sphere", 4, "0.1", "EM-stable"), ("sphere", 4, "0.5", "unstable"),
    ("plane", 6, 0, "EM-stable"), ("plane", 8, 0, "unstable"),
    ("plane", 7, 0, "spectrally-stable-only"),           # b_7 = 0: a_3 = 0 exactly
    ("sinh", 10, "-1/7", "spectrally-stable-only"),       # b_10 = -1/7: a_5 = 0 exactly
    ("mt", 3, 1, "spectrally-stable-only"),               # equator: a_1 = 0
])
def test_verdicts(model, N, value, verdict):
    s = stability(model, N, value)
    assert s.verdict == verdict and s.status == hv.Status.PROVED
    assert [k.m for k in s.blocks] == list(range(1, N // 2 + 1))


def test_degenerate_point_is_explained():
    s = stability("sinh", 10, "-1/7")
    assert "a_m = 0 exactly for m = [5]" in s.basis and s.blocks[4].a == 0 and s.blocks[4].b == 25


def test_blocks_and_eigenvalues():
    s = stability("sinh", 9, Fraction(-1, 20))
    k = s.blocks[3]                                       # m = 4, the critical mode
    assert k.a == Fraction(2 * 8 * (1 + Fraction(1, 400)), 1) / (1 + Fraction(-1, 20)) ** 2 - 20 and k.a < 0
    sig = k.eigenvalues
    assert abs(sig[0].imag) < 1e-30 and sig[0].real > 0 and abs(sig[0] + sig[1]) < 1e-12   # real pair +-sigma
    stable = stability("sinh", 9, "-0.3")
    assert all(abs(e.real) < 1e-30 for k in stable.blocks for e in k.eigenvalues)


def test_float_value_is_exact_binary():
    s = stability("sinh", 9, 0.1 - 0.4)
    assert s.value == Fraction(0.1 - 0.4)


@pytest.mark.parametrize("args", [("sinh", 9, "0.2"), ("sphere", 4, "0"), ("mt", 5, "1.5"), ("mt", 5, "x"),
                                  ("mt", 1, "0.1")])
def test_invalid(args):
    with pytest.raises(hv.InvalidInput):
        stability(*args)


@pytest.mark.parametrize("model,N,value,verdict", [
    ("mt", 4, 1, "unstable"),          # equator: a_2 = 3 - 4 < 0
    ("sphere", 4, "1", "unstable"),
    ("sinh", 9, "-0.999999", "EM-stable"),
])
def test_review_focus_domain_edges(model, N, value, verdict):
    # Review Focus 4: values on a domain edge get exact verdicts
    assert stability(model, N, value).verdict == verdict


def test_review_focus_float_near_threshold_is_decided_exactly():
    # Review Focus 5: float(b_8) is not b_8; its exact binary value decides the verdict, never "degenerate"
    import sympy as sp
    from havelock.models.mt import bn_exact
    b = float(sp.N(bn_exact(8), 30))
    s = stability("sinh", 8, b)
    exact_below = sp.Rational(*Fraction(b).as_integer_ratio()) < bn_exact(8)
    assert s.verdict == ("EM-stable" if exact_below else "unstable")
