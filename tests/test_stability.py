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
def test_domain_edges(model, N, value, verdict):
    # values on a domain edge get exact verdicts
    assert stability(model, N, value).verdict == verdict


@pytest.mark.parametrize("model,N,verdict", [("mt", 3, "spectrally-stable-only"), ("sphere", 2, "spectrally-stable-only"),
                                             ("mt", 4, "unstable"), ("sphere", 7, "unstable")])
def test_equator_is_not_labelled_proved(model, N, verdict):
    # lam_r2 = 1 lies outside the range |b| < 1 of the paper's Propositions 3 and 4
    s = stability(model, N, 1)
    assert s.verdict == verdict and s.status == hv.Status.CERTIFIED
    assert "equator" in s.basis and "not a paper theorem" in s.basis
    assert "see havelock.nonlinear" not in s.basis and s.blocks[0].a == 0
    assert stability(model, N, "0.999").status == hv.Status.PROVED


def test_float_near_threshold_is_decided_exactly():
    # float(b_8) is not b_8; its exact binary value decides the verdict, never "degenerate"
    import sympy as sp
    from havelock.models.mt import bn_exact
    b = float(sp.N(bn_exact(8), 30))
    s = stability("sinh", 8, b)
    exact_below = sp.Rational(*Fraction(b).as_integer_ratio()) < bn_exact(8)
    assert s.verdict == ("EM-stable" if exact_below else "unstable")


def _annotate_cases():
    from havelock.status import ModeBlock
    F = Fraction
    e = (0j, 0j)
    return {
        "b0": [ModeBlock(m=2, a=F(3), b=F(0), g=F(0), eigenvalues=e)],
        "neg": [ModeBlock(m=2, a=F(-2), b=F(-3), g=F(1), eigenvalues=e)],
    }


def test_b_zero_with_a_nonzero_is_explained():
    from havelock.stability import _annotate
    out = _annotate("BASIS", _annotate_cases()["b0"], "spectrally-stable-only")
    assert out.startswith("BASIS") and "b_m = 0 exactly for m = [2]" in out and "a_m != 0" in out
    assert "Krein window" not in out and "a_m = 0 exactly" not in out


def test_negative_definite_block_is_explained():
    from havelock.stability import _annotate
    out = _annotate("BASIS", _annotate_cases()["neg"], "spectrally-stable-only")
    assert "negative definite for m = [2]" in out and "not an energy minimum" in out
    assert "Krein window" not in out and "b_m = 0" not in out


def test_unstable_verdict_gets_no_annotation():
    from havelock.stability import _annotate
    assert _annotate("BASIS", _annotate_cases()["neg"], "unstable") == "BASIS"


def test_mt_views_report_r():
    import mpmath as mp
    s = stability("sinh", 9, "-0.3")
    assert "r" in s.to_dict() and stability("plane", 7, 0).r == 1
    with mp.workdps(50):
        assert abs(s.r - 2 * mp.atanh(mp.sqrt(mp.mpf(3) / 10))) < mp.mpf("1e-40")
        assert abs(stability("sphere", 4, "1/4").r - 2 * mp.atan(mp.mpf(1) / 2)) < mp.mpf("1e-40")
