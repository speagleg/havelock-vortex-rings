import mpmath as mp
import pytest
import sympy as sp

import havelock as hv
from havelock.threshold import threshold

# paper Table c4, column b_n (10 decimals), n = 4..25
B_TABLE = {4: "0.2679491924", 5: "0.1715728753", 6: "0.0557280900", 7: "0", 8: "-0.0627460668",
           9: "-0.1010205144", 10: "-0.1428571429", 11: "-0.1715728753", 12: "-0.2024108640",
           13: "-0.2251482266", 14: "-0.2492775954", 15: "-0.2679491924", 16: "-0.2875950927",
           17: "-0.3033370453", 18: "-0.3197938296", 19: "-0.3333333333", 20: "-0.3474164357",
           21: "-0.3592455180", 22: "-0.3714998241", 23: "-0.3819660113", 24: "-0.3927725152",
           25: "-0.4021298312"}


@pytest.mark.parametrize("n", sorted(B_TABLE))
def test_mt_threshold_matches_paper_table(n):
    t = threshold("mt", n)
    assert isinstance(t, hv.Threshold) and t.status == hv.Status.PROVED
    assert t.mode == n // 2 and t.stable_side == "below" and t.fallback_from is None
    assert abs(t.value - mp.mpf(B_TABLE[n])) < mp.mpf("5e-11")
    assert "eq. (2)" in t.basis and "Havelock.Conj13" in t.basis


def test_exact_forms_and_geometry():
    t4 = threshold("sphere", 4)
    assert sp.simplify(t4.exact - (2 - sp.sqrt(3))) == 0 and t4.geometry == "sphere"
    assert t4.xi is None and t4.cosh_r is None
    with mp.workdps(50):
        assert abs(t4.r - mp.acos(1 / mp.sqrt(3))) < mp.mpf(10) ** -45
    t7 = threshold("plane", 7)
    assert t7.value == 0 and t7.geometry == "plane" and t7.r == 1
    assert t7.stable_side is None and "degenerate view" in t7.basis and threshold("mt", 7).stable_side == "below"
    t8 = threshold("sinh", 8)
    assert sp.simplify(t8.exact + 8 - 3 * sp.sqrt(7)) == 0 and t8.geometry == "H2"
    with mp.workdps(50):
        assert abs(t8.cosh_r ** 2 - mp.mpf(9) / 7) < mp.mpf(10) ** -45 and "q(N)" in t8.basis
    assert abs(t8.r - mp.mpf("0.511875")) < mp.mpf("5e-7")


@pytest.mark.parametrize("model,N,verdict", [
    ("sphere", 9, "unstable at every radius"), ("sphere", 7, "unstable at every radius"),
    ("plane", 6, "stable at every radius"), ("plane", 8, "unstable at every radius"),
    ("sinh", 7, "stable at every radius"), ("sinh", 5, "stable at every radius"),
    ("sinh", 3, "stable at every radius"), ("plane", 2, "stable at every radius"),
])
def test_no_threshold_views(model, N, verdict):
    t = threshold(model, N)
    assert isinstance(t, hv.NoThreshold) and t.verdict == verdict and t.status == hv.Status.PROVED


@pytest.mark.parametrize("model,N", [("mt", 3), ("mt", 2), ("sphere", 3)])
def test_equator_cases_raise(model, N):
    with pytest.raises(hv.OutOfValidityRange, match="equator"):
        threshold(model, N)


@pytest.mark.parametrize("args,kw", [(("mt", 1), {}), (("nope", 8), {}), (("mt", 8.0), {}),
                                     (("mt", 8), {"precision": 10}), (("mt", 8), {"strict": "yes"})])
def test_invalid_input(args, kw):
    with pytest.raises(hv.InvalidInput):
        threshold(*args, **kw)


def test_value_has_50_digits_and_json():
    t = threshold("sinh", 15)
    with mp.workdps(60):
        assert mp.nstr(t.value, 45) == mp.nstr(-(2 - mp.sqrt(3)), 45)
    d = t.to_dict()
    assert d["type"] == "Threshold" and d["value"].startswith("-0.26794919243112270647255365849")


def test_large_n_stays_exact():
    # huge N in the MT family is closed-form, exact and fast
    t = threshold("sinh", 10001)
    assert t.status == hv.Status.PROVED and t.mode == 5000
    with mp.workdps(50):
        assert abs(t.cosh_r ** 2 - mp.mpf(10001 - 3) / 4) < mp.mpf(10) ** -40
