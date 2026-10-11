from fractions import Fraction

import flint
import pytest

from havelock.certify import (IsolatedRoot, canonical_factors, certify_no_threshold, certify_threshold,
                              mod_p_witness, real_root_balls, sign_at, squarefree, strip_roots, sturm_count)
from havelock.models import get_model

F = Fraction
x = flint.fmpz_poly([0, 1])


def test_sturm_counts():
    p = (x**2 - 2) * (x - 1) ** 3 * (3 * x - 1)           # roots -sqrt2, 1/3, 1 (triple), sqrt2
    assert sturm_count(p, F(-2), F(2)) == 4
    assert sturm_count(p, F(0), F(1)) == 2               # (0, 1] contains 1/3 and 1
    assert sturm_count(p, F(1), F(2)) == 1               # (1, 2] contains sqrt2 only
    assert sturm_count(p, F(1, 3), F(1, 2)) == 0          # (1/3, 1/2]: 1/3 excluded
    assert sturm_count(flint.fmpz_poly([5]), F(0), F(1)) == 0


def test_sign_squarefree_strip():
    p = (x - 1) ** 2 * (x + 2)
    assert sign_at(p, F(0)) == 1 and sign_at(p, F(1)) == 0 and sign_at(p, F(-3)) == -1
    assert squarefree(p) == (x - 1) * (x + 2)
    assert strip_roots(p * x**3, (F(0), F(1))) == x + 2


def test_real_root_balls_and_factors():
    p = (x**2 - 2) * (2 * x - 1)
    balls = real_root_balls(p, F(0), F(2), 128)
    assert len(balls) == 2
    (L1, U1), (L2, U2) = balls
    assert L1 <= F(1, 2) <= U1 and L2 ** 2 < 2 < U2 ** 2 and U2 - L2 < F(1, 10**30)
    assert canonical_factors(-p) == [2 * x - 1, x**2 - 2]


def test_mod_p_witness():
    assert mod_p_witness(x**2 - 2) == {"p": 3, "degrees": [2]}          # x^2 - 2 irreducible mod 3
    assert mod_p_witness(x**4 + 1, max_prime=200) is None              # reducible mod every prime


def test_isolated_root_approx():
    r = IsolatedRoot(poly=(-2, 0, 1), lo=F(141421356237, 10**11), hi=F(141421356238, 10**11))
    assert abs(float(r.approx(20)) - 2**0.5) < 1e-10


def test_certify_threshold_mt_family():
    dom = get_model("mt").family_domain
    c = certify_threshold("mt", 8, dom, "below")
    assert c.mode == 4 and c.sturm == 1 and c.endpoint_signs == (1, -1)
    L, U = c.interval
    q8 = x**2 + 16 * x + 1                                 # b_8 = -(8 - 3 sqrt7) is its root near -0.0627
    assert sign_at(q8, L) * sign_at(q8, U) < 0 and -0.06275 < float(L) < -0.06274 and U - L < F(1, 10**60)
    assert c.factor == q8
    assert all(v == {"A": 0, "B": 0} for v in c.stable_side["root_counts"].values())
    assert c.mode_signs[4]["above"]["A"] == -1 and c.mode_signs[4]["below"]["A"] == 1
    c7 = certify_threshold("mt", 7, dom, "below")
    assert c7.factor == x and c7.mode == 3


def test_certify_no_threshold():
    dom = get_model("mt").family_domain
    assert certify_no_threshold("mt", 8, (F(-1), F(-1, 2), False, False)) is True    # stable for b < b_8
    assert certify_no_threshold("mt", 8, dom) is False


def test_isolated_root_str():
    r = IsolatedRoot(poly=(-2, 0, 1), lo=F(1414, 1000), hi=F(1415, 1000))
    assert str(r) == "the root in [1.4139999999999999, 1.415] (width 1.0e-03) of the degree-2 polynomial [-2, 0, 1] (constant term first)"
