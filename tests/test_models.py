from fractions import Fraction

import mpmath as mp
import numpy as np
import pytest
import sympy as sp

import havelock as hv
from havelock.models import check_N, get_model, to_fraction
from havelock.models.mt import bn_exact, c4_closed, sinh_q


def test_registry_and_errors():
    assert get_model("MT").name == "mt"
    with pytest.raises(hv.InvalidInput, match="unknown model"):
        get_model("vortex")
    for bad in (1, 0, -3, 2.0, True, "8"):
        with pytest.raises(hv.InvalidInput):
            check_N(bad)
    assert check_N(np.int64(9)) == 9


@pytest.mark.parametrize("value,expected", [
    (3, Fraction(3)), ("1/7", Fraction(1, 7)), ("-0.3", Fraction(-3, 10)), (0.5, Fraction(1, 2)),
    (sp.Rational(2, 9), Fraction(2, 9)), (mp.mpf("0.25"), Fraction(1, 4)), (Fraction(5, 3), Fraction(5, 3)),
])
def test_to_fraction(value, expected):
    assert to_fraction(value) == expected


@pytest.mark.parametrize("bad", [float("nan"), float("inf"), "abc", True, None, [0.1], sp.sqrt(2)])
def test_to_fraction_rejects(bad):
    with pytest.raises(hv.InvalidInput):
        to_fraction(bad)


def test_float_is_converted_exactly():
    assert to_fraction(0.1) == Fraction(3602879701896397, 36028797018963968)


@pytest.mark.parametrize("name,ok,bad", [
    ("mt", ["-0.99", "0", "1"], ["-1", "1.01"]),
    ("sphere", ["0.01", "1"], ["0", "-0.1", "2"]),
    ("plane", ["0"], ["0.1", "-0.1"]),
    ("sinh", ["-0.5", "-0.999"], ["0", "-1", "0.2"]),
])
def test_domains(name, ok, bad):
    m = get_model(name)
    for v in ok:
        m.check_value(Fraction(v))
    for v in bad:
        with pytest.raises(hv.InvalidInput, match="outside the domain"):
            m.check_value(Fraction(v))


def test_bn_closed_forms_match_paper():
    # paper: b_4 = 2 - sqrt3, b_6 = 9 - 4 sqrt5, b_7 = 0, b_8 = -(8 - 3 sqrt7); Lean: b_5 = 3 - 2 sqrt2, b_10 = -1/7
    s = sp.sqrt
    for n, e in [(4, 2 - s(3)), (5, 3 - 2 * s(2)), (6, 9 - 4 * s(5)), (7, 0), (8, -(8 - 3 * s(7))),
                 (10, sp.Rational(-1, 7)), (2, 1), (3, 1)]:
        assert sp.simplify(bn_exact(n) - e) == 0, n
    for n in range(4, 40):
        if n % 2 == 0:
            paper = (4 * (n - 2) * s(n - 1) - n**2) / sp.Integer(n * n - 8 * n + 8)
        elif n != 7:
            paper = (4 * s(n - 3) - (n + 1)) / sp.Integer(n - 7)
        else:
            continue
        assert abs(sp.N(bn_exact(n) - paper, 60)) < sp.Float(10) ** -55, n


def test_bn_solves_eq_bn_and_table_decimals():
    table = {4: "0.2679491924", 5: "0.1715728753", 6: "0.0557280900", 8: "-0.0627460668", 9: "-0.1010205144",
             10: "-0.1428571429", 13: "-0.2251482266", 19: "-0.3333333333", 25: "-0.4021298312"}
    for n in range(4, 26):
        b = sp.N(bn_exact(n), 60)
        lhs = (1 + b**2) / (1 + b) ** 2
        rhs = sp.Rational(n * n // 4, 2 * (n - 1))
        assert abs(lhs - rhs) < sp.Float(10) ** -50
        if n in table:
            assert f"{float(b):.10f}" == table[n]


def test_sinh_q_and_bn_agree():
    for N in range(8, 40):
        q = sinh_q(N)
        xi = (sp.sqrt(sp.Rational(q.numerator, q.denominator)) - 1) / (sp.sqrt(sp.Rational(q.numerator, q.denominator)) + 1)
        assert abs(sp.N(bn_exact(N) + xi, 60)) < sp.Float(10) ** -55, N
    assert [sinh_q(N) for N in (8, 9, 10, 11, 15)] == [Fraction(9, 7), Fraction(3, 2), Fraction(16, 9), 2, 3]


def test_c4_closed_paper_values():
    assert [c4_closed(n) for n in (4, 6, 8, 10, 12)] == [Fraction(1, 4), Fraction(23, 30), Fraction(95, 56),
                                                          Fraction(49, 15), Fraction(755, 132)]
    assert c4_closed(7) == Fraction(9, 14) and c4_closed(3) == 0


def test_conversions():
    m = get_model("mt")
    with mp.workdps(30):
        assert mp.almosteq(m.radius(2 - mp.sqrt(3)), mp.acos(1 / mp.sqrt(3)), 1e-25)  # critical colatitude, n = 4
        assert m.radius(0) == 1
        b8 = -(8 - 3 * mp.sqrt(7))
        assert mp.almosteq(m.radius(b8), mp.mpf("0.511875"), 1e-6)  # sinh note table, r* for N = 8
        assert mp.almosteq(m.cosh_r(b8) ** 2, mp.mpf(9) / 7, 1e-25)
        assert mp.almosteq(m.xi(b8), 8 - 3 * mp.sqrt(7), 1e-25)
        assert m.xi(mp.mpf("0.2")) is None and m.cosh_r(mp.mpf(0)) is None
        assert m.kappa0(0) == 2


def test_cartesian_kernels_consistent():
    m = get_model("mt")
    N, v = 6, 0.3
    _, _, H, J, _ = m.numpy_kernel(N, v)
    z = np.exp(2j * np.pi * np.arange(N) / N) * (1 + 0.01 * np.arange(N))
    with mp.workdps(30):
        zm = [mp.mpc(complex(w)) for w in z]
        assert abs(float(m.H_mp(zm, mp.mpf(v))) - H(z)) < 1e-12
        assert np.allclose([float(t) for t in m.J_mp(zm, mp.mpf(v))], J(z), atol=1e-12)
