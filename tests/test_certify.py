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
    assert str(r) == ("the root ~ 1.414 (isolated in an exact rational interval of width 1.0e-03) of the degree-2 "
                      "polynomial [-2, 0, 1] (constant term first)")


def test_sturm_count_raises_on_inverted_interval():
    """Test that sturm_count raises ValueError if lo > hi. (M5)"""
    p = x**2 - 2
    with pytest.raises(ValueError, match="lo.*hi"):
        sturm_count(p, F(2), F(0))


def test_certify_threshold_interval_clear_records_checks():
    """Test that certify_threshold records the interval_clear field with roots checked. (M6.1 + I1)"""
    dom = get_model("mt").family_domain
    c = certify_threshold("mt", 8, dom, "below")
    assert c.interval_clear is not None
    # All recorded blocks (A_m for m != mode, and B_m for all m) should have 0 roots in the interval
    assert all(count == 0 for count in c.interval_clear.values())
    # Verify that we checked the expected blocks
    assert "B_1" in c.interval_clear and "B_2" in c.interval_clear
    assert "B_3" in c.interval_clear and "B_4" in c.interval_clear


def test_certify_no_threshold_closed_endpoint_mt_n2():
    """Test that certify_no_threshold correctly fails at MT N=2 with a closed endpoint at b=1. (M6.3 + I2)"""
    # At N=2, A_1(1) = 0, so certify_no_threshold should fail if we include the endpoint
    dom_closed_at_hi = (F(-1), F(1), False, True)  # (-1, 1]
    result = certify_no_threshold("mt", 2, dom_closed_at_hi)
    assert result is False, "At N=2, b=1 is a root of A_1, so certify_no_threshold must return False"


def test_certify_threshold_stable_side_above():
    """Test certify_threshold with stable_side 'above'. (M6.2)"""
    from havelock.errors import CertificationFailed
    # MT N=8 is stable below, so asking for stable_side "above" should raise
    # because the sign pattern won't match (all A_m, B_m > 0 above the threshold)
    dom = get_model("mt").family_domain
    with pytest.raises(CertificationFailed, match="sign pattern"):
        certify_threshold("mt", 8, dom, "above")


def test_certify_threshold_ladder_exhaustion(monkeypatch):
    """Test that CertificationFailed message includes the ladder tried when exhausted. (M6.2)"""
    from havelock.errors import CertificationFailed
    from havelock import certify as certify_module

    # Monkeypatch real_root_balls to always raise _Unresolved to force ladder exhaustion
    original_real_root_balls = certify_module.real_root_balls

    def fake_real_root_balls(p, lo, hi, bits):
        raise certify_module._Unresolved()

    monkeypatch.setattr(certify_module, "real_root_balls", fake_real_root_balls)

    dom = get_model("mt").family_domain
    # This must raise CertificationFailed with the ladder in the message
    with pytest.raises(CertificationFailed) as exc_info:
        certify_threshold("mt", 8, dom, "below")

    error_msg = str(exc_info.value)
    # Verify the ladder is mentioned in the error message
    assert "ladder" in error_msg.lower()
    # Verify the precision values from PRECISION_LADDER appear in the message
    assert "128" in error_msg or "256" in error_msg or "bits" in error_msg


def test_certify_threshold_another_block_root_inside_raises(monkeypatch):
    """Test that CertificationFailed is raised when another B_1 shares the threshold root. (M6.3 raise)"""
    from havelock.errors import CertificationFailed
    from havelock import certify as certify_module
    from havelock.blocks import Numerators

    # Create a fake family with A = 1 - 3x (root at 1/3) and B_1 = (3x - 1)(x^2 + 1) (shares root)
    def fake_numerators(family, N, m):
        if family == "fake" and m == 1:
            A = 1 - 3 * x
            B = (3 * x - 1) * (x**2 + 1)
            G = flint.fmpz_poly([0])  # no coupling
            D = flint.fmpz_poly([0])
            return Numerators(A=A, B=B, G=G, D=D, den_a=flint.fmpz_poly([1]), den_b=flint.fmpz_poly([1]))
        return None

    monkeypatch.setattr(certify_module, "numerators", fake_numerators)

    # Use domain (0, 1) open, stable_side "below" (where A > 0)
    dom = (F(0), F(1), False, False)
    # This must raise because B_1 has a root at 1/3 inside the threshold interval
    with pytest.raises(CertificationFailed, match="B_1"):
        certify_threshold("fake", 2, dom, "below")


def test_certify_threshold_another_block_no_root_passes(monkeypatch):
    """Test that certify_threshold succeeds and records interval_clear when B_1 has no shared roots. (M6.3 pass)"""
    from havelock import certify as certify_module
    from havelock.blocks import Numerators

    # Create a fake family with A = 1 - 3x (root at 1/3) and B_1 = x^2 + 1 (no real roots)
    def fake_numerators(family, N, m):
        if family == "fake" and m == 1:
            A = 1 - 3 * x
            B = x**2 + 1
            G = flint.fmpz_poly([0])  # no coupling
            D = flint.fmpz_poly([0])
            return Numerators(A=A, B=B, G=G, D=D, den_a=flint.fmpz_poly([1]), den_b=flint.fmpz_poly([1]))
        return None

    monkeypatch.setattr(certify_module, "numerators", fake_numerators)

    # Use domain (0, 1) open, stable_side "below" (where A > 0)
    dom = (F(0), F(1), False, False)
    # This should succeed
    c = certify_threshold("fake", 2, dom, "below")
    assert c is not None
    # interval_clear should record exactly {"B_1": 0}
    assert c.interval_clear == {"B_1": 0}
    # The interval must bracket 1/3 (this proves the monkeypatch took effect)
    L, U = c.interval
    assert L < F(1, 3) <= U
