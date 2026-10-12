import copy
import re
from fractions import Fraction
from pathlib import Path

import mpmath as mp
import pytest
import sympy as sp

import havelock as hv
from havelock.threshold import _tanh_numerical, threshold

REPO = Path(__file__).resolve().parents[1]

# tanh/note_tanh.tex Table 1: xi* (25 digits), r*, cosh r*, degree and leading coefficient of the minimal polynomial;
# mod-p witnesses: note Section irred ("p = 7, 13, 107, 17, 11, 101, 3, 29, 163 for N = 8..16").
TABLE = {
    8: ("0.0673562230040812622684374", "0.5312104560", "1.1444414784", 10, 1, 7),
    9: ("0.1128710632366027404396940", "0.6990695396", "1.2544637167", 16, 4, 13),
    10: ("0.1661142415262391180817586", "0.8653901256", "1.3984100696", 12, 7, 107),
    11: ("0.2046011768826692906355358", "0.9752480466", "1.5144618547", 20, 10, 17),
    12: ("0.2474177768040677091739650", "1.0917203545", "1.6575169309", 16, 7, 11),
    13: ("0.2798613886439066870299965", "1.1775646302", "1.7772431147", 24, 18, 101),
    14: ("0.3149389631155066233935054", "1.2691494599", "1.9194478920", 18, 23, 3),
    15: ("0.3424585835121702566010226", "1.3406681986", "2.0416335000", 28, 28, 29),
    16: ("0.3716818139060454566367287", "1.4167509599", "2.1831006077", 22, 17, 163),
}


@pytest.mark.parametrize("N", sorted(TABLE))
def test_note_table_1(N):
    xs, r, ch, deg, lead, p = TABLE[N]
    t = threshold("tanh", N)
    assert isinstance(t, hv.Threshold) and t.status == hv.Status.CERTIFIED and t.fallback_from is None
    assert t.mode == N // 2 and t.stable_side == "above" and t.geometry == "H2"
    assert mp.nstr(t.value, 40).startswith(xs)                      # the table truncates to 25 digits
    assert abs(t.r - mp.mpf(r)) < mp.mpf("6e-11") and abs(t.cosh_r - mp.mpf(ch)) < mp.mpf("6e-11")
    f = t.exact.poly
    assert len(f) - 1 == deg and f[-1] == lead and list(f) == list(reversed(f))       # palindromic
    assert t.exact.lo < Fraction(xs) + Fraction(1, 10**25) and t.exact.hi > Fraction(xs) - Fraction(1, 10**25)
    assert t._certification.irreducibility_witness == {"p": p, "degrees": [deg]}
    assert f"cert_N{N}" in t.basis
    ic = t._certification.interval_clear                 # review fix: no B_m, and no A_m (m != mode), in the interval
    assert ic and all(v == 0 for v in ic.values())
    assert set(ic) == {f"B_{m}" for m in range(1, N // 2 + 1)} | {f"A_{m}" for m in range(1, N // 2 + 1) if m != t.mode}


@pytest.mark.parametrize("N", range(2, 8))
def test_small_rings_stable_everywhere(N):
    t = threshold("tanh", N)
    assert isinstance(t, hv.NoThreshold) and t.verdict == "stable at every radius"
    assert t.status == hv.Status.CERTIFIED and f"stable_N{N}" in t.basis


def test_reproduces_certify_v2_out():
    # research reproduction: minimal polynomials and 50-digit xi* recorded in tanh/certify_v2.out
    text = (REPO / "tanh/certify_v2.out").read_text()
    blocks = re.split(r"^N=(\d+) critical mode", text, flags=re.M)
    x = sp.Symbol("x")
    seen = 0
    for i in range(1, len(blocks), 2):
        N, body = int(blocks[i]), blocks[i + 1]
        xs = re.search(r"xi\*  = ([0-9.]+)", body).group(1)
        mp_ = re.search(r"minpoly\(xi\*\) = (.+)", body).group(1)
        coeffs = [int(c) for c in reversed(sp.Poly(sp.sympify(mp_), x).all_coeffs())]
        t = threshold("tanh", N)
        assert list(t.exact.poly) == coeffs
        assert mp.nstr(t.value, 50, strip_zeros=False)[:48] == xs[:48]
        seen += 1
    assert seen == 9


def test_size_cap_strict_and_fallback():
    assert hv.limits()["tanh_exact_max_N"] == 64
    with pytest.raises(hv.CertificationFailed, match="size cap"):
        threshold("tanh", 65)
    t = threshold("tanh", 65, strict=False)
    assert t.status == hv.Status.NUMERICAL and t.fallback_from == hv.Status.CERTIFIED and t.exact is None
    assert t.mode == 32 and mp.mpf("0.79") < t.value < mp.mpf("0.80")


def test_numerical_fallback_agrees_with_certified():
    v, m = _tanh_numerical(16)
    t = threshold("tanh", 16)
    assert m == 8 and abs(v - t.value) < mp.mpf("1e-40")


def test_precision_argument():
    t = threshold("tanh", 9, precision=512)
    assert t._certification.bits == 512


@pytest.mark.slow
def test_size_cap_n64_certifies():
    t = threshold("tanh", 64)
    assert t.status == hv.Status.CERTIFIED and t.mode == 32 and len(t.exact.poly) - 1 == 94


@pytest.mark.parametrize("N,deg", [(17, 32), (18, 24), (20, 28), (24, 34)])
def test_review_focus_beyond_the_table(N, deg):
    # Review Focus 2: N above the paper's table but under the cap: mode floor(N/2), degree pattern (eq. deg:
    # 2N-2 odd, 3N/2-2 for N = 0 mod 4, 3N/2-3 for N = 2 mod 4), palindromic minimal polynomial
    t = threshold("tanh", N)
    assert t.status == hv.Status.CERTIFIED and t.mode == N // 2
    assert len(t.exact.poly) - 1 == deg and list(t.exact.poly) == list(reversed(t.exact.poly))


# Controller ruling C2: certificates exist only for certified results.
def _fallback_threshold():
    return hv.Threshold(status=hv.Status.NUMERICAL, fallback_from=hv.Status.CERTIFIED, basis="x",
                        provenance=hv.status.make_provenance("threshold", {}, 0.0), model="tanh", N=65,
                        value=mp.mpf("0.79"), exact=None, mode=32, geometry="H2", stable_side="above")


def test_certificate_refuses_numerical_threshold_fast():
    with pytest.raises(hv.CertificationFailed, match="no certificate for a NUMERICAL result"):
        hv.certificate(_fallback_threshold())


def test_certificate_refuses_fallback_from_even_if_status_certified():
    t = _fallback_threshold()
    t.status = hv.Status.CERTIFIED
    with pytest.raises(hv.CertificationFailed, match="no certificate for a NUMERICAL result"):
        hv.certificate(t)


@pytest.mark.slow
def test_certificate_refuses_real_fallback_above_cap():
    t = threshold("tanh", 65, strict=False)
    with pytest.raises(hv.CertificationFailed, match="no certificate for a NUMERICAL result"):
        hv.certificate(t)


KREIN = "the energy-momentum threshold xi_EM lies strictly above xi* (exact D_"


def test_basis_states_spectral_threshold_and_krein_window():
    odd, even = threshold("tanh", 9), threshold("tanh", 10)
    assert "spectral threshold xi* (note_tanh Theorem 7" in odd.basis and KREIN in odd.basis
    assert "spectral threshold xi*" in even.basis and "Krein" not in even.basis
    assert "Theorem 7(b)" in odd.basis and "Thm 6" not in odd.basis


def test_small_ring_basis_cites_theorem_7a():
    b = threshold("tanh", 5).basis
    assert "Theorem 7(a)" in b and "Thm 6" not in b


def test_fallback_reason_does_not_advise_strict_false():
    t = threshold("tanh", 65, strict=False)
    assert "size cap" in t.basis and "use strict=False" not in t.basis


@pytest.mark.parametrize("strict", [True, False])
def test_small_ring_fails_closed_when_checks_disagree(monkeypatch, strict):
    import sys
    th = sys.modules["havelock.threshold"]
    monkeypatch.setattr(th, "certify_no_threshold", lambda *a, **k: False)
    monkeypatch.setattr(th, "certify_threshold", lambda *a, **k: None)
    if strict:
        with pytest.raises(hv.CertificationFailed):
            th.threshold("tanh", 5)
    else:
        # N = 5 has no sign change for the NUMERICAL fallback to find, so it raises too; it never returns CERTIFIED
        with pytest.raises(hv.CertificationFailed, match="no sign change"):
            th.threshold("tanh", 5, strict=False)


def test_review_focus_n33():
    t = threshold("tanh", 33)
    assert t.status == hv.Status.CERTIFIED and t.mode == 16
    assert len(t.exact.poly) - 1 == 64       # 2N - 2 for odd N (eq. deg)


# ---- final review F2: no certified verdict when certify_threshold returns None, for any N ----
@pytest.mark.parametrize("strict", [True, False])
def test_n8_fails_closed_when_no_sign_change_is_reported(monkeypatch, strict):
    import sys
    th = sys.modules["havelock.threshold"]
    monkeypatch.setattr(th, "certify_threshold", lambda *a, **k: None)
    if strict:
        with pytest.raises(hv.CertificationFailed, match="no certified verdict"):
            th.threshold("tanh", 8)
    else:
        t = th.threshold("tanh", 8, strict=False)
        assert isinstance(t, hv.Threshold) and t.status == hv.Status.NUMERICAL
        assert t.fallback_from == hv.Status.CERTIFIED and "no certified verdict" in t.basis


# ---- final review F3: Krein sentence only when D_crit < 0 is certified; citation scoped to N <= 16 ----
def test_n17_basis_states_the_krein_gap_and_the_n_above_16_citation():
    t = threshold("tanh", 17)
    c = t._certification
    assert c.mode_signs[c.mode]["above"]["D"] == -1
    assert KREIN in t.basis and "Krein window" in t.basis
    assert "as in note_tanh Theorem 7 (proved there for N <= 16); certified here at run time" in t.basis


def test_n20_basis_carries_the_n_above_16_citation():
    t = threshold("tanh", 20)
    assert "proved there for N <= 16" in t.basis and "Krein" not in t.basis and "Lean-checked" not in t.basis


def test_n9_basis_cites_theorem_7_directly():
    assert "proved there for N <= 16" not in threshold("tanh", 9).basis


def test_krein_sentence_omitted_when_d_is_not_certified_negative(monkeypatch):
    import dataclasses
    import sys
    th = sys.modules["havelock.threshold"]
    real = th.certify_threshold

    def fake(*a, **k):
        c = real(*a, **k)
        ms = copy.deepcopy(c.mode_signs)
        ms[c.mode]["above"]["D"] = 1
        return dataclasses.replace(c, mode_signs=ms)

    monkeypatch.setattr(th, "certify_threshold", fake)
    t = th.threshold("tanh", 9)
    assert t.status == hv.Status.CERTIFIED and "Krein" not in t.basis and "xi_EM" not in t.basis
