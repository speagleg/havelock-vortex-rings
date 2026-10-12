import copy
import json
from fractions import Fraction

import pytest

import havelock as hv
from havelock.certificates import SCHEMA, Certificate, certificate
from havelock.threshold import threshold


@pytest.mark.parametrize("model,N", [("mt", 4), ("sphere", 5), ("plane", 7), ("sinh", 8), ("mt", 9), ("sinh", 25)])
def test_roundtrip(model, N):
    c = certificate(threshold(model, N))
    d = json.loads(c.to_json())
    assert d["schema"] == SCHEMA and d["model"] == model and d["N"] == N and d["mode"] == N // 2
    assert d["sturm"] == 1 and d["stable_side"]["side"] == "below"
    rep = Certificate.verify(c.to_json())
    assert rep.ok, rep.failures
    assert rep.checks > 10


def test_sinh8_contents():
    d = certificate(threshold("sinh", 8)).data
    assert d["polynomial"] == [-2, -32, -2]                      # 14(1+b^2) - 16(1+b)^2 = -2(b^2 + 16b + 1)
    assert d["irreducibility_witness"] == {"p": 5, "degrees": [2]}    # b^2 + 16b + 1 = (b - 1)^2 mod 3
    lo, hi = (Fraction(s) for s in d["interval"])
    assert Fraction(-6274606680623, 10**14) < lo < hi < Fraction(-6274606680622, 10**14)   # b_8 = -0.06274606680622822...


def test_interval_clear_is_recorded():
    d = certificate(threshold("sinh", 8)).data
    ic = d["interval_clear"]
    assert set(ic) == {"B_1", "B_2", "B_3", "B_4", "A_1", "A_2", "A_3"}      # mode is 4, so no A_4
    assert all(v == 0 for v in ic.values())


def _corrupt(d, how):
    d = copy.deepcopy(d)
    if how == "shift endpoint":
        lo, hi = (Fraction(s) for s in d["interval"])
        d["interval"] = [str(hi + (hi - lo)), str(hi + 2 * (hi - lo))]
    elif how == "change coefficient":
        d["polynomial"][0] += 1
    elif how == "flip mode sign":
        d["mode_signs"]["2"]["below"]["A"] *= -1
    elif how == "wrong sturm":
        d["sturm"] = 2
    elif how == "bad witness":
        d["irreducibility_witness"] = {"p": 2, "degrees": [2]}
    elif how == "wrong schema":
        d["schema"] = "other/1"
    elif how == "missing field":
        del d["mode_signs"]
    return d


_REASON = {"shift endpoint": "Sturm count", "change coefficient": "not the numerator", "flip mode sign": "mode_signs",
           "wrong sturm": "Sturm", "bad witness": "not irreducible", "wrong schema": "schema",
           "missing field": "malformed"}


@pytest.mark.parametrize("how", ["shift endpoint", "change coefficient", "flip mode sign", "wrong sturm",
                                 "bad witness", "wrong schema", "missing field"])
def test_corrupted_certificates_fail(how):
    d = certificate(threshold("sinh", 8)).data
    rep = Certificate.verify(_corrupt(d, how))
    assert not rep.ok and rep.failures
    assert any(_REASON[how] in f for f in rep.failures), rep.failures


def test_verify_from_path(tmp_path):
    p = tmp_path / "c.json"
    p.write_text(certificate(threshold("mt", 5)).to_json())
    assert Certificate.verify(p).ok and Certificate.verify(str(p)).ok


def test_certificate_requires_threshold_and_lean_is_tanh_only(tmp_path):
    with pytest.raises(hv.InvalidInput):
        certificate(threshold("sinh", 5))
    with pytest.raises(hv.InvalidInput, match="tanh"):
        certificate(threshold("mt", 8)).to_lean(tmp_path / "x.lean")


def test_large_n_certificate():
    # an MT certificate for large N is a quadratic and verifies quickly
    c = certificate(threshold("sinh", 501))
    assert len(c.data["polynomial"]) == 3 and Certificate.verify(c).ok


@pytest.mark.parametrize("bad", ['[1, 2]', '{"schema": "havelock-threshold-certificate/1", "model": "sinh", "N": 8}',
                                 '{"schema": "havelock-thr'])
def test_malformed_json_is_invalid_not_a_crash(bad):
    # hand-edited or truncated certificates give an INVALID report, never a traceback
    rep = Certificate.verify(bad)
    assert not rep.ok and rep.failures


def test_garbage_interval():
    d = copy.deepcopy(certificate(threshold("sinh", 8)).data)
    d["interval"] = ["abc", "1/0"]
    rep = Certificate.verify(d)
    assert not rep.ok and "malformed" in rep.failures[-1]


# ---- wrong-typed or malformed fields give INVALID, never a traceback ----
def _mutate(field, value):
    d = copy.deepcopy(certificate(threshold("sinh", 8)).data)
    d[field] = value
    return d


@pytest.mark.parametrize("field,value", [
    ("stable_side", []), ("stable_side", "x"), ("irreducibility_witness", "x"), ("irreducibility_witness", []),
    ("irreducibility_witness", {"p": "5", "degrees": [2]}), ("irreducibility_witness", {"p": 10**400, "degrees": [2]}),
    ("polynomial", [-2, float("inf"), -2]), ("polynomial", [-2, -32.5, -2]), ("polynomial", [-2, True, -2]),
    ("polynomial", "abc"), ("polynomial", 5), ("polynomial", []),
    ("interval", "ab"), ("interval", [1, 2]), ("interval", None), ("interval", ["1/2"]),
    ("points", 7), ("points", [None, None]), ("mode_signs", []), ("mode_signs", {"1": []}),
    ("mode_signs", {"1": {"below": 3}}), ("endpoint_signs", 3), ("endpoint_signs", None),
    ("factor_index", "0"), ("factor_index", [0]), ("sturm", "1"), ("sturm", None),
    ("N", "8"), ("N", 8.0), ("N", True), ("N", 10**30), ("mode", [4]), ("model", ["sinh"]), ("model", None),
    ("interval_clear", []), ("interval_clear", "x"), ("interval_clear", {"B_1": "0"}),
])
def test_wrong_typed_fields_are_invalid_not_a_crash(field, value):
    rep = Certificate.verify(_mutate(field, value))
    assert not rep.ok and rep.failures


def test_json_infinity_in_polynomial_is_invalid():
    d = certificate(threshold("sinh", 8)).data
    text = json.dumps(d).replace("[-2, -32, -2]", "[-2, Infinity, -2]")
    assert "Infinity" in text
    rep = Certificate.verify(text)
    assert not rep.ok and rep.failures


def test_verify_dict_rejects_non_dict():
    from havelock.certificates import verify_dict
    assert not verify_dict([1, 2]).ok and not verify_dict(None).ok


# ---- interval_clear is re-checked independently ----
def test_interval_clear_deleted_is_invalid():
    d = certificate(threshold("sinh", 8)).data
    del d["interval_clear"]
    rep = Certificate.verify(d)
    assert not rep.ok and any("interval_clear" in f for f in rep.failures)


def test_interval_clear_altered_is_invalid():
    d = copy.deepcopy(certificate(threshold("sinh", 8)).data)
    d["interval_clear"]["A_1"] = 1
    rep = Certificate.verify(d)
    assert not rep.ok and any("interval_clear" in f for f in rep.failures)
    d = copy.deepcopy(certificate(threshold("sinh", 8)).data)
    del d["interval_clear"]["B_2"]
    rep = Certificate.verify(d)
    assert not rep.ok and any("interval_clear" in f for f in rep.failures)


def test_interval_containing_another_blocks_root_is_invalid():
    # widen the interval so that it swallows a root of A_3 (the critical mode is 4): the recorded zero counts are
    # then false, and the verifier's own Sturm counts must say so
    d = copy.deepcopy(certificate(threshold("sinh", 8)).data)
    lo, hi = (Fraction(s) for s in d["interval"])
    d["interval"] = [str(lo), "-1/100000"]
    rep = Certificate.verify(d)
    assert not rep.ok
    assert any("interval_clear" in f and "A_" in f for f in rep.failures), rep.failures


# ---- closed domain endpoints on the stable side ----
def test_closed_stable_endpoint_is_checked(monkeypatch):
    from havelock.certificates import verify as V
    d = certificate(threshold("sinh", 8)).data
    base = Certificate.verify(d)
    assert base.ok
    monkeypatch.setitem(V._CLOSED, "mt", (True, True))        # stable side is "below": the lower end becomes closed
    closed = Certificate.verify(d)
    assert closed.ok and closed.checks > base.checks

    real = V._REBUILD["mt"]

    def bad(N, m):
        A, B, D, G = real(N, m)
        return (A if m != 1 else V._add(A, [-10**6]), B, D, G)

    monkeypatch.setitem(V._REBUILD, "mt", bad)
    rep = Certificate.verify(d)
    assert not rep.ok and any("closed" in f for f in rep.failures)


# ---- size limits ----
def test_certificate_refuses_n_above_verify_max():
    from havelock.limits import VERIFY_MAX_N, limits
    assert VERIFY_MAX_N == 100000 and limits()["certificate_max_N"] == VERIFY_MAX_N
    with pytest.raises(hv.CertificationFailed, match="limited to N <= 100000"):
        certificate(threshold("sinh", VERIFY_MAX_N + 2))


def test_verifier_literals_mirror_limits():
    from havelock.certificates import verify as V
    from havelock.limits import CERT_MAX_DIGITS, TANH_EXACT_MAX_N, VERIFY_MAX_N, WITNESS_MAX_PRIME, limits
    assert V._MAX_N == VERIFY_MAX_N and V._MAX_PRIME == WITNESS_MAX_PRIME
    assert V._MAX_DIGITS == CERT_MAX_DIGITS < 4300 and limits()["certificate_max_digits"] == CERT_MAX_DIGITS
    assert V._TANH_MAX_N == TANH_EXACT_MAX_N     # raising the tanh cap must update the verifier mirror


def test_json_parse_limits_are_invalid_not_a_crash():
    big = '{"schema": 1' + "0" * 5000 + "}"
    assert not Certificate.verify(big).ok
    deep = "[" * 100000
    rep = Certificate.verify(deep)
    assert not rep.ok and rep.failures


def test_strict_types_in_mode_signs_counts_degrees():
    base = certificate(threshold("sinh", 8)).data
    d = copy.deepcopy(base); d["mode_signs"]["1"]["below"]["A"] = True
    assert not Certificate.verify(d).ok
    d = copy.deepcopy(base); d["mode_signs"]["1"]["below"]["B"] = 1.0
    assert not Certificate.verify(d).ok
    d = copy.deepcopy(base); d["stable_side"]["root_counts"]["1"]["A"] = False
    assert not Certificate.verify(d).ok
    d = copy.deepcopy(base); d["stable_side"]["root_counts"]["1"]["B"] = 0.0
    assert not Certificate.verify(d).ok
    d = copy.deepcopy(base); d["irreducibility_witness"]["degrees"] = [2.0]
    assert not Certificate.verify(d).ok


# ---- the model's view must contain the certified root ----
@pytest.mark.parametrize("view", ["sphere", "plane"])
def test_relabelled_view_not_containing_the_threshold_is_invalid(view):
    d = copy.deepcopy(certificate(threshold("mt", 9)).data)       # b_9 < 0: threshold("sphere"|"plane", 9) has none
    assert Certificate.verify(d).ok
    d["model"] = view
    rep = Certificate.verify(d)
    assert not rep.ok and any(f"{view!r} view" in f for f in rep.failures), rep.failures


def test_relabelled_sphere_threshold_to_sinh_or_plane_is_invalid():
    d = copy.deepcopy(certificate(threshold("mt", 4)).data)       # b_4 > 0
    for view in ("sinh", "plane"):
        e = dict(d, model=view)
        rep = Certificate.verify(e)
        assert not rep.ok and any(f"{view!r} view" in f for f in rep.failures), rep.failures
    assert Certificate.verify(dict(d, model="sphere")).ok


def test_plane_view_requires_the_root_zero():
    d = certificate(threshold("plane", 7)).data
    assert Certificate.verify(d).ok and Certificate.verify(dict(d, model="mt")).ok
    assert not Certificate.verify(dict(d, model="sinh")).ok and not Certificate.verify(dict(d, model="sphere")).ok


# ---- the verifier caps tanh N ----
def test_verifier_refuses_tanh_n_above_cap_quickly():
    import time
    d = copy.deepcopy(certificate(threshold("tanh", 8)).data)
    d["N"] = 100000
    t0 = time.perf_counter()
    rep = Certificate.verify(d)
    assert time.perf_counter() - t0 < 5
    assert not rep.ok and any("limited to N <= 64" in f for f in rep.failures), rep.failures
    d["N"] = 65
    assert any("limited to N <= 64" in f for f in Certificate.verify(d).failures)


# ---- duplicate keys and fields outside the schema ----
def test_duplicate_json_keys_are_invalid():
    text = certificate(threshold("sinh", 8)).to_json()
    dup = text.replace('"model": "sinh",', '"model": "tanh", "model": "sinh",', 1)
    assert dup != text and json.loads(dup)["model"] == "sinh"          # plain json.loads keeps the last value
    rep = Certificate.verify(dup)
    assert not rep.ok and any("duplicate key 'model'" in f for f in rep.failures), rep.failures


def test_duplicate_keys_in_a_nested_object_are_invalid(tmp_path):
    text = certificate(threshold("sinh", 8)).to_json()
    dup = text.replace('"side": "below"', '"side": "above", "side": "below"', 1)
    assert dup != text
    p = tmp_path / "dup.json"
    p.write_text(dup)
    rep = Certificate.verify(p)
    assert not rep.ok and any("duplicate key 'side'" in f for f in rep.failures), rep.failures


@pytest.mark.parametrize("path", [(), ("stable_side",), ("irreducibility_witness",), ("mode_signs", "1")])
def test_unknown_fields_are_invalid(path):
    d = copy.deepcopy(certificate(threshold("sinh", 8)).data)
    assert Certificate.verify(d).ok
    target = d
    for k in path:
        target = target[k]
    target["extra"] = 1
    rep = Certificate.verify(d)
    assert not rep.ok and any("extra" in f or "fields" in f for f in rep.failures), rep.failures


def test_unknown_mode_keys_in_root_counts_are_invalid():
    d = copy.deepcopy(certificate(threshold("sinh", 8)).data)
    d["stable_side"]["root_counts"]["9"] = {"A": 0, "B": 0}
    assert not Certificate.verify(d).ok


def test_provenance_and_witness_are_optional():
    d = copy.deepcopy(certificate(threshold("sinh", 8)).data)
    del d["provenance"]
    d["irreducibility_witness"] = None
    assert Certificate.verify(d).ok


# ---- rational size: a clear message below Python's int/str limit, and huge endpoints rejected quickly ----
def test_rational_digit_cap_has_a_clear_message():
    from havelock.limits import CERT_MAX_DIGITS
    d = copy.deepcopy(certificate(threshold("sinh", 8)).data)
    d["interval"] = ["-1/" + "1" * (CERT_MAX_DIGITS + 1), d["interval"][1]]
    rep = Certificate.verify(d)
    assert not rep.ok and any(f"more than {CERT_MAX_DIGITS} digits" in f for f in rep.failures), rep.failures
    d["interval"][0] = "-1/" + "1" * 4400                       # beyond Python's 4300-digit limit: same message
    rep = Certificate.verify(d)
    assert not rep.ok and any(f"more than {CERT_MAX_DIGITS} digits" in f for f in rep.failures), rep.failures


def test_padded_but_equal_endpoints_still_verify():
    d = copy.deepcopy(certificate(threshold("sinh", 8)).data)
    lo, hi = (Fraction(s) for s in d["interval"])
    pad = 10 ** 2000
    d["interval"] = [f"{lo.numerator * pad}/{lo.denominator * pad}", f"{hi.numerator * pad}/{hi.denominator * pad}"]
    assert Certificate.verify(d).ok


def _huge(q, digits):
    den = 10 ** digits + 1
    return str(Fraction(int(q * den), den))


def test_huge_wrong_endpoints_are_rejected_quickly():
    import time
    from havelock.limits import CERT_MAX_DIGITS
    d = copy.deepcopy(certificate(threshold("tanh", 12)).data)
    lo, hi = (Fraction(s) for s in d["interval"])
    w = hi - lo
    d["interval"] = [_huge(lo - 1000 * w, CERT_MAX_DIGITS - 1), _huge(hi + 1000 * w, CERT_MAX_DIGITS - 1)]
    t0 = time.perf_counter()
    rep = Certificate.verify(d)
    assert time.perf_counter() - t0 < 20
    assert not rep.ok


def test_huge_but_correct_enclosure_verifies_with_the_exact_counts():
    # endpoints with ~2000 digits that still enclose the threshold tightly: every recomputed count is exact
    import time
    d = copy.deepcopy(certificate(threshold("tanh", 12)).data)
    lo, hi = (Fraction(s) for s in d["interval"])
    d["interval"] = [_huge(lo, 2000), _huge(hi, 2000)]
    lo2, hi2 = (Fraction(s) for s in d["interval"])
    assert lo2 <= lo and len(d["interval"][0]) > 3000
    d["interval"][1] = str(hi2 + Fraction(1, 10 ** 2000 + 1))     # round the upper end outwards
    d["stable_side"]["interval"][0] = d["interval"][1]
    t0 = time.perf_counter()
    rep = Certificate.verify(d)
    assert rep.ok, rep.failures
    assert time.perf_counter() - t0 < 60


def test_certificate_refuses_endpoints_above_the_digit_cap():
    t = threshold("tanh", 8, precision=9000)            # endpoints of ~18000 bits: more than 2600 digits
    with pytest.raises(hv.CertificationFailed, match="more than 2600 digits"):
        certificate(t)
