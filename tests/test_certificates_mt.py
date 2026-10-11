import copy
import json
from fractions import Fraction

import pytest

import havelock as hv
from havelock.certificates import SCHEMA, Certificate, certificate
from havelock.threshold import threshold


@pytest.mark.parametrize("model,N", [("mt", 4), ("plane", 7), ("sinh", 8), ("mt", 9), ("sinh", 25)])
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


@pytest.mark.parametrize("how", ["shift endpoint", "change coefficient", "flip mode sign", "wrong sturm",
                                 "bad witness", "wrong schema", "missing field"])
def test_corrupted_certificates_fail(how):
    d = certificate(threshold("sinh", 8)).data
    rep = Certificate.verify(_corrupt(d, how))
    assert not rep.ok and rep.failures


def test_verify_from_path(tmp_path):
    p = tmp_path / "c.json"
    p.write_text(certificate(threshold("mt", 5)).to_json())
    assert Certificate.verify(p).ok and Certificate.verify(str(p)).ok


def test_certificate_requires_threshold_and_lean_is_tanh_only(tmp_path):
    with pytest.raises(hv.InvalidInput):
        certificate(threshold("sinh", 5))
    with pytest.raises(hv.InvalidInput, match="tanh"):
        certificate(threshold("mt", 8)).to_lean(tmp_path / "x.lean")


def test_review_focus_large_n_certificate():
    # Review Focus 1: an MT certificate for large N is a quadratic and verifies quickly
    c = certificate(threshold("sinh", 501))
    assert len(c.data["polynomial"]) == 3 and Certificate.verify(c).ok


@pytest.mark.parametrize("bad", ['[1, 2]', '{"schema": "havelock-threshold-certificate/1", "model": "sinh", "N": 8}',
                                 '{"schema": "havelock-thr'])
def test_review_focus_malformed_json_is_invalid_not_a_crash(bad):
    # Review Focus 3: hand-edited or truncated certificates give an INVALID report, never a traceback
    rep = Certificate.verify(bad)
    assert not rep.ok and rep.failures


def test_review_focus_garbage_interval():
    d = copy.deepcopy(certificate(threshold("sinh", 8)).data)
    d["interval"] = ["abc", "1/0"]
    rep = Certificate.verify(d)
    assert not rep.ok and "malformed" in rep.failures[-1]


# ---- Ruling C1: wrong-typed or malformed fields give INVALID, never a traceback ----
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


# ---- Ruling T4-I1: interval_clear is re-checked independently ----
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


# ---- Ruling T4-I2: closed domain endpoints on the stable side ----
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
