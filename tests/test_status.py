import json
import time
from fractions import Fraction

import mpmath as mp
import numpy as np
import pytest

import havelock as hv
from havelock.status import Status, jsonable, make_provenance


def test_version_and_errors():
    assert hv.__version__ == "0.1.0"
    assert issubclass(hv.InvalidInput, hv.HavelockError) and issubclass(hv.InvalidInput, ValueError)
    assert issubclass(hv.OutOfValidityRange, hv.HavelockError) and issubclass(hv.OutOfValidityRange, ValueError)
    assert issubclass(hv.CertificationFailed, hv.HavelockError)


def test_status_values():
    assert [s.value for s in Status] == ["PROVED", "CERTIFIED", "NUMERICAL"]
    assert str(Status.PROVED) == "PROVED"


def test_jsonable_keeps_exact_values():
    with mp.workdps(30):
        x = mp.mpf(1) / 3
    out = jsonable({"q": Fraction(1, 7), "x": x, "z": 1 + 2j, "a": np.array([1.0, 2.0]), "s": Status.CERTIFIED})
    assert out["q"] == "1/7"
    assert out["x"].startswith("0.33333333333333333333333333")
    assert out["z"] == [1.0, 2.0] and out["a"] == [1.0, 2.0] and out["s"] == "CERTIFIED"
    json.dumps(out)


def test_result_to_json_roundtrip():
    prov = make_provenance("test", {"N": 3}, time.perf_counter())
    r = hv.NoThreshold(status=Status.PROVED, basis="b", provenance=prov, model="sinh", N=3,
                       verdict="stable at every radius", geometry="H2")
    d = json.loads(r.to_json())
    assert d["type"] == "NoThreshold" and d["status"] == "PROVED" and d["fallback_from"] is None
    assert d["provenance"]["version"] == "0.1.0" and "python-flint" in d["provenance"]["libraries"]


def test_limits_is_a_copy():
    lim = hv.limits()
    assert lim["precision_ladder_bits"] == [256, 512, 1024, 2048, 4096]
    lim["work_dps"] = -1
    assert hv.limits()["work_dps"] == 50


def test_threshold_private_field_not_serialised():
    prov = make_provenance("t", {}, time.perf_counter())
    t = hv.Threshold(status=Status.PROVED, basis="", provenance=prov, model="mt", N=8, value=mp.mpf(0),
                     exact=0, mode=4, geometry="H2", stable_side="below", _certification=object())
    assert "_certification" not in t.to_dict()
    with pytest.raises(TypeError):
        hv.Threshold(Status.PROVED, "", prov)  # keyword-only
