import copy
import hashlib
import shutil
import subprocess
from fractions import Fraction
from pathlib import Path

import mpmath as mp
import pytest

import havelock as hv
from havelock.certificates import Certificate, certificate
from havelock.nonlinear import nonlinear
from havelock.simulate import simulate
from havelock.threshold import threshold

REPO = Path(__file__).resolve().parents[1]


def test_tanh_nonlinear_n8():
    # verified while planning: 186.6363411877127469 (dps 40) and 186.6363411877127482 (dps 60, two eps grids)
    r = nonlinear("tanh", 8)
    assert r.status == hv.Status.NUMERICAL and r.verdict == "stable" and r.method == "level-set engine"
    assert abs(r.c4 / mp.mpf("186.63634118771274") - 1) < mp.mpf("1e-12")


@pytest.mark.slow
def test_tanh_nonlinear_n10():
    # verified while planning: 58.0912759473557533698 (dps 60)
    assert abs(nonlinear("tanh", 10).c4 / mp.mpf("58.09127594735575") - 1) < mp.mpf("1e-12")


@pytest.mark.parametrize("N", [9, 5])
def test_tanh_nonlinear_out_of_validity(N):
    with pytest.raises(hv.OutOfValidityRange):
        nonlinear("tanh", N)


def test_tanh_simulation():
    assert simulate("tanh", 8, "0.3").growth < 1.5
    s = simulate("tanh", 8, "0.05")
    assert s.growth > 100 and s.energy_drift < s.drift_tolerance and s.momentum_drift < s.drift_tolerance


@pytest.mark.parametrize("N", [8, 9, 12, 16])
def test_tanh_certificate_roundtrip(N):
    c = certificate(threshold("tanh", N))
    assert c.data["stable_side"]["side"] == "above" and "D" in c.data["mode_signs"]["1"]["below"]
    rep = Certificate.verify(c.to_json())
    assert rep.ok, rep.failures


@pytest.mark.parametrize("N", [17, 18, 20, 24])
def test_review_focus_beyond_the_table_certificates_verify(N):
    # Review Focus 2: certificates beyond the paper's table verify too
    assert Certificate.verify(certificate(threshold("tanh", N))).ok


@pytest.mark.slow
def test_size_cap_n64_certificate_verifies():
    # benchmark while planning: verification of the N = 64 certificate takes 19-21 s
    assert Certificate.verify(certificate(threshold("tanh", 64))).ok


@pytest.mark.parametrize("how", ["shift endpoint", "change coefficient", "flip mode sign", "flip D sign"])
def test_tanh_corrupted_certificates_fail(how):
    d = copy.deepcopy(certificate(threshold("tanh", 9)).data)
    if how == "shift endpoint":
        lo, hi = (Fraction(s) for s in d["interval"])
        d["interval"] = [str(lo - 3 * (hi - lo)), str(lo - 2 * (hi - lo))]
    elif how == "change coefficient":
        d["polynomial"][5] -= 1
    elif how == "flip mode sign":
        d["mode_signs"]["4"]["above"]["A"] *= -1
    else:
        d["mode_signs"]["3"]["above"]["D"] *= -1
    assert not Certificate.verify(d).ok


# sha256 of the emitted N = 8 module, recorded while planning
LEAN_N8_SHA256 = "db9ab9284d7f9594f4685d267b617366ad82caa61ea95930aa038860d1d4cd1f"


def test_lean_emitter_n8_golden(tmp_path):
    p = certificate(threshold("tanh", 8)).to_lean(tmp_path / "CertN8.lean")
    text = p.read_text(encoding="utf-8")
    assert text.startswith("import Havelock.Tanh.General\n") and "theorem cert_N8" in text
    assert hashlib.sha256(text.encode("utf-8")).hexdigest() == LEAN_N8_SHA256
    committed = REPO / "lean/Havelock/Tanh/CertN8.lean"
    if committed.exists():
        marker = "/-! ### N = 8 -/\n"
        assert text.split(marker, 1)[1] == committed.read_text(encoding="utf-8").split(marker, 1)[1]


@pytest.mark.parametrize("N", range(9, 17))
def test_lean_emitter_runs_for_table(N, tmp_path):
    p = certificate(threshold("tanh", N)).to_lean(tmp_path / f"C{N}.lean")
    assert f"theorem cert_N{N}" in p.read_text(encoding="utf-8")


@pytest.mark.slow
@pytest.mark.skipif(shutil.which("lake") is None, reason="needs a Lean toolchain (lake) and a built lean/ project")
def test_lean_emitted_n8_compiles(tmp_path):
    p = certificate(threshold("tanh", 8)).to_lean(tmp_path / "ExportN8.lean")
    r = subprocess.run(["lake", "env", "lean", str(p)], cwd=REPO / "lean", capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr


def test_tanh_certificate_carries_interval_clear_and_verify_recomputes_it():
    d = copy.deepcopy(certificate(threshold("tanh", 9)).data)
    assert d["interval_clear"] and all(v == 0 for v in d["interval_clear"].values())
    assert Certificate.verify(d).ok
    d["interval_clear"]["A_1"] = 1
    rep = Certificate.verify(d)
    assert not rep.ok and any("interval_clear" in f for f in rep.failures)
    d2 = copy.deepcopy(certificate(threshold("tanh", 9)).data)
    del d2["interval_clear"]["B_1"]
    assert not Certificate.verify(d2).ok


def test_tanh_nonlinear_odd_n_message_explains_the_krein_window():
    with pytest.raises(hv.OutOfValidityRange, match="indefinite") as e:
        nonlinear("tanh", 9)
    assert "Krein window" in str(e.value) and "linearly stable" in str(e.value)


def test_lean_export_refuses_mt_cleanly(tmp_path):
    with pytest.raises(hv.InvalidInput, match="tanh model only"):
        certificate(threshold("mt", 5)).to_lean(tmp_path / "x.lean")


# ---- final review F4: Lean export range is tanh N = 8..16 ----
def test_lean_export_refuses_n17_with_a_clear_message(tmp_path):
    with pytest.raises(hv.CertificationFailed, match=r"tanh N = 8\.\.16 only"):
        certificate(threshold("tanh", 17)).to_lean(tmp_path / "C17.lean")
    assert not (tmp_path / "C17.lean").exists()


# ---- final review F7: r on Stability / Nonlinear / Simulation (tanh) ----
def test_tanh_results_report_r():
    with mp.workdps(50):
        r = 2 * mp.atanh(mp.sqrt(mp.mpf(3) / 10))
        s = hv.stability("tanh", 8, "0.3")
        assert abs(s.r - r) < mp.mpf("1e-40") and "r" in s.to_dict()
        sim = simulate("tanh", 8, "0.3", t_end=1.0)
        assert abs(sim.r - r) < mp.mpf("1e-40") and "r" in sim.to_dict()
    nl = nonlinear("tanh", 8)
    assert nl.r == threshold("tanh", 8).r and abs(nl.r - mp.mpf("0.5312104560")) < mp.mpf("1e-9")
