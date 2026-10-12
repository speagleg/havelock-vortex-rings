import json
import subprocess
import sys

import pytest

from havelock.cli import main


def run(capsys, *argv):
    rc = main(list(argv))
    out, err = capsys.readouterr()
    return rc, out, err


@pytest.mark.parametrize("argv", [
    ["threshold", "sinh", "12"], ["threshold", "sphere", "9"], ["stability", "sinh", "9", "-0.3"],
    ["nonlinear", "mt", "8"], ["simulate", "sinh", "9", "-0.3", "--t-end", "5"], ["limits"],
])
def test_json_output_parses(capsys, argv):
    rc, out, _ = run(capsys, *argv, "--json")
    assert rc == 0
    json.loads(out)


def test_human_output(capsys):
    rc, out, _ = run(capsys, "threshold", "mt", "8")
    assert rc == 0 and "critical mode 4" in out and "PROVED" in out
    rc, out, _ = run(capsys, "nonlinear", "mt", "7")
    assert "c4 = 9/14" in out


def test_cert_and_verify(capsys, tmp_path):
    p = tmp_path / "c.json"
    rc, out, _ = run(capsys, "cert", "sinh", "8", "--json", str(p))
    assert rc == 0 and p.exists()
    rc, out, _ = run(capsys, "verify", str(p))
    assert rc == 0 and out.startswith("VALID")
    d = json.loads(p.read_text())
    d["polynomial"][0] += 1
    p.write_text(json.dumps(d))
    rc, out, _ = run(capsys, "verify", str(p), "--json")
    assert rc == 1 and json.loads(out)["ok"] is False


def test_errors_exit_2(capsys):
    rc, _, err = run(capsys, "threshold", "mt", "3")
    assert rc == 2 and "OutOfValidityRange" in err
    rc, _, err = run(capsys, "stability", "sinh", "9", "0.5")
    assert rc == 2 and "InvalidInput" in err
    rc, _, err = run(capsys, "verify", "no-such-file.json")
    assert rc == 2 and "FileNotFoundError" in err


def test_verify_non_utf8_exits_2(capsys, tmp_path):
    p = tmp_path / "bad.json"
    p.write_bytes(b'{"schema": "\xff\xfe\x80"}')
    rc, out, err = run(capsys, "verify", str(p))
    assert rc == 2 and "UnicodeDecodeError" in err and "Traceback" not in err


def test_negative_rational_via_value_option(capsys):
    rc, out, _ = run(capsys, "stability", "sinh", "9", "--value=-1/7", "--json")
    assert rc == 0 and json.loads(out)["value"]
    rc, out, _ = run(capsys, "simulate", "sinh", "9", "--value=-1/7", "--t-end", "2", "--json")
    assert rc == 0 and json.loads(out)["value"]


def test_value_positional_xor_option(capsys):
    with pytest.raises(SystemExit) as e:
        main(["stability", "sinh", "9"])
    assert e.value.code == 2
    with pytest.raises(SystemExit) as e:
        main(["stability", "sinh", "9", "-0.3", "--value=-0.3"])
    assert e.value.code == 2
    capsys.readouterr()


def test_console_script_module():
    r = subprocess.run([sys.executable, "-m", "havelock", "threshold", "sinh", "8", "--json"],
                       capture_output=True, text=True)
    assert r.returncode == 0 and json.loads(r.stdout)["status"] == "PROVED"


def test_negative_rational_positional_value(capsys):
    rc, out, _ = run(capsys, "stability", "sinh", "9", "-1/7", "--json")
    assert rc == 0 and json.loads(out)["verdict"] == "EM-stable"


def test_threshold_line_names_the_parameter_and_stable_direction(capsys):
    rc, out, _ = run(capsys, "threshold", "sinh", "12")
    assert rc == 0 and "threshold lam_r2 = -0.20241086" in out and "stable for lam_r2 < threshold" in out
    assert "thm:" not in out and "prop:" not in out


@pytest.mark.parametrize("argv", [
    ["stability", "sinh", "9", "-1e-3"], ["stability", "sinh", "9", "--value", "-1e-3"],
    ["stability", "sinh", "9", "-1E-3"], ["stability", "sinh", "9", "-.3"], ["stability", "sinh", "9", "-3e-1"],
    ["stability", "sinh", "9", "--value", "-1/7"], ["stability", "sinh", "9", "--value=-1e-3"],
])
def test_negative_scientific_and_option_values(capsys, argv):
    rc, out, _ = run(capsys, *argv, "--json")
    assert rc == 0 and json.loads(out)["verdict"] in ("EM-stable", "unstable")


def test_negative_value_after_another_option_is_not_taken_as_the_value(capsys):
    rc, out, err = run(capsys, "simulate", "sinh", "9", "-0.3", "--t-end", "-1e-3")
    assert rc == 2 and "t_end" in err


def test_negative_scientific_value_is_used_exactly(capsys):
    rc, out, _ = run(capsys, "stability", "sinh", "9", "-1e-3", "--json")
    assert rc == 0 and json.loads(out)["value"] == "-1/1000"


def test_tanh_threshold_line_says_spectral_and_krein_for_odd_n(capsys):
    rc, out, _ = run(capsys, "threshold", "tanh", "9")
    assert rc == 0 and "(spectrally stable for xi > threshold; Krein window just above it)" in out
    rc, out, _ = run(capsys, "threshold", "tanh", "8")
    assert rc == 0 and "(stable for xi > threshold)" in out


def test_plane_threshold_line_is_flagged_degenerate(capsys):
    rc, out, _ = run(capsys, "threshold", "plane", "7")
    first = out.splitlines()[0]
    assert rc == 0 and "degenerate one-point view" in first and "stable for lam_r2 <" not in first
