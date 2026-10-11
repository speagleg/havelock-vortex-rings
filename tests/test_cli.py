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
