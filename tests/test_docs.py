import shlex
from pathlib import Path

import pytest

from havelock.cli import main

REPO = Path(__file__).resolve().parents[1]


def _commands():
    text = (REPO / "docs/quickstart.md").read_text()
    return [shlex.split(line[len("$ havelock "):]) for line in text.splitlines() if line.startswith("$ havelock ")]


def test_docs_exist():
    for name in ("quickstart.md", "models.md", "certificates.md"):
        assert (REPO / "docs" / name).read_text().strip()
    readme = (REPO / "README.md").read_text(encoding="utf-8")
    assert "## Using the `havelock` Python package" in readme and "is not built yet" not in readme
    assert "is planned" not in readme and "planned Python package" not in readme


@pytest.mark.parametrize("argv", _commands())
def test_quickstart_commands_run(argv, tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    if argv[0] == "verify":
        assert main(["cert", "sinh", "8", "--json", "cert.json"]) == 0
    assert main(argv) == 0
