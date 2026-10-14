"""elastemp must never depend on, or change, the process working directory."""

import re
from pathlib import Path

import pytest

from elastemp.base.strain_analysis import get_curvature
from elastemp.cli import main
from test_fit_points import CURVATURE_GPA, STRAINS, _make_deformation

SRC = Path(__file__).resolve().parents[1] / "src" / "elastemp"


def test_no_chdir_or_shell_calls_in_package():
    offenders = []
    for path in SRC.rglob("*.py"):
        for n, line in enumerate(path.read_text().splitlines(), start=1):
            if re.search(r"os\.chdir|os\.system|Path\.cwd\(\)|os\.getcwd\(\)", line):
                offenders.append(f"{path.relative_to(SRC)}:{n}: {line.strip()}")
    # resolve_workspace() in base/utils.py is the only place the current directory is used,
    # as the default when no project folder is given.
    assert len(offenders) == 1 and offenders[0].startswith("base/utils.py"), offenders


def test_no_pyplot_global_state():
    for path in SRC.rglob("*.py"):
        assert "matplotlib.pyplot" not in path.read_text(), path


def test_curvature_with_explicit_folder_leaves_cwd_untouched(tmp_path, monkeypatch):
    project = tmp_path / "project"
    project.mkdir()
    deformation = _make_deformation(project, unconverged={4})
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    monkeypatch.chdir(elsewhere)

    curvature = get_curvature(STRAINS, deformation_folder=deformation, legacy_oszicar=True)

    assert curvature == pytest.approx(CURVATURE_GPA, rel=1e-6)
    assert (deformation / "energies.txt").exists()
    assert (deformation / "curvature_constant").exists()
    assert list(elsewhere.iterdir()) == []


def test_cli_workdir_must_exist(tmp_path):
    with pytest.raises(SystemExit) as exc:
        main(["--operations", "get_elastic", "--workdir", str(tmp_path / "nope")])
    assert exc.value.code == 2


def test_cli_workdir_errors_name_the_project_folder(tmp_path, monkeypatch, caplog):
    project = tmp_path / "project"
    project.mkdir()
    monkeypatch.chdir(tmp_path)
    assert main(["--operations", "get_elastic", "--workdir", str(project), "--log-file", ""]) == 1
    assert str(project / "strains.txt") in caplog.text
