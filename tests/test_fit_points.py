"""Dropping unconverged strain points and refitting."""

import numpy as np
import pandas as pd
import pytest
from pymatgen.core import Lattice, Structure

from elastemp.base.strain_analysis import MIN_FIT_POINTS, get_curvature, select_fit_points
from elastemp.base.utils import EV_PER_A3_TO_GPA
from elastemp.exceptions import ConvergenceError

STRAINS = [-0.05, -0.04, -0.03, -0.02, -0.01, 0.0, 0.01, 0.02, 0.03, 0.04, 0.05]
CURVATURE_GPA = 400.0


def test_min_points_is_seven():
    assert MIN_FIT_POINTS == 7


def test_all_converged_keeps_everything():
    assert select_fit_points(STRAINS, [True] * 11, "d") == list(range(11))


def test_drops_unconverged_points():
    converged = [True] * 11
    converged[3] = converged[8] = False
    assert select_fit_points(STRAINS, converged, "d") == [0, 1, 2, 4, 5, 6, 7, 9, 10]


def test_exactly_seven_left_is_allowed():
    converged = [True] * 7 + [False] * 4  # keeps -0.05 .. 0.01
    assert len(select_fit_points(STRAINS, converged, "d")) == 7


def test_six_left_raises():
    converged = [True] * 6 + [False] * 5
    with pytest.raises(ConvergenceError, match="only 6 of 11"):
        select_fit_points(STRAINS, converged, "d")


def test_one_sided_strains_raise():
    strains = [-0.05, -0.04, -0.03, -0.02, -0.01, 0.0, 0.0, 0.01]
    converged = [True] * 7 + [False]
    with pytest.raises(ConvergenceError, match="compression\\s+and tension|both compression"):
        select_fit_points(strains, converged, "d")


def _make_deformation(tmp_path, unconverged=()):
    """A fake deformation folder whose energies lie exactly on a parabola."""
    struct = Structure(Lattice.cubic(3.57), ["C"], [[0, 0, 0]])
    struct.to(filename=str(tmp_path / "POSCAR"), fmt="poscar")
    deformation = tmp_path / "deformation-1"
    deformation.mkdir()
    V0 = struct.volume
    for i, e in enumerate(STRAINS, start=1):
        folder = deformation / f"strain-{i}"
        folder.mkdir()
        struct.to(filename=str(folder / "POSCAR"), fmt="poscar")
        energy = -10.0 + CURVATURE_GPA * V0 / EV_PER_A3_TO_GPA * e**2
        if i in unconverged:
            last = "  60 RMM:  -0.10E+02   0.1E-03  0.1E-04   100   0.2E-02\n"
        else:
            last = f"   1 F= {energy:.12E} E0= {energy:.12E}  d E =0.0E+00\n"
        (folder / "OSZICAR").write_text("header\n" + last)
    return deformation


def test_get_curvature_refits_without_unconverged_point(tmp_path, monkeypatch):
    deformation = _make_deformation(tmp_path, unconverged={4})
    monkeypatch.chdir(deformation)
    curvature = get_curvature(STRAINS, legacy_oszicar=True)
    assert curvature == pytest.approx(CURVATURE_GPA, rel=1e-6)

    table = pd.read_csv("energies.txt", sep="\t")
    assert list(table.columns) == ["energies", "strains", "converged"]
    assert table["converged"].tolist() == [i != 3 for i in range(11)]
    assert np.isnan(table["energies"][3])


def test_get_curvature_raises_when_too_few_points(tmp_path, monkeypatch):
    deformation = _make_deformation(tmp_path, unconverged={1, 2, 3, 4, 5})
    monkeypatch.chdir(deformation)
    with pytest.raises(ConvergenceError, match="deformation-1: only 6 of 11"):
        get_curvature(STRAINS, legacy_oszicar=True)
