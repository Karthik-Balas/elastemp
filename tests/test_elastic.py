"""Physics self-consistency tests for the elastic-constant assembly (no VASP/phonopy needed)."""

import numpy as np
import pytest

from elastemp.base.symmetry import NUM_CONSTANTS
from elastemp.elastic.constants import stiffness_from_curvatures
from elastemp.elastic.deformations import STRAIN_PATTERNS
from elastemp.elastic.response import compute_response_functions
from elastemp.elastic.stability import evaluate_stability


def _sym(c):
    return np.triu(c) + np.triu(c, 1).T


# Physically reasonable stiffness matrices (GPa) obeying each Laue class's constraints.
def _reference_stiffness(symmetry):
    c = np.zeros((6, 6))
    c11, c12, c13, c33, c44, c14, c15, c16 = 300.0, 120.0, 90.0, 280.0, 80.0, 15.0, 10.0, 12.0
    if symmetry == "cubic":
        c[0, 0] = c[1, 1] = c[2, 2] = c11
        c[0, 1] = c[0, 2] = c[1, 2] = c12
        c[3, 3] = c[4, 4] = c[5, 5] = c44
    elif symmetry in ("hexagonal", "trigonal1", "trigonal2"):
        c[0, 0] = c[1, 1] = c11
        c[0, 1] = c12
        c[0, 2] = c[1, 2] = c13
        c[2, 2] = c33
        c[3, 3] = c[4, 4] = c44
        c[5, 5] = (c11 - c12) / 2
        if symmetry != "hexagonal":
            c[0, 3], c[1, 3], c[4, 5] = c14, -c14, c14
        if symmetry == "trigonal2":
            c[0, 4], c[1, 4], c[3, 5] = -c15, c15, c15
    elif symmetry in ("tetragonal1", "tetragonal2"):
        c[0, 0] = c[1, 1] = c11
        c[0, 1] = c12
        c[0, 2] = c[1, 2] = c13
        c[2, 2] = c33
        c[3, 3] = c[4, 4] = c44
        c[5, 5] = 100.0
        if symmetry == "tetragonal2":
            c[0, 5], c[1, 5] = c16, -c16
    else:
        rng = np.random.default_rng(0)
        c = np.diag([300.0, 280.0, 260.0, 90.0, 80.0, 70.0])
        nonzero = {
            "orthorhombic": [(0, 1), (0, 2), (1, 2)],
            "monoclinic": [(0, 1), (0, 2), (0, 4), (1, 2), (1, 4), (2, 4), (3, 5)],
            "triclinic": [(i, j) for i in range(6) for j in range(i + 1, 6)],
        }[symmetry]
        for i, j in nonzero:
            c[i, j] = rng.uniform(5, 40)
    return _sym(c)


@pytest.mark.parametrize("symmetry", sorted(NUM_CONSTANTS))
def test_one_pattern_per_constant(symmetry):
    assert len(STRAIN_PATTERNS[symmetry]) == NUM_CONSTANTS[symmetry]


@pytest.mark.parametrize("symmetry", sorted(NUM_CONSTANTS))
def test_patterns_and_assembly_recover_stiffness(symmetry):
    """E(e) = V/2 e^T C e for each pattern; the fitted curvature is (1/2) p^T C p.

    Feeding those exact curvatures back must reproduce C.
    """
    c_ref = _reference_stiffness(symmetry)
    curvatures = [0.5 * np.asarray(p) @ c_ref @ np.asarray(p) for p in STRAIN_PATTERNS[symmetry]]
    c = stiffness_from_curvatures(symmetry, curvatures)
    np.testing.assert_allclose(c, c_ref, atol=1e-9)


def test_voigt_reuss_hill_cubic():
    c = _reference_stiffness("cubic")
    resp = compute_response_functions(c, "cubic")
    Eh, Gh, hc, nu, pugh = resp.get_moduli(resp.get_compliance_matrix(False), False)
    c11, c12, c44 = c[0, 0], c[0, 1], c[3, 3]
    gv = (c11 - c12 + 3 * c44) / 5
    gr = 5 * (c11 - c12) * c44 / (4 * c44 + 3 * (c11 - c12))
    assert Gh == pytest.approx((gv + gr) / 2, abs=1e-3)


def test_voigt_equals_reuss_for_isotropic_solid():
    # Isotropic: C44 = (C11 - C12) / 2, so Voigt and Reuss shear moduli coincide.
    c = np.zeros((6, 6))
    c[:3, :3] = 120.0
    np.fill_diagonal(c, [300.0, 300.0, 300.0, 90.0, 90.0, 90.0])
    resp = compute_response_functions(c, "cubic")
    _, Gh, _, _, _ = resp.get_moduli(resp.get_compliance_matrix(False), False)
    assert Gh == pytest.approx(90.0, abs=1e-3)


@pytest.mark.parametrize("symmetry", sorted(NUM_CONSTANTS))
def test_reference_matrices_are_stable(symmetry):
    results = evaluate_stability(_reference_stiffness(symmetry), symmetry)
    assert all(met for _, met in results), results


def test_unstable_cubic_detected():
    c = _reference_stiffness("cubic")
    c[3, 3] = c[4, 4] = c[5, 5] = -10.0
    results = dict(evaluate_stability(c, "cubic"))
    assert results["C44 > 0"] is False
