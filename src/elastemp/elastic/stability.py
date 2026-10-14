"""Born mechanical-stability criteria (Mouhat & Coudert, PRB 90, 224104, 2014)."""

import logging

import numpy as np

from elastemp.base.utils import resolve_workspace, results_dir

logger = logging.getLogger(__name__)


def _C(c, i, j):
    """1-based Voigt accessor: _C(c, 1, 2) == C12."""
    return c[i - 1][j - 1]


def _minor_det(c, idx):
    """Determinant of the principal sub-matrix of C with 1-based Voigt indices ``idx``."""
    idx0 = [i - 1 for i in idx]
    return np.linalg.det(np.asarray(c, dtype=float)[np.ix_(idx0, idx0)])


def _positive_definite(c):
    return bool(np.all(np.linalg.eigvalsh(np.asarray(c, dtype=float)) > 0))


_HEX_TET_COMMON = [
    ("C11 > |C12|", lambda c: _C(c, 1, 1) > abs(_C(c, 1, 2))),
    ("2*C13^2 < C33*(C11+C12)", lambda c: 2 * _C(c, 1, 3) ** 2 < _C(c, 3, 3) * (_C(c, 1, 1) + _C(c, 1, 2))),
    ("C44 > 0", lambda c: _C(c, 4, 4) > 0),
]

_TRIGONAL_COMMON = [
    ("C11 > |C12|", lambda c: _C(c, 1, 1) > abs(_C(c, 1, 2))),
    ("C44 > 0", lambda c: _C(c, 4, 4) > 0),
    ("C13^2 < 1/2*C33*(C11+C12)", lambda c: _C(c, 1, 3) ** 2 < 0.5 * _C(c, 3, 3) * (_C(c, 1, 1) + _C(c, 1, 2))),
]

# symmetry -> list of (description, predicate on the 6x6 stiffness matrix)
CRITERIA = {
    "cubic": [
        ("C11 - C12 > 0", lambda c: _C(c, 1, 1) > _C(c, 1, 2)),
        ("C11 + 2*C12 > 0", lambda c: _C(c, 1, 1) + 2 * _C(c, 1, 2) > 0),
        ("C44 > 0", lambda c: _C(c, 4, 4) > 0),
    ],
    "hexagonal": _HEX_TET_COMMON + [("C66 > 0", lambda c: _C(c, 6, 6) > 0)],
    "tetragonal1": _HEX_TET_COMMON + [("C66 > 0", lambda c: _C(c, 6, 6) > 0)],
    "tetragonal2": _HEX_TET_COMMON
    + [
        ("2*C16^2 < C66*(C11-C12)", lambda c: 2 * _C(c, 1, 6) ** 2 < _C(c, 6, 6) * (_C(c, 1, 1) - _C(c, 1, 2))),
    ],
    "trigonal1": _TRIGONAL_COMMON
    + [
        ("C14^2 < 1/2*C44*(C11-C12)", lambda c: _C(c, 1, 4) ** 2 < 0.5 * _C(c, 4, 4) * (_C(c, 1, 1) - _C(c, 1, 2))),
    ],
    "trigonal2": _TRIGONAL_COMMON
    + [
        (
            "C14^2 + C15^2 < 1/2*C44*(C11-C12)",
            lambda c: _C(c, 1, 4) ** 2 + _C(c, 1, 5) ** 2 < 0.5 * _C(c, 4, 4) * (_C(c, 1, 1) - _C(c, 1, 2)),
        ),
    ],
    "orthorhombic": [
        ("C11 > 0", lambda c: _C(c, 1, 1) > 0),
        ("C11*C22 > C12^2", lambda c: _C(c, 1, 1) * _C(c, 2, 2) > _C(c, 1, 2) ** 2),
        ("det[C11..C33 block] > 0", lambda c: _minor_det(c, [1, 2, 3]) > 0),
        ("C44 > 0", lambda c: _C(c, 4, 4) > 0),
        ("C55 > 0", lambda c: _C(c, 5, 5) > 0),
        ("C66 > 0", lambda c: _C(c, 6, 6) > 0),
    ],
    "monoclinic": [
        ("C11 > 0", lambda c: _C(c, 1, 1) > 0),
        ("C22 > 0", lambda c: _C(c, 2, 2) > 0),
        ("C33 > 0", lambda c: _C(c, 3, 3) > 0),
        ("C44 > 0", lambda c: _C(c, 4, 4) > 0),
        ("C55 > 0", lambda c: _C(c, 5, 5) > 0),
        ("C66 > 0", lambda c: _C(c, 6, 6) > 0),
        (
            "C11+C22+C33+2*(C12+C13+C23) > 0",
            lambda c: _C(c, 1, 1) + _C(c, 2, 2) + _C(c, 3, 3) + 2 * (_C(c, 1, 2) + _C(c, 1, 3) + _C(c, 2, 3)) > 0,
        ),
        ("C33*C55 - C35^2 > 0", lambda c: _C(c, 3, 3) * _C(c, 5, 5) - _C(c, 3, 5) ** 2 > 0),
        ("C44*C66 - C46^2 > 0", lambda c: _C(c, 4, 4) * _C(c, 6, 6) - _C(c, 4, 6) ** 2 > 0),
        ("C22 + C33 - 2*C23 > 0", lambda c: _C(c, 2, 2) + _C(c, 3, 3) - 2 * _C(c, 2, 3) > 0),
        (
            "C22*(C33*C55-C35^2) + 2*C23*C25*C35 - C23^2*C55 - C25^2*C33 > 0  (det of {2,3,5} block)",
            lambda c: _minor_det(c, [2, 3, 5]) > 0,
        ),
        ("det of {1,2,3,5} block > 0", lambda c: _minor_det(c, [1, 2, 3, 5]) > 0),
    ],
    "triclinic": [
        ("C positive definite (all eigenvalues > 0)", _positive_definite),
    ],
}

GENERAL_CRITERION = ("Born criterion: C positive definite (all eigenvalues > 0)", _positive_definite)


def evaluate_stability(c, symmetry):
    """Evaluate the stability criteria for ``symmetry``.

    :returns: list of (description, met) tuples. The last entry is always the general
              positive-definiteness criterion, which is necessary and sufficient for any symmetry.
    """
    criteria = list(CRITERIA[symmetry])
    if criteria[-1][1] is not _positive_definite:
        criteria.append(GENERAL_CRITERION)
    return [(desc, bool(pred(c))) for desc, pred in criteria]


class stability_check:
    """Check mechanical stability and write results_dir/mechanical_stability_criteria.txt."""

    def __init__(self, c, symmetry, workspace=None):
        """
        :param c: 6x6 stiffness matrix (GPa)
        :param symmetry: symmetry label of the structure
        """
        self.c = c
        self.symmetry = symmetry.strip()
        self.workspace = resolve_workspace(workspace)
        self.results = self.check_stability()

    @property
    def is_stable(self):
        return all(met for _, met in self.results)

    def check_stability(self):
        """Evaluate criteria, write them to file, and return list of (description, met)."""
        results = evaluate_stability(self.c, self.symmetry)
        with open(results_dir(self.workspace) / "mechanical_stability_criteria.txt", "w") as f:
            f.write("**" * 30 + "\n")
            f.write(f"Detected symmetry of structure = {self.symmetry} \n\n")
            f.write(f"Criteria for mechanical stability of {self.symmetry} structures : \n\n")
            for n, (desc, _) in enumerate(results, start=1):
                f.write(f"condition {n} = {desc}\n\n")
            f.write("---" * 30 + "\n")
            for n, (_, met) in enumerate(results, start=1):
                f.write("condition {} constraint:\t {}\n".format(n, "met" if met else "NOT met"))
            f.write("**" * 30 + "\n")

        failed = [desc for desc, met in results if not met]
        if failed:
            logger.warning("Mechanical stability criteria NOT met: %s", "; ".join(failed))
        else:
            logger.info("All %d mechanical stability criteria met", len(results))
        return results
