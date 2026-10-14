"""Symmetry-dependent strain deformations used to extract zero-temperature elastic constants.

For each crystal class, ``STRAIN_PATTERNS`` lists one Voigt strain pattern per
independent elastic constant. A pattern is scaled by each strain value ``e`` to
give the applied strain vector. The comment next to each pattern is the
combination of C_ij that the energy-strain curvature yields (see
``elastemp.elastic.constants``). A final hydrostatic pattern is always added for the
bulk-modulus fit and written to ``deformation-bulk``.
"""

import logging
import shutil

from pymatgen.analysis.elasticity.strain import Strain

from elastemp.base.strain_analysis import get_strain_matrix
from elastemp.base.symmetry import NUM_CONSTANTS
from elastemp.base.utils import deformation_dir, resolve_workspace, strain_dir

logger = logging.getLogger(__name__)

BULK_PATTERN = (1, 1, 1, 0, 0, 0)

STRAIN_PATTERNS = {
    "cubic": (
        (1, 1, 0, 0, 0, 0),  # C11 + C12
        (1, 1, 1, 0, 0, 0),  # 3/2 (C11 + 2 C12)
        (0, 0, 0, 1, 1, 1),  # 3/2 C44
    ),
    "hexagonal": (
        (1, 1, 0, 0, 0, 0),  # C11 + C12
        (0, 0, 0, 0, 0, 1),  # 1/4 (C11 - C12)
        (0, 0, 1, 0, 0, 0),  # 1/2 C33
        (0, 0, 0, 1, 1, 0),  # C44
        (1, 1, 1, 0, 0, 0),  # C11 + C12 + 2 C13 + C33/2
    ),
    "trigonal1": (
        (1, 1, 0, 0, 0, 0),  # C11 + C12
        (0, 0, 0, 0, 0, 1),  # 1/4 (C11 - C12)
        (0, 0, 1, 0, 0, 0),  # 1/2 C33
        (0, 0, 0, 1, 1, 0),  # C44
        (1, 1, 1, 0, 0, 0),  # C11 + C12 + 2 C13 + C33/2
        (0, 0, 0, 0, 1, 1),  # C11/4 - C12/4 + C14 + C44/2
    ),
    "trigonal2": (
        (1, 1, 0, 0, 0, 0),  # C11 + C12
        (0, 0, 0, 0, 0, 1),  # 1/4 (C11 - C12)
        (0, 0, 1, 0, 0, 0),  # 1/2 C33
        (0, 0, 0, 1, 1, 0),  # C44
        (1, 1, 1, 0, 0, 0),  # C11 + C12 + 2 C13 + C33/2
        (0, 0, 0, 0, 1, 1),  # C11/4 - C12/4 + C14 + C44/2
        (0, 0, 0, 1, 0, 1),  # C11/4 - C12/4 - C15 + C44/2
    ),
    "tetragonal1": (
        (1, 1, 0, 0, 0, 0),  # C11 + C12
        (0, 0, 0, 0, 0, 1),  # C66/2
        (0, 0, 1, 0, 0, 0),  # C33/2
        (0, 0, 0, 1, 1, 0),  # C44
        (1, 1, 1, 0, 0, 0),  # C11 + C12 + 2 C13 + C33/2
        (0, 1, 1, 0, 0, 0),  # C11/2 + C13 + C33/2
    ),
    "tetragonal2": (
        (1, 1, 0, 0, 0, 0),  # C11 + C12
        (0, 0, 0, 0, 0, 1),  # C66/2
        (0, 0, 1, 0, 0, 0),  # C33/2
        (0, 0, 0, 1, 1, 0),  # C44
        (1, 1, 1, 0, 0, 0),  # C11 + C12 + 2 C13 + C33/2
        (0, 1, 1, 0, 0, 0),  # C11/2 + C13 + C33/2
        (1, 0, 0, 0, 0, 1),  # C11/2 + C16 + C66/2
    ),
    "orthorhombic": (
        (1, 0, 0, 0, 0, 0),  # C11/2
        (0, 1, 0, 0, 0, 0),  # C22/2
        (0, 0, 1, 0, 0, 0),  # C33/2
        (0, 0, 0, 1, 0, 0),  # C44/2
        (0, 0, 0, 0, 1, 0),  # C55/2
        (0, 0, 0, 0, 0, 1),  # C66/2
        (1, 1, 0, 0, 0, 0),  # C11/2 + C12 + C22/2
        (1, 0, 1, 0, 0, 0),  # C11/2 + C13 + C33/2
        (0, 1, 1, 0, 0, 0),  # C22/2 + C23 + C33/2
    ),
    "monoclinic": (
        (1, 0, 0, 0, 0, 0),  # C11/2
        (0, 1, 0, 0, 0, 0),  # C22/2
        (0, 0, 1, 0, 0, 0),  # C33/2
        (0, 0, 0, 1, 0, 0),  # C44/2
        (0, 0, 0, 0, 1, 0),  # C55/2
        (0, 0, 0, 0, 0, 1),  # C66/2
        (1, 1, 0, 0, 0, 0),  # C11/2 + C12 + C22/2
        (1, 0, 1, 0, 0, 0),  # C11/2 + C13 + C33/2
        (1, 0, 0, 0, 1, 0),  # C11/2 + C15 + C55/2
        (0, 1, 1, 0, 0, 0),  # C22/2 + C23 + C33/2
        (0, 1, 0, 0, 1, 0),  # C22/2 + C25 + C55/2
        (0, 0, 1, 0, 1, 0),  # C33/2 + C35 + C55/2
        (0, 0, 0, 1, 0, 1),  # C44/2 + C46 + C66/2
    ),
    "triclinic": (
        (1, 0, 0, 0, 0, 0),  # C11/2
        (0, 1, 0, 0, 0, 0),  # C22/2
        (0, 0, 1, 0, 0, 0),  # C33/2
        (0, 0, 0, 1, 0, 0),  # C44/2
        (0, 0, 0, 0, 1, 0),  # C55/2
        (0, 0, 0, 0, 0, 1),  # C66/2
        (1, 1, 0, 0, 0, 0),  # C11/2 + C12 + C22/2
        (1, 0, 1, 0, 0, 0),  # C11/2 + C13 + C33/2
        (1, 0, 0, 1, 0, 0),  # C11/2 + C14 + C44/2
        (1, 0, 0, 0, 1, 0),  # C11/2 + C15 + C55/2
        (1, 0, 0, 0, 0, 1),  # C11/2 + C16 + C66/2
        (0, 1, 1, 0, 0, 0),  # C22/2 + C23 + C33/2
        (0, 1, 0, 1, 0, 0),  # C22/2 + C24 + C44/2
        (0, 1, 0, 0, 1, 0),  # C22/2 + C25 + C55/2
        (0, 1, 0, 0, 0, 1),  # C22/2 + C26 + C66/2
        (0, 0, 1, 1, 0, 0),  # C33/2 + C34 + C44/2
        (0, 0, 1, 0, 1, 0),  # C33/2 + C35 + C55/2
        (0, 0, 1, 0, 0, 1),  # C33/2 + C36 + C66/2
        (0, 0, 0, 1, 1, 0),  # C44/2 + C45 + C55/2
        (0, 0, 0, 1, 0, 1),  # C44/2 + C46 + C66/2
        (0, 0, 0, 0, 1, 1),  # C55/2 + C56 + C66/2
    ),
}

# Sanity check: one pattern per independent constant.
assert all(len(STRAIN_PATTERNS[s]) == n for s, n in NUM_CONSTANTS.items())


def get_strain_vector(symmetry, k, strain_val):
    """Voigt strain vector for deformation ``k`` (1-based; ``num_constants + 1`` is bulk)."""
    patterns = STRAIN_PATTERNS[symmetry] + (BULK_PATTERN,)
    return [strain_val * x for x in patterns[k - 1]]


def _write_strained_inputs(struct, pattern, strain_values, folder, workspace):
    """Create ``folder/strain-i`` for every strain value with a strained POSCAR and VASP inputs.

    INCAR, POTCAR and KPOINTS are copied from ``workspace``.
    """
    if folder.exists():
        shutil.rmtree(folder)
    folder.mkdir()
    for i, e in enumerate(strain_values, start=1):
        target = folder / strain_dir(i)
        target.mkdir()
        strain = Strain(get_strain_matrix([e * x for x in pattern]))
        new_struct = strain.get_deformation_matrix().apply_to_structure(struct)
        new_struct.to(fmt="poscar", filename=str(target / "POSCAR"))
        for name in ("INCAR", "POTCAR", "KPOINTS"):
            shutil.copy(workspace / name, target / name)


def make_elastic_deformations(struct, symmetry, strain_values, workspace=None):
    """Create ``deformation-1 .. deformation-n`` and ``deformation-bulk`` in ``workspace``.

    Existing deformation folders are deleted and recreated.
    """
    workspace = resolve_workspace(workspace)
    patterns = STRAIN_PATTERNS[symmetry]
    logger.info("Making %d deformations + bulk for %s symmetry", len(patterns), symmetry)
    for k, pattern in enumerate(patterns, start=1):
        _write_strained_inputs(struct, pattern, strain_values, workspace / deformation_dir(k), workspace)
    _write_strained_inputs(struct, BULK_PATTERN, strain_values, workspace / deformation_dir("bulk"), workspace)
    logger.info("Finished making deformations")
