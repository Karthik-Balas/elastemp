"""Steps 1-2: set up and analyse the zero-temperature elastic-constant calculations."""

import logging

from pymatgen.core import Structure

from elastemp.base.makefiles import make_incar_static, make_kpt_static
from elastemp.base.strain_analysis import MIN_FIT_POINTS
from elastemp.base.symmetry import get_num_constants, get_symmetry
from elastemp.base.utils import (
    STRAINS_FILE,
    SYMMETRY_FILE,
    read_strains,
    read_symmetry,
    require_file,
    resolve_workspace,
    results_dir,
    write_text,
)
from elastemp.elastic.constants import compute_elastic_constants
from elastemp.elastic.deformations import make_elastic_deformations
from elastemp.elastic.response import compute_response_functions
from elastemp.elastic.stability import stability_check
from elastemp.exceptions import ElastempError, MissingInputError

logger = logging.getLogger(__name__)

# 11 points: one or two failed runs can be dropped and still leave the 7-point minimum.
DEFAULT_STRAINS = [-0.05, -0.04, -0.03, -0.02, -0.01, 0, 0.01, 0.02, 0.03, 0.04, 0.05]


def make_strains(workspace=None):
    """Write the default strains.txt into ``workspace`` and return the strain list."""
    write_text(resolve_workspace(workspace) / STRAINS_FILE, "".join(f"{s}\n" for s in DEFAULT_STRAINS))
    return list(DEFAULT_STRAINS)


def load_or_create_strains(workspace=None):
    """Read strains.txt (validating it), or write and return the defaults if it is missing/unreadable."""
    try:
        strain_values = read_strains(workspace)
    except (MissingInputError, ValueError, IndexError):
        logger.warning(
            "Unable to read %s (missing or corrupt); using default strains %s",
            STRAINS_FILE,
            DEFAULT_STRAINS,
        )
        return make_strains(workspace)
    if len(strain_values) < MIN_FIT_POINTS:
        raise ElastempError(
            f"{STRAINS_FILE} has {len(strain_values)} strain values; at least {MIN_FIT_POINTS} are required "
            "because get_elastic needs that many converged points per deformation. "
            "Use more, so that one or two failed runs can be dropped."
        )
    if len(strain_values) == MIN_FIT_POINTS:
        logger.warning(
            "%s has exactly %d strains: if any strain point fails to converge, get_elastic will stop. "
            "Consider adding a few more.",
            STRAINS_FILE,
            MIN_FIT_POINTS,
        )
    if min(strain_values) > 0:
        logger.warning("Some strain values should be negative (ideally about half of them)")
    return strain_values


def make_elastic_wrapper(symprec=0.05, angle_tolerance=10, workspace=None):
    """Step 1: detect symmetry and create the deformation folders.

    Writes default INCAR/KPOINTS/strains.txt if they are missing.

    :param symprec: tolerance in lattice parameters for symmetry detection
    :param angle_tolerance: tolerance in lattice angles for symmetry detection
    :param workspace: project folder containing POSCAR and POTCAR (default: current directory)
    """
    ws = resolve_workspace(workspace)
    require_file(ws / "POSCAR")
    require_file(ws / "POTCAR")
    struct = Structure.from_file(ws / "POSCAR")
    logger.info("Read POSCAR (%s, %d atoms)", struct.composition.reduced_formula, len(struct))

    strain_values = load_or_create_strains(ws)
    if not (ws / "INCAR").exists():
        make_incar_static(ws)
    if not (ws / "KPOINTS").exists():
        make_kpt_static(ws)
    make_elastic_deformation(struct, strain_values, symprec, angle_tolerance, ws)


def write_symmetry(struct, symprec, angle_tolerance, workspace=None):
    """Detect the symmetry of ``struct`` and write it to the ``symmetry`` file."""
    symmetry = get_symmetry(struct, symprec, angle_tolerance)
    write_text(resolve_workspace(workspace) / SYMMETRY_FILE, symmetry)
    return symmetry


def make_elastic_deformation(struct, strain_values, symprec, angle_tolerance, workspace=None):
    """Create symmetry-dependent deformation folders (reuses an existing ``symmetry`` file)."""
    ws = resolve_workspace(workspace)
    if not (ws / SYMMETRY_FILE).exists():
        write_symmetry(struct, symprec, angle_tolerance, ws)
    symmetry = read_symmetry(ws)
    logger.info("Symmetry: %s (%d independent elastic constants)", symmetry, get_num_constants(symmetry))
    make_elastic_deformations(struct, symmetry, strain_values, ws)
    logger.info("Run VASP in every deformation-*/strain-* folder, then run get_elastic")


def get_elastic_wrapper(make_plots, workspace=None):
    """Step 2: compute zero-temperature elastic constants from finished VASP runs."""
    get_elastic(read_strains(workspace), make_plots, workspace)


def get_elastic(strain_values, make_plots, workspace=None):
    """Fit energy-strain curves, then write stiffness, compliance, moduli and stability files.

    :returns: (stiffness matrix, stability_check)
    """
    ws = resolve_workspace(workspace)
    symmetry = read_symmetry(ws)
    num_constants = get_num_constants(symmetry)
    e = compute_elastic_constants(num_constants, strain_values, symmetry, workspace=ws)
    e.get_bulk_curvatures()
    curvatures = e.get_curvatures()

    logger.info("Computing elastic constants")
    c = e.get_elastic_constants(curvatures, write_elastic=True)
    e.write_zero_elastic_constants(c)
    if make_plots:
        e.get_plots()

    logger.info("Computing response functions")
    response = compute_response_functions(c, symmetry, workspace=ws)
    s = response.get_compliance_matrix(write_compliance=True)
    response.get_moduli(s, write_moduli=True)

    stability = stability_check(c, symmetry, workspace=ws)
    logger.info("Results written to %s", results_dir(ws))
    return c, stability
