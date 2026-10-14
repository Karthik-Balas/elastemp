"""Set up and post-process phonopy (DFPT) calculations for each strained structure."""

import logging
import shutil

from elastemp.base.symmetry import get_num_constants
from elastemp.base.utils import deformation_dir, require_file, resolve_workspace, run_phonopy, strain_dir
from elastemp.exceptions import ConvergenceError
from elastemp.validation.convergence import validate_vasp_run

logger = logging.getLogger(__name__)

PHONON_DIR = "Phonon_calculation"
HIGHTEMP_DIR = "Hightemp_static"

# Project-level input file -> name it gets inside each Phonon_calculation folder
_PHONON_INPUTS = (
    ("KPOINTS_dynamic", "KPOINTS"),
    ("INCAR_dynamic", "INCAR"),
    ("INCAR", "INCAR_static"),
    ("POTCAR", "POTCAR"),
)


class dynamic_deformations:
    """Create phonon supercells for every strain folder of a project."""

    def __init__(self, symmetry, strain_values, dim, struct, workspace=None):
        """
        :param symmetry: symmetry label of the structure
        :param strain_values: list of strain values
        :param dim: supercell dimension, e.g. [2, 2, 2]
        :param struct: reference pymatgen Structure
        :param workspace: project folder containing deformation-* and the VASP inputs
                          (default: current directory)
        """
        self.num_constants = get_num_constants(symmetry)
        self.strain_values = strain_values
        self.dim = dim
        self.struct = struct
        self.workspace = resolve_workspace(workspace)

    def _make_phonon_folder(self, strain_folder, source_subdir):
        """Create ``strain_folder/Phonon_calculation`` and run ``phonopy -d`` in it.

        POSCAR and ZPE.txt come from ``strain_folder/source_subdir``; KPOINTS_dynamic,
        INCAR_dynamic, INCAR and POTCAR come from the project folder. Every input is checked
        before an existing Phonon_calculation folder is deleted.
        """
        source = strain_folder / source_subdir
        zpe, poscar = source / "ZPE.txt", source / "POSCAR"
        require_file(
            zpe,
            f"{zpe} not found: the static calculation has not converged, or the matching "
            "get_* step (get_elastic / get_temp_static) has not been run",
        )
        require_file(poscar)
        for name, _ in _PHONON_INPUTS:
            require_file(self.workspace / name)

        phonon = strain_folder / PHONON_DIR
        if phonon.exists():
            shutil.rmtree(phonon)
        phonon.mkdir()
        shutil.copy(zpe, phonon / "ZPE.txt")
        shutil.copy(poscar, phonon / "POSCAR")

        run_phonopy("-d", "--dim={} {} {}".format(*self.dim), cwd=phonon)
        shutil.copy(phonon / "POSCAR", phonon / "POSCAR-unitcell")
        shutil.move(phonon / "SPOSCAR", phonon / "POSCAR")

        for name, destination in _PHONON_INPUTS:
            shutil.copy(self.workspace / name, phonon / destination)

    def make_phonon_deformation(self, strain_folder):
        """Phonon folder from the 0 K static calculation in ``strain_folder``."""
        self._make_phonon_folder(resolve_workspace(strain_folder), ".")

    def make_phonon_deformation_T(self, strain_folder):
        """Phonon folder from the volume-scaled (tmax) static calculation in ``strain_folder/Hightemp_static``."""
        self._make_phonon_folder(resolve_workspace(strain_folder), HIGHTEMP_DIR)

    def make_bulk_deformation(self):
        """Set up DFPT calculations for every strain of ``deformation-bulk``."""
        for i in range(1, len(self.strain_values) + 1):
            folder = self.workspace / deformation_dir("bulk") / strain_dir(i)
            logger.info("Setting up phonon calculation in %s", folder)
            self.make_phonon_deformation(folder)

    def make_deformation(self):
        """Set up DFPT calculations for every strain of every ``deformation-k``."""
        for k in range(1, self.num_constants + 1):
            for i in range(1, len(self.strain_values) + 1):
                folder = self.workspace / deformation_dir(k) / strain_dir(i)
                logger.info("Setting up phonon calculation in %s", folder)
                self.make_phonon_deformation_T(folder)


def extract_force_constants(phonon_folder=None):
    """Run ``phonopy --fc vasprun.xml`` in ``phonon_folder`` (default: current directory).

    Raises ConvergenceError if the DFPT run's SCF did not converge or vasprun.xml is unreadable.
    """
    folder = resolve_workspace(phonon_folder)
    require_file(folder / "vasprun.xml")
    result = validate_vasp_run(folder)
    if not result.converged:
        raise ConvergenceError(f"DFPT run cannot be used: {result.reason}")
    logger.info("Extracting force constants in %s", folder)
    run_phonopy("--fc", "vasprun.xml", cwd=folder)


class dynamic_response:
    """Backwards-compatible wrapper: instantiating it extracts force constants in ``phonon_folder``."""

    def __init__(self, phonon_folder=None):
        self.phonon_folder = phonon_folder
        self.extract_phonon()

    def extract_phonon(self):
        extract_force_constants(self.phonon_folder)
