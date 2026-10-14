"""Shared helpers: project-folder paths, input-file readers, and the phonopy runner.

Every function that touches the file system takes an explicit ``workspace`` (the project
folder containing POSCAR, strains.txt, deformation-*, results_dir, ...). Nothing in elastemp
changes the process working directory, so independent projects can be processed side by side.
``workspace=None`` means the current directory.
"""

import functools
import logging
import subprocess
from pathlib import Path

import numpy as np

from elastemp.exceptions import ExternalToolError, MissingInputError

logger = logging.getLogger(__name__)

# Temperature step (K) written to mesh.conf. The temperature grid of every
# downstream table (energy_temp.dat, volume_temp.txt, *_temp_*.csv) follows it.
TSTEP = 2

# Unit conversion: eV/Angstrom^3 -> GPa
EV_PER_A3_TO_GPA = 160.217

STRAINS_FILE = "strains.txt"
DIM_FILE = "dim.txt"
TMAX_FILE = "tmax.txt"
SYMMETRY_FILE = "symmetry"
RESULTS_DIR = "results_dir"


def resolve_workspace(workspace=None):
    """Absolute Path of ``workspace``; the current directory if ``None``."""
    return Path(workspace).expanduser().resolve() if workspace is not None else Path.cwd()


def deformation_dir(index):
    """Folder name for the ``index``-th deformation (1-based) or ``'bulk'``."""
    return f"deformation-{index}"


def strain_dir(index):
    """Folder name for the ``index``-th strain (1-based)."""
    return f"strain-{index}"


def read_strains(workspace=None):
    """Read strains.txt from ``workspace`` and return a list of floats."""
    path = resolve_workspace(workspace) / STRAINS_FILE
    try:
        with open(path) as f:
            return [float(line.split()[0]) for line in f if line.strip()]
    except FileNotFoundError as exc:
        raise MissingInputError(f"{path} not found") from exc


def read_dim(workspace=None):
    """Read the supercell dimension (three integers) from dim.txt in ``workspace``."""
    path = resolve_workspace(workspace) / DIM_FILE
    try:
        with open(path) as f:
            dim = [int(x) for x in f.readline().split()]
    except FileNotFoundError as exc:
        raise MissingInputError(f"{path} not found") from exc
    if len(dim) != 3:
        raise ValueError(f"{path} must contain three integers, got {dim}")
    return dim


def read_tmax(workspace=None):
    """Read the maximum temperature (K) from tmax.txt in ``workspace``."""
    path = resolve_workspace(workspace) / TMAX_FILE
    try:
        with open(path) as f:
            return int(f.readline().strip())
    except FileNotFoundError as exc:
        raise MissingInputError(f"{path} not found") from exc


def read_symmetry(workspace=None):
    """Read the crystal-symmetry label written by ``make_elastic`` in ``workspace``."""
    path = resolve_workspace(workspace) / SYMMETRY_FILE
    try:
        with open(path) as f:
            return f.readline().strip()
    except FileNotFoundError as exc:
        raise MissingInputError(f"{path} not found. Run the make_elastic step first.") from exc


def results_dir(workspace=None):
    """``workspace/results_dir``, created if needed."""
    path = resolve_workspace(workspace) / RESULTS_DIR
    path.mkdir(parents=True, exist_ok=True)
    return path


def write_text(path, text):
    """Write ``text`` to ``path`` (overwrites)."""
    with open(path, "w") as f:
        f.write(text)


def format_matrix(m, decimals):
    """Tab-separated text block of a 6x6 matrix, in the layout used by results_dir/*_matrix.txt."""
    rows = ["\t".join(str(np.round(m[i][j], decimals)) for j in range(6)) + "\t" for i in range(6)]
    return " " + "\n".join(rows) + "\n"


def require_file(path, message=None):
    """Raise MissingInputError if ``path`` does not exist."""
    if not Path(path).exists():
        raise MissingInputError(message or f"{Path(path).resolve()} not found")


@functools.cache
def _check_phonopy_version():
    try:
        import phonopy
    except ImportError as exc:
        raise ExternalToolError("phonopy is not installed") from exc
    major = int(phonopy.__version__.split(".")[0])
    if major >= 4:
        raise ExternalToolError(
            f"phonopy {phonopy.__version__} detected. elastemp uses the phonopy v2/v3 command-line interface; "
            "install 'phonopy<4'."
        )
    return phonopy.__version__


def run_phonopy(*args, cwd=None):
    """Run the phonopy CLI in an explicit subprocess directory.

    Output is captured and logged at DEBUG level. A non-zero exit status raises
    ExternalToolError instead of being silently ignored.
    """
    _check_phonopy_version()
    workspace = resolve_workspace(cwd)
    cmd = ["phonopy", *args]
    logger.debug("Running %s in %s", " ".join(cmd), workspace)
    try:
        result = subprocess.run(cmd, cwd=workspace, capture_output=True, text=True)
    except FileNotFoundError as exc:
        raise ExternalToolError("phonopy executable not found on PATH") from exc
    if result.stdout:
        logger.debug(result.stdout)
    if result.returncode != 0:
        raise ExternalToolError(
            f"'{' '.join(cmd)}' failed in {workspace} (exit {result.returncode}):\n{result.stderr.strip()}"
        )
    return result
