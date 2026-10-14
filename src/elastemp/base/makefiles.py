"""Default VASP input files, written only when the user has not supplied their own."""

import logging

from elastemp.base.utils import resolve_workspace, write_text

logger = logging.getLogger(__name__)

DEFAULT_KPOINT_MESH = (12, 12, 12)

INCAR_STATIC = (
    "PREC = Accurate\nEDIFF = 1e-6\nISMEAR = -5\nSIGMA = 0.060\nPOTIM = 0.020\n"
    "ISTART = 0\nNCORE=4\nLCHARG = FALSE\nLWAVE = FALSE\n"
)

INCAR_DYNAMIC = (
    "PREC=Accurate\nIBRION=8\nEDIFF=1e-8\nIALGO=38\nISMEAR=0\nSIGMA=0.06\n"
    "LREAL=.FALSE.\nADDGRID=.TRUE.\nLWAVE=.FALSE.\nLCHARG=.FALSE.\n"
)


def _gamma_kpoints(mesh=DEFAULT_KPOINT_MESH):
    kx, ky, kz = mesh
    return f"Automatic\n0\nGamma\n{kx:2d} {ky:2d} {kz:2d}\n 0. 0. 0.\n"


def _write_default(workspace, name, text, note=""):
    path = resolve_workspace(workspace) / name
    logger.warning("%s not found; writing a generic %s.%s", name, path, note)
    write_text(path, text)


def make_kpt_static(workspace=None):
    """Write a default Gamma-centred KPOINTS file for static calculations into ``workspace``."""
    _write_default(workspace, "KPOINTS", _gamma_kpoints(), " Converge this for production runs.")


def make_kpt_dynamic(workspace=None):
    """Write a default Gamma-centred KPOINTS_dynamic file for phonon (DFPT) calculations."""
    _write_default(workspace, "KPOINTS_dynamic", _gamma_kpoints(), " Converge this for production runs.")


def make_incar_static(workspace=None):
    """Write a default INCAR for static calculations into ``workspace``."""
    _write_default(workspace, "INCAR", INCAR_STATIC)


def make_incar_dynamic(workspace=None):
    """Write a default INCAR_dynamic for DFPT phonon calculations into ``workspace``."""
    _write_default(workspace, "INCAR_dynamic", INCAR_DYNAMIC)
