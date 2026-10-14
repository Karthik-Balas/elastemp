"""Energy-strain analysis: parsing VASP output, fitting parabolas/EOS, plotting."""

import logging
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from pymatgen.analysis.eos import EOS
from pymatgen.core import Structure

from elastemp.base.utils import EV_PER_A3_TO_GPA, resolve_workspace, strain_dir, write_text
from elastemp.exceptions import ConvergenceError, MissingInputError

logger = logging.getLogger(__name__)

# Minimum number of converged strain points needed to fit an energy-strain curve or EOS.
MIN_FIT_POINTS = 7

# Returned by check_convergence in place of the last OSZICAR line when the file is unreadable.
UNCONVERGED_PLACEHOLDER = 1e4


def get_strain_matrix(vector):
    """Convert a 6-component Voigt strain vector to a symmetric 3x3 strain matrix."""
    e1, e2, e3, e4, e5, e6 = vector
    return np.array(
        [[e1, e6 / 2, e5 / 2], [e6 / 2, e2, e4 / 2], [e5 / 2, e4 / 2, e3]],
        dtype=float,
    )


def check_convergence(fil, *, require_ionic=False, legacy_oszicar=False):
    """Check whether the VASP run that wrote ``fil`` (its OSZICAR) converged.

    By default this validates ``vasprun.xml`` in the same folder (see
    ``elastemp.validation.convergence``). ``legacy_oszicar=True`` falls back to the old test
    for an ``F=`` line in OSZICAR, which only shows that an ionic step finished, not that the
    SCF converged. Use it only for old data without vasprun.xml.

    :returns: (converged, detail). ``detail`` is the failure reason (XML mode) or the last
              OSZICAR line (legacy mode).
    """
    from elastemp.validation.convergence import validate_vasp_run

    fil = Path(fil)
    if not legacy_oszicar:
        result = validate_vasp_run(fil.parent, require_ionic=require_ionic)
        if not result.converged:
            logger.warning("%s", result.reason)
        return result.converged, result.reason or "converged (validated with vasprun.xml)"

    try:
        lines = fil.read_text().splitlines()
    except OSError:
        logger.warning("Could not read %s", fil)
        return False, [UNCONVERGED_PLACEHOLDER] * 3
    if not lines:
        logger.warning("%s is empty", fil)
        return False, [UNCONVERGED_PLACEHOLDER] * 3
    if "F= " not in lines[-1]:
        logger.warning("%s has no final F= line", fil)
        return False, lines[-1]
    return True, lines[-1]


def fit_parabola(energies, strain_array, V0):
    """Fit E(strain) to a parabola.

    :returns: (curvature constant in GPa, polynomial coefficients [a, b, c])
    """
    p = np.polyfit(strain_array, energies, 2)
    curvature_const = np.round(p[0] * EV_PER_A3_TO_GPA / V0, 4)
    return curvature_const, p


def fit_parabola_bulk(energies, volumes, calc_bulk_parabola=True):
    """Fit E(V) to the Murnaghan equation of state.

    :param calc_bulk_parabola: if True, return the bulk modulus (GPa);
                               otherwise return the equilibrium volume.
    """
    eos_fit = EOS(eos_name="murnaghan").fit(volumes, energies).results
    if calc_bulk_parabola:
        return np.round(eos_fit["b0"] * EV_PER_A3_TO_GPA, 4)
    return eos_fit["v0"]


def get_energy_volume(dir_path):
    """Read the final energy from ``dir_path/OSZICAR`` and the volume from ``dir_path/POSCAR``.

    Side effect: writes the energy per atom to ``dir_path/ZPE.txt``; later phonon
    steps read it from there.

    :returns: (total energy in eV, volume in A^3)
    """
    dir_path = Path(dir_path)
    struct = Structure.from_file(dir_path / "POSCAR")
    with open(dir_path / "OSZICAR") as g:
        energy = float(g.readlines()[-1].split()[2])
    write_text(dir_path / "ZPE.txt", str(energy / len(struct.sites)))
    return energy, struct.volume


def select_fit_points(strain_array, converged, label, min_points=MIN_FIT_POINTS):
    """Return the indices of converged strain points, checking there are enough to fit.

    Unconverged points are dropped from the fit and logged. Raises ConvergenceError if fewer
    than ``min_points`` remain, or if the remaining strains are not on both sides of zero.

    :param strain_array: all strain values
    :param converged: one bool per strain value
    :param label: folder name used in messages
    """
    keep = [i for i, ok in enumerate(converged) if ok]
    dropped = [round(float(strain_array[i]), 6) for i, ok in enumerate(converged) if not ok]
    if dropped:
        logger.warning(
            "%s: dropping %d unconverged strain point(s) %s from the fit; fitting %d of %d points",
            label,
            len(dropped),
            dropped,
            len(keep),
            len(strain_array),
        )
    kept_strains = [round(float(strain_array[i]), 6) for i in keep]
    if len(keep) < min_points:
        raise ConvergenceError(
            f"{label}: only {len(keep)} of {len(strain_array)} strain points converged; "
            f"at least {min_points} are needed for a reliable fit. Unconverged strains: {dropped}"
        )
    if not (any(e < 0 for e in kept_strains) and any(e > 0 for e in kept_strains)):
        raise ConvergenceError(
            f"{label}: the converged strain points {kept_strains} do not span both compression "
            "and tension; the curvature cannot be fitted reliably"
        )
    return keep


def get_curvature(
    strain_array, calc_bulk=False, min_points=MIN_FIT_POINTS, deformation_folder=None, *, legacy_oszicar=False
):
    """Collect the energies of every strain folder in a deformation folder and fit them.

    Unconverged strain points are left out of the fit (see ``select_fit_points``).
    Writes to ``deformation_folder``: energies.txt (every point, with a ``converged`` column and
    NaN energy/volume for unconverged points), plus B.txt for the bulk fit or
    curvature_constant and Fitting_func otherwise. The reference POSCAR is read from the
    parent (project) folder.

    :param calc_bulk: if True, fit an EOS and return the bulk modulus;
                      otherwise fit a parabola and return the curvature constant.
    :param min_points: minimum number of converged strain points required
    :param deformation_folder: the ``deformation-*`` folder (default: current directory)
    :param legacy_oszicar: judge convergence from OSZICAR instead of vasprun.xml (see ``check_convergence``)
    """
    folder = resolve_workspace(deformation_folder)
    label = folder.name
    energies, volumes, converged = [], [], []
    for i in range(1, len(strain_array) + 1):
        strain_folder = folder / strain_dir(i)
        if not strain_folder.exists():
            raise MissingInputError(f"{strain_folder} does not exist")
        ok, _ = check_convergence(strain_folder / "OSZICAR", legacy_oszicar=legacy_oszicar)
        energy, volume = get_energy_volume(strain_folder) if ok else (np.nan, np.nan)
        energies.append(energy)
        volumes.append(volume)
        converged.append(ok)

    if calc_bulk:
        table = {"energies": energies, "volumes": volumes, "strains": strain_array, "converged": converged}
    else:
        table = {"energies": energies, "strains": strain_array, "converged": converged}
    pd.DataFrame(table).to_csv(folder / "energies.txt", sep="\t", index=False)

    keep = select_fit_points(strain_array, converged, label, min_points)
    fit_energies = [energies[i] for i in keep]

    if calc_bulk:
        bulk_mod = fit_parabola_bulk(fit_energies, [volumes[i] for i in keep], calc_bulk_parabola=True)
        write_text(folder / "B.txt", f"B = {bulk_mod} GPa\n")
        return bulk_mod

    V0 = Structure.from_file(folder.parent / "POSCAR").volume
    curvature_const, p = fit_parabola(fit_energies, [strain_array[i] for i in keep], V0)
    write_text(folder / "curvature_constant", f"curvature_constant = {curvature_const} GPa\n")
    write_text(folder / "Fitting_func", f"{p[0]},{p[1]},{p[2]}")
    return curvature_const


def plot_parabola(strains, energies, title, per_atom=True, folder=None):
    """Plot energy vs strain with its quadratic fit and save as ``folder/<title>.png``."""
    fit = np.polyfit(strains, energies, 2)
    st_vals = np.linspace(min(strains), max(strains), 101)
    fig = Figure()
    ax = fig.subplots()
    ax.plot(strains, energies, "ro", label="energies")
    ax.plot(st_vals, np.polyval(fit, st_vals), "r--", label="fit")
    ax.set_title(title)
    ax.set_xlabel("strain")
    ax.set_ylabel("Energy (eV/atom)" if per_atom else "Energy (eV)")
    ax.legend()
    fig.savefig(resolve_workspace(folder) / f"{title}.png")
