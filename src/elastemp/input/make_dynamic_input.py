"""Steps 3-8: phonon set-up/extraction for the bulk and per-deformation calculations.

Every step takes ``workspace``, the project folder (default: current directory).
"""

import logging
import shutil

import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from pymatgen.core import Structure

from elastemp.base.makefiles import make_incar_dynamic, make_kpt_dynamic
from elastemp.base.strain_analysis import check_convergence, fit_parabola_bulk, get_energy_volume, plot_parabola
from elastemp.base.symmetry import get_num_constants
from elastemp.base.utils import (
    DIM_FILE,
    EV_PER_A3_TO_GPA,
    TMAX_FILE,
    deformation_dir,
    read_dim,
    read_strains,
    read_symmetry,
    read_tmax,
    require_file,
    resolve_workspace,
    results_dir,
    strain_dir,
    write_text,
)
from elastemp.dynamic.deformations import HIGHTEMP_DIR, PHONON_DIR, dynamic_deformations, extract_force_constants
from elastemp.exceptions import ConvergenceError, MissingInputError
from elastemp.thermal.response import VOLUME_TEMP_FILE, phonon_thermal

logger = logging.getLogger(__name__)

DEFAULT_DIM = [2, 2, 2]
DEFAULT_TMAX = 1000


def make_dim(workspace=None):
    """Write the default supercell dimension file dim.txt (2 2 2)."""
    write_text(resolve_workspace(workspace) / DIM_FILE, "{} {} {}".format(*DEFAULT_DIM))


def make_tmax(workspace=None):
    """Write the default tmax.txt (1000 K)."""
    write_text(resolve_workspace(workspace) / TMAX_FILE, f"{DEFAULT_TMAX}\n")


def load_or_create_dim(workspace=None):
    """Read dim.txt, or write and return the default 2 2 2 if it is missing or unreadable."""
    try:
        dim = read_dim(workspace)
        logger.info("Supercell dimension: %s", dim)
        return dim
    except (MissingInputError, ValueError):
        logger.warning("Unable to read %s; using default supercell dimension %s", DIM_FILE, DEFAULT_DIM)
        make_dim(workspace)
        return list(DEFAULT_DIM)


def load_or_create_tmax(workspace=None):
    """Read tmax.txt, or write and return the default 1000 K if it is missing or unreadable."""
    try:
        return read_tmax(workspace)
    except (MissingInputError, ValueError):
        logger.warning("Unable to read %s; using default tmax %d K", TMAX_FILE, DEFAULT_TMAX)
        make_tmax(workspace)
        return DEFAULT_TMAX


def _phonon_dir(workspace, deformation, i):
    return workspace / deformation_dir(deformation) / strain_dir(i) / PHONON_DIR


def make_dynamic_wrapper(calc_bulk_dynamic, workspace=None):
    """Steps 3 and 7: create phonon supercells.

    :param calc_bulk_dynamic: True for deformation-bulk (step 3), False for every deformation (step 7)
    """
    ws = resolve_workspace(workspace)
    require_file(ws / "POTCAR")
    struct = Structure.from_file(ws / "POSCAR")
    strain_values = read_strains(ws)
    if not (ws / "KPOINTS_dynamic").exists():
        make_kpt_dynamic(ws)
    if not (ws / "INCAR_dynamic").exists():
        make_incar_dynamic(ws)
    dim = load_or_create_dim(ws)

    if calc_bulk_dynamic:
        make_bulk_dynamic(struct, strain_values, dim, ws)
    else:
        make_dynamic(struct, strain_values, dim, ws)


def make_bulk_dynamic(struct, strain_values, dim, workspace=None):
    """Step 3: create phonon supercells for each strain of deformation-bulk."""
    ws = resolve_workspace(workspace)
    dynamic_deformations(read_symmetry(ws), strain_values, dim, struct, ws).make_bulk_deformation()
    logger.info("Run VASP DFPT in deformation-bulk/strain-*/%s, then run get_bulk_dynamic", PHONON_DIR)


def get_dynamic_wrapper(calc_bulk_dynamic, workspace=None):
    """Steps 4 and 8: extract force constants (and, for the bulk, thermal expansion)."""
    ws = resolve_workspace(workspace)
    strain_values = read_strains(ws)
    if calc_bulk_dynamic:
        get_bulk_dynamic(strain_values, ws)
    else:
        get_dynamic(strain_values, ws)


def _merge_energy_temp(folders, column, label):
    """Merge one column of energy_temp.dat from each folder into a single table keyed on temperature.

    :param column: 1 for Helmholtz energy, 2 for Cv
    """
    merged = None
    for i, folder in enumerate(folders, start=1):
        df = pd.read_csv(
            folder / "energy_temp.dat",
            usecols=[0, column],
            names=["temperature", f"strain-{i}_{label}"],
            skiprows=1,
            sep="\t",
        )
        merged = df if merged is None else merged.merge(df, on="temperature")
    return merged


def get_bulk_dynamic(strain_values, workspace=None):
    """Step 4: QHA on deformation-bulk -> equilibrium volume, thermal expansion and phi(T).

    Writes results_dir/volume_temp.txt and results_dir/Thermal_expansion_temperature.png.
    """
    ws = resolve_workspace(workspace)
    bulk = ws / deformation_dir("bulk")
    n = len(strain_values)

    # Check the static volumes first, before spending time on phonopy.
    vols = pd.read_csv(bulk / "energies.txt", sep="\t")["volumes"]
    if vols.isna().any():
        missing = [strain_values[i] for i in np.flatnonzero(vols.isna().values)]
        raise ConvergenceError(
            f"deformation-bulk: static calculations for strains {missing} did not converge, "
            "so their volumes are unknown. Fix them and rerun get_elastic before get_bulk_dynamic."
        )

    folders = [_phonon_dir(ws, "bulk", j) for j in range(1, n + 1)]
    for folder in folders:
        extract_force_constants(folder)

    tmax = load_or_create_tmax(ws)
    dim = read_dim(ws)
    for folder in folders:
        thermal = phonon_thermal(dim, tmax, phonon_folder=folder)
        thermal.run_phonopy_thermal()
        thermal.extract_thermal_yaml()

    df_en_temp = _merge_energy_temp(folders, 1, "energy")
    df_cv_temp = _merge_energy_temp(folders, 2, "cv")
    df_en_temp.to_csv(bulk / "energy_strain_temperature.csv", sep="\t", index=False)
    df_cv_temp.to_csv(bulk / "Cv_strain_temperature.csv", sep="\t", index=False)
    plot_parabola(
        strain_values, df_en_temp.iloc[-1].values[1:], "Deformation-bulk-curve-tmax", per_atom=True, folder=bulk
    )

    struct = Structure.from_file(ws / "POSCAR")
    natoms = len(struct.sites)
    V0 = struct.volume

    temps, vols_T, lin_exp = [], [], []
    lines = [
        "temperature \t Volume_T\tVolume_0\tVolume_thermal_expansion\tLinear_thermal_expansion"
        "\tCv_calc\tLinear_expansion_factor\tphi\n"
    ]
    for i in range(len(df_en_temp)):
        energies = df_en_temp.iloc[i].values[1:]
        cvs = df_cv_temp.iloc[i].values[1:]
        temperature = df_en_temp.iloc[i].iloc[0]

        V_T = fit_parabola_bulk(energies, vols, calc_bulk_parabola=False)
        cv_fit = np.polyfit(strain_values, cvs, 2)
        vols_T.append(V_T)
        temps.append(temperature)

        vol_thermal_exp = 0
        if i != 0:
            dV = (vols_T[i] - vols_T[i - 1]) / vols_T[i - 1]
            vol_thermal_exp = dV / (temps[i] - temps[i - 1])
        lin_thermal_exp = vol_thermal_exp / 3
        lin_exp.append(lin_thermal_exp * 1e6)

        eq_strain = lin_thermal_exp * temperature
        cv_calc = np.polyval(cv_fit, eq_strain)
        lin_expansion_factor = 1 + eq_strain
        phi = 0
        if i != 0 and cv_calc != 0:
            phi = (1e5 / EV_PER_A3_TO_GPA) * (temperature * (V_T / natoms) / cv_calc) * lin_thermal_exp**2
        lines.append(
            f"{temperature}\t{V_T}\t{V0}\t{vol_thermal_exp}\t{lin_thermal_exp}"
            f"\t{cv_calc}\t{lin_expansion_factor}\t{phi}\n"
        )

    out = results_dir(ws)
    write_text(ws / VOLUME_TEMP_FILE, "".join(lines))

    fig = Figure()
    ax = fig.subplots()
    ax.plot(temps, lin_exp, "r")
    ax.set_title("Linear_thermal_expansion vs temperature")
    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel(" Linear_thermal_expansion * 1e6 (/K)")
    fig.savefig(out / "Thermal_expansion_temperature.png")
    logger.info(
        "Wrote %s; linear expansion factor at %s K = %s", ws / VOLUME_TEMP_FILE, temps[-1], lin_expansion_factor
    )


def _read_lin_expansion_factor(workspace):
    with open(workspace / VOLUME_TEMP_FILE) as f:
        return float(f.readlines()[-1].split()[6])


def make_T_Zpe(workspace=None):
    """Step 5: create Hightemp_static folders with POSCARs scaled to the tmax lattice parameter."""
    ws = resolve_workspace(workspace)
    require_file(ws / VOLUME_TEMP_FILE, f"{ws / VOLUME_TEMP_FILE} not found. Run get_bulk_dynamic first.")
    lin_expansion_factor = _read_lin_expansion_factor(ws)
    logger.info("Linear expansion factor at tmax: %s", lin_expansion_factor)
    num_constants = get_num_constants(read_symmetry(ws))
    strain_values = read_strains(ws)

    for i in range(1, num_constants + 1):
        for j in range(1, len(strain_values) + 1):
            folder = ws / deformation_dir(i) / strain_dir(j)
            high = folder / HIGHTEMP_DIR
            high.mkdir(parents=True, exist_ok=True)
            for name in ("INCAR", "KPOINTS", "POTCAR", "POSCAR"):
                shutil.copy(folder / name, high / name)
            poscar = high / "POSCAR"
            data = poscar.read_text().splitlines(keepends=True)
            data[1] = f"{lin_expansion_factor}\n"  # POSCAR line 2 = lattice scaling factor
            poscar.write_text("".join(data))
    logger.info("Created %s folders. Run VASP in them, then run get_temp_static", HIGHTEMP_DIR)


def get_T_Zpe(workspace=None):
    """Step 6: read the tmax-volume static energies and write ZPE.txt (energy/atom) in each folder."""
    ws = resolve_workspace(workspace)
    num_constants = get_num_constants(read_symmetry(ws))
    strain_values = read_strains(ws)
    for i in range(1, num_constants + 1):
        for j in range(1, len(strain_values) + 1):
            folder = ws / deformation_dir(i) / strain_dir(j) / HIGHTEMP_DIR
            oszicar = folder / "OSZICAR"
            require_file(oszicar)
            converged, detail = check_convergence(oszicar)
            if not converged:
                raise ConvergenceError(f"{HIGHTEMP_DIR} run cannot be used: {detail}")
            get_energy_volume(folder)
    logger.info("Wrote ZPE at tmax volume. Next: make_dynamic")


def make_dynamic(struct, strain_values, dim, workspace=None):
    """Step 7: create phonon supercells for every deformation/strain from Hightemp_static."""
    ws = resolve_workspace(workspace)
    dynamic_deformations(read_symmetry(ws), strain_values, dim, struct, ws).make_deformation()
    logger.info("Run VASP DFPT in deformation-*/strain-*/%s, then run get_dynamic", PHONON_DIR)


def get_dynamic(strain_values, workspace=None):
    """Step 8: extract force constants for every deformation/strain."""
    ws = resolve_workspace(workspace)
    num_constants = get_num_constants(read_symmetry(ws))
    for i in range(1, num_constants + 1):
        for j in range(1, len(strain_values) + 1):
            extract_force_constants(_phonon_dir(ws, i, j))
