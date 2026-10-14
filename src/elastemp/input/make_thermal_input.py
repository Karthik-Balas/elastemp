"""Step 9: temperature-dependent elastic constants and moduli."""

import logging

import pandas as pd
from matplotlib.figure import Figure
from pymatgen.core import Structure

from elastemp.base.strain_analysis import plot_parabola
from elastemp.base.symmetry import get_num_constants, get_symmetry
from elastemp.base.utils import (
    deformation_dir,
    read_dim,
    read_strains,
    read_symmetry,
    read_tmax,
    resolve_workspace,
    results_dir,
    strain_dir,
)
from elastemp.dynamic.deformations import PHONON_DIR
from elastemp.exceptions import MissingInputError
from elastemp.thermal.response import constants_thermal, phonon_thermal

logger = logging.getLogger(__name__)


def _symmetry_or_detect(workspace, struct):
    try:
        return read_symmetry(workspace)
    except MissingInputError:
        logger.warning("symmetry file missing; detecting symmetry with default tolerances")
        return get_symmetry(struct, symprec=0.05, angle_tolerance=10)


def get_thermal_wrapper(workspace=None):
    """Run phonopy QHA in every deformation folder and tabulate energy vs strain vs T."""
    ws = resolve_workspace(workspace)
    get_thermal(read_tmax(ws), read_dim(ws), read_strains(ws), ws)


def get_thermal(tmax, dim, strain_values, workspace=None):
    """Write deformation-*/energy_strain_temperature.csv and a tmax energy-strain plot per deformation."""
    ws = resolve_workspace(workspace)
    symmetry = _symmetry_or_detect(ws, Structure.from_file(ws / "POSCAR"))
    num_constants = get_num_constants(symmetry)

    for i in range(1, num_constants + 1):
        deformation = ws / deformation_dir(i)
        logger.info("Thermal properties for %s", deformation)
        merged = None
        for j in range(1, len(strain_values) + 1):
            folder = deformation / strain_dir(j) / PHONON_DIR
            thermal = phonon_thermal(dim, tmax, phonon_folder=folder)
            thermal.run_phonopy_thermal()
            thermal.extract_thermal_yaml()
            df_en = pd.read_csv(
                folder / "energy_temp.dat",
                usecols=[0, 1],
                names=["temperature", f"strain-{j}_energy"],
                skiprows=1,
                sep="\t",
            )
            merged = df_en if merged is None else merged.merge(df_en, on="temperature")

        plot_parabola(
            strain_values,
            merged.iloc[-1].values[1:],
            f"Deformation-{i}-strain-curve_tmax",
            per_atom=True,
            folder=deformation,
        )
        merged.to_csv(deformation / "energy_strain_temperature.csv", sep="\t", index=False)


def get_thermal_constants(workspace=None):
    """Compute isothermal/adiabatic elastic constants and moduli for 0..tmax."""
    ws = resolve_workspace(workspace)
    struct = Structure.from_file(ws / "POSCAR")
    symmetry = _symmetry_or_detect(ws, struct)
    V0 = struct.volume / len(struct.sites)
    constants_thermal(
        get_num_constants(symmetry), read_strains(ws), symmetry, read_tmax(ws), V0, workspace=ws
    ).extract_thermal_constants()


def _plot_vs_temperature(T, isothermal, adiabatic, ylabel, filename):
    fig = Figure()
    ax = fig.subplots()
    ax.plot(T, isothermal, "b--", label="Isothermal")
    ax.plot(T, adiabatic, "r", label="Adiabatic")
    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel(ylabel)
    ax.legend()
    fig.savefig(filename)


def get_elastic_temp_plots(workspace=None):
    """Plot isothermal and adiabatic moduli and C_ij versus temperature into results_dir/moduli_plots."""
    out = results_dir(workspace)
    df_mod_t = pd.read_csv(out / "moduli_temp_isothermal.csv", sep="\t")
    df_mod_s = pd.read_csv(out / "moduli_temp_adiabatic.csv", sep="\t")
    df_ct = pd.read_csv(out / "elastic_constants_temp_isothermal.csv")
    df_cs = pd.read_csv(out / "elastic_constants_temp_adiabatic.csv")

    plot_dir = out / "moduli_plots"
    plot_dir.mkdir(exist_ok=True)

    moduli = {
        "K": "Bulk_modulus (GPa)",
        "Eh": "Elastic_modulus (GPa)",
        "Gh": "Shear_modulus (GPa)",
        "hc": "Hardness (GPa)",
        "nu": "poisson_ratio",
    }
    logger.info("Plotting moduli vs temperature")
    for col, label in moduli.items():
        _plot_vs_temperature(
            df_mod_t["T"].values,
            df_mod_t[col].values,
            df_mod_s[col].values,
            label,
            plot_dir / f"{label}_temperature_isothermal_adiabatic.png",
        )

    logger.info("Plotting elastic constants vs temperature")
    for i in range(1, 7):
        for j in range(i, 7):
            name = f"c{i}{j}"
            _plot_vs_temperature(
                df_ct["temp"].values,
                df_ct[name].values,
                df_cs[name].values,
                name,
                plot_dir / f"{name}_temperature_isothermal_adiabatic.png",
            )
    logger.info("Plots written to %s", plot_dir)
