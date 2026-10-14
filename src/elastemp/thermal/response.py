"""Quasi-harmonic thermal properties and temperature-dependent elastic constants."""

import logging

import numpy as np
import pandas as pd
import yaml
from pymatgen.core import Structure

from elastemp.base.strain_analysis import fit_parabola
from elastemp.base.utils import (
    RESULTS_DIR,
    TSTEP,
    deformation_dir,
    resolve_workspace,
    results_dir,
    run_phonopy,
    write_text,
)
from elastemp.elastic.constants import compute_elastic_constants
from elastemp.elastic.response import compute_response_functions
from elastemp.exceptions import MissingInputError

logger = logging.getLogger(__name__)

# phonopy reports the free energy in kJ/mol per unit cell. The original workflow converts
# with a factor of 0.01; the exact value is 1/96.485 = 0.010364 eV per kJ/mol.
# Kept at 0.01 so results match the published version -- see README "Known issues".
KJMOL_TO_EV_APPROX = 0.01

VOLUME_TEMP_FILE = f"{RESULTS_DIR}/volume_temp.txt"  # relative to the project folder


class phonon_thermal:
    """Run phonopy thermal-property calculations in one Phonon_calculation folder."""

    def __init__(self, dim, tmax, phonon_folder=None):
        """
        :param dim: supercell dimension
        :param tmax: maximum temperature (K)
        :param phonon_folder: the Phonon_calculation folder to work in (default: current directory)
        """
        self.dim = dim
        self.tmax = tmax
        self.folder = resolve_workspace(phonon_folder)

    def run_phonopy_thermal(self):
        """Write mesh.conf and run ``phonopy -c POSCAR-unitcell -t mesh.conf`` in the phonon folder."""
        self.make_mesh_conf()
        logger.info("Running phonopy thermal properties in %s", self.folder)
        run_phonopy("-c", "POSCAR-unitcell", "-t", "mesh.conf", cwd=self.folder)

    def make_mesh_conf(self):
        """Write the phonopy settings file mesh.conf."""
        write_text(
            self.folder / "mesh.conf",
            "DIM={} {} {}\nMP= 48 48 48\nTMAX={}\nTSTEP={}\nFORCE_CONSTANTS=READ".format(*self.dim, self.tmax, TSTEP),
        )

    def extract_thermal_yaml(self):
        """Convert thermal_properties.yaml into energy_temp.dat (T, F + E0 per atom, Cv per atom)."""
        yaml_path = self.folder / "thermal_properties.yaml"
        try:
            with open(yaml_path) as f:
                data = yaml.load(f, Loader=getattr(yaml, "CSafeLoader", yaml.SafeLoader))
        except OSError as exc:
            raise MissingInputError(f"Unable to read {yaml_path}") from exc
        with open(self.folder / "ZPE.txt") as f:
            zpe = float(f.readline())
        natoms = len(Structure.from_file(self.folder / "POSCAR-unitcell").sites)

        rows = data["thermal_properties"]
        df = pd.DataFrame(
            {
                "temperature": [r["temperature"] for r in rows],
                "helm_energy": [r["free_energy"] * KJMOL_TO_EV_APPROX / natoms + zpe for r in rows],
                "Cv": [r["heat_capacity"] / natoms for r in rows],
            }
        )
        df.to_csv(self.folder / "energy_temp.dat", index=False, sep="\t")


def _average_bulk_modulus(row):
    """(<C11> + 2<C12>)/3 from a [T, c11, c12, ..., c66] row, averaging over the three axes."""
    c11_avg = (row[1] + row[8] + row[15]) / 3.0
    c12_avg = (row[2] + row[3] + row[9]) / 3.0
    return (c11_avg + 2 * c12_avg) / 3.0


class constants_thermal:
    """Extract elastic constants at different temperatures."""

    def __init__(self, num_constants, strain_values, symmetry, tmax, V0, workspace=None):
        """
        :param num_constants: number of independent elastic constants
        :param strain_values: list of strain values
        :param symmetry: symmetry label of the structure
        :param tmax: maximum temperature (K)
        :param V0: volume per atom of the reference structure
        :param workspace: project folder (default: current directory)
        """
        self.num_constants = num_constants
        self.strain_values = strain_values
        self.symmetry = symmetry.strip()
        self.tmax = tmax
        self.V0 = V0
        self.workspace = resolve_workspace(workspace)

    def extract_thermal_constants(self):
        """Fit energy-strain curves at tmax to get isothermal C_ij(tmax), then interpolate and convert."""
        logger.info("Extracting elastic constants vs temperature")
        obj_elastic = compute_elastic_constants(
            self.num_constants, self.strain_values, self.symmetry, workspace=self.workspace
        )

        with open(self.workspace / VOLUME_TEMP_FILE) as f:
            V_T = float(f.readlines()[-1].split()[1])
        natoms = len(Structure.from_file(self.workspace / "POSCAR").sites)
        V_Tpa = V_T / natoms
        c0_file = self.workspace / RESULTS_DIR / "elastic_constants_rawdata_0K.txt"
        c0 = pd.read_csv(c0_file, sep="\t", header=None)[1].values

        curvatures_t = []
        for j in range(1, self.num_constants + 1):
            df = pd.read_csv(self.workspace / deformation_dir(j) / "energy_strain_temperature.csv", sep="\t")
            energies_t = list(df.iloc[-1].values[1 : len(self.strain_values) + 1])
            curvature_const, _ = fit_parabola(energies_t, self.strain_values, V_Tpa)
            curvatures_t.append(curvature_const)

        ct = obj_elastic.get_elastic_constants(curvatures_t, write_elastic=False)
        self.convert_isothermal_to_adiabatic(list(np.asarray(ct).ravel()), c0)

    def convert_isothermal_to_adiabatic(self, ct_array, c0):
        """Linearly interpolate C_ij(T) between 0 K and tmax and convert isothermal -> adiabatic.

        Writes elastic_constants_temp_{isothermal,adiabatic}.csv, moduli_temp_{isothermal,adiabatic}.csv
        and lambda_vals.txt to results_dir.

        :param ct_array: 36 isothermal elastic constants at tmax (row-major)
        :param c0: 36 elastic constants at 0 K (row-major)
        """
        a = [1, 1, 1, 0, 0, 0]
        if self.symmetry == "monoclinic":
            a[4] = 1 / 2
        elif self.symmetry == "triclinic":
            a[3] = a[4] = a[5] = 1 / 2
        elif self.symmetry not in (
            "cubic",
            "trigonal1",
            "trigonal2",
            "hexagonal",
            "orthorhombic",
            "tetragonal1",
            "tetragonal2",
        ):
            raise ValueError(f"Unknown symmetry '{self.symmetry}'")

        # Slope dC_ij/dT of the linear interpolation; zero where C_ij(0 K) is zero by symmetry.
        l_arr = [0 if c0[i] == 0 else (ct_array[i] - c0[i]) / self.tmax for i in range(36)]
        write_text(results_dir(self.workspace) / "lambda_vals.txt", "".join(f"{v}\n" for v in l_arr))

        # phi(T) is the last column of volume_temp.txt; row k corresponds to T = k * TSTEP.
        with open(self.workspace / VOLUME_TEMP_FILE) as f:
            phi_values = [float(line.split()[-1]) for line in f.readlines()[1:]]

        num_T_points = 1 + int(self.tmax / TSTEP)
        ct_rows, cs_rows, mod_t_rows, mod_s_rows = [], [], [], []

        for k in range(num_T_points):
            temp = float(k * TSTEP)
            ct_row = [temp] + [c0[j] + l_arr[j] * temp for j in range(36)]
            c_temp_t = np.reshape(ct_row[1:], (6, 6))

            m = c_temp_t @ np.array(a, dtype=float)
            c_temp_s = c_temp_t + phi_values[k] * np.outer(m, m)
            cs_row = [temp] + list(c_temp_s.ravel())

            resp_t = compute_response_functions(c_temp_t, self.symmetry)
            s_t = resp_t.get_compliance_matrix(write_compliance=False)
            Eh, Gh, hc, nu, _ = resp_t.get_moduli(s_t, write_moduli=False)
            mod_t_rows.append([temp, _average_bulk_modulus(ct_row), Eh, Gh, hc, nu])

            resp_s = compute_response_functions(c_temp_s, self.symmetry)
            s_s = resp_s.get_compliance_matrix(write_compliance=False)
            Eh, Gh, hc, nu, _ = resp_s.get_moduli(s_s, write_moduli=False)
            mod_s_rows.append([temp, _average_bulk_modulus(cs_row), Eh, Gh, hc, nu])

            ct_rows.append(ct_row)
            cs_rows.append(cs_row)

        elastic_cols = ["temp"] + [f"c{i + 1}{j + 1}" for i in range(6) for j in range(6)]
        moduli_cols = ["T", "K", "Eh", "Gh", "hc", "nu"]

        pd.DataFrame(ct_rows, columns=elastic_cols).to_csv(
            self.workspace / RESULTS_DIR / "elastic_constants_temp_isothermal.csv", index=False
        )
        pd.DataFrame(cs_rows, columns=elastic_cols).to_csv(
            self.workspace / RESULTS_DIR / "elastic_constants_temp_adiabatic.csv", index=False
        )
        pd.DataFrame(mod_t_rows, columns=moduli_cols).to_csv(
            self.workspace / RESULTS_DIR / "moduli_temp_isothermal.csv", sep="\t", index=False
        )
        pd.DataFrame(mod_s_rows, columns=moduli_cols).to_csv(
            self.workspace / RESULTS_DIR / "moduli_temp_adiabatic.csv", sep="\t", index=False
        )
        logger.info("Wrote temperature-dependent elastic constants and moduli to %s", results_dir(self.workspace))
