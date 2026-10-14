"""Zero-temperature elastic constants from energy-strain curvatures."""

import logging

import numpy as np
import pandas as pd

from elastemp.base.strain_analysis import get_curvature, plot_parabola
from elastemp.base.utils import deformation_dir, format_matrix, resolve_workspace, results_dir

logger = logging.getLogger(__name__)


def _symmetrize_upper(c):
    """Copy the upper triangle of a 6x6 matrix into its lower triangle."""
    return np.triu(c) + np.triu(c, 1).T


def stiffness_from_curvatures(symmetry, curvatures):
    """Assemble the 6x6 stiffness matrix (GPa) from the curvature of each deformation.

    Only the upper triangle is set explicitly; it is then mirrored.
    """
    sym = symmetry.lower()
    k = curvatures
    c = np.zeros((6, 6))

    if sym == "cubic":
        c[0][1] = (2 / 3) * k[1] - k[0]
        c[0][0] = k[0] - c[0][1]
        c[3][3] = (2 / 3) * k[2]
        c[4][4] = c[5][5] = c[3][3]
        c[0][2] = c[1][2] = c[0][1]
        c[1][1] = c[2][2] = c[0][0]

    elif sym in ("hexagonal", "trigonal1", "trigonal2"):
        c[0][1] = (k[0] - 4 * k[1]) / 2
        c[0][0] = k[0] - c[0][1]
        c[2][2] = 2 * k[2]
        c[3][3] = k[3]
        c[0][2] = (k[4] - c[0][0] - c[0][1] - (c[2][2] / 2)) / 2
        c[1][1] = c[0][0]
        c[1][2] = c[0][2]
        c[4][4] = c[3][3]
        c[5][5] = (c[0][0] - c[0][1]) / 2
        if sym in ("trigonal1", "trigonal2"):
            c[0][3] = k[5] - (c[0][0] / 4) + (c[0][1] / 4) - (c[3][3] / 2)
            c[1][3] = -c[0][3]
            c[4][5] = c[0][3]
        if sym == "trigonal2":
            c[0][4] = -1 * (k[6] - (c[0][0] / 4) + (c[0][1] / 4) - (c[3][3] / 2))
            c[1][4] = -c[0][4]
            c[3][5] = -c[0][4]

    elif sym in ("tetragonal1", "tetragonal2"):
        c[5][5] = k[1] * 2
        c[2][2] = k[2] * 2
        c[3][3] = k[3]
        c[0][2] = (k[4] - k[0] - 0.5 * c[2][2]) / 2
        c[0][0] = (k[5] - c[0][2] - 0.5 * c[2][2]) * 2
        c[0][1] = k[0] - c[0][0]
        c[1][1] = c[0][0]
        c[1][2] = c[0][2]
        c[4][4] = c[3][3]
        if sym == "tetragonal2":
            c[0][5] = k[6] - 0.5 * (c[0][0] + c[5][5])
            c[1][5] = -c[0][5]

    elif sym in ("orthorhombic", "monoclinic", "triclinic"):
        # Diagonal terms come from single-component strains (curvature = Cii/2) ...
        for i in range(6):
            c[i][i] = k[i] * 2
        # ... and each off-diagonal Cij from a two-component strain (Cii/2 + Cij + Cjj/2).
        pairs = {
            "orthorhombic": [(0, 1), (0, 2), (1, 2)],
            "monoclinic": [(0, 1), (0, 2), (0, 4), (1, 2), (1, 4), (2, 4), (3, 5)],
            "triclinic": [(i, j) for i in range(6) for j in range(i + 1, 6)],
        }[sym]
        for n, (i, j) in enumerate(pairs, start=6):
            c[i][j] = k[n] - 0.5 * (c[i][i] + c[j][j])

    else:
        raise ValueError(f"Unknown symmetry '{symmetry}'")

    return _symmetrize_upper(c)


class compute_elastic_constants:
    """Compute zero-temperature elastic constants from energy-strain data on disk."""

    def __init__(self, num_constants, strain_array, symmetry, workspace=None):
        """
        :param num_constants: number of independent elastic constants
        :param strain_array: list of strain values
        :param symmetry: symmetry label of the structure
        :param workspace: project folder containing deformation-* (default: current directory)
        """
        self.num_constants = num_constants
        self.strain_array = strain_array
        self.symmetry = symmetry.strip()
        self.workspace = resolve_workspace(workspace)

    def get_curvatures(self):
        """Fit every ``deformation-k`` folder and return the list of curvatures (GPa)."""
        curvatures = []
        for i in range(1, self.num_constants + 1):
            folder = self.workspace / deformation_dir(i)
            curvatures.append(get_curvature(self.strain_array, calc_bulk=False, deformation_folder=folder))
        logger.debug("Curvatures: %s", curvatures)
        return curvatures

    def get_bulk_curvatures(self):
        """Fit the EOS of ``deformation-bulk``; writes B.txt and energies.txt there.

        :returns: bulk modulus (GPa)
        """
        folder = self.workspace / deformation_dir("bulk")
        bulk_mod = get_curvature(self.strain_array, calc_bulk=True, deformation_folder=folder)
        logger.info("Bulk modulus from EOS fit: %s GPa", bulk_mod)
        return bulk_mod

    def get_elastic_constants(self, curvatures, write_elastic):
        """Build the stiffness matrix from curvatures.

        :param write_elastic: if True, round to 3 decimals and write results_dir/stiffness_matrix.txt
        :returns: 6x6 stiffness matrix (GPa)
        """
        c = stiffness_from_curvatures(self.symmetry, curvatures)
        if write_elastic:
            c = np.round(c.astype(float), 3)
            with open(results_dir(self.workspace) / "stiffness_matrix.txt", "w") as f:
                f.write("**" * 30 + "\n")
                f.write("STIFFNESS MATRIX\n")
                f.write("**" * 30 + "\n")
                f.write(format_matrix(c, 3))
        return c

    def write_zero_elastic_constants(self, c):
        """Write the 0 K elastic constants, one ``c_ij=<tab>value`` per line, to results_dir."""
        with open(results_dir(self.workspace) / "elastic_constants_rawdata_0K.txt", "w") as f:
            for i in range(6):
                for j in range(6):
                    f.write(f"c_{i + 1}{j + 1}=\t{c[i][j]}\n")

    def get_plots(self):
        """Save an energy-strain plot in every deformation folder."""
        for k in range(1, self.num_constants + 1):
            folder = self.workspace / deformation_dir(k)
            df = pd.read_csv(folder / "energies.txt", sep="\t")
            if "converged" in df:
                df = df[df["converged"]]
            plot_parabola(
                df["strains"].values,
                df["energies"].values,
                f"Deformation-{k}-strain-curve",
                per_atom=False,
                folder=folder,
            )
