"""Polycrystalline moduli (Voigt-Reuss-Hill), Poisson/Pugh ratios and Chen hardness."""

import logging

import numpy as np

from elastemp.base.utils import format_matrix, resolve_workspace, results_dir
from elastemp.exceptions import ElastempError

logger = logging.getLogger(__name__)

POISSON_BRITTLE_BELOW = 0.28
PUGH_BRITTLE_ABOVE = 0.5


class compute_response_functions:
    """Compute response functions (E, G, K, H) from a stiffness matrix."""

    def __init__(self, c, symmetry, workspace=None):
        """
        :param c: 6x6 stiffness matrix (GPa)
        :param symmetry: symmetry label of the structure
        """
        self.c = c
        self.symmetry = symmetry
        self.workspace = resolve_workspace(workspace)

    def get_moduli(self, s, write_moduli):
        """Voigt/Reuss/Hill moduli, Poisson and Pugh ratios, and Chen hardness.

        :param s: 6x6 compliance matrix (1/GPa)
        :param write_moduli: if True, write results_dir/moduli.txt
        :returns: (Eh, Gh, hc, nu, pugh)
        """
        c = self.c
        Kv = np.round(((c[0][0] + c[1][1] + c[2][2]) + 2 * (c[0][1] + c[1][2] + c[2][0])) / 9.0, 4)
        Gv = np.round(
            ((c[0][0] + c[1][1] + c[2][2]) - (c[0][1] + c[1][2] + c[2][0]) + 3 * (c[3][3] + c[4][4] + c[5][5])) / 15.0,
            4,
        )
        l1 = (s[0][0] + s[1][1] + s[2][2]) + 2 * (s[0][1] + s[1][2] + s[2][0])
        Kr = np.round(1.0 / l1, 4)
        l2 = 4 * (s[0][0] + s[1][1] + s[2][2]) - 4 * (s[0][1] + s[1][2] + s[2][0]) + 3 * (s[3][3] + s[4][4] + s[5][5])
        Gr = np.round(15.0 / l2, 4)
        Kh = np.round((Kv + Kr) / 2, 4)
        Gh = np.round((Gv + Gr) / 2, 4)

        Ev = np.round((9 * Kv * Gv) / (3 * Kv + Gv), 4)
        Er = np.round((9 * Kr * Gr) / (3 * Kr + Gr), 4)
        Eh = np.round((Ev + Er) / 2, 4)

        nu = np.round((3 * Kh - 2 * Gh) / (6 * Kh + Gh), 4)
        ductility_variable = "Brittle" if nu < POISSON_BRITTLE_BELOW else "Ductile"

        pugh = np.round(Gh / Kh, 4)
        pugh_variable = "Brittle" if pugh > PUGH_BRITTLE_ABOVE else "Ductile"

        # Chen's model: H = 2 (k^2 G)^0.585 - 3, with k = G/K (Pugh's ratio)
        base = Gh * pugh**2
        if base > 0:
            hc = np.round(2 * base**0.585 - 3, 4)
        else:
            logger.warning("Cannot compute Chen hardness (G*k^2 = %s <= 0); setting it to 0", base)
            hc = 0

        if write_moduli:
            with open(results_dir(self.workspace) / "moduli.txt", "w") as f:
                f.write("**" * 30 + "\n")
                f.write(f"Voigt bulk modulus ={Kv} GPa\n ")
                f.write(f"Voigt shear modulus ={Gv} GPa\n ")
                f.write(f"Voigt elastic modulus ={Ev} GPa\n\n ")
                f.write(f"Reuss bulk modulus ={Kr} GPa\n ")
                f.write(f"Reuss shear modulus ={Gr} GPa\n ")
                f.write(f"Reuss elastic modulus ={Er} GPa\n\n ")
                f.write(f"Hill bulk modulus = {Kh} GPa\n")
                f.write(f"Hill shear modulus = {Gh} GPa\n")
                f.write(f"Hill elastic modulus = {Eh} GPa\n\n")
                f.write("----" * 30 + "\n")
                f.write(
                    f"poisson ratio = {nu}.\t. Classified as {ductility_variable} (poisson ratio less than 0.28 "
                    "classified as brittle, else ductile)\n"
                )
                f.write(
                    f"pugh's ratio = {pugh}.\t. Classified as {pugh_variable} (pugh's ratio greater than 0.5 "
                    "classified as brittle, else ductile)\n"
                )
                f.write(f"Hardness from Chen's model = {hc} GPa\n")
                f.write("**" * 30 + "\n")
            logger.info("Hill moduli: K=%s G=%s E=%s GPa, nu=%s, H=%s GPa", Kh, Gh, Eh, nu, hc)
        return Eh, Gh, hc, nu, pugh

    def get_compliance_matrix(self, write_compliance):
        """Invert the stiffness matrix.

        :param write_compliance: if True, write results_dir/compliance_matrix.txt
        :returns: 6x6 compliance matrix (1/GPa)
        """
        if np.linalg.det(self.c) == 0:
            raise ElastempError("Stiffness matrix is singular; cannot compute compliance matrix")
        s = np.linalg.inv(self.c)
        if write_compliance:
            with open(results_dir(self.workspace) / "compliance_matrix.txt", "w") as f:
                f.write("**" * 30 + "\n")
                f.write("COMPLIANCE MATRIX\n")
                f.write("**" * 30 + "\n")
                f.write(format_matrix(s, 4))
        return s
