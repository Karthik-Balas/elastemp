"""Crystal-symmetry detection and the number of independent elastic constants per class."""

from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

# Laue-class label -> number of independent elastic constants
NUM_CONSTANTS = {
    "cubic": 3,
    "hexagonal": 5,
    "trigonal1": 6,
    "tetragonal1": 6,
    "trigonal2": 7,
    "tetragonal2": 7,
    "orthorhombic": 9,
    "monoclinic": 13,
    "triclinic": 21,
}

# (first space-group number, label), in descending order
_SPACEGROUP_RANGES = (
    (195, "cubic"),
    (168, "hexagonal"),
    (149, "trigonal1"),
    (143, "trigonal2"),
    (89, "tetragonal1"),
    (75, "tetragonal2"),
    (16, "orthorhombic"),
    (3, "monoclinic"),
)


def get_symmetry(struct, symprec, angle_tolerance):
    """Return the Laue-class label of ``struct``.

    :param struct: pymatgen Structure
    :param symprec: tolerance in lattice parameters for symmetry detection
    :param angle_tolerance: tolerance in lattice angles for symmetry detection
    :returns: one of the keys of ``NUM_CONSTANTS``
    """
    num = SpacegroupAnalyzer(struct, symprec=symprec, angle_tolerance=angle_tolerance).get_space_group_number()
    for first, label in _SPACEGROUP_RANGES:
        if num >= first:
            return label
    return "triclinic"


def get_num_constants(symmetry):
    """Number of independent elastic constants for a symmetry label (0 if unknown)."""
    return NUM_CONSTANTS.get(symmetry.strip(), 0)
