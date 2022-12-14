
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer as sga


def get_symmetry(struct,symprec,angle_tolerance):
    """ Detects symmetry and writes to file.
    :param struct          : Structure
           symprec         : Tolerance in lattice parameters for symmetry detection.
           angle_tolerance : Tolerance in lattice angles for symmetry detection.
    :type  struct          : Pymatgen structure object
           symprec         : float
           angle_tolerance : float
    :returns symmetry file
    :rtype plaintext file
    """
    num = sga(struct,symprec=symprec,angle_tolerance=angle_tolerance).get_space_group_number()
    if num in range(195,231):
        return 'cubic'
    elif num in range(168,195):
        return 'hexagonal'
    elif num in range(149,168):
        return 'trigonal1'
    elif num in range(143,149):
        return 'trigonal2'
    elif num in range(89,143):
        return 'tetragonal1'
    elif num in range(75,89):
        return 'tetragonal2'
    elif num in range(16,75):
        return 'orthorhombic'
    elif num in range(3,16):
        return 'monoclinic'
    else:
        return 'triclinic'

def get_num_constants(symmetry):
    """ Detects number of independent elastic constants based on symmetry and writes to file.
    :param symmetry          : Symmtery of structure
    :type  plaintext
    :returns num_constants
    :rtype int
    """
    if symmetry=='cubic':
        num_constants = 3
    elif symmetry=='hexagonal':
        num_constants = 5
    elif symmetry=='trigonal1' or symmetry=='tetragonal1':
        num_constants= 6
    elif symmetry=='trigonal2'or symmetry== 'tetragonal2':
        num_constants=7
    elif symmetry=='orthorhombic':
        num_constants=9
    elif symmetry=='monoclinic':
        num_constants=13
    elif symmetry=='triclinic':
        num_constants=21
    else:
        num_constants=0
    return num_constants

