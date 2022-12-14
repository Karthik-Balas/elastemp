from pymatgen.core import Structure
import os,pandas as pd, sys
from elastemp.base.symmetry import get_symmetry,get_num_constants
from elastemp.base.makefiles import make_incar_static, make_kpt_static
from elastemp.elastic.constants import compute_elastic_constants
from elastemp.elastic.response import compute_response_functions
from elastemp.elastic.stability import stability_check
from elastemp.elastic.deformations import cubic_def,hexagonal_def, trigonal1_def, trigonal2_def
from elastemp.elastic.deformations import tetragonal1_def, tetragonal2_def, orthorhombic_def, monoclinic_def, triclinic_def



def make_strains():
    """ Returns a default strain vector if none is present.
    :param None
    :returns Strain vector, strains.txt
    :rtype: list, .txt file
    """
    strains = [-0.05,-0.03,-0.01,0,0.01,0.03,0.05]
    f = open('strains.txt','w')
    for i in range(len(strains)):
        f.write(str(strains[i])+'\n')
    f.close()
    return strains


def make_elastic_wrapper(symprec=0.05,angle_tolerance=10):
    """ Wrapper function to preprocess data to call make_elastic_deformation function.
        Creates a default INCAR and KPOINTS file if not provided by user.
    
    :param symprec         : Tolerance parameter in lattice parameters for symmetry detection
           angle_tolerance : Tolerance parameter in lattice angles for symmetry detection.
    :type  symprec         : float
           angle_tolerance : float
    :returns: INCAR (if absent) : VASP file needed to run vasp calculation.
             KPOINTS (if absent): VASP file needed to run vasp calculation.
    :rtype :  plaintext file
    """
    try:
        struct = Structure.from_file('POSCAR')
        print('Input POSCAR file read successfully')
        try:
            strain_values = list(pd.read_csv('strains.txt',header=None)[0].values)
            if len(strain_values)<6:
                print('Warning: Atleast 6 strain values are needed for a good fit!. ')
                
            if min(strain_values)>0:
                print('Warning: Atleast some strain values needs to be less than zero. Ideally half the values. ')
                
            
        except:
            print('Unable to read strain values from Input folder. File is either  missing or corrputed')
            print('Default values of -0.05 to 0.05 are added')
            strain_values=make_strains()
        if not os.path.exists('INCAR'):
            make_incar_static()
        if not os.path.exists('KPOINTS'):
            make_kpt_static()
    except:
        print('Unable to read POSCAR file. Program will stop')
        sys.exit()
    make_elastic_deformation(struct,strain_values,symprec,angle_tolerance)

def write_symmetry(struct,symprec,angle_tolerance):
    """ Returns a symmetry file to store the symmetry (one of 9 classes) of the structure.
    :param struct :          Structure
           symprec :         Tolerance parameter in lattice parameters for symmetry detection
           angle_tolerance : Tolerance parameter in lattice angles for symmetry detection.
    :type  struct :          Pymatgen Structure object
           symprec:          float
           angle_tolerance : float
    :returns symmetry
    :rtype : plaintext file
    """
    symmetry = get_symmetry(struct,symprec,angle_tolerance)
    f = open('symmetry','w')
    f.write(str(symmetry))
    f.close()

def make_elastic_deformation(struct, strain_values,symprec,angle_tolerance):
    """ Calls functions for each symmetry to create the deformations. 
    :param struct          : Structure
           strain_values   : strain values needed to fit the energy-strain curve. 
           symprec         : Tolerance parameter in lattice parameters for symmetry detection.
           angle_tolerance : Tolerance parameter in lattice angles for symmetry detection. 
    :type  struct          : Pymatgen structure object
           strain_values   : list
           symprec         : float
           angle_tolerance : float
    :returns Calls function to save symmetry and function to create symmetry dependent strain
             deformations with the strain values
    """
    if not os.path.exists('symmetry'):
        write_symmetry(struct,symprec,angle_tolerance)
    symmetry = open('symmetry').readlines()[0]
    print('Spacegroup detected : {}'.format(symmetry))
    
    if symmetry == 'cubic':
        print('Making 3 deformations')
        cubic_def(struct,symmetry,strain_values)
        print('Finished making 3 deformations')
    elif symmetry == 'hexagonal':
        print('Making 5 deformations')
        hexagonal_def(struct,symmetry, strain_values)
    elif symmetry == 'trigonal1':
        print('Making 6 deformations')
        trigonal1_def(struct,symmetry, strain_values)
    elif symmetry == 'trigonal2':
        print('Making 7 deformations')
        trigonal2_def(struct,symmetry, strain_values)
    elif symmetry == 'tetragonal1':
        print('Making 6 deformations')
        tetragonal1_def(struct,symmetry, strain_values)
    elif symmetry == 'tetragonal2':
        print('Making 7 deformations')
        tetragonal2_def(struct,symmetry, strain_values)
    elif symmetry == 'orthorhombic':
        print('Making 9 deformations')
        orthorhombic_def(struct,symmetry, strain_values)
    elif symmetry == 'monoclinic':
        print('Making 13 deformations')
        monoclinic_def(struct,symmetry, strain_values)
    elif symmetry == 'triclinic':
        print('Making 21 deformations')
        triclinic_def(struct,symmetry, strain_values)

def get_elastic_wrapper(make_plots):
    """ Wrapper function to preprocess data to call get_elastic function. 
    :param make_plots  : Boolean value to decide if plotting needs to be done. 
    :type  make_plots  : Bool (True/False)
    :returns Calls get_elastic function 
    """
    strain_values = list(pd.read_csv('strains.txt',header=None)[0].values)
    get_elastic(strain_values, make_plots)

def get_elastic(strain_values,make_plots):
    """ Calls function to get the elastic constants from energy-strain curves.
    :param strain_values : list of strain values
           make_plots    : Boolean value to decide if plotting needs to be done.
    :type  strain_values : list
           make_plots    : Bool (True/False)
    :returns Calls function to compute elastic constants, elastic moduli, check mechanical stability.
    """
    symmetry = open('symmetry').readlines()[0]
    num_constants = get_num_constants(symmetry)
    e = compute_elastic_constants(num_constants,strain_values,symmetry)
    curvatures = e.get_curvatures()
    print('Getting elastic constants')
    c = e.get_elastic_constants(curvatures,write_elastic=True)
    e.write_zero_elastic_constants(c)
    if make_plots==True:
        e.get_plots()
    print('Calculating response functions')
    obj = compute_response_functions(c,symmetry)
    s = obj.get_compliance_matrix(write_compliance=True)
    Eh,Gh, hc, nu, pugh = obj.get_moduli(s,write_moduli=True)
    obj = stability_check(c,symmetry)

