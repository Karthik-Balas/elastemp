from pymatgen.core import Structure
import numpy as np,os, pandas as pd,sys, shutil
import matplotlib.pyplot as plt
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer as sga
from elastemp.base.strain_analysis import get_strain_matrix as smv,get_curvature,check_convergence,get_energy_volume
from elastemp.base.symmetry import get_symmetry,get_num_constants
from elastemp.base.makefiles import make_incar_static, make_kpt_static
from elastemp.base.makefiles import make_incar_dynamic, make_kpt_dynamic
from elastemp.elastic.constants import compute_elastic_constants
from elastemp.elastic.response import compute_response_functions
from elastemp.elastic.stability import stability_check
from elastemp.elastic.deformations import cubic_def,hexagonal_def, trigonal1_def, trigonal2_def
from elastemp.elastic.deformations import tetragonal1_def, tetragonal2_def, orthorhombic_def, monoclinic_def, triclinic_def
from elastemp.dynamic.deformations import dynamic_deformations
from elastemp.dynamic.response import get_dynamic_properties
from elastemp.thermal.response import phonon_thermal,constants_thermal
#from elastemp.thermal.response import phonon_thermal


def make_strains():
    strains = [-0.05,-0.03,-0.01,0,0.01,0.03,0.05]
    f = open('strains.txt','w')
    for i in range(len(strains)):
        f.write(str(strains[i])+'\n')
    f.close()
    return strains

def make_dim():
    f = open('dim.txt','w')
    f.write('{} {} {}'.format(2,2,2))
    f.close()
def make_tmax():
    f = open('tmax.txt','w')
    f.write('1000\n')
    f.close()

def make_elastic_wrapper(symprec=0.05,angle_tolerance=10):
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
    symmetry = get_symmetry(struct,symprec,angle_tolerance)
    f = open('symmetry','w')
    f.write(str(symmetry))
    f.close()

def make_elastic_deformation(struct, strain_values,symprec,angle_tolerance):
    if not os.path.exists('symmetry'):
        write_symmetry(struct,symprec,angle_tolerance)
    symmetry = open('symmetry').readlines()[0]
    print('Spacegroup detected : {}'.format(symmetry))
    
    if symmetry == 'cubic':
        print('Making 3 deformations')
        c = cubic_def(struct,symmetry,strain_values)
        print('Finished making 3 deformations')
    elif symmetry == 'hexagonal':
        print('Making 5 deformations')
        c = hexagonal_def(struct,symmetry, strain_values)
    elif symmetry == 'trigonal1':
        print('Making 6 deformations')
        c = trigonal1_def(struct,symmetry, strain_values)
    elif symmetry == 'trigonal2':
        print('Making 7 deformations')
        c = trigonal2_def(struct,symmetry, strain_values)
    elif symmetry == 'tetragonal1':
        print('Making 6 deformations')
        c = tetragonal1_def(struct,symmetry, strain_values)
    elif symmetry == 'tetragonal2':
        print('Making 7 deformations')
        c = tetragonal2_def(struct,symmetry, strain_values)
    elif symmetry == 'orthorhombic':
        print('Making 9 deformations')
        c = orthorhombic_def(struct,symmetry, strain_values)
    elif symmetry == 'monoclinic':
        print('Making 13 deformations')
        c = monoclinic_def(struct,symmetry, strain_values)
    elif symmetry == 'triclinic':
        print('Making 21 deformations')
        c = triclinic_def(struct,symmetry, strain_values)

def get_elastic_wrapper(make_plots):
    strain_values = list(pd.read_csv('strains.txt',header=None)[0].values)
    get_elastic(strain_values, make_plots)

def get_elastic(strain_values,make_plots):
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

