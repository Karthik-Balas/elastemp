import numpy as np
import numpy.polynomial.polynomial as poly
import pandas as pd, os
from pymatgen.core import Structure
from pymatgen.analysis.eos import EOS
import shutil
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer as sga

def get_strain_matrix(vector):
    eps1 = vector[0]
    eps2 = vector[1]
    eps3 = vector[2]
    eps4 = vector[3]
    eps5 = vector[4]
    eps6 = vector[5]
    strain_matrix = np.array([[eps1,eps6/2, eps5/2],[eps6/2,eps2,eps4/2],[eps5/2,eps4/2,eps3]],dtype=object)
    return strain_matrix


def check_convergence(fil):
    try:
        g = open(fil,'r')
        l = [line for line in g]
        if l[-1].find('F= ')==-1:
            print('{} not converged!'.format(dir_path))
            return False, l[-1]
        return True, l[-1]
    except:
        l = [1e4,1e4,1e4]
        return False,l

def fit_parabola(energies,strain_array,V0):
    p = np.polyfit(strain_array,energies,2)
    curvature_const = np.round((p[0]*160.217/V0),4)
    return curvature_const,p

def fit_parabola_temp(energies,strain_array,emin,V0):
    p = np.polyfit(strain_array,energies,2)
    curvature_const = np.round((p[0]*160.217/V0),4)
    return curvature_const,p

def fit_parabola_bulk(energies,volumes):
    eb = EOS(eos_name='murnaghan')
    eb_fit = eb.fit(volumes,energies).results
    b0 = eb_fit['b0']*160.217
    bulk_mod = np.round(b0,4)
    f = open('B.txt','w')
    f.write('B = {} GPa\n'.format(bulk_mod))
    f.close()
    return bulk_mod

def get_energy_volume(dir_path):
    g = open('{}/OSZICAR'.format(dir_path),'r')
    g1 = Structure.from_file('{}/POSCAR'.format(dir_path))
    natoms = len(g1.sites) 
    l = [line for line in g]
    energy = float(l[-1].split()[2])
    zpe = energy/natoms
    f = open('{}/ZPE.txt'.format(dir_path),'w')
    f.write(str(zpe))
    f.close()
    vol = g1.volume
    return energy, vol


def get_curvature(strain_array,calc_bulk=False):
    energies = []
    volumes=[]
    dir_list = ['strain-{}'.format(i) for i in range(1,len(strain_array)+1)]
    V0 = Structure.from_file('../POSCAR').volume
    for i in dir_list:
        if os.path.exists(i):
            convergence_val,l = check_convergence('{}/OSZICAR'.format(i))
            if convergence_val:
                energy,volume = get_energy_volume(i)
                energies.append(energy)
                volumes.append(volume)
            else:
                print(os.getcwd())
                print('File not found or not converged. Adding arbitrarily high energy value and volume')
                energies.append(1e4)
                volumes.append(1e4)
        else:
            print('{} folder does not exist'.format(i))

    if calc_bulk == False:
        df = pd.DataFrame()
        df['energies']=energies
        df['strains']=strain_array
        df.to_csv('energies.txt',sep='\t',index=False)
        curvature_const,p = fit_parabola(energies,strain_array,V0)
        f = open('curvature_constant','w')
        f.write('curvature_constant = {} GPa\n'.format(curvature_const))
        f.close()
        g = open('Fitting_func','w')
        g.write('{},{},{}'.format(p[0],p[1],p[2]))
        g.close()
        return curvature_const
    
    elif calc_bulk==True:
        df = pd.DataFrame()
        df['energies'] = energies
        df['volumes'] = volumes
        df.to_csv('energies.txt',sep='\t',index=False)
        bulk_mod = fit_parabola_bulk(energies,volumes)

        return bulk_mod

