import numpy as np
import pandas as pd, os
import matplotlib.pyplot as plt
from pymatgen.core import Structure
from pymatgen.analysis.eos import EOS


def get_strain_matrix(vector):
    """ returns 3x3 strain matrix from 6x1 strain vector
    :param vector : 6 element strain vector in a list format
    :type  vector : list
    :returns 3x3 strain matrix
    :rtype matrix
    """
    eps1 = vector[0]
    eps2 = vector[1]
    eps3 = vector[2]
    eps4 = vector[3]
    eps5 = vector[4]
    eps6 = vector[5]
    strain_matrix = np.array([[eps1,eps6/2, eps5/2],[eps6/2,eps2,eps4/2],[eps5/2,eps4/2,eps3]],dtype=object)
    return strain_matrix


def check_convergence(fil):
    """ checks convergence of OSZICAR file and returns converged last line. 
        If not converged, returns arbitrary 1e4 value.
    :param fil :file to check for convergence
    :type  fil : plaintext
    :returns Boolean value (True for converged, False for not), last line with energy values. 
    :rtype   Bool (True/False), list
    """
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
    """ fits parabola with energies and strains to return fitting parameters.
    :param energies    : energies to fit parabola
           strain_array: strain values to fit parabola
           V0          : volume of reference structure.
    :type  energies    : list
           strain_array: list
           V0          : float
    :returns curvature and fitting parameters of the fitted parabola
    :rtype   float,list
    """
    p = np.polyfit(strain_array,energies,2)
    curvature_const = np.round((p[0]*160.217/V0),4)
    return curvature_const,p

def fit_parabola_bulk(energies,volumes,calc_bulk_parabola=True):
    """ fits equation of state using Birch-Murnaghan eos to get equilibrium volume, bulk modulus
    :param energies    : energies to fit parabola
           volumes     : volumes to fit to EOS
    :type  energies    : list
           volumes     : list
    :returns bulk modulues
    :rtype   float
    """
 
    eb = EOS(eos_name='murnaghan')
    eb_fit = eb.fit(volumes,energies).results
    b0 = eb_fit['b0']*160.217
    v = eb_fit['v0']
    if calc_bulk_parabola:
        bulk_mod = np.round(b0,4)
        f = open('B.txt','w')
        f.write('B = {} GPa\n'.format(bulk_mod))
        f.close()
        return bulk_mod
    else:
        return v

def get_energy_volume(dir_path):
    """ Returns enegy/atom for chosen structure - OSZICAR file
    :param dir_path    : directory path with OSZICAR file
    :type  dir_path    : string
    :returns energy/atom, vol of structure
    :rtype   float,float
    """

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
    """ Returns curvature constant after fitting parabola. 
        Returns bulk mod in case of bulk modulus fitting.
    :param strain_array    : list of strain values
           calc_bulk       : Boolean value determining if fit is for bulk modulus
    :type  strain_array    : list
           calc_bulk       : Bool (True/False)
    :returns curvature constant/bulk modulus
    :rtype   float
    """

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
        bulk_mod = fit_parabola_bulk(energies,volumes,calc_bulk_parabola=True)

        return bulk_mod

def plot_parabola(strains,energies,title,per_atom=True):
    """ plots parabolic function of strains and energies.
    :param strains  : strain values to fit parabola
           energies : energies to fit parabola
           title    : title to be used in plot
           per_atom : Boolean value to determine ylabel (eV or eV/atom)
    :type  strains  : list
           energies : list
           title    : string
           per_atom : Bool(True/False)
    :returns plot of energy-strains.
    :rtype   png figure.
    """
    fit = np.polyfit(strains,energies,2)
    st_min = min(strains)
    st_max = max(strains)
    st_vals = [st_min + (st_max-st_min)*i/100 for i in range(101)]
    e_fit = [fit[0]*st_vals[k]**2 + fit[1]*st_vals[k] + fit[2] for k in range(len(st_vals))]
    plt.figure()
    plt.plot(strains,energies,'ro',label='energies')
    plt.plot(st_vals,e_fit,'r--',label='fit')
    plt.title(title)
    plt.xlabel('strain')
    if per_atom==True:
        plt.ylabel('Energy (eV/atom)')
    elif per_atom == False:
        plt.ylabel('Energy (eV)')
    plt.legend()
    plt.savefig(title+'.png')
    plt.close()
