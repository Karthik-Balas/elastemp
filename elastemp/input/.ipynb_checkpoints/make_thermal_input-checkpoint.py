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

def get_thermal_wrapper():
    tmax = int(open('tmax.txt').readlines()[0])

    f = open('dim.txt').readlines()[0]
    dim = [int(f.split()[i]) for i in range(len(f.split()))]
    strain_values = list(pd.read_csv('strains.txt',header=None)[0].values)
    get_thermal(tmax,dim,strain_values)

def get_thermal(tmax,dim,strain_values):
    struct = Structure.from_file('POSCAR')
    try:
        symmetry = open('symmetry').readlines()[0]
    except:
        symmetry = get_symmetry(struct,symprec=0.05,angle_tolerance=10)
    num_constants = get_num_constants(symmetry)
    num_strains = len(strain_values)
    
    e = phonon_thermal(dim,tmax)
    
    for i in range(1,num_constants+1):
       
        print('curr_dir :',os.getcwd())
        os.chdir('deformation-{}'.format(i))
      
        for j in range(1,num_strains+1):
            os.chdir('strain-{}/Phonon_calculation'.format(j))   
            e.run_phonopy_thermal()
            e.extract_thermal_yaml()
            os.chdir('../..')
            df2 = pd.DataFrame()
       
        for i in range(1,num_strains+1):
            df_en = pd.read_csv('strain-{}/Phonon_calculation/energy_temp.dat'.format(i),usecols=[0,1],names=['temperature','strain-{}_energy'.format(i)],skiprows=1,sep='\t')
            if len(df2)==0:
                df2 = df_en.copy()
            else:
                df2 = df2.merge(df_en,on='temperature')
                
        energies_temp = df2.iloc[-1][1:].values     
        fit = np.polyfit(strain_values,energies_temp,2)
        st_min = min(strain_values)
        st_max = max(strain_values)
        st_vals = [st_min + (st_max-st_min)*i/100 for i in range(101)]
        e_fit = [fit[0]*st_vals[k]**2 + fit[1]*st_vals[k] + fit[2] for k in range(len(st_vals))]
        #e_fit = [fit[0]*strain_values[k]**2+fit[1]*strain_values[k]+fit[2] for k in range(len(strain_values))]
        fig = plt.figure()
        plt.plot(strain_values,energies_temp,'ro',label='energies')
        #plt.plot(strain_values,e_fit,'r--',label='fit')
        plt.plot(st_vals,e_fit,'r--',label='fit')
        plt.title('Energy-strain curve at tmax')
        plt.xlabel('strain')
        plt.legend()
        plt.ylabel('Energy (eV/atom)')
        plt.savefig('energy-strain_tmax.png')
    
        df2.to_csv('energy_strain_temperature.csv',sep='\t',index=False)
        os.chdir('..')

def get_thermal_constants():
    print(os.getcwd())
    struct = Structure.from_file('POSCAR')
    try:
        symmetry = open('symmetry').readlines()[0]
    except:
        symmetry = get_symmetry(struct,symprec=0.05,angle_tolerance=10)
    num_constants = get_num_constants(symmetry)
    tmax = int(open('tmax.txt').readlines()[0]) 
    natoms = len(struct.sites)
    V0 = struct.volume/natoms
    strain_values = list(pd.read_csv('strains.txt',header=None)[0].values)
  
    e = constants_thermal(num_constants,strain_values,symmetry,tmax,V0)
    e.extract_thermal_constants() 

def get_elastic_temp_plots():
    os.chdir('results_dir')
    df_mod_t = pd.read_csv('moduli_temp_isothermal.csv',sep='\t')
    df_mod_s = pd.read_csv('moduli_temp_adiabatic.csv',sep='\t')
    df_ct = pd.read_csv('elastic_constants_temp_isothermal.csv')
    df_cs = pd.read_csv('elastic_constants_temp_adiabatic.csv')
    if not os.path.exists('moduli_plots'):
        os.mkdir('moduli_plots')
        
    os.chdir('moduli_plots')
    moduli_names = ['K','Eh','Gh','hc','nu']
    moduli_label_names= ['Bulk_modulus (GPa)','Elastic_modulus (GPa)','Shear_modulus (GPa)','Hardness (GPa)','poisson_ratio']
 
  
    for a,b in zip(moduli_names,moduli_label_names):
        fig = plt.figure()
        plt.plot(df_mod_t['T'].values,df_mod_t[a].values,'b--',label='Isothermal')
        plt.plot(df_mod_t['T'].values,df_mod_s[a].values,'r',label='Adiabatic')
        plt.xlabel('Temperature (K)')
        plt.ylabel(b)
        plt.legend()
        plt.savefig('{}_temperature_isothermal_adiabatic.png'.format(b))

    print('Plotting elastic_constants vs temperature') 
    constant_num = [str(i) for i in range(1,37)]
    constant_names = []
    for i in range(1,7):
        for j in range(i,7):
            constant_names.append('c{}{}'.format(i,j))

   
    
    for a,b in zip(constant_num, constant_names):
        fig = plt.figure()
        plt.plot(df_ct['temp'].values,df_ct[b].values,'b--',label='Isothermal')
        plt.plot(df_ct['temp'].values,df_cs[b].values,'r',label='Adiabatic')
        plt.xlabel('Temperature (K)')
        plt.legend()
        plt.ylabel(b)
        plt.savefig('{}_temperature_isothermal_adiabatic.png'.format(b))

    os.chdir('../..')
    
