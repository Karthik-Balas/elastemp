from pymatgen.core import Structure
import os, pandas as pd, matplotlib.pyplot as plt
from elastemp.base.strain_analysis import get_curvature,check_convergence,get_energy_volume,plot_parabola
from elastemp.base.symmetry import get_symmetry,get_num_constants
from elastemp.thermal.response import phonon_thermal,constants_thermal

def get_thermal_wrapper():
    """ Wrapper to preprocess data and call get_thermal function.
    :param None
    :returns None. calls get_thermal function
    """
    tmax = int(open('tmax.txt').readlines()[0])
    f = open('dim.txt').readlines()[0]
    dim = [int(f.split()[i]) for i in range(len(f.split()))]
    strain_values = list(pd.read_csv('strains.txt',header=None)[0].values)
    get_thermal(tmax,dim,strain_values)

def get_thermal(tmax,dim,strain_values):
    """ Function which calls routines from thermal.response class and process data to write energy-strain values at each temperature.
    :param tmax          : Maximum temperature at which elastic constants are calculated.
           dim           : Supercell dimension for phonon calculations
           strain_values : list of strain values
    :type  tmax          : float
           dim           : list
           strain_values : list
    :returns energy_strain_temperature.csv : file storing energy-strain values at each temperature
    :rtype   energy_strain_temperature.csv : csv file
    """
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
       
        for k in range(1,num_strains+1):
            df_en = pd.read_csv('strain-{}/Phonon_calculation/energy_temp.dat'.format(k),usecols=[0,1],names=['temperature','strain-{}_energy'.format(k)],skiprows=1,sep='\t')
            if len(df2)==0:
                df2 = df_en.copy()
            else:
                df2 = df2.merge(df_en,on='temperature')
                
        energies_temp = df2.iloc[-1][1:].values     
        plot_parabola(strain_values,energies_temp,'Deformation-{}-strain-curve_tmax'.format(i),per_atom=True)
    
        df2.to_csv('energy_strain_temperature.csv',sep='\t',index=False)
        os.chdir('..')

def get_thermal_constants():
    """ Function which calls routines from thermal-class to get elastic constants at different temperatures.
    :param None
    :returns None.calls functions to get thermal elastic constants.
    :rtype   None
    """
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
    """ Function to read adiabatic & isothermal elastic constants & moduli and plots them. 
    :param None
    :returns elastic constants and moduli (adiabatic/isothermal) plots as a function of temperature.
    :rtype .png files
    """
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
        plt.figure()
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
        plt.figure()
        plt.plot(df_ct['temp'].values,df_ct[b].values,'b--',label='Isothermal')
        plt.plot(df_ct['temp'].values,df_cs[b].values,'r',label='Adiabatic')
        plt.xlabel('Temperature (K)')
        plt.legend()
        plt.ylabel(b)
        plt.savefig('{}_temperature_isothermal_adiabatic.png'.format(b))

    os.chdir('../..')
    
