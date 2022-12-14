from pymatgen.core import Structure
import numpy as np,os, pandas as pd,sys, shutil
import matplotlib.pyplot as plt
from elastemp.base.strain_analysis import get_curvature,check_convergence,get_energy_volume,fit_parabola_bulk,plot_parabola
from elastemp.base.symmetry import get_symmetry,get_num_constants
from elastemp.base.makefiles import make_incar_dynamic, make_kpt_dynamic
from elastemp.dynamic.deformations import dynamic_deformations,dynamic_response
from elastemp.thermal.response import phonon_thermal,constants_thermal


def make_dim():
    """ Returns a default file dim.txt to be used in supercell dimensions for phonon calculations.
    :param None
    :returns: dim.txt with default value [2,2,2]
    :rtype : .txt file
    """
    f = open('dim.txt','w')
    f.write('{} {} {}'.format(2,2,2))
    f.close()


def make_tmax():
    """ Returns a default file tmax.txt which stores the maximum temperature at which elastic constants are calculated.
    :param None
    :returns: tmax.txt with default value 1000K
    :rtype : .txt file
    """
    f = open('tmax.txt','w')
    f.write('1000\n')
    f.close()

def make_dynamic_wrapper(calc_bulk_dynamic):
    """ Wrapper for make_dynamic function to do preprocesing. 
    :param calc_bulk_dynamic : Boolean value to determine if dynamic calculation is for deformation-bulk or other modes.
    :type  calc_bulk_dynamic : Bool (True/False)
    :returns None. calls make_bulk_dynamic or make_dynamic functions. 
    """
    struct = Structure.from_file('POSCAR')
    strain_values = list(pd.read_csv('strains.txt',header=None)[0].values)
    
    if os.path.exists('KPOINTS_dynamic'):
        pass
    else:
        print('KPOINTS_dynamic is missing. Making a generic KPOINTS_dynamic file')
        make_kpt_dynamic()
    
    if os.path.exists('INCAR_dynamic'):
        pass
    else:
        print('INCAR_dynamic is missing. Making a generic INCAR_dynamic file')
        make_incar_dynamic()

    if not os.path.exists('POTCAR'):
        print('POTCAR file is missing! Cannot proceed with calculation.')
        sys.exit()

    try:
        f = open('dim.txt').readlines()[0]
        dim = [int(f.split()[i]) for i in range(len(f.split()))]
        print('Read supercell dimension sucessfully!')
    except:
        print('Unable to read dim file from Input_files. Adding a default supercell dimension of 2 2 2')
        dim = [2,2,2]
        make_dim()
    
    if calc_bulk_dynamic:
        make_bulk_dynamic(struct,strain_values,dim)
    else:
        make_dynamic(struct,strain_values,dim)    

def make_bulk_dynamic(struct,strain_values,dim):
    """ Function to call dynamic_deformations routine from dynamic class to make volumetric strain deformations.
    :param struct       : Structure
           strain_values: list of strain values
           dim          : dimension of supercell
    :type  struct       : Pymatgen structure object
           strain_values: list
           dim          : list
    :returns None. calls routine to make bulk deformations.
    """
    print(os.getcwd())
    symmetry = open('symmetry').readlines()[0]
    e = dynamic_deformations(symmetry,strain_values, dim, struct)
    e.make_bulk_deformation()


def get_dynamic_wrapper(calc_bulk_dynamic):
    """ Wrapper to call function to get either bulk or other deformation modes.
    :param calc_bulk_dynamic : Boolean value to call either bulk or other modes.
    :type  calc_bulk_dynamic : Bool (True/False)
    :returns  None. calls get_bulk_dynamic or get_dynamic functions.
    """
    strain_values = list(pd.read_csv('strains.txt',header=None)[0].values)
    if calc_bulk_dynamic==True:
        get_bulk_dynamic(strain_values)
    else:
        get_dynamic(strain_values)

def get_bulk_dynamic(strain_values):
    """ Function to get thermal properties like thermal expansion coeff, Vols at different temp
    :param strain_values : list of strain values
    :type  strain_values : list
    :return volume_temp.txt : file which stores vols at different temp, thermal expansion coeff, Cv and phi values
            Themal_expansion_temperature.png :plot of thermal expansion at different temperatures
    :rtype  volume_temp.txt : .txt file
            Thermal_expansion_temperature.png : .png file
    """
    os.chdir('deformation-bulk')
    df_en_temp = pd.DataFrame()
    df_cv_temp = pd.DataFrame()
    
    for j in range(1,len(strain_values)+1):
        os.chdir('strain-{}/Phonon_calculation'.format(j))
        print('Curr_dir :', os.getcwd())
        e = dynamic_response()
        os.chdir('../..')
    os.chdir('..')

    try:
        tmax = int(open('tmax.txt').readlines()[0])
    except:
        print('Unable to read max temperature from tmax.txt. Setting default value of 1000K')
        tmax=1000
        make_tmax()

    f = open('dim.txt').readlines()[0]
    dim = [int(f.split()[i]) for i in range(len(f.split()))]
    e = phonon_thermal(dim,tmax)
    
    os.chdir('deformation-bulk')

    for j in range(1,len(strain_values)+1):
        os.chdir('strain-{}/Phonon_calculation'.format(j))
        print('Curr_dir: ',os.getcwd())
        e.run_phonopy_thermal()
        e.extract_thermal_yaml()
        os.chdir('../..')
    
    for i in range(1,len(strain_values)+1):
        df_en = pd.read_csv('strain-{}/Phonon_calculation/energy_temp.dat'.format(i),usecols=[0,1],names=['temperature','strain-{}_energy'.format(i)],skiprows=1,sep='\t')
        df_cv = pd.read_csv('strain-{}/Phonon_calculation/energy_temp.dat'.format(i),usecols=[0,2],names=['temperature','strain-{}_cv'.format(i)],skiprows=1,sep='\t')
        
        if len(df_en_temp)==0:
            df_en_temp = df_en.copy()
        else:
            df_en_temp = df_en_temp.merge(df_en,on='temperature')
        
        if len(df_cv_temp)==0:
            df_cv_temp = df_cv.copy()
        else:
            df_cv_temp = df_cv_temp.merge(df_cv,on='temperature')
            
        
        df_en_temp.to_csv('energy_strain_temperature.csv',sep='\t',index=False)
        df_cv_temp.to_csv('Cv_strain_temperature.csv',sep='\t',index=False)
    
    
    energies_temp = df_en_temp.iloc[-1][1:].values
    plot_parabola(strain_values,energies_temp,'Deformation-bulk-curve-tmax',per_atom=True) 
    
    df3 = pd.read_csv('energies.txt',sep='\t')
    vols = df3['volumes']
    
    temp = [] ;     e_min = [] ;vols_T=[] ; lin_exp=[]
    os.chdir('..')
    natoms = len(Structure.from_file('POSCAR').sites)
    V0 = Structure.from_file('POSCAR').volume
    
    f = open('results_dir/volume_temp.txt','w')
    f.write('temperature \t Volume_T\tVolume_0\tVolume_thermal_expansion\tLinear_thermal_expansion\tCv_calc\tLinear_expansion_factor\tphi\n')
    
   
    
    for i in range(len(df_en_temp)):
        energies = df_en_temp.iloc[i][1:].values
        cvs = df_cv_temp.iloc[i][1:].values
        V_T1 = fit_parabola_bulk(energies,vols,calc_bulk_parabola=False)
        cv_fit = np.polyfit(strain_values,cvs,2) 
        temperature = df_en_temp.iloc[i][0]

        vols_T.append(V_T1)
        temp.append(temperature)
        vol_thermal_exp=0
        vpa_T = V_T1/natoms

        if i!=0:
            dV = (vols_T[i] - vols_T[i-1])/vols_T[i-1]
            vol_thermal_exp = dV/(temp[i]-temp[i-1])
            
        lin_thermal_exp=vol_thermal_exp/3
        lin_exp.append(lin_thermal_exp*1e6)
        
        eq_strain= lin_thermal_exp*temp[i]
        cv_calc = cv_fit[0]*eq_strain**2+cv_fit[1]*eq_strain+cv_fit[2]
        lin_expansion_factor = 1+eq_strain
        phi=0
        if i!=0 and cv_calc!=0:
            phi = (1e5/160.217)*(temperature*vpa_T/cv_calc)*lin_thermal_exp**2
        f.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(temperature,V_T1,V0,vol_thermal_exp,lin_thermal_exp,cv_calc,lin_expansion_factor,phi))
    f.close()
    
    fig = plt.figure()
    plt.plot(temp,lin_exp,'r')
    plt.title('Linear_thermal_expansion vs temperature')
    plt.xlabel('Temperature (K)')
    plt.ylabel(' Linear_thermal_expansion * 1e6 (/K)')
    plt.savefig('results_dir/Thermal_expansion_temperature.png')
    
    

def make_T_Zpe():
    """ function to make files to run enthalpy calculations for structures with lattice parameters corresponding to tmax.
    :param None.
    :returns saves strained POSCAR files corresponding to each deformation/strain.
    :rtype POSCAR plaintext file
    """
    if not os.path.exists('results_dir/volume_temp.txt'):
        print('volume_temp file does not exist. Program will stop!')
        sys.exit()   
    print('Getting linear expansion factor at tmax')
    symmetry = open('symmetry').readlines()[0]
    lin_expansion_factor = float(open('results_dir/volume_temp.txt').readlines()[-1].split()[6])
    print('linear_exp',lin_expansion_factor)
    num_constants = get_num_constants(symmetry)
    strain_values = list(pd.read_csv('strains.txt',header=None)[0].values)

    for i in range(1,num_constants+1):
        for j in range(1,len(strain_values)+1):
            os.chdir('deformation-{}/strain-{}'.format(i,j))
            if not os.path.exists('Hightemp_static'):
                os.mkdir('Hightemp_static')
            os.chdir('Hightemp_static')
            shutil.copy('../INCAR','INCAR')
            shutil.copy('../KPOINTS','KPOINTS')
            shutil.copy('../POTCAR','POTCAR')
            shutil.copy('../POSCAR','POSCAR')
            
            with open('POSCAR', 'r') as file:
                data = file.readlines()  
                data[1] = "{}\n".format(lin_expansion_factor)
            with open('POSCAR', 'w') as file:
                file.writelines(data)
            os.chdir('../../..')
    print('Created hightemp_static dir and corresponding files to obtain the ZPE at tmax')


def get_T_Zpe():
    """ function to return zero temp enthalpies of structures corresponding to lattice parameters at tmax
    :param   None.Reads symmetry, num_constants, strain_values from stored files.
    :returns None. saves ZPE.txt files (enthalpy/atom) in plaintext.
    """
    symmetry = open('symmetry').readlines()[0]
    num_constants = get_num_constants(symmetry)
    strain_values = list(pd.read_csv('strains.txt',header=None)[0].values)
    for i in range(1,num_constants+1):
        for j in range(1,len(strain_values)+1):
            os.chdir('deformation-{}/strain-{}/Hightemp_static'.format(i,j))
            if not os.path.exists('OSZICAR'):
                print('OSZICAR file does not exist. Program will stop!')
                sys.exit()
            convergence_val,l = check_convergence('OSZICAR')
            if convergence_val:
                energy,volume = get_energy_volume('../Hightemp_static')
            os.chdir('../../..')
    print('Calculated ZPE at tmax. Please proceed to do phonon calculations for each deformation/strain')

def make_dynamic(struct,strain_values,dim):
    """ function calling routine from dynamic class to make deformations for phonon calculations.
    :param struct        : Structure
           strain_values : list of strain values
           dim           : supercell dimension
    :type  struct        : Pymatgen structure object
           strain_values : list
           dim           : list
    :returns None
    """
    symmetry = open('symmetry').readlines()[0]
    e = dynamic_deformations(symmetry,strain_values,dim,struct)
    e.make_deformation()

def get_dynamic(strain_values):
    """ function to call routine from dynamic class to get dynamic properties. 
    :param strain_values : list of strain values
    :type  strain_values : list
    :returns None
    """
    symmetry = open('symmetry').readlines()[0]
    num_constants = get_num_constants(symmetry)

    for i in range(1,num_constants+1):
        os.chdir('deformation-{}'.format(i))
        for j in range(1,len(strain_values)+1):
            os.chdir('strain-{}/Phonon_calculation'.format(j))
            dynamic_response()
            os.chdir('../..')
        os.chdir('..')
