import os,yaml,numpy as np,sys, pandas as pd
from elastemp.base.strain_analysis  import fit_parabola,fit_parabola_bulk
from elastemp.elastic.constants import compute_elastic_constants 
from elastemp.elastic.response import compute_response_functions 
from pymatgen.core import Structure

class phonon_thermal:
    """ 
    Phonon class consisting of functions to run phonopy phonon calculations.
    """
    def __init__(self,dim,tmax):
        """ Constructor function for phonon_thermal class
        :param dim  : dimension of supercell
               tmax : Maximum temperature at which elastic constants calculations are run.
        """
        self.dim=dim
        self.tmax = tmax

    def run_phonopy_thermal(self):
        """ Function to run thermal properties calculations via phonopy.
        """
        self.make_mesh_conf()
        print(os.getcwd())
        os.system('phonopy -c POSCAR-unitcell -t mesh.conf')

    def make_mesh_conf(self):
        """ Function to create mesh file needed for phonopy calculations.
        :param None
        :returns mesh.conf  : phonopy settings file
        :rtype   mesh.conf  : plaintext file
        """
        dim1,dim2,dim3 = self.dim
        f = open('mesh.conf','w')
        mesh_string = 'DIM={} {} {}\nMP= 48 48 48\nTMAX={}\nTSTEP=2\n'.format(dim1,dim2,dim3,self.tmax)
        f.write(mesh_string)
        f.write('FORCE_CONSTANTS=READ')
        

    def extract_thermal_yaml(self):
        """ Function to extract thermodynamic properties from thermal_properties.yaml file.
        :param None
        :returns energy_temp file containing thermodynamic properties
        :rtype   energy_temp :dat file
        """
        fil = 'thermal_properties.yaml'
        try:
            data = yaml.load(open(fil),Loader=yaml.FullLoader)
        except:
            print('Unable to read thermal_properties file. Program will stop!')
            sys.exit()
        zpe = float(open('ZPE.txt').readlines()[0])
        
        struct = Structure.from_file('POSCAR-unitcell')
        natoms = len(struct.sites)
        df = pd.DataFrame()
        temp_array = []; energy_array = []; cv_array=[]
        data_th = data['thermal_properties']
        
        for m in range(len(data_th)):
            temp = data_th[m]['temperature']
            en1 = (data_th[m]['free_energy']*0.01)/(natoms)
            cv = (data_th[m]['heat_capacity'])/(natoms)
            en1 = en1
            en = en1+zpe
            temp_array.append(temp)
            cv_array.append(cv)
            energy_array.append(en)
            
        df['temperature']=temp_array
        df['helm_energy']=energy_array
        df['Cv']=cv_array
        df.to_csv('energy_temp.dat',index=False,sep='\t')



class constants_thermal:
    """ class which extract elastic constants at different temperatures.
    """
    def __init__(self,num_constants,strain_values,symmetry,tmax,V0):
        """ Constructor function for constant thermal class
        :param num_constants : number of independent constants for the structure based on its symmetry
               strain_values : list of strain values
               symmetry      : symmetry of the structure
               tmax          : maximum temperature at which elastic constant calculations are run
               V0            : volume of reference structure
        """
        self.num_constants= num_constants
        self.strain_values = strain_values
        self.symmetry = symmetry
        self.tmax = tmax
        self.V0 = V0

    def extract_thermal_constants(self):
        
        """ function to extract isothermal elastic constants from parabolic fits of energy,strains at temperatures.
        :param None
        :returns None. calls the class function to convert isothermal elastic constants to adiabatic elastic constants.
        """

        print('Extracting elastic constants vs temperature')
        obj_elastic = compute_elastic_constants(self.num_constants,self.strain_values,self.symmetry)
        
        f = open('results_dir/volume_temp.txt','r')
        V_T=float(f.readlines()[-1].split()[1])
        natoms = len(Structure.from_file('POSCAR').sites)
        V_Tpa = V_T/natoms
        df = pd.read_csv('results_dir/elastic_constants_rawdata_0K.txt',sep='\t',header=None)

        c0 = df[1].values
        curvatures_t = []
        for j in range(1,self.num_constants+1):
            df = pd.read_csv('deformation-{}/energy_strain_temperature.csv'.format(j),sep='\t')
            energies_t = [df.iloc[-1][1:].values[k] for k in range(len(self.strain_values))]
            
            curvature_const,p = fit_parabola(energies_t,self.strain_values,V_Tpa)
            curvatures_t.append(curvature_const)

        ct = obj_elastic.get_elastic_constants(curvatures_t,write_elastic=False)
        ct_array = []
        for i in range(6):
            for j in range(6):
                val = ct[i][j]
                ct_array.append(val)
        
        
        self.convert_isothermal_to_adiabatic(ct_array,c0)
       


    
    def convert_isothermal_to_adiabatic(self,ct_array,c0):
        
        """ Function to convert isothermal elastic constants to adiabatic. Calculates isothermal and adiabatic moduli. 
        :param ct_array : array of isothermal elastic constants.
               c0       : list of elastic constants at 0 K
        :returns isothermal & adiabatic elastic constants and moduli in pandas dataframes and saved as csv files.
        :rtype   .csv file
        """
        
        symmetry = open('symmetry').readlines()[0]
        num_T_points = 1+int(self.tmax/2)
        df_ct = pd.DataFrame()
        df_cs = pd.DataFrame()
        df_mod_t = pd.DataFrame()
        df_mod_s = pd.DataFrame()
    
        a = [1,1,1,0,0,0]
        if symmetry in ['cubic','trigonal1','trigonal2','hexagonal','orthorhombic','tetragonal1','tetragonal2']:
             pass
        elif symmetry == 'monoclinic':
             a[4]= 1/2; 
        elif symmetry == 'triclinic':
            a[3]=a[4]=a[5]=1/2
        else: 
            print('symmetry found not in known ones. Program will end')
            sys.exit()
        

        
        l_arr = []
        for i in range(36):
            if c0[i]==0:
                l_arr.append(0)
            else:
                val = (ct_array[i]-c0[i])/self.tmax
                l_arr.append(val)

        f = open('results_dir/lambda_vals.txt','w')
        for i in range(len(l_arr)):
            f.write(str(l_arr[i])+'\n')
        f.close()

        for k in range(num_T_points):
            temp = k*2
            ct_row= [temp] 
            cs_row = [temp]
            cs = np.zeros((6,6))
            m = np.zeros(6)

            for j in range(36):
                val = c0[j]+l_arr[j]*temp
                ct_row.append(val)
           
            c_temp_t = np.reshape(ct_row[1:], (-1, 6))
            
          
            for i in range(6):
                for t in range(6):
                    m[i] +=c_temp_t[i][t]*a[t]
            
           
            phi_T = float(open('results_dir/volume_temp.txt','r').readlines()[k+1].split()[-1])
            
           
            for i in range(6):
                for j in range(6):
                  
                    cs[i][j] = c_temp_t[i][j] + phi_T*m[i]*m[j]
                    cs_row.append(cs[i][j])
            c_temp_s = np.reshape(cs_row[1:], (-1, 6))
            c11_ct_avg = (ct_row[1]+ct_row[8]+ct_row[15])/3.0
            c12_ct_avg = (ct_row[2]+ct_row[3]+ct_row[9])/3.0
            bulk_mod_ct = (c11_ct_avg+2*c12_ct_avg)/3.0
            
            c11_cs_avg = (cs_row[1]+cs_row[8]+cs_row[15])/3.0
            c12_cs_avg = (cs_row[2]+cs_row[3]+cs_row[9])/3.0
            bulk_mod_cs = (c11_cs_avg+2*c12_cs_avg)/3.0
             
            obj_response_t = compute_response_functions(c_temp_t,self.symmetry)
            s_t = obj_response_t.get_compliance_matrix(write_compliance=False)       
            
            
            obj_response_s = compute_response_functions(c_temp_s,self.symmetry)
            s_s = obj_response_s.get_compliance_matrix(write_compliance=False)   
            
            Eh,Gh, hc, nu, pugh = obj_response_t.get_moduli(s_t,write_moduli=False)
            moduli_row_t = [temp,bulk_mod_ct, Eh, Gh, hc, nu]
            
            df = pd.DataFrame(moduli_row_t)
            df_mod_t = pd.concat([df_mod_t,df.T],ignore_index=True)
            
            Eh,Gh, hc, nu, pugh = obj_response_s.get_moduli(s_s,write_moduli=False)
            moduli_row_s = [temp,bulk_mod_cs, Eh, Gh, hc, nu]
            
            
        
            df = pd.DataFrame(moduli_row_s)
            df_mod_s = pd.concat([df_mod_s,df.T],ignore_index=True)
            
            df = pd.DataFrame(ct_row)
            df_ct = pd.concat([df_ct,df.T],ignore_index=True)
            
            df = pd.DataFrame(cs_row)
            df_cs = pd.concat([df_cs,df.T],ignore_index=True)
                        
            
        elastic_col_names = ['temp']
        
        for i in range(6):
            for j in range(6):
                elastic_col_names.append('c{}{}'.format(i+1,j+1))
        
        df_ct.columns = [elastic_col_names]
        df_cs.columns = [elastic_col_names]
        
        
        df_ct.to_csv('results_dir/elastic_constants_temp_isothermal.csv',index=False)      
        df_cs.to_csv('results_dir/elastic_constants_temp_adiabatic.csv',index=False)
        
        
        df_mod_t.columns = ['T','K','Eh','Gh','hc','nu']
        df_mod_s.columns = ['T','K','Eh','Gh','hc','nu']
        
        
        df_mod_t.to_csv('results_dir/moduli_temp_isothermal.csv',sep='\t',index=False)  
        df_mod_s.to_csv('results_dir/moduli_temp_adiabatic.csv',sep='\t',index=False) 

        
        

    
