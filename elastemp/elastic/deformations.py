from pymatgen.core import Structure
import numpy as np
import os,shutil
from pymatgen.analysis.elasticity.strain import Deformation,Strain
from pymatgen.analysis.elasticity.strain import convert_strain_to_deformation as csd 
from elastemp.base.symmetry import get_symmetry
from elastemp.base.strain_analysis import get_strain_matrix as smv
from elastemp.base.symmetry import get_num_constants 

class parent_def:
    def __init__(self,struct,symmetry,strain_values):
        self.struct = struct
        self.symmetry = symmetry
        self.strain_values = strain_values
        self.num_constants = get_num_constants(self.symmetry) 
        self.make_folders()
        
    
    def make_folders(self):
        for k in range(1,self.num_constants+2):
            if os.path.exists('deformation-{}'.format(k)):
                shutil.rmtree('deformation-{}'.format(k))
                os.mkdir('deformation-{}'.format(k))
            else:
                os.mkdir('deformation-{}'.format(k))
            for i in range(1,len(self.strain_values)+1):
                try:
                    os.mkdir('deformation-{}/strain-{}'.format(k,i))
                except:
                    print('Folders cannot be made or already exists!')
    

class cubic_def(parent_def):
    def __init__(self,struct,symmetry,strain_values):
        super().__init__(struct,symmetry,strain_values)
        self.struct = struct
        self.strain_values = strain_values
        self.symmetry = symmetry
        self.num_constants = get_num_constants(self.symmetry)
        self.make_deformation_folders()
        self.rename_bulk_folder()
    

    def get_strain_array(self,strain_val,k):
        '''
        1: calculates C11+C12; 2: calculates 3/2(C11+2C12)
        3: calculates 3/2 (C44) 4: calculats B
        '''
        e = strain_val
        strain_dict={'1':[e,e,0,0,0,0],'2':[e,e,e,0,0,0],'3':[0,0,0,e,e,e],'4':[e,e,e,0,0,0]}
        strain_vector = strain_dict[str(k)]
        return strain_vector
   

    def make_deformation_folders(self):
        for k in range(1,self.num_constants+2):
            for i in range(1,len(self.strain_values)+1):
                e = self.strain_values[i-1]
                strain_vector = self.get_strain_array(e,k)
                strain_arr = smv(strain_vector)
                strain = Strain(strain_arr)
                deformation = strain.get_deformation_matrix()
                new_struct = deformation.apply_to_structure(self.struct)
                new_struct.to(fmt='poscar',filename='deformation-{}/strain-{}/POSCAR'.format(k,i))
                shutil.copy('INCAR','deformation-{}/strain-{}/INCAR'.format(k,i))   
                shutil.copy('POTCAR','deformation-{}/strain-{}/POTCAR'.format(k,i))   
                shutil.copy('KPOINTS','deformation-{}/strain-{}/KPOINTS'.format(k,i))   
    
    def rename_bulk_folder(self):
        if os.path.exists('deformation-bulk'):
            shutil.rmtree('deformation-bulk')
        shutil.move('deformation-{}'.format(self.num_constants+1),'deformation-bulk')

#-------------------------------------------------------------------------#
class hexagonal_def(parent_def):
    def __init__(self,struct,symmetry,strain_values):
        super().__init__(struct,symmetry,strain_values)
        self.struct = struct
        self.strain_values = strain_values
        self.symmetry = symmetry
        self.num_constants = get_num_constants(self.symmetry)
        self.make_deformation_folders()
        self.rename_bulk_folder()


    def get_strain_array(self,strain_val,k):
        '''
        1: calculates C11+C12 ; 2: calculates 1/4 (C11-C12) ; 3 : calculates 1/2 C33
        4: calculates (C44) ; 5: calculates C11 + C12 +2C13+C33/2 ; 6: calculates B
        '''
        e = strain_val
        strain_dict={'1':[e,e,0,0,0,0],'2':[0,0,0,0,0,e],'3':[0,0,e,0,0,0],'4':[0,0,0,e,e,0],'5':[e,e,e,0,0,0], '6':[e,e,e,0,0,0]}
        strain_vector = strain_dict[str(k)]
        return strain_vector

    def make_deformation_folders(self):
        for k in range(1,self.num_constants+2):
            for i in range(1,len(self.strain_values)+1):
                e = self.strain_values[i-1]
                strain_vector = self.get_strain_array(e,k)
                strain_arr = smv(strain_vector)
                strain = Strain(strain_arr)
                deformation = strain.get_deformation_matrix()
                new_struct = deformation.apply_to_structure(self.struct)
                new_struct.to(fmt='poscar',filename='deformation-{}/strain-{}/POSCAR'.format(k,i))
                shutil.copy('INCAR','deformation-{}/strain-{}/INCAR'.format(k,i))
                shutil.copy('POTCAR','deformation-{}/strain-{}/POTCAR'.format(k,i))
                shutil.copy('KPOINTS','deformation-{}/strain-{}/KPOINTS'.format(k,i))


    def rename_bulk_folder(self):
        if os.path.exists('deformation-bulk'):
            shutil.rmtree('deformation-bulk')
        shutil.move('deformation-{}'.format(self.num_constants+1),'deformation-bulk')


#------------------------------------------------------------------------#

class trigonal1_def(parent_def):
    def __init__(self,struct,symmetry,strain_values):
        super().__init__(struct,symmetry,strain_values)
        self.struct = struct
        self.strain_values = strain_values
        self.symmetry = symmetry
        self.num_constants = get_num_constants(self.symmetry)
        self.make_deformation_folders()
        self.rename_bulk_folder()


    def get_strain_array(self,strain_val,k):
        '''
        1: calculates C11+C12 ; 2: calculates 1/4 (C11-C12) ; 3 : calculates 1/2 C33
        4: calculates (C44) ; 5: calculates C11 + C12 +2C13+C33/2 ; 6: calculates C11/4 - C12/4 + C14+ C44/2
        7: calculates B
        '''
        e = strain_val
        strain_dict={'1':[e,e,0,0,0,0],'2':[0,0,0,0,0,e],'3':[0,0,e,0,0,0],'4':[0,0,0,e,e,0],'5':[e,e,e,0,0,0],'6':[0,0,0,0,e,e],'7':[e,e,e,0,0,0]}
        strain_vector = strain_dict[str(k)]
        return strain_vector

    def make_deformation_folders(self):
        for k in range(1,self.num_constants+2):
            for i in range(1,len(self.strain_values)+1):
                e = self.strain_values[i-1]
                strain_vector = self.get_strain_array(e,k)
                strain_arr = smv(strain_vector)
                strain = Strain(strain_arr)
                deformation = strain.get_deformation_matrix()
                new_struct = deformation.apply_to_structure(self.struct)
                new_struct.to(fmt='poscar',filename='deformation-{}/strain-{}/POSCAR'.format(k,i))
                shutil.copy('INCAR','deformation-{}/strain-{}/INCAR'.format(k,i))
                shutil.copy('POTCAR','deformation-{}/strain-{}/POTCAR'.format(k,i))
                shutil.copy('KPOINTS','deformation-{}/strain-{}/KPOINTS'.format(k,i))

    def rename_bulk_folder(self):
        if os.path.exists('deformation-bulk'):
            shutil.rmtree('deformation-bulk')
        shutil.move('deformation-{}'.format(self.num_constants+1),'deformation-bulk')

#--------------------------------------------------------------------------#


class trigonal2_def(parent_def):
    def __init__(self,struct,symmetry,strain_values):
        super().__init__(struct,symmetry,strain_values)
        self.struct = struct
        self.strain_values = strain_values
        self.symmetry = symmetry
        self.num_constants = get_num_constants(self.symmetry)
        self.make_deformation_folders()
        self.rename_bulk_folder()

    def get_strain_array(self,strain_val,k):
        '''
        1: calculates C11+C12 ; 2: calculates 1/4 (C11-C12) ; 3 : calculates 1/2 C33 ; 4: calculates (C44) 
        5: calculates C11 + C12 +2C13+C33/2 ; 6: calculates C11/4 - C12/4 + C14+ C44/2 ; 7: calculates C11/4 - C12/4 - C15 + C44/2 ; 8: calculates B
        '''
        e = strain_val
        strain_dict={'1':[e,e,0,0,0,0],'2':[0,0,0,0,0,e],'3':[0,0,e,0,0,0],'4':[0,0,0,e,e,0],'5':[e,e,e,0,0,0],'6':[0,0,0,0,e,e],'7':[0,0,0,e,0,e],'8':[e,e,e,0,0,0]}
        strain_vector = strain_dict[str(k)]
        return strain_vector

    def make_deformation_folders(self):
        for k in range(1,self.num_constants+2):
            for i in range(1,len(self.strain_values)+1):
                e = self.strain_values[i-1]
                strain_vector = self.get_strain_array(e,k)
                strain_arr = smv(strain_vector)
                strain = Strain(strain_arr)
                deformation = strain.get_deformation_matrix()
                new_struct = deformation.apply_to_structure(self.struct)
                new_struct.to(fmt='poscar',filename='deformation-{}/strain-{}/POSCAR'.format(k,i))
                shutil.copy('INCAR','deformation-{}/strain-{}/INCAR'.format(k,i))
                shutil.copy('POTCAR','deformation-{}/strain-{}/POTCAR'.format(k,i))
                shutil.copy('KPOINTS','deformation-{}/strain-{}/KPOINTS'.format(k,i))

    def rename_bulk_folder(self):
        if os.path.exists('deformation-bulk'):
            shutil.rmtree('deformation-bulk')
        shutil.move('deformation-{}'.format(self.num_constants+1),'deformation-bulk')
#----------------------------------------------------------------------------#


class tetragonal1_def(parent_def):
    def __init__(self,struct,symmetry,strain_values):
        super().__init__(struct,symmetry,strain_values)
        self.struct = struct
        self.strain_values = strain_values
        self.symmetry = symmetry
        self.num_constants = get_num_constants(self.symmetry)
        self.make_deformation_folders()
        self.rename_bulk_folder()

    def get_strain_array(self,strain_val,k):
        '''
        1: calculates C11+C12 ; 2: calculates C66/2 ; 3 : calculates 1/2 C33 ; 4: calculates (C44) 
        5: calculates C11 + C12 +2C13+C33/2 ; 6: calculates C11/2 + C13 + C33/2 7: calculates B
        '''
        e = strain_val
        strain_dict={'1':[e,e,0,0,0,0],'2':[0,0,0,0,0,e],'3':[0,0,e,0,0,0],'4':[0,0,0,e,e,0],'5':[e,e,e,0,0,0],'6':[0,e,e,0,0,0],'7':[e,e,e,0,0,0]}
        strain_vector = strain_dict[str(k)]
        return strain_vector

    def make_deformation_folders(self):
        for k in range(1,self.num_constants+2):
            for i in range(1,len(self.strain_values)+1):
                e = self.strain_values[i-1]
                strain_vector = self.get_strain_array(e,k)
                strain_arr = smv(strain_vector)
                strain = Strain(strain_arr)
                deformation = strain.get_deformation_matrix()
                new_struct = deformation.apply_to_structure(self.struct)
                new_struct.to(fmt='poscar',filename='deformation-{}/strain-{}/POSCAR'.format(k,i))
                shutil.copy('INCAR','deformation-{}/strain-{}/INCAR'.format(k,i))
                shutil.copy('POTCAR','deformation-{}/strain-{}/POTCAR'.format(k,i))
                shutil.copy('KPOINTS','deformation-{}/strain-{}/KPOINTS'.format(k,i))

    def rename_bulk_folder(self):
        if os.path.exists('deformation-bulk'):
            shutil.rmtree('deformation-bulk')
        shutil.move('deformation-{}'.format(self.num_constants+1),'deformation-bulk')

#----------------------------------------------------------------------------#


class tetragonal2_def(parent_def):
    def __init__(self,struct,symmetry,strain_values):
        super().__init__(struct,symmetry,strain_values)
        self.struct = struct
        self.strain_values = strain_values
        self.symmetry = symmetry
        self.num_constants = get_num_constants(self.symmetry)
        self.make_deformation_folders()
        self.rename_bulk_folder()


    def get_strain_array(self,strain_val,k):
        '''
        1: calculates C11+C12 ; 2: calculates C66/2 ; 3 : calculates 1/2 C33 ; 4: calculates (C44) 
        5: calculates C11 + C12 +2C13+C33/2 ; 6: calculates C11/2 + C13 + C33/2 ; 7: calculates C11/2+C16+C66/2 8: calculates B
        '''
        e = strain_val
        strain_dict={'1':[e,e,0,0,0,0],'2':[0,0,0,0,0,e],'3':[0,0,e,0,0,0],'4':[0,0,0,e,e,0],'5':[e,e,e,0,0,0],'6':[0,e,e,0,0,0], '7':[e,0,0,0,0,e],'8':[e,e,e,0,0,0]}
        strain_vector = strain_dict[str(k)]
        return strain_vector

    def make_deformation_folders(self):
        for k in range(1,self.num_constants+2):
            for i in range(1,len(self.strain_values)+1):
                e = self.strain_values[i-1]
                strain_vector = self.get_strain_array(e,k)
                strain_arr = smv(strain_vector)
                strain = Strain(strain_arr)
                deformation = strain.get_deformation_matrix()
                new_struct = deformation.apply_to_structure(self.struct)
                new_struct.to(fmt='poscar',filename='deformation-{}/strain-{}/POSCAR'.format(k,i))
                shutil.copy('INCAR','deformation-{}/strain-{}/INCAR'.format(k,i))
                shutil.copy('POTCAR','deformation-{}/strain-{}/POTCAR'.format(k,i))
                shutil.copy('KPOINTS','deformation-{}/strain-{}/KPOINTS'.format(k,i))

    def rename_bulk_folder(self):
        if os.path.exists('deformation-bulk'):
            shutil.rmtree('deformation-bulk')
        shutil.move('deformation-{}'.format(self.num_constants+1),'deformation-bulk')
        
#----------------------------------------------------------------------------#

class orthorhombic_def(parent_def):
    def __init__(self,struct,symmetry,strain_values):
        super().__init__(struct,symmetry,strain_values)
        self.struct = struct
        self.strain_values = strain_values
        self.symmetry = symmetry
        self.num_constants = get_num_constants(self.symmetry)
        self.make_deformation_folders()
        self.rename_bulk_folder()

    
    def get_strain_array(self,strain_val,k):
        '''
        1: 1/2 C11 ; 2: 1/2 C22 ; 3: 1/2 C33 ; 4: 1/2 C44 ; 5: 1/2 C55 6: 1/2 C66
        7: C11/2 + C12 +C22/2 ; 8: C11/2 +C13 +C33/2 ; 9: C22/2 + C23 + C33/2 10: calculates B
        
        '''
        e = strain_val
        strain_dict={'1':[e,0,0,0,0,0],'2':[0,e,0,0,0,0],'3':[0,0,e,0,0,0],'4':[0,0,0,e,0,0],'5':[0,0,0,0,e,0],
                     '6':[0,0,0,0,0,e],'7':[e,e,0,0,0,0],'8':[e,0,e,0,0,0],'9':[0,e,e,0,0,0],'10':[e,e,e,0,0,0]}
        strain_vector = strain_dict[str(k)]
        return strain_vector

    def make_deformation_folders(self):
        for k in range(1,self.num_constants+2):
            for i in range(1,len(self.strain_values)+1):
                e = self.strain_values[i-1]
                strain_vector = self.get_strain_array(e,k)
                strain_arr = smv(strain_vector)
                strain = Strain(strain_arr)
                deformation = strain.get_deformation_matrix()
                new_struct = deformation.apply_to_structure(self.struct)
                new_struct.to(fmt='poscar',filename='deformation-{}/strain-{}/POSCAR'.format(k,i))
                shutil.copy('INCAR','deformation-{}/strain-{}/INCAR'.format(k,i))
                shutil.copy('POTCAR','deformation-{}/strain-{}/POTCAR'.format(k,i))
                shutil.copy('KPOINTS','deformation-{}/strain-{}/KPOINTS'.format(k,i))

    def rename_bulk_folder(self):
        if os.path.exists('deformation-bulk'):
            shutil.rmtree('deformation-bulk')
        shutil.move('deformation-{}'.format(self.num_constants+1),'deformation-bulk')
#---------------------------------------------------------------------------#

class monoclinic_def(parent_def):
    def __init__(self,struct,symmetry,strain_values):
        super().__init__(struct,symmetry,strain_values)
        self.struct = struct
        self.strain_values = strain_values
        self.symmetry = symmetry
        self.num_constants = get_num_constants(self.symmetry)
        self.make_deformation_folders()
        self.rename_bulk_folder()

    def get_strain_array(self,strain_val,k):
        '''
        1: 1/2 C11 ; 2: 1/2 C22 ; 3: 1/2 C33 ; 4: 1/2 C44 ; 5: 1/2 C55 6: 1/2 C66
        7: C11/2 + C12 +C22/2 ; 8: C11/2 +C13 +C33/2 ; 9: C11/2 + C15 + C55/2 ; 10: C22/2 + C23 + C33/2
        11: C22/2 + C25 + C55/2 ; 12: C33/2 + C35 + C55/2 ; 13: C44/2 + C46 + C66/2 ; 14: calculates B
        
        '''
        e = strain_val
        strain_dict={'1':[e,0,0,0,0,0],'2':[0,e,0,0,0,0],'3':[0,0,e,0,0,0],'4':[0,0,0,e,0,0],'5':[0,0,0,0,e,0],'6':[0,0,0,0,0,e],'7':[e,e,0,0,0,0],
                     '8':[e,0,e,0,0,0],'9':[e,0,0,0,e,0],'10':[0,e,e,0,0,0],'11':[0,e,0,0,e,0], '12':[0,0,e,0,e,0],'13':[0,0,0,e,0,e],'14':[e,e,e,0,0,0]}
        strain_vector = strain_dict[str(k)]
        return strain_vector

    def make_deformation_folders(self):
        for k in range(1,self.num_constants+2):
            for i in range(1,len(self.strain_values)+1):
                e = self.strain_values[i-1]
                strain_vector = self.get_strain_array(e,k)
                strain_arr = smv(strain_vector)
                strain = Strain(strain_arr)
                deformation = strain.get_deformation_matrix()
                new_struct = deformation.apply_to_structure(self.struct)
                new_struct.to(fmt='poscar',filename='deformation-{}/strain-{}/POSCAR'.format(k,i))
                shutil.copy('INCAR','deformation-{}/strain-{}/INCAR'.format(k,i))
                shutil.copy('POTCAR','deformation-{}/strain-{}/POTCAR'.format(k,i))
                shutil.copy('KPOINTS','deformation-{}/strain-{}/KPOINTS'.format(k,i))

    def rename_bulk_folder(self):
        if os.path.exists('deformation-bulk'):
            shutil.rmtree('deformation-bulk')
        shutil.move('deformation-{}'.format(self.num_constants+1),'deformation-bulk')
#-------------------------------------------------------------------------------#
class triclinic_def(parent_def):
    def __init__(self,struct,symmetry,strain_values):
        super().__init__(struct,symmetry,strain_values)
        self.struct = struct
        self.strain_values = strain_values
        self.symmetry = symmetry
        self.num_constants = get_num_constants(self.symmetry)
        self.make_deformation_folders()
        self.rename_bulk_folder() 

    def get_strain_array(self,strain_val,k):
        '''
        1: 1/2 C11 ; 2: 1/2 C22 ; 3: 1/2 C33 ; 4: 1/2 C44 ; 5: 1/2 C55 6: 1/2 C66 ;7: C11/2 + C12 +C22/2 ; 8: C11/2 +C13 +C33/2 ; 9: C11/2+C14+C44/2 ; 10: C11/2 +C15 + C55/2
        11: C11/2+C16+C66/2 ; 12: C22/2+C23+C33/2 ; 13: C22/2+C24+C44/2 ; 14: C22/2 + C25 + C55/2 ; 15: C22/2+C26+C66/2 ; 16: C33/2+ C34+ C44/2 ; 17: C33/2+C35+C55/2 ; 
        18: C33/2 + C36 + C66/2 ; 19: C44/2 + C45 + C55/2 ; 20: C44/2 + C46 + C66/2 ; 21: C55/2 + C56 + C66/2 ; 22: calculates B        
        '''
        e = strain_val
        strain_dict={'1':[e,0,0,0,0,0],'2':[0,e,0,0,0,0],'3':[0,0,e,0,0,0],'4':[0,0,0,e,0,0],'5':[0,0,0,0,e,0],'6':[0,0,0,0,0,e],'7':[e,e,0,0,0,0],
                     '8':[e,0,e,0,0,0],'9':[e,0,0,e,0,0], '10':[e,0,0,0,e,0], '11':[e,0,0,0,0,e] ,'12':[0,e,e,0,0,0], '13':[0,e,0,e,0,0], '14':[0,e,0,0,e,0],
                     '15':[0,e,0,0,0,e], '16':[0,0,e,e,0,0],'17':[0,0,e,0,e,0], '18':[0,0,e,0,0,e],'19':[0,0,0,e,e,0],'20':[0,0,0,e,0,e],'21':[0,0,0,0,e,e],
                     '22':[e,e,e,0,0,0]
                }
        strain_vector = strain_dict[str(k)]
        return strain_vector

    def make_deformation_folders(self):
        for k in range(1,self.num_constants+2):
            for i in range(1,len(self.strain_values)+1):
                e = self.strain_values[i-1]
                strain_vector = self.get_strain_array(e,k)
                strain_arr = smv(strain_vector)
                strain = Strain(strain_arr)
                deformation = strain.get_deformation_matrix()
                new_struct = deformation.apply_to_structure(self.struct)
                new_struct.to(fmt='poscar',filename='deformation-{}/strain-{}/POSCAR'.format(k,i))
                shutil.copy('INCAR','deformation-{}/strain-{}/INCAR'.format(k,i))
                shutil.copy('POTCAR','deformation-{}/strain-{}/POTCAR'.format(k,i))
                shutil.copy('KPOINTS','deformation-{}/strain-{}/KPOINTS'.format(k,i))

    def rename_bulk_folder(self):
        if os.path.exists('deformation-bulk'):
            shutil.rmtree('deformation-bulk')
        shutil.move('deformation-{}'.format(self.num_constants+1),'deformation-bulk')
#--------------------------------------------------------------------------#


