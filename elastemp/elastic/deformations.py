import os,shutil
from pymatgen.analysis.elasticity.strain import Deformation,Strain
from elastemp.base.strain_analysis import get_strain_matrix as smv
from elastemp.base.symmetry import get_num_constants 

class parent_def:
    """
    Class parent_def which is inherited by the subsequent classes to create folders.
    """
    def __init__(self,struct,symmetry,strain_values):
        """ Constructor function for parent_def class
        :param struct        : Reference structure
               symmetry      : symmetry of the reference structure
               strain_values : list of strain values
        :type  struct        : Pymatgen structure object
               symmetry      : plaintext file
               strain_values : list
        :returns None.calls make_folders function
        """
        self.struct = struct
        self.symmetry = symmetry
        self.strain_values = strain_values
        self.num_constants = get_num_constants(self.symmetry) 
        self.make_base_folders()
        
    
    def make_base_folders(self):
        """ Function to make folders in which vasp calculations are run.
        :param   None
        :returns creates num_constants + 1 folders and subfolders corresponding to strain values to run vasp calculations.
        :rtype   folders
        """
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

class make_folders:
    """
    Class to make the different deformation folders, strain sub folders and get the required vasp files.
    """
    def __init__(self,struct,strain_values,strain_vector):
        """
        Constructor function for make_folders class
        :param struct        : Structure of object
               strain_values : list of strain values
               strain_vector : 6 element strain vector for each strain value
        :type  struct        : Pymatgen structure object
               strain_values : list
               strain_vector : list
        :returns calls make_deformation_folders and rename_bulk_folder functions
        """
        self.struct   = struct
        self.strain_values = strain_values
        self.strain_vector = strain_vector
        self.make_deformation_folders()
        self.rename_bulk_folder()
    
    def make_deformation_folders(self):
        """ Function to make folders and subfolders for different deformation-strains and moves INCAR,strained POSCARS, POTCAR and KPOINTS file
                 in each folder/sub folder
        :param None
        returns folders
        """
        for k in range(1,self.num_constants+2):
            for i in range(1,len(self.strain_values)+1):
                e = self.strain_values[i-1]
                strain_arr = smv(self.strain_vector)
                strain = Strain(strain_arr)
                deformation = strain.get_deformation_matrix()
                new_struct = deformation.apply_to_structure(self.struct)
                new_struct.to(fmt='poscar',filename='deformation-{}/strain-{}/POSCAR'.format(k,i))
                shutil.copy('INCAR','deformation-{}/strain-{}/INCAR'.format(k,i))
                shutil.copy('POTCAR','deformation-{}/strain-{}/POTCAR'.format(k,i))
                shutil.copy('KPOINTS','deformation-{}/strain-{}/KPOINTS'.format(k,i))

    def rename_bulk_folder(self):
        """ Function to rename folder as deformation-bulk
        :param None
        returns None
        """
        if os.path.exists('deformation-bulk'):
            shutil.rmtree('deformation-bulk')
        shutil.move('deformation-{}'.format(self.num_constants+1),'deformation-bulk')

    
#--------------------------------------------------------------------#
class cubic_def(parent_def):
    """
    Class to create deformations for cubic structures
    """
    def __init__(self,struct,symmetry,strain_values):
        super().__init__(struct,symmetry,strain_values)
        """
        Constructor function which inherits from parent_def parent class
        :param struct   : Reference structure
               symmetry : Symmetry of reference structure
               strain_values : list of strain values
        :type  struct   : Pymatgen structure object
               strain_values : list
               symmetry  : plaintext file
        """
        self.struct = struct
        self.strain_values = strain_values
        self.symmetry = symmetry
        self.num_constants = get_num_constants(self.symmetry)
    

    def get_strain_array(self,strain_val,k):
        '''
        Function to create unique tailored strain vectors for each symmetry. Calls make folders function from make_folders class

        1: calculates C11+C12; 2: calculates 3/2(C11+2C12)
        3: calculates 3/2 (C44) 4: calculats B
        '''
        e = strain_val
        strain_dict={'1':[e,e,0,0,0,0],'2':[e,e,e,0,0,0],'3':[0,0,0,e,e,e],'4':[e,e,e,0,0,0]}
        strain_vector = strain_dict[str(k)]
        make_folders(self.struct,self.strain_values,strain_vector)
   
#-------------------------------------------------------------------------#
class hexagonal_def(parent_def):
    """
    Class to create deformations for hexagonal structures
    """
    def __init__(self,struct,symmetry,strain_values):
        super().__init__(struct,symmetry,strain_values)
        self.struct = struct
        self.strain_values = strain_values
        self.symmetry = symmetry
        self.num_constants = get_num_constants(self.symmetry)


    def get_strain_array(self,strain_val,k):
        '''
        1: calculates C11+C12 ; 2: calculates 1/4 (C11-C12) ; 3 : calculates 1/2 C33
        4: calculates (C44) ; 5: calculates C11 + C12 +2C13+C33/2 ; 6: calculates B
        '''
        e = strain_val
        strain_dict={'1':[e,e,0,0,0,0],'2':[0,0,0,0,0,e],'3':[0,0,e,0,0,0],'4':[0,0,0,e,e,0],'5':[e,e,e,0,0,0], '6':[e,e,e,0,0,0]}
        strain_vector = strain_dict[str(k)]
        make_folders(self.struct,self.strain_values,strain_vector)

#------------------------------------------------------------------------#

class trigonal1_def(parent_def):
    """
    Class to create deformations for trigonal1 structures
    """
    def __init__(self,struct,symmetry,strain_values):
        super().__init__(struct,symmetry,strain_values)
        self.struct = struct
        self.strain_values = strain_values
        self.symmetry = symmetry
        self.num_constants = get_num_constants(self.symmetry)


    def get_strain_array(self,strain_val,k):
        '''
        1: calculates C11+C12 ; 2: calculates 1/4 (C11-C12) ; 3 : calculates 1/2 C33
        4: calculates (C44) ; 5: calculates C11 + C12 +2C13+C33/2 ; 6: calculates C11/4 - C12/4 + C14+ C44/2
        7: calculates B
        '''
        e = strain_val
        strain_dict={'1':[e,e,0,0,0,0],'2':[0,0,0,0,0,e],'3':[0,0,e,0,0,0],'4':[0,0,0,e,e,0],'5':[e,e,e,0,0,0],'6':[0,0,0,0,e,e],'7':[e,e,e,0,0,0]}
        strain_vector = strain_dict[str(k)]
        make_folders(self.struct,self.strain_values,strain_vector)

#--------------------------------------------------------------------------#


class trigonal2_def(parent_def):
    """
    Class to create deformations for trigonal2 structures
    """
    def __init__(self,struct,symmetry,strain_values):
        super().__init__(struct,symmetry,strain_values)
        self.struct = struct
        self.strain_values = strain_values
        self.symmetry = symmetry
        self.num_constants = get_num_constants(self.symmetry)

    def get_strain_array(self,strain_val,k):
        '''
        1: calculates C11+C12 ; 2: calculates 1/4 (C11-C12) ; 3 : calculates 1/2 C33 ; 4: calculates (C44) 
        5: calculates C11 + C12 +2C13+C33/2 ; 6: calculates C11/4 - C12/4 + C14+ C44/2 ; 7: calculates C11/4 - C12/4 - C15 + C44/2 ; 8: calculates B
        '''
        e = strain_val
        strain_dict={'1':[e,e,0,0,0,0],'2':[0,0,0,0,0,e],'3':[0,0,e,0,0,0],'4':[0,0,0,e,e,0],'5':[e,e,e,0,0,0],'6':[0,0,0,0,e,e],'7':[0,0,0,e,0,e],'8':[e,e,e,0,0,0]}
        strain_vector = strain_dict[str(k)]
        make_folders(self.struct,self.strain_values,strain_vector)
#----------------------------------------------------------------------------#


class tetragonal1_def(parent_def):
    """
    Class to create deformations for tetragonal1 structures
    """
    def __init__(self,struct,symmetry,strain_values):
        super().__init__(struct,symmetry,strain_values)
        self.struct = struct
        self.strain_values = strain_values
        self.symmetry = symmetry
        self.num_constants = get_num_constants(self.symmetry)

    def get_strain_array(self,strain_val,k):
        '''
        1: calculates C11+C12 ; 2: calculates C66/2 ; 3 : calculates 1/2 C33 ; 4: calculates (C44) 
        5: calculates C11 + C12 +2C13+C33/2 ; 6: calculates C11/2 + C13 + C33/2 7: calculates B
        '''
        e = strain_val
        strain_dict={'1':[e,e,0,0,0,0],'2':[0,0,0,0,0,e],'3':[0,0,e,0,0,0],'4':[0,0,0,e,e,0],'5':[e,e,e,0,0,0],'6':[0,e,e,0,0,0],'7':[e,e,e,0,0,0]}
        strain_vector = strain_dict[str(k)]
        make_folders(self.struct,self.strain_values,strain_vector)


#----------------------------------------------------------------------------#


class tetragonal2_def(parent_def):
    """
    Class to create deformations for tetragonal2 structures
    """
    def __init__(self,struct,symmetry,strain_values):
        super().__init__(struct,symmetry,strain_values)
        self.struct = struct
        self.strain_values = strain_values
        self.symmetry = symmetry
        self.num_constants = get_num_constants(self.symmetry)


    def get_strain_array(self,strain_val,k):
        '''
        1: calculates C11+C12 ; 2: calculates C66/2 ; 3 : calculates 1/2 C33 ; 4: calculates (C44) 
        5: calculates C11 + C12 +2C13+C33/2 ; 6: calculates C11/2 + C13 + C33/2 ; 7: calculates C11/2+C16+C66/2 8: calculates B
        '''
        e = strain_val
        strain_dict={'1':[e,e,0,0,0,0],'2':[0,0,0,0,0,e],'3':[0,0,e,0,0,0],'4':[0,0,0,e,e,0],'5':[e,e,e,0,0,0],'6':[0,e,e,0,0,0], '7':[e,0,0,0,0,e],'8':[e,e,e,0,0,0]}
        strain_vector = strain_dict[str(k)]
        make_folders(self.struct,self.strain_values,strain_vector)
        
#----------------------------------------------------------------------------#

class orthorhombic_def(parent_def):
    """
    Class to create deformations for orthorhombic structures
    """
    def __init__(self,struct,symmetry,strain_values):
        super().__init__(struct,symmetry,strain_values)
        self.struct = struct
        self.strain_values = strain_values
        self.symmetry = symmetry
        self.num_constants = get_num_constants(self.symmetry)

    
    def get_strain_array(self,strain_val,k):
        '''
        1: 1/2 C11 ; 2: 1/2 C22 ; 3: 1/2 C33 ; 4: 1/2 C44 ; 5: 1/2 C55 6: 1/2 C66
        7: C11/2 + C12 +C22/2 ; 8: C11/2 +C13 +C33/2 ; 9: C22/2 + C23 + C33/2 10: calculates B
        
        '''
        e = strain_val
        strain_dict={'1':[e,0,0,0,0,0],'2':[0,e,0,0,0,0],'3':[0,0,e,0,0,0],'4':[0,0,0,e,0,0],'5':[0,0,0,0,e,0],
                     '6':[0,0,0,0,0,e],'7':[e,e,0,0,0,0],'8':[e,0,e,0,0,0],'9':[0,e,e,0,0,0],'10':[e,e,e,0,0,0]}
        strain_vector = strain_dict[str(k)]
        make_folders(self.struct,self.strain_values,strain_vector)

#---------------------------------------------------------------------------#

class monoclinic_def(parent_def):
    """
    Class to create deformations for monoclinic structures
    """
    def __init__(self,struct,symmetry,strain_values):
        super().__init__(struct,symmetry,strain_values)
        self.struct = struct
        self.strain_values = strain_values
        self.symmetry = symmetry
        self.num_constants = get_num_constants(self.symmetry)

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
        make_folders(self.struct,self.strain_values,strain_vector)
#-------------------------------------------------------------------------------#
class triclinic_def(parent_def):
    """
    Class to create deformations for triclinic structures
    """
    def __init__(self,struct,symmetry,strain_values):
        super().__init__(struct,symmetry,strain_values)
        self.struct = struct
        self.strain_values = strain_values
        self.symmetry = symmetry
        self.num_constants = get_num_constants(self.symmetry)

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
        make_folders(self.struct,self.strain_values,strain_vector)
#--------------------------------------------------------------------------#


