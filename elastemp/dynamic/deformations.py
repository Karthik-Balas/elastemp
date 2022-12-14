import os,shutil
import sys
from elastemp.base.symmetry import get_num_constants

class dynamic_deformations:
    """ Class to create deformations for phonon calculations
    """
    def __init__(self,symmetry,strain_values,dim,struct):
        """ Constructor function for dynamic deformations class
        :param symmetry      : symmetry of the structure
               strain_values : list of strain values
               dim           : dimension of the supercell
               struct        : Reference structure
        :type  symmetry      : plain text file
               strain_values : list
               dim           : list
               struct        : Pymatgen structure object
        """
        self.num_constants = get_num_constants(symmetry)
        self.strain_values = strain_values
        self.dim = dim
        self.struct = struct
    
    def make_phonon_deformation(self):
        """ Function to make deformations for phonon calculations using Phonopy for bulk deformation
        :param None
        :returns None. Saves supercells of structures for phonon calculations in POSCAR format
        """
        dim0,dim1,dim2 = self.dim
        if os.path.exists('Phonon_calculation'):
            shutil.rmtree('Phonon_calculation')
        os.mkdir('Phonon_calculation')
        os.chdir('Phonon_calculation')
        if not os.path.exists('../ZPE.txt'):
            print('ZPE file does not exist!. Looks like static calculation is not converged.! Dynamic calculation will not go forward')
            sys.exit()
        else:
            shutil.copy('../ZPE.txt','ZPE.txt')
        if not os.path.exists('../POSCAR'):
            print('POSCAR file does not exist!. Calculation will not go forward')
            sys.exit()
        else:
            shutil.copy('../POSCAR','POSCAR')
        
        os.system("phonopy -d --dim='{} {} {}'".format(dim0,dim1,dim2))
                                    
        shutil.copy('POSCAR','POSCAR-unitcell')
        shutil.move('SPOSCAR','POSCAR')


    def make_phonon_deformation_T(self):
        """ Function to make deformations for phonon calculations using Phonopy for different deformations
        :param None
        :returns None. Saves supercells of structures for phonon calculations in POSCAR format
        """
        dim0,dim1,dim2 = self.dim
        if os.path.exists('Phonon_calculation'):
            shutil.rmtree('Phonon_calculation')
        os.mkdir('Phonon_calculation')
        os.chdir('Phonon_calculation')



        if not os.path.exists('../Hightemp_static/ZPE.txt'):
            print('ZPE file does not exist!. Looks like High temperature static calculation is not converged.! Dynamic calculation will not go forward')
            sys.exit()
        else:
            shutil.copy('../Hightemp_static/ZPE.txt','ZPE.txt')
        
        if not os.path.exists('../Hightemp_static/POSCAR'):
            print('Volume corrected POSCAR file does not exist in Hightemp_static folder!. Calculation will not go forward')
            sys.exit()
        else:
            shutil.copy('../Hightemp_static/POSCAR','POSCAR')

        os.system("phonopy -d --dim='{} {} {}'".format(dim0,dim1,dim2))
                                    
        shutil.copy('POSCAR','POSCAR-unitcell')
        shutil.move('SPOSCAR','POSCAR')

    def make_bulk_deformation(self):
        """ 
        Function to get files needed to run vasp DFPT calculations for bulk deformation.
        """
        os.chdir('deformation-bulk')
        strain_array = self.strain_values
        for i in range(1,len(strain_array)+1):
            os.chdir('strain-{}'.format(i))
            print(os.getcwd())
            self.make_phonon_deformation()
            shutil.copy('../../../KPOINTS_dynamic','KPOINTS')
            shutil.copy('../../../INCAR_dynamic','INCAR')
            shutil.copy('../../../INCAR','INCAR_static')
            shutil.copy('../../../POTCAR','POTCAR')
            os.chdir('../..')
        os.chdir('..')
        

    def make_deformation(self):
        """ 
        Function to get files needed to run vasp DFPT calculations for the different deformations.
        """
        for k in range(1,self.num_constants+1):
            print(os.getcwd())
            os.chdir('deformation-{}'.format(k))
            for i in range(1,len(self.strain_values)+1):
                os.chdir('strain-{}'.format(i))
                self.make_phonon_deformation_T()
                print(os.getcwd())
                shutil.copy('../../../KPOINTS_dynamic','KPOINTS')
                shutil.copy('../../../INCAR_dynamic','INCAR')
                shutil.copy('../../../INCAR','INCAR_static')
                shutil.copy('../../../POTCAR','POTCAR')
                os.chdir('../..')
            os.chdir('..')
     
    
class dynamic_response:
    """
    Class to get force constants from phonon calculations.
    """
    def __init__(self):
        """ Constructor function which calls the extract phonon function.
        """
        self.extract_phonon()

    def extract_phonon(self):
        """ Function which extracts force constants from vasprun.xml files.
        """
        print('Extracting Force constants:\n')
        if not os.path.exists('vasprun.xml'):
            print('Current directory: ',os.getcwd())
            print('Vasprun.xml file is missing!. Program will stop')
            sys.exit()
        os.system('phonopy --fc vasprun.xml')




    
