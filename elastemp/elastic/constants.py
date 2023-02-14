import numpy as np,os
import pandas as pd
from elastemp.base.strain_analysis import get_strain_matrix, get_curvature,plot_parabola


class compute_elastic_constants:
    """
    Class to compute elastic constants from energy-strain values.
    """
    def __init__(self,num_constants,strain_array,symmetry):
        """ Constructor function for compute_elastic_constants class.
        :param num_constants : number of independent elastic constants 
               strain_array  : list of strain values
               symmetry      : symmetry of the structure
        :type  num_constants : int 
               strain_array  : list 
               symmetry      : plaintext file
        """       
        self.num_constants = num_constants
        self.strain_array = strain_array
        self.symmetry = symmetry
        self.get_bulk_curvatures()

    def get_curvatures(self):
        """ Function to get curvatures of parabolic fitting for each deformation.
        :param None
        :return curvatures of parabolic fit.
        :rtype  list
        """
        num_constants = self.num_constants
        curvatures = []
        for i in range(1,num_constants+1):
            os.chdir('deformation-{}'.format(i))
            curvature = get_curvature(self.strain_array,calc_bulk=False)
            curvatures.append(curvature)
            os.chdir('..')
        return curvatures

    def get_bulk_curvatures(self):
        """ Function to get curvatures of parabolic fitting for bulk deformation
        :param None
        :returns None. get_curvature function writes bulk modulus to file when calc_bulk is True.
        """
        os.chdir('deformation-bulk')
        mod = get_curvature(self.strain_array,calc_bulk=True)
        os.chdir('..')


    def get_elastic_constants(self,curvatures,write_elastic):
        """ Function to get elastic constants from curvatures 
        :param  curvatures    : list of curvatures obtained from fitting energy-strain curves for each deformation
                write_elastic : Boolean value (True/False) to determine if elastic constants needs to be written to file
        :type   curvatures    : list
                write_elastic : Bool (True/False)
        :returns stiffness matrix which is written to the stiffness_matrix.txt file
        :rtype   stiffness_matrix.txt : .txt file
        """
        c = np.zeros((6,6))

        def mirror(c):
            for i in range(6):
                for j in range(6):
                    if i!=j:
                        c[j][i] = c[i][j]
            return c
        
        if self.symmetry.lower()=='cubic':
            c[0][1] = (2/3)*curvatures[1] - curvatures[0]
            c[0][0] = curvatures[0] - c[0][1]
            c[3][3] = (2/3)*curvatures[2]  
            c[4][4] = c[3][3]
            c[5][5] = c[3][3]
            c[0][2] = c[0][1]
            c[1][2] = c[0][1]
            c[1][1] = c[0][0]
            c[2][2] = c[0][0]
            c=mirror(c)
        
        elif self.symmetry.lower()=='hexagonal':
            c[0][1] = (curvatures[0]-4*curvatures[1])/2
            c[0][0] = curvatures[0]-c[0][1]
            c[2][2] = 2*curvatures[2]
            c[3][3] = curvatures[3]
            c[0][2] = (curvatures[4]-c[0][0]-c[0][1]-(c[2][2]/2))/2
            c[1][1] = c[0][0]
            c[4][4] = c[3][3]
            c[5][5] = (c[0][0]-c[0][1])/2
            c[1][2] = c[0][2]
            c = mirror(c)
        
        elif self.symmetry.lower()=='trigonal1':
            c[0][1] = (curvatures[0]-4*curvatures[1])/2
            c[0][0] = curvatures[0]-c[0][1]
            c[2][2] = 2*curvatures[2]
            c[3][3] = curvatures[3]
            c[0][2] = (curvatures[4]-c[0][0]-c[0][1]-(c[2][2]/2))/2
            c[0][3] = curvatures[5]-(c[0][0]/4)+(c[0][1]/4)-(c[3][3]/2)
            c[1][1] = c[0][0]
            c[1][2] = c[0][2]
            c[1][3] = -c[0][3]
            c[4][5] = c[0][3]
            c[4][4] = c[3][3]
            c[5][5] = (c[0][0]-c[0][1])/2
            c = mirror(c)
        
        elif self.symmetry.lower()=='trigonal2':
            c[0][1] = (curvatures[0]-4*curvatures[1])/2
            c[0][0] = curvatures[0]-c[0][1]
            c[2][2] = 2*curvatures[2]
            c[3][3] = curvatures[3]
            c[0][2] = (curvatures[4]-c[0][0]-c[0][1]-(c[2][2]/2))/2
            c[0][3] = curvatures[5]-(c[0][0]/4)+(c[0][1]/4)-(c[3][3]/2)
            c[0][4] = -1*(curvatures[6]-(c[0][0]/4)+(c[0][1]/4)-(c[3][3]/2))
            c[1][1] = c[0][0]
            c[1][2] = c[0][2]
            c[1][3] = -c[0][3]
            c[1][4] = -c[0][4]
            c[3][5] = -c[0][4]
            c[4][4] = c[3][3]
            c[4][5] = c[0][3]
            c[5][5] = (c[0][0]-c[0][1])/2
            c = mirror(c)
        
        elif self.symmetry.lower()=='tetragonal1':
            c[5][5] = curvatures[1]*2
            c[2][2] = curvatures[2]*2
            c[3][3] = curvatures[3]
            c[0][2] = (curvatures[4]-curvatures[0]-0.5*c[2][2])/2
            c[0][0] = (curvatures[5]-c[0][2]-0.5*c[2][2])*2
            c[0][1] = curvatures[0]-c[0][0]
            c[1][1] = c[0][0]
            c[1][2] = c[0][2]
            c[4][4] = c[3][3]
            c = mirror(c)
        elif self.symmetry.lower()=='tetragonal2':
            c[5][5] = curvatures[1]*2
            c[2][2] = curvatures[2]*2
            c[3][3] = curvatures[3]
            c[0][2] = (curvatures[4]-curvatures[0]-0.5*c[2][2])/2
            c[0][0] = (curvatures[5]-c[0][2]-0.5*c[2][2])*2
            c[0][1] = curvatures[0]-c[0][0]
            c[0][5] = curvatures[6]-0.5*(c[0][0]+c[5][5])
            c[1][1] = c[0][0]
            c[1][2] = c[0][2]
            c[1][5] = -c[0][5]
            c[4][4] = c[3][3]
            c = mirror(c)
        elif self.symmetry.lower()=='orthorhombic':
            c[0][0] = curvatures[0]*2
            c[1][1] = curvatures[1]*2
            c[2][2] = curvatures[2]*2
            c[3][3] = curvatures[3]*2
            c[4][4] = curvatures[4]*2
            c[5][5] = curvatures[5]*2
            c[0][1] = curvatures[6] - 0.5*(c[0][0]+c[1][1])
            c[0][2] = curvatures[7] - 0.5*(c[0][0]+c[2][2])
            c[1][2] = curvatures[8] - 0.5*(c[1][1]+c[2][2])
            c = mirror(c)
        elif self.symmetry.lower()=='monoclinic':
            c[0][0] = curvatures[0]*2
            c[1][1] = curvatures[1]*2
            c[2][2] = curvatures[2]*2
            c[3][3] = curvatures[3]*2
            c[4][4] = curvatures[4]*2
            c[5][5] = curvatures[5]*2
            c[0][1] = curvatures[6] - 0.5*(c[0][0]+c[1][1])
            c[0][2] = curvatures[7] - 0.5*(c[0][0]+c[2][2])
            c[0][4] = curvatures[8] - 0.5*(c[0][0]+c[4][4])
            c[1][2] = curvatures[9] - 0.5*(c[1][1]+c[2][2])
            c[1][4] = curvatures[10] - 0.5*(c[1][1]+c[4][4])
            c[2][4] = curvatures[11] - 0.5*(c[2][2]+c[4][4])
            c[3][5] = curvatures[12] - 0.5*(c[3][3]+c[5][5])
            c = mirror(c)
        elif self.symmetry.lower()=='triclinic':
            c[0][0] = curvatures[0]*2
            c[1][1] = curvatures[1]*2
            c[2][2] = curvatures[2]*2
            c[3][3] = curvatures[3]*2
            c[4][4] = curvatures[4]*2
            c[5][5] = curvatures[5]*2
            c[0][1] = curvatures[6] - 0.5*(c[0][0]+c[1][1])
            c[0][2] = curvatures[7] - 0.5*(c[0][0]+c[2][2])
            c[0][3] = curvatures[8] - 0.5*(c[0][0]+c[3][3])
            c[0][4] = curvatures[9] - 0.5*(c[0][0]+c[4][4])
            c[0][5] = curvatures[10] - 0.5*(c[0][0]+c[5][5])
            c[1][2] = curvatures[11] - 0.5*(c[1][1]+c[2][2])
            c[1][3] = curvatures[12] - 0.5*(c[1][1]+c[3][3])
            c[1][4] = curvatures[13] - 0.5*(c[1][1]+c[4][4])
            c[1][5] = curvatures[14] - 0.5*(c[1][1]+c[5][5])
            c[2][3] = curvatures[15] - 0.5*(c[2][2]+c[3][3])
            c[2][4] = curvatures[16] - 0.5*(c[2][2]+c[4][4])
            c[2][5] = curvatures[17] - 0.5*(c[2][2]+c[5][5])
            c[3][4] = curvatures[18] - 0.5*(c[3][3]+c[4][4])
            c[3][5] = curvatures[19] - 0.5*(c[3][3]+c[5][5])
            c[4][5] = curvatures[20] - 0.5*(c[4][4]+c[5][5])
            c = mirror(c) 
        
        if write_elastic ==True:
            if not os.path.exists('results_dir'):
                os.mkdir('results_dir')
            os.chdir('results_dir')

            elast_const_str = " "
            for i in range(6):
                for j in range(6):
                    c[i][j]=np.round(float(c[i][j]),3)
                    c[j][i] = c[i][j]
                    elast_const_str +=str(c[i][j])+'\t'
                elast_const_str +='\n'

            f = open('stiffness_matrix.txt','w')
            f.write('**'*30+'\n')
            f.write('STIFFNESS MATRIX\n')
            f.write('**'*30+'\n')
            f.write(elast_const_str)
            f.close()

            os.chdir('..') 
        return c

    def write_zero_elastic_constants(self,c):
        """ Function to write zero temperature elastic constants to file
        :param  c   : stiffness matrix
        :type   c   : matrix
        :returns elastic_constants_rawdata_0K.txt : file which dumps elastic constants at 0 K to file.
        :rtype   elastic_constants_rawdata_0K.txt : .txt file
        """

        f= open('results_dir/elastic_constants_rawdata_0K.txt','w')
        for i in range(6):
            for j in range(6):
                f.write('c_{}{}='.format(i+1,j+1)+'\t'+str(c[i][j])+'\n')
        f.close()


    def get_plots(self):
        """ Function to plot energy-strain curves for each deformation.
        :param None
        :return Deformation-strain-curve : plot of energy-strain for each deformation
        :rtype  Deformation-strain-curve.png : .png file
        """
        for k in range(1,self.num_constants+1):
            os.chdir('deformation-{}'.format(k))
            df = pd.read_csv('energies.txt',sep='\t')
            energies = df['energies'].values
            strains = df['strains'].values
            
            plot_parabola(strains,energies,'Deformation-{}-strain-curve'.format(k),per_atom=False)
            os.chdir('..')

