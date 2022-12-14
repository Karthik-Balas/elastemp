from pymatgen.core import Structure
import numpy as np,os
import pandas as pd
import matplotlib.pyplot as plt
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer as sga
from elastemp.base.symmetry import get_symmetry
from elastemp.base.strain_analysis import get_strain_matrix, get_curvature


class compute_elastic_constants:
    def __init__(self,num_constants,strain_array,symmetry):
        self.num_constants = num_constants
        self.strain_array = strain_array
        self.symmetry = symmetry
        self.get_bulk_curvatures()

    def get_curvatures(self):
        num_constants = self.num_constants
        curvatures = []
        for i in range(1,num_constants+1):
            os.chdir('deformation-{}'.format(i))
            curvature = get_curvature(self.strain_array,calc_bulk=False)
            curvatures.append(curvature)
            os.chdir('..')
        return curvatures

    def get_bulk_curvatures(self):
        os.chdir('deformation-bulk')
        mod = get_curvature(self.strain_array,calc_bulk=True)
        os.chdir('..')


    def get_elastic_constants(self,curvatures,write_elastic):
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
            c[0][2] = (curvatures[4]-c[0][0]-c[0][1]-(c[2][2]/2))
            c[0][3] = curvatures[5]-(c[0][0]/4)+(c[0][1]/4)-(c[3][3]/2)
            c[1][1] = c[0][0]
            c[1][2] = c[0][2]
            c[1][3] = -c[0][3]
            c[4][5] = c[0][3]
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
        f= open('results_dir/elastic_constants_rawdata_0K.txt','w')
        for i in range(6):
            for j in range(6):
                f.write('c_{}{}='.format(i+1,j+1)+'\t'+str(c[i][j])+'\n')
        f.close()


    def get_plots(self):
        for k in range(1,self.num_constants+1):
            os.chdir('deformation-{}'.format(k))
            df = pd.read_csv('energies.txt',sep='\t')
            energies = df['energies'].values
            strains = df['strains'].values
                    
            fit = np.polyfit(strains,energies,2)
            
            st_min = min(strains)
            st_max = max(strains)
            st_vals = [st_min + (st_max-st_min)*i/100 for i in range(101)]
            e_fit = [fit[0]*st_vals[k]**2 + fit[1]*st_vals[k] + fit[2] for k in range(len(st_vals))]
            fig = plt.figure()
            plt.plot(strains,energies,'ro',label='energies')
            plt.plot(st_vals,e_fit,'r--',label='fit')            
            plt.legend()
            plt.xlabel('strain')
            plt.ylabel('Energy (eV)')
            plt.title('Deformation-{}-strain-curve'.format(k))
            plt.savefig('deformation-{}-strain-curve'.format(k))
            plt.close()
            os.chdir('..')

