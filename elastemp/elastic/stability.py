import numpy as np,os


class stability_check:
    """ 
    Class to check mechanical stability of materials 
    """
    def __init__(self,c,symmetry):
        """ Constructor function for stability check class
        :param c        : elastic constants matrix
               symmetry : symmetry of structure
        :type  c        : matrix
               symmetry : plaintext file
        :returns None. calls check_stability function
        """
        self.c = c
        self.symmetry = symmetry
        self.check_stability()

    def check_stability(self):
        """
        Function to check mechanical stability of materials.
        :param None
        :returns mechanical_stability_criteria.txt : text file which mentions stability criteria and if they have been met.
        :rtype   mechanical_stability_criteria.txt : .txt file
        """
        print('Checking mechanical stability. Results will be in results_dir')
        c = self.c
        symmetry = self.symmetry
        if not os.path.exists('results_dir'):
            os.mkdir('results_dir')
        os.chdir('results_dir')
        f = open('mechanical_stability_criteria.txt','w')
        f.write('**'*30+'\n')
        f.write('Detected symmetry of structure = {} \n\n'.format(symmetry))
        f.write('Criteria for mechanical stability of {} structures : \n\n'.format(symmetry))
        
        if symmetry=='cubic':
            condition1 = 'c11-c12>0'
            condition2 = 'c11+2c12>0'
            condition3 = 'c44>0'
            f.write('conditon 1 = {}\n\n'.format(condition1))
            f.write('conditon 2 = {}\n\n'.format(condition2))
            f.write('conditon 3 = {}\n\n'.format(condition3))
            f.write('---'*30+'\n')
            if c[0][0]>c[0][1]:
                f.write('condition 1 constraint:\t  met\n')
            else:
                f.write('condition 1 constraint: \t NOT met \n')
            if c[0][0]+2*c[0][1]>0:
                f.write('condition 2 constraint:\t met\n')
            else:
                f.write('condition 2 constraint:\t NOT met \n')
            if c[3][3]>0:
                f.write('condition 3 constraint:\t met \n')
            else:
                f.write('condition 3 constraint:\t NOT met\n')

        elif symmetry=='hexagonal' or symmetry=='tetragonal1' or symmetry=='tetragonal2':
            condition1 = 'c11>|c12|'
            condition2 = '2C13^2 < C33*(C11+C12)'
            condition3 = 'C44>0'
            if symmetry=='tetragonal2':
                condition4 =='2C16^2 < C66*(C11-C12)'
            else:
                condition4 = 'C66>0'
            f.write('condition 1 = {}\n\n'.format(condition1))
            f.write('condition 2 = {}\n\n'.format(condition2))
            f.write('condition 3 = {}\n\n'.format(condition3))
            f.write('condition 4 = {}\n\n'.format(condition4))
            f.write('---'*30+'\n')
            if c[0][0]>abs(c[0][1]):
                f.write('condition 1 constraint :\t met\n')
            else:
                f.write('condition 1 constraint :\t NOT met\n')
            if 2*c[0][2]**2 < c[2][2]*(c[0][0]+c[0][1]):
                f.write('condition 2 constraint : \t met\n')
            else:
                f.write('condition 2 constraint: \t NOT met\n')
            if c[3][3]>0:
                f.write('condition 3 constraint :\t met \n')
            else:
                f.write('condition 3 constraint :\t NOT met\n')
            if symmetry =='tetragonal2':
                if 2*c[0][5]**2 < c[5][5]*(c[0][0]-c[0][1]):
                    f.write('condition 4 constraint :\t met\n')
                else:
                    f.write('condition 4 constraint :\t NOT met\n')
            else:
                if c[5][5]>0:
                    f.write('condition 4 constraint :\t met \n')
                else:
                    f.write('condition 4 constraint:\t NOT met\n')

        elif symmetry=='trigonal1' or symmetry=='trigonal2':
            
            condition1 = 'c11>|c12|'
            condition2 = 'c44>0'
            condition3 = 'c13^2 < 1/2*(c33*(c11+c12))'
            if symmetry=='trigonal1':
                condition4 = 'C14^2  <1/2*(c44*(c11-c12))'
            elif symmetry=='trigonal2':
                condition4 = 'C14^2 +C15^2 <1/2*(c44*(c11-c12))'                        
            f.write('condition 1 = {}\n\n'.format(condition1))
            f.write('condition 2 = {}\n\n'.format(condition2))
            f.write('condition 3 = {}\n\n'.format(condition3))
            f.write('condition 4 = {}\n\n'.format(condition4))

            if c[0][0]>abs(c[0][1]):
                f.write('condition 1 constraint :\t met \n')
            else:
                f.write('condition 1 constraint :\t NOT met\n')
            if c[3][3]>0:
                f.write('condition 2  constraint:\t met \n')
            else:
                f.write('condition 2  constraint:\t NOT met\n')
            if symmetry=='trigonal1':
                if c[0][3]**2 < 0.5* (c[3][3]*(c[0][0]-c[0][1])):
                    f.write('condition 3 constraint:\t met\n')
                else:
                    f.write('condition 3 constraint:\t NOT met\n')
            else:
                if c[0][3]**2 + c[0][4]**2 < 0.5*(c[3][3]*(c[0][0]-c[0][1])):
                    f.write('condition 4 constraint:\t met\n')
                else:
                    f.write('condition 4 constraint:\t NOT met\n')
       
        elif  symmetry=='orthorhombic':
            condition1 = 'c11>0'
            condition2 = 'c11*c22>c12^2'
            condition3 = 'c11*c22*c33+2*c12*c13*c23-c11*c23^2-c22*c13^2-c33*c12^2>0'
            condition4 = 'c44>0'
            condition5 = 'c55>0'
            condition6 = 'c66>0'
            
            f.write('condition 1 = {}\n\n'.format(condition1))
            f.write('condition 2 = {}\n\n'.format(condition2))
            f.write('condition 3 = {}\n\n'.format(condition3))
            f.write('condition 4 = {}\n\n'.format(condition4))
            f.write('condition 5 = {}\n\n'.format(condition5))
            f.write('condition 6 = {}\n\n'.format(condition6))
            
            if c[0][0]>0:
                f.write('condition 1 constraint :\t met \n')
            else:
                f.write('condition 1 constraint: \t NOT met\n')
            if c[0][0]*c[1][1]>c[0][1]**2:
                f.write('condition 2 constraint: \t met\n')
            else:
                f.write('condition 2 constraint: \t NOT met\n')
            if c[0][0]*c[1][1]*c[2][2]+2*c[0][1]*c[0][2]*c[1][2]-c[0][0]*c[1][2]**2-c[1][1]*c[0][2]**2-c[2][2]*c[0][1]**2>0:
                f.write('condition 3 constraint:\t met\n')
            else:
                f.write('condition 3 constraint:\t NOT met\n')
            if c[3][3]>0:
                f.write('condition 4 constraint :\t met \n')
            else:
                f.write('condition 4 constraint:\t NOT met\n')
            if c[4][4]>0:
                f.write('condition 5 constraint: \t met \n')
            else:
                f.write('condition 5 constraint:\t NOT met\n')
            if c[5][5]>0:
                f.write('condition 6 constraint:\t met \n')
            else:
                f.write('condition 6  constraint:\t NOT met\n')

        elif symmetry=='monoclinic':
            condition1 = 'c11 >0'
            condition2 = 'c22 >0'
            condition3 = 'c33 >0'
            condition4 = 'c44 >0'
            condition5 = 'c55 >0'
            condition6 = 'c66 >0'
            condition7 = 'c11+c22+c33+2*(c12+c13+c23) >0'
            condition8 = 'c33*c55 - c35^2 >0'
            condition9 = 'c44*c66-c46^2 >0'
            condition10 = 'c22+c33-2c23 >0'
            condition11 = 'c22*(c33*c55-c35^2)+2(c23*c25*c35 - c23^2*c55 - c25^2*c33) >0'
            condition12 = '2*(c15*c25*(c33c12 - c13c23)+c15*c35*(c22c13-c12c23)+c25*c35(c11c23 - c12c13)) - (c15^2(c22c33 - c23^2)+c25^2(c11c33 - c13^2) + c35^2(c11c22-c12^2))+c55*(c11c22c33-c11c23^2-c22c13^2 - c33c12^2+2c12c13c23) >0'
            

            f.write('condition 1 = {}\n\n'.format(condition1))
            f.write('condition 2 = {}\n\n'.format(condition2))
            f.write('condition 3 = {}\n\n'.format(condition3))
            f.write('condition 4 = {}\n\n'.format(condition4))
            f.write('condition 5 = {}\n\n'.format(condition5))
            f.write('condition 6 = {}\n\n'.format(condition6))
            f.write('condition 7 = {}\n\n'.format(condition7))
            f.write('condition 8 = {}\n\n'.format(condition8))
            f.write('condition 9 = {}\n\n'.format(condition9))
            f.write('condition 10 = {}\n\n'.format(condition10))
            f.write('condition 11 = {}\n\n'.format(condition11))
            f.write('condition 12 = {}\n\n'.format(condition12))
            
            if c[0][0]>0:
                f.write('condition 1 constraint :\t met \n')
            else:
                f.write('condition 1 constraint: \t NOT met\n')
            if c[1][1]>0:
                f.write('condition 2 constraint :\t met \n')
            else:
                f.write('condition 2 constraint: \t NOT met\n')
            if c[2][2]>0:
                f.write('condition 3 constraint :\t met \n')
            else:
                f.write('condition 3 constraint: \t NOT met\n')
            if c[3][3]>0:
                f.write('condition 4 constraint :\t met \n')
            else:
                f.write('condition 4 constraint: \t NOT met\n')
            if c[4][4]>0:
                f.write('condition 5 constraint :\t met \n')
            else:
                f.write('condition 5 constraint: \t NOT met\n')
            if c[5][5]>0:
                f.write('condition 6 constraint :\t met \n')
            else:
                f.write('condition 6 constraint: \t NOT met\n')
            if c[0][0]+c[1][1]+c[2][2]+2*(c[0][1]+c[0][2]+c[1][2])>0:
                f.write('condition 7 constraint :\t met \n')
            else:
                f.write('condition 7 constraint: \t NOT met\n')
            if c[2][2]*c[4][4]-c[2][4]**2>0:
                f.write('condition 8 constraint :\t met \n')
            else:
                f.write('condition 8 constraint: \t NOT met\n')
            if c[3][3]*c[5][5]-c[3][5]**2>0:
                f.write('condition 9 constraint :\t met \n')
            else:
                f.write('condition 9 constraint: \t NOT met\n')
            if c[1][1]+c[2][2]-2*c[1][2]>0:
                f.write('condition 10 constraint :\t met \n')
            else:
                f.write('condition 10 constraint: \t NOT met\n')
            if c[1][1]*(c[2][2]*c[4][4]-c[2][4]**2)+2*(c[1][2]*c[1][4]*c[2][4]-c[1][2]**2*c[4][4]-c[1][4]**2*c[2][2])>0:
                f.write('condition 11 constraint :\t met \n')
            else:
                f.write('condition 11 constraint: \t NOT met\n')
            g1 = c[0][4]*c[1][4]*(c[2][2]*c[0][1]-c[0][2]*c[1][2])
            g2 = c[0][4]*c[2][4]*(c[1][1]*c[0][2]-c[0][1]*c[1][2])
            g3 = c[1][4]*c[2][4]*(c[0][0]*c[1][2]-c[0][1]*c[0][2])
            h1 = c[0][4]**2*(c[1][1]*c[2][2]-c[1][2]**2)
            h2 = c[1][4]**2*(c[0][0]*c[2][2]-c[0][2]**2)
            h3 = c[2][4]**2*(c[0][0]*c[1][1]-c[0][1]**2)
            k1 = c[0][0]*c[1][1]*c[2][2] - c[0][0]*c[1][2]*c[1][2] - c[1][1]*c[0][2]*c[0][2] - c[2][2]*c[0][1]*c[0][1] + 2*c[0][1]*c[0][2]*c[1][2]
            if 2*(g1+g2+g3)-(h1+h2+h3)+c[4][4]*k1>0:
                f.write('condition 12 constraint :\t met \n')
            else:
                f.write('condition 12 constraint: \t NOT met\n')
        
        elif symmetry=='triclinic':
            condition = 'Stability criteria for triclinic is too complex to solve. Determinant greater than zero is reasonably sufficient'
            
            f.write('condition  = {}\n\n'.format(condition))
            determinant = np.linalg.det(c)
            if determinant >0 :
                f.write(' Basic stability criteria for triclinic structure:\t met\n')
            else:
                f.write('Basic stability criteria for triclinic structure:\t NOT met\n')



        f.write('**'*30 +'\n')
        f.close()
        os.chdir('..')




