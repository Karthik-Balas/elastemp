import numpy as np,os

class compute_response_functions:
    """
    Class to compute the response functions (E,G,K,H) from elastic constants
    """
    def __init__(self,c,symmetry):
        """Constructor function for compute_response_functions class.
        :param c        : elastic constant matrix
               symmetry : Symmetry of structure
        :type  c        : matrix
               symmetry : Plaintext file
        """
        self.c = c
        self.symmetry = symmetry
        
    def get_moduli(self,s,write_moduli):
        """
        Function to get moduli from compliance matrix.
        :param s            : compliance matrix
               write_moduli : Boolean value to see if moduli needs to be written to file. 
        :type  s            : matrix
               write_moduli : Bool(True/False)
        :returns Voigt, Reuss and Hills averages for moduli, pugh and poisson ratio, Chen's hardness
                 Information is written in moduli.txt file
        :rtype   moduli.txt : .txt file
        """
        c = self.c
        Kv = np.round(((c[0][0]+c[1][1]+c[2][2])+2*(c[0][1]+c[1][2]+c[2][0]))/9.0,4)
        Gv = np.round(((c[0][0]+c[0][1]+c[2][2])-(c[0][1]+c[1][2]+c[2][0])+4*(c[3][3]+c[4][4]+c[5][5]))/15.0,4)
        l1 = (s[0][0]+s[1][1]+s[2][2])+2*(s[0][1]+s[1][2]+s[2][0])
        Kr = np.round(1.0/l1,4)
        l2 = 4*(s[0][0]+s[1][1]+s[2][2])-4*(s[0][1]+s[1][2]+s[2][0])+3*(s[3][3]+s[4][4]+s[5][5])
        Gr = np.round(15.0/l2,4)
        Kh = np.round((Kv+Kr)/2,4)
        Gh = np.round((Gv+Gr)/2,4)

        Ev = np.round((9*Kv*Gv)/(3*Kv+Gv),4)
        Er = np.round((9*Kr*Gr)/(3*Kr+Gr),4)
        Eh = np.round((Ev+Er)/2,4)

        nu = np.round((3*Kh-2*Gh)/(6*Kh+Gh),4)
        ductility_variable='Ductile'
        if nu<0.28:
            ductility_variable='Brittle'

        pugh = np.round(Gh/Kh,4)
        pugh_variable = 'Ductile'
        if pugh>0.5:
            pugh_variable = 'Brittle'
        

        
        try:
            hc = np.round(2*(Gh*pugh**2)**0.585 - 3,4)
        except:
            print("Unable to compute hardness via chen's model. Could be because of negative moduli. setting it to zero")
            hc = 0
        if write_moduli==True:
            os.chdir('results_dir')
            f = open('moduli.txt','w')
            f.write('**'*30+'\n')
            f.write('Voigt bulk modulus ={} GPa\n '.format(Kv))
            f.write('Voigt shear modulus ={} GPa\n '.format(Gv))
            f.write('Voigt elastic modulus ={} GPa\n\n '.format(Ev))
            f.write('Reuss bulk modulus ={} GPa\n '.format(Kr))
            f.write('Reuss shear modulus ={} GPa\n '.format(Gr))
            f.write('Reuss elastic modulus ={} GPa\n\n '.format(Er))
            f.write('Hill bulk modulus = {} GPa\n'.format(Kh))
            f.write('Hill shear modulus = {} GPa\n'.format(Gh))
            f.write('Hill elastic modulus = {} GPa\n\n'.format(Eh))
            f.write('----'*30+'\n')
            f.write('poisson ratio = {}.\t. Classified as {} (poisson ratio less than 0.28 classified as brittle, else ductile)\n'.format(nu,ductility_variable))
            f.write("pugh's ratio = {}.\t. Classified as {} (pugh's ratio greater than 0.5 classified as brittle, else ductile)\n".format(pugh,pugh_variable))
            f.write("Hardness from Chen's model = {} GPa\n".format(hc))
            f.write('**'*30+'\n')
            f.close()
            os.chdir('..')
        return Eh, Gh, hc, nu, pugh



    def get_compliance_matrix(self,write_compliance):
        """ Function to get compliance matrix from stiffness matrix.
        :param write_compliance : Boolean value to determine if compliance matrix needs to be written to file.
        :type  write_compliance : Bool (True/False)
        :returns complaince matrix : compliance matrix is written in compliance_matrix.txt
        :rtype   compliance_matrix.txt : .txt file
        """

        det = np.linalg.det(self.c)
        if det==0:
            print('Cannot compute compliance!. Determinant is zero')
        else:
            s = np.linalg.inv(self.c)
        compliance_const_str = " "
        for i in range(6):
            for j in range(6):
                sij = np.round(s[i][j],4)
                compliance_const_str +=str(sij)+'\t'
            compliance_const_str +='\n'
        if write_compliance==True:
            if not os.path.exists('results_dir'):
                os.mkdir('results_dir')
            os.chdir('results_dir')
            f = open('compliance_matrix.txt','w')
            f.write('**'*30+'\n')
            f.write('COMPLIANCE MATRIX\n')
            f.write('**'*30+'\n')
            f.write(compliance_const_str)
            f.close()
            os.chdir('..')
        return s 
