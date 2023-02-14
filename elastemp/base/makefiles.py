
def make_kpt_dynamic():
    """ Generates a k point file for dynamic calculations if absent
    :param None
    :returns KPOINTS file
    :rtype plaintext file
    """
    kpt_string="Automatic\n0\nGamma\n%2i %2i %2i\n 0. 0. 0.\n" %(12,12,12)
    with open('KPOINTS_dynamic', 'w') as f:
        f.write(kpt_string) 
    f.close()

def make_incar_static():
    """ Generates a INCAR file for static calculations if absent
    :param None
    :returns INCAR file
    :rtype plaintext file
    """
    f = open('INCAR','w')
    f.write('PREC = Accurate\nEDIFF = 1e-6\nISMEAR = -5\nSIGMA = 0.060\nPOTIM = 0.020\nISTART = 0\nNCORE=4\nLCHARG = FALSE\nLWAVE = FALSE\n')
    f.close()

def make_incar_dynamic():
    """ Generates a INCAR file for dynamic calculations if absent
    :param None
    :returns INCAR file
    :rtype plaintext file
    """
    f = open('INCAR_dynamic','w')
    f.write('PREC=Accurate\nIBRION=8\nEDIFF=1e-8\nIALGO=38\nISMEAR=0\nSIGMA=0.06\nLREAL=.FALSE.\nADDGRID=.TRUE.\nLWAVE=.FALSE.\nLCHARG=.FALSE.\n')
    f.close()


def make_kpt_static():
    """ Generates a k point file for static calculations if absent
    :param None
    :returns KPOINTS file
    :rtype plaintext file
    """
    kpt_string="Automatic\n0\nGamma\n%2i %2i %2i\n 0. 0. 0.\n" %(14,14,14)
    f = open('KPOINTS','w')
    f.write(kpt_string)
    f.close()

