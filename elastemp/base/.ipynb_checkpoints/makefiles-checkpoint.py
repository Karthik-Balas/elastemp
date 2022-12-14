import numpy as np
import pandas as pd, os
from pymatgen.core import Structure
from pymatgen.analysis.eos import EOS
import shutil
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer as sga


def make_kpt_dynamic():
    kpt_string="kpt_density=%2i\n0\nGamma\n%2i %2i %2i\n 0. 0. 0.\n" %(3000,12,12,12)
    with open('KPOINTS_dynamic', 'w') as f:
        f.write(kpt_string) 
    f.close()

def make_incar_static():
    f = open('INCAR','w')
    f.write('PREC = Accurate\nEDIFF = 1e-6\nISMEAR = -5\nSIGMA = 0.060\nPOTIM = 0.020\nISTART = 0\nNCORE=4\nLCHARG = FALSE\nLWAVE = FALSE\n')
    f.close()

def make_incar_dynamic():
    f = open('INCAR_dynamic','w')
    f.write('PREC=Accurate\nIBRION=8\nEDIFF=1e-8\nIALGO=38\nISMEAR=0\nSIGMA=0.06\nLREAL=.FALSE.\nADDGRID=.TRUE.\nLWAVE=.FALSE.\nLCHARG=.FALSE.\n')
    f.close()


def make_kpt_static():
     kpt_string="kpt_density=%2i\n0\nGamma\n%2i %2i %2i\n 0. 0. 0.\n" %(1000,12,12,12)
     f = open('KPOINTS','w')
     f.write(kpt_string)
     f.close()

