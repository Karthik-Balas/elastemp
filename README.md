# Elastemp

Elastemp is a quasi harmonic approximation to estimate temperature dependent elastic constants. It integrates VASP qm simulation package with Phonopy to set up the calculations and extract information. To use this this automated workflow,clone this repo to your directory. Please have both the elastemp folder and elastic_input.py in the same folder where you want to run calculations. We are working on making this a pip package in subsequent versions. 

# Files needed 

To run vasp calculations, user needs to input
1. POSCAR (mandatory). 
2. POTCAR (mandatory).
3. KPOINTS (optional) - a default KPOINTS file will be provided. But for best results, users need to tune it depending on the structure.
4. INCAR   (optional) - a default INCAR file will be provided. But for best results, users need to tune it depending on the structure.
5. strains.txt (optional) - a default strains.txt file will be provided with strain values -0.05, -0.03, -0.01, 0, 0.01, 0.03, 0.05. 
6. dim.txt (optional) - a default dim.txt file with default values 2 2 2 will be provided. This file contains the dimensions needed to make supercells                           for phonon calculations. For best results, users need to tune it depending on the structure.
7. tmax.txt (optional) - a default tmax.txt file with defualt value of 1000 will be provided. This file contains the maximum temperature at which elastic                           constant calculations are to be performed.

# Workflow

<p align="center">
<img src="https://user-images.githubusercontent.com/120595580/207714085-196181cf-5f77-46b5-9c53-b64c73da9e68.png" width="600" height="600">
</p>

The workflow consists of three modules -  Getting zero temp elastic constants, Getting thermal expansion coeff and Getting elastic constants at temp T. 

### Steps before starting calculations.
1. clone repo and elastemp_input.py to your directory. Have POSCAR and POTCAR in the same directory. INCAR and KPOINTS are optional. A default one will     be used otherwise.
2. A file strains.txt containing strain values as a column. However a default set of strains will be used otherwise and they are accurate for most cases.

## Getting zero temp. elastic constants (Steps 1 & 2)

### Step 1: This creates the deformation folders to run vasp calculations

python command : python3 elastic_input.py --operations make_elastic

1. This will recognize the symmetry of the structure and make deformation-1, deformation-2,.... deformation-n and deformation-bulk folders. n is the      number of independent elastic constants.
2. User needs to run vasp in these folders.

### Step 2: This gets the results from vasp calculations for zero temp.elastic constants

python command : python3 elastic_input.py --operations get_elastic

1. Elastemp extracts enthalpies from folders and outputs results to results_dir. The files in results_dir are<br/>
   a. In each deformation folder, a plot of energy vs strain, fitting function and curvature constant for fitting.\
   b. In results_dir the following files will be created\
      (i)   stiffness_matrix.txt - Stiffness matrix will be created\
      (ii)  compliance_matrix.txt - Compliance matrix will be created\
      (iii) moduli.txt            - A text file with different moduli (K, E, G, H, nu, pugh) values created\
      (iv)  stability_criteria.txt - A file which checks the mechanical stability criteria of materials. 
 
