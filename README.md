# Elastemp

Elastemp is a quasi harmonic approximation to estimate temperature dependent elastic constants. It integrates VASP qm simulation package with Phonopy to set up the calculations and extract information. To use this this automated workflow,clone this repo to your directory. Please have both the elastemp folder and elastic_input.py in the same folder where you want to run calculations. We are working on making this a pip package in subsequent versions. 

# Files needed 

To run vasp calculations, user needs to input
1. POSCAR (mandatory). 
2. POTCAR (mandatory).
3. KPOINTS (optional) - a default KPOINTS file will be provided. But for best results, users need to tune it depending on the structure.
4. INCAR   (optional) - a default INCAR file will be provided. But for best results, users need to tune it depending on the structure.
5. strain_values.txt (optional) - a default strain_values.txt file will be provided with strain values -0.05, -0.03, -0.01, 0, 0.01, 0.03, 0.05. 
6. dim.txt (optional) - a default dim.txt file with default values 2 2 2 will be provided. This file contains the dimensions needed to make supercells for                         phonon calculations. For best results, users need to tune it depending on the structure.
7. tmax.txt (optional) - a default tmax.txt file with defualt value of 1000 will be provided. This file contains the maximum temperature at which elastic                           constant calculations are to be performed.

# Workflow


<img src="https://user-images.githubusercontent.com/120595580/207714085-196181cf-5f77-46b5-9c53-b64c73da9e68.png" width="600" height="600">
