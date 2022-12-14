# Elastemp

Elastemp is a quasi harmonic approximation to estimate temperature dependent elastic constants. It integrates VASP qm simulation package with Phonopy to set up the calculations and extract information. To use this this automated workflow,clone this repo to your directory. Please have both the elastemp folder and elastic_input.py in the same folder where you want to run calculations. We are working on making this a pip package in subsequent versions. 

# Files needed 

To run vasp calculations, user needs to input
1. POSCAR (mandatory). 
2. POTCAR (mandatory).
3. KPOINTS (optional) - a default KPOINTS file will be provided. But for best results, users need to tune it depending on the structure.
4. INCAR   (optional) - a default INCAR file will be provided. But for best results, users need to tune it depending on the structure.
