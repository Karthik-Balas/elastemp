# Elastemp

Elastemp is a quasi harmonic approximation to estimate temperature dependent elastic constants. It integrates VASP qm simulation package with Phonopy to set up the calculations and extract information.

# Installation

From source:

```
git clone https://github.com/Karthik-Balas/elastemp.git
cd elastemp
pip install .                 # or `pip install -e ".[dev]"` for development
```

Or with conda: `conda env create -f environment.yml && conda activate Elastemp-env`.

Requires Python 3.10 or newer. Dependencies (installed automatically): numpy, pandas, matplotlib, pyyaml, pymatgen and phonopy. phonopy must be below version 4, because its v4 command-line interface is not supported yet.

Run every step below from the folder that holds your POSCAR/POTCAR:

```
elastemp --operations <step>          # or: python -m elastemp --operations <step>
```

The old `python3 elastic_input.py --operations <step>` still works: the script is in `scripts/` and needs elastemp to be installed.

To run a step on a project in another folder, add `--workdir <project folder>`. Nothing is read from or written to the folder you run the command from.

# The elastemp agent (optional)

An assistant you can ask about a project in plain language. It is built with LangGraph and uses Claude. It can tell you where a project stands, run the next analysis step, explain why a VASP run failed, and summarize the results. It never runs or submits VASP jobs; you still do that.

```
pip install "elastemp[agent]"
cp .env.example .env          # then put your Claude API key in .env
elastemp-agent --workdir runs/TiN
```

You can also ask a single question: `elastemp-agent --workdir runs/TiN --ask "What should I run next?"`.

**API key:** the agent reads `ANTHROPIC_API_KEY` from the environment or from a `.env` file in the current or project folder. `.env` is in `.gitignore`. Never commit a real key; `.env.example` holds empty placeholders only. Set `ELASTEMP_AGENT_MODEL` in `.env` (or pass `--model`) to choose the Claude model.

**What it can do** (its four tools):

| Tool | Purpose |
|---|---|
| `project_status` | Which steps are done, which VASP runs are finished, failed or missing, and what comes next |
| `run_step` | Run one workflow step. A step that would delete finished VASP runs (for example rerunning `make_bulk_dynamic`) is refused unless you explicitly confirm |
| `inspect_run` | Convergence verdict and OSZICAR/OUTCAR evidence for one run |
| `read_results` | Elastic constants, moduli, stability, temperature dependence, consistency checks and flags |

**Number check:** every number in the agent's answer is compared, in code, with the tool results. If it states a value that came from nowhere (an "experimental" value from memory, say), it is asked to rewrite the answer once. If it still does, the answer is shown with an `[Unverified: ...]` note.

The physics package never imports the agent, so `pip install elastemp` without `[agent]` installs no LLM libraries.

# Using elastemp from Python

Every step takes the project folder explicitly, and elastemp never changes the process's working directory. Independent projects can therefore be processed side by side, for example in threads:

```python
from concurrent.futures import ThreadPoolExecutor
from elastemp.cli import run_operation

def analyse(project):
    for step in ["get_elastic", "get_bulk_dynamic", "get_temp_static", "get_dynamic", "get_thermal"]:
        run_operation(step, workspace=project)

with ThreadPoolExecutor() as pool:
    list(pool.map(analyse, ["runs/diamond", "runs/TiN"]))
```

Two steps on the *same* project folder must still run one after the other, because each step reads the files the previous one wrote.

# Repository layout

```
src/elastemp/          the package
  cli.py               command-line entry point
  base/                symmetry, strain fitting, shared I/O helpers
  elastic/             deformations, C_ij assembly, moduli, stability
  dynamic/             phonon supercell set-up and force constants
  thermal/             QHA thermal properties and C_ij(T)
  input/               the workflow steps called by the CLI
tests/                 pytest suite (no VASP or phonopy needed)
examples/diamond/      a complete worked example
scripts/               legacy elastic_input.py
```

# Logging

Progress is reported through Python's `logging` module rather than `print`.

* `--log-level {DEBUG,INFO,WARNING,ERROR}` sets the console verbosity (default `INFO`).
* A full DEBUG log, including phonopy's own output, is written to `elastemp.log` in the working directory. Use `--log-file <path>` to change it or `--log-file ''` to disable it.
* Failures (missing POTCAR/OSZICAR/vasprun.xml, unconverged runs, phonopy errors) end the run with a non-zero exit status and an `ERROR:` message.

When elastemp is used as a library, it only emits log records; configure logging in your own application (e.g. `logging.basicConfig(level=logging.INFO)`).

# Files needed 

To run vasp calculations, user needs to input
1. POSCAR (mandatory). 
2. POTCAR (mandatory).
3. KPOINTS (optional) - a default KPOINTS file will be provided. But for best results, users need to tune it depending on the structure.
4. INCAR   (optional) - a default INCAR file will be provided. But for best results, users need to tune it depending on the structure.
5. strains.txt (optional) - a default strains.txt file will be provided with 11 strain values from -0.05 to 0.05 in steps of 0.01. 
6. dim.txt (optional) - a default dim.txt file with default values 2 2 2 will be provided. This file contains the dimensions needed to make supercells                           for phonon calculations. For best results, users need to tune it depending on the structure.
7. tmax.txt (optional) - a default tmax.txt file with defualt value of 1000 will be provided. This file contains the maximum temperature at which elastic                           constant calculations are to be performed.
8. KPOINTS_dynamic (optional) - a default KPOINTS file for dynamic calculations will be provided. But for best results, users need to tune it depending                                   on the structure. 
9. INCAR_dynamic (optional) - a default INCAR file for dynamic calculations will be provided. But for best results, users need to tune it depending                                   on the structure. 
10. mesh.conf           - This is automatically created by Elastemp and needed for phonopy calculations.

# Workflow

<p align="center">
<img src="https://user-images.githubusercontent.com/120595580/207714085-196181cf-5f77-46b5-9c53-b64c73da9e68.png" width="600" height="600">
</p>

The workflow consists of three modules -  Getting zero temp elastic constants, Getting thermal expansion coeff and Getting elastic constants at temp T. 

### Steps before starting calculations.
1. Install elastemp (see above). Have POSCAR and POTCAR in the directory where you run the calculations. INCAR and KPOINTS are optional. A default one will be used otherwise.
2. A file strains.txt containing strain values as a column. However a default set of strains will be used otherwise and they are accurate for most cases.

## Getting zero temp. elastic constants (Steps 1 & 2)

### Step 1: This creates the deformation folders to run vasp calculations

command : elastemp --operations make_elastic

1. This will recognize the symmetry of the structure and make deformation-1, deformation-2,.... deformation-n and deformation-bulk folders. n is the      number of independent elastic constants.
2. User needs to run vasp in these folders.

### Step 2: This gets the results from vasp calculations for zero temp.elastic constants

command : elastemp --operations get_elastic

1. Elastemp extracts enthalpies from folders and outputs results to results_dir. The files in results_dir are<br/>
   a. In each deformation folder, a plot of energy vs strain, fitting function and curvature constant for fitting.\
   b. In results_dir the following files will be created\
      (i)   stiffness_matrix.txt - Stiffness matrix will be created\
      (ii)  compliance_matrix.txt - Compliance matrix will be created\
      (iii) moduli.txt            - A text file with different moduli (K, E, G, H, nu, pugh) values created\
      (iv)  stability_criteria.txt - A file which checks the mechanical stability criteria of materials. 
      
 ## Getting thermal expansion coefficient (Steps 3 & 4)
 
 ### Step 3: This makes the supercell for phonon calculation for volumetric calculations to get thermal expansion coeff.
 
 command : elastemp --operations make_bulk_dynamic
 
 1. This creates files in the deformation bulk folder using phonopy. 
 2. Users needs to run vasp dfpt calculations in deformation-bulk/strain/Phonon_calculation folders. Currently, we support only DFPT calculations. 
 
 ### Step 4: This extracts results from phonon calculations in deformation-bulk folder. 
 
 command : elastemp --operations get_bulk_dynamic
 
 1. The force constants are extracted by phonopy and QHA calculations are done to get Helmholtz energies. Birch-murnaghan equation of state is fit to get equilibrium volumes at different temps.
 2. The following files are written to results_dir.<br/>
    a. volume_temp.txt - file containing equilibrium volumes, thermal expansion coefficient at all temperatures.\
    b. Plot of thermal expansion coefficient vs temperature.
 3. The scaling factor of lattice parameters corresponding to the structure at tmax is calculated. 
 
 ## Getting elastic constants at high temperature (Steps 5-9). 
  
  This has 2 submodules - getting zero temp enthalpies of structure with lattice parameters corresponding to tmax, getting elastic constants at tmax and interpoltaing them to all temperatures t<tmax. 
  
  ### Step 5: This makes the folders for running vasp calculations to get zero temp enthalpies. 
  
  command : elastemp --operations make_temp_static
  
  1. This creates subfolders in deformation_i (deformation-1,..deformation-n) folders.
  2. Users needs to run vasp calculations in these sub folders. 
  
  ### Step 6: Getting zero temp. enthalpies
  
  command : elastemp --operations get_temp_static
  
  1. This extracts the zero temp. enthalpies and stores them in ZPE.txt files as enthalpy/atom. 
  
  ### Step 7: This creates the folders for running phonon calculations in each deformation folder. 
  
  command : elastemp --operations make_dynamic
  
  1. This creates the subfolders (Phonon_calculation) in each deformation folder and copies the ZPE.txt file to be saved here.
  2. User needs to run DFPT calculations similar to the deformation-bulk folder.

   ### Step 8: This extracts the force constants after phonon calculations are run in each deformation folder.
   
   command : elastemp --operations get_dynamic
   
   1. This extracts the force constants from phonon calculations. 

   ### Step 9: Extracting the elastic constants as a function of temperature. 
   
   command : elastemp --operations get_thermal
   
   1. This runs phonon QHA calculations in each phonon calculation folder, fits them to a polynomial of strain to get elastic constants at tmax.
   2. Linear interpolation is done using elastic constants at 0K and tmax to get elastic constants, and subsequently moduli, at all temperatues in between. 
   3. In each deformation folder, a plot of energy-strain at tmax is created.
   4. In results_dir, moduli_plots folder which has plots of Cij, E, G, H, nu, K as a function of temperature are created. Csv files with the adiabatic and isothermal elastic constants at each temperature is also stored. 

# Tests

```
pip install -e ".[test]"
pytest
```

To build the distribution files:

```
pip install build twine
python -m build          # writes dist/elastemp-<version>.tar.gz and .whl
twine check dist/*
```

The tests check, for every crystal class, that the strain patterns and the C_ij assembly reproduce a known stiffness matrix, plus the Voigt-Reuss-Hill averages and the stability criteria. They do not need VASP or phonopy.

# Limitations :
1. Isotropic linear expansion coefficient assumption. 
2. At higher temperatures, for 'softer' materials, fitting can be bad. Users need to experiment with supercell size and Kpoint convergence. 
3. Convergence is judged from each run's `vasprun.xml`, not from OSZICAR. VASP writes an `F=` line to OSZICAR even when the SCF stopped at NELM without converging, so OSZICAR alone can't tell a converged run from an unconverged one. A run counts as converged if its final ionic step needed fewer SCF steps than NELM (the same rule pymatgen uses). A missing or truncated `vasprun.xml`, which is what a job killed at walltime leaves behind, counts as not converged. So keep `vasprun.xml` for every static and DFPT run.
4. Strain points whose VASP run did not converge are left out of the energy-strain (or EOS) fit, and a warning names them. If fewer than 7 converged points remain, or the remaining strains are not on both sides of zero, the step stops with an error. Leaving out a point is not entirely neutral: real energy-strain curves have small cubic terms, so an uneven set of strains shifts the fitted curvature slightly (about 1% in C11 for the diamond example when one point is removed). Rerunning the failed calculation is still best.
5. The phonopy free energy (kJ/mol) is converted to eV with a factor of 0.01. The exact factor is 1/96.485 = 0.010364. It is kept at 0.01 so results match the published version.

See [CHANGELOG.md](CHANGELOG.md) for release notes.
