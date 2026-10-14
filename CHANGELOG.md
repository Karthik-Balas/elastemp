# Changelog

All notable changes to this project are documented here.

## [0.2.1]

* Default strains.txt now has 11 values (-0.05 to 0.05 in steps of 0.01) instead of 7, so a failed VASP run can be dropped without falling below the 7-point minimum.
* Agent answers are plain text for the terminal; any Markdown the model still uses is stripped.
* Agent status lists every unfinished run as a compact per-deformation summary (e.g. "deformation-2: strain-3, strain-5") instead of a truncated list.

## [0.2.0]

* `print` replaced by `logging`; `sys.exit` replaced by exceptions (`elastemp.exceptions`), so failures now give a non-zero exit status.
* `os.system` replaced by `subprocess`, so phonopy failures are no longer silently ignored.
* The working directory is always restored after each step, even if the step fails.
* The nine per-symmetry deformation classes are replaced by one table of strain patterns (`elastemp/elastic/deformations.py`). The stability criteria are also table-driven.
* Installable package with an `elastemp` console command; `elastic_input.py` still works.
* Works with current numpy, pandas 2 and pymatgen. Strain files are now parsed exactly; pandas' parser could be off in the last digit.
* **Bug fixes that change results:**
  * Voigt shear modulus. It now uses the standard G_V = [(C11+C22+C33) - (C12+C23+C13) + 3(C44+C55+C66)]/15. The previous code used C12 in place of C22 and 4 in place of 3, which overestimated G_V (by about 10% for diamond). This also changed G_H, E, Poisson's ratio, Pugh's ratio and hardness.
  * Monoclinic stability condition 11 now uses the correct determinant.
  * Trigonal stability: the C13 condition is now evaluated. Previously it was listed but never checked.
  * Tetragonal-II stability: a typo raised an error.
  * Triclinic stability: now tests positive-definiteness. det(C) > 0 is not sufficient. Every symmetry also reports this general Born criterion.
* Other fixes: an undefined variable in the convergence check; a crash on a singular stiffness matrix; matplotlib figures were never closed; the bulk EOS fit was re-run as a side effect of every elastic-constant calculation.

# Citing our work

If you use our work, please cite the following paper in your work. 

"Karthik Balasubramanian, Sukriti Manna, Subramanian K.R.S. Sankaranarayanan,
Elastemp — A workflow to compute the quasi-harmonic temperature dependent elastic constants of materials,
Computational Materials Science,Volume 226, 2023"

### Unconverged strain points

* Convergence is checked from `vasprun.xml` (`elastemp.validation.convergence.validate_vasp_run`) instead of the OSZICAR `F=` line, which VASP also writes when the SCF hits NELM without converging. This applies to static strain runs, `Hightemp_static` runs and DFPT phonon runs. A missing or truncated `vasprun.xml` counts as not converged. `check_convergence(..., legacy_oszicar=True)` keeps the old test for data without `vasprun.xml`.
* The check reads only the elements it needs instead of using pymatgen's full `Vasprun` parser. That parser fails on valid DFPT output where Fortran drops the `E` from 3-digit exponents (e.g. `3.37026300-199`). That affected 17 of the 44 phonon runs in the diamond example. On the 104 runs pymatgen can parse, both give identical results.

* An unconverged strain point is now dropped from the fit, with a warning, instead of being given a placeholder energy of 1e4 eV. The placeholder made that deformation's fit, and every C_ij built from it, meaningless while the run still reported success.
* The fit needs at least 7 converged points, on both sides of zero strain. Otherwise the step stops with an error. `make_elastic` now refuses a `strains.txt` with fewer than 7 values, and warns if it has exactly 7, because then any single failure stops `get_elastic`.
* `energies.txt` has a new `converged` column. Unconverged points have empty (NaN) energies.

### Project folder instead of current directory

* Every workflow function takes a `workspace` (project folder) argument; `None` means the current directory. elastemp no longer calls `os.chdir` anywhere, and phonopy runs with `subprocess.run(cwd=...)`.
* New CLI option `--workdir`. The log file is written inside the project folder.
* Plots use matplotlib's `Figure` API instead of `pyplot`, which keeps global state.
* Readers take the project folder: `read_strains(workspace)` and so on, replacing `read_strains(path)`. `get_curvature` and `extract_force_constants` take the deformation folder and the phonon folder respectively.
* B.txt is written to `deformation-bulk/`, as before.
* Phonon set-up checks every input file before deleting an existing `Phonon_calculation` folder.

### Agent (optional, `elastemp[agent]`)

* New `elastemp-agent` command: a LangGraph agent using Claude, with four tools (`project_status`, `run_step`, `inspect_run`, `read_results`).
* `run_step` refuses steps that would delete finished VASP runs unless the user confirms.
* A verification node checks every number in the agent's answer against the tool results; an invented value triggers one rewrite, then an `[Unverified]` note.
* The API key comes from `ANTHROPIC_API_KEY` or a local `.env` file (template: `.env.example`); no key is stored in the code.

### Packaging

* `src/` layout with a PEP 621 `pyproject.toml`; the version is defined once, in `elastemp.__version__`.
* Console script `elastemp` and `python -m elastemp`. The old `elastic_input.py` now lives in `scripts/`.
* Optional extras: `pip install "elastemp[test]"` and `pip install "elastemp[dev]"`.
* GitHub Actions CI runs the tests and builds the package.

## [0.1.0] - 2023

* Initial release accompanying the Computational Materials Science paper.
