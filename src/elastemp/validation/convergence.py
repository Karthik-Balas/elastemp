"""Validate VASP completion and convergence using vasprun.xml.

OSZICAR's ``F=`` line only shows that an ionic step finished: VASP writes it even when the
SCF loop stopped at NELM without converging. vasprun.xml records every SCF step and the NELM
setting, so the electronic loop can be checked properly.

The rules follow pymatgen's ``Vasprun.converged_electronic`` / ``converged_ionic``, but only
the few elements they need are read. pymatgen's full parser is not used because it converts
every number in the file and fails on valid VASP output where Fortran drops the ``E`` from
3-digit exponents (e.g. ``3.37026300-199``), which happens in DFPT runs.

A missing, truncated (job killed before VASP closed the file) or malformed vasprun.xml is
never treated as a successful run.
"""

import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path

DEFAULT_NELM = 60


@dataclass(frozen=True)
class ConvergenceResult:
    electronic: bool
    ionic: bool
    converged: bool
    energy: float | None  # final free energy TOTEN (e_fr_energy, = OSZICAR F) in eV, if the SCF converged
    reason: str | None


def _value(node):
    return (node.text or "").strip() if node is not None else None


def _as_bool(text):
    return text is not None and text.strip().upper().strip(".") in {"T", "TRUE"}


def _as_int(text, default):
    try:
        return int(float(text))
    except (TypeError, ValueError):
        return default


def _as_float(text, default):
    try:
        return float(text)
    except (TypeError, ValueError):
        return default


def _parameter(parameters, name, section=None):
    """First value of ``name`` in <parameters>, preferring the separator called ``section``."""
    if parameters is None:
        return None
    if section is not None:
        for separator in parameters.iter("separator"):
            if separator.get("name") == section:
                node = separator.find(f".//*[@name='{name}']")
                if node is not None:
                    return _value(node)
    return _value(parameters.find(f".//*[@name='{name}']"))


def _electronic_converged(incar, nelm, scsteps):
    if _as_bool(incar.get("ML_LMLFF")):
        return True  # machine-learned force-field MD has no SCF steps
    if (incar.get("ALGO") or "").lower() == "exact" and _as_int(incar.get("NELM"), None) == 1:
        return True
    if _as_bool(incar.get("LEPSILON")):
        # Linear-response steps appended after the SCF carry only these three energies.
        response_only = {"e_wo_entrp", "e_fr_energy", "e_0_energy"}
        idx = 1
        while idx < len(scsteps) and {i.get("name") for i in scsteps[idx].iter("i")} == response_only:
            idx += 1
        return idx + 1 != nelm
    return len(scsteps) < nelm


def _ionic_converged(parameters, n_ionic_steps):
    nsw = _as_int(_parameter(parameters, "NSW", "ionic"), 0)
    ibrion = _as_int(_parameter(parameters, "IBRION", "ionic"), -1 if nsw in (-1, 0) else 0)
    if ibrion == 0:  # molecular dynamics: "converged" once all NSW steps ran
        return nsw <= 1 or n_ionic_steps == nsw
    ediffg = _as_float(_parameter(parameters, "EDIFFG", "ionic"), 1.0)
    if ibrion in (1, 2) and ediffg == 0:
        return nsw <= 1 or n_ionic_steps == nsw
    return nsw <= 1 or n_ionic_steps < nsw


def validate_vasp_run(run_folder: str | Path, *, require_ionic: bool = False) -> ConvergenceResult:
    """Inspect the finished VASP run in ``run_folder`` using its vasprun.xml.

    Call this after VASP has exited, not while vasprun.xml is still being written.

    :param require_ionic: also require ionic (relaxation) convergence
    """
    path = Path(run_folder) / "vasprun.xml"
    if not path.is_file():
        return ConvergenceResult(False, False, False, None, f"Missing {path}")
    try:
        root = ET.parse(path).getroot()
    except (ET.ParseError, OSError) as exc:
        return ConvergenceResult(False, False, False, None, f"Truncated or unreadable {path} ({exc})")

    incar_node = root.find("incar")
    incar = {node.get("name"): _value(node) for node in incar_node} if incar_node is not None else {}
    parameters = root.find("parameters")
    nelm = _as_int(_parameter(parameters, "NELM", "electronic"), DEFAULT_NELM)
    calculations = root.findall("calculation")
    if not calculations:
        return ConvergenceResult(False, False, False, None, f"No completed ionic step in {path}")

    scsteps = calculations[-1].findall("scstep")
    electronic = _electronic_converged(incar, nelm, scsteps)
    ionic = _ionic_converged(parameters, len(calculations))
    converged = electronic and (ionic or not require_ionic)

    energy = None
    if electronic:
        energy = _as_float(_value(calculations[-1].find("energy/i[@name='e_fr_energy']")), None)

    if converged:
        reason = None
    elif not electronic:
        reason = f"Electronic SCF did not converge ({len(scsteps)} steps, NELM = {nelm}) in {path}"
    else:
        reason = f"Ionic relaxation did not converge ({len(calculations)} ionic steps) in {path}"
    return ConvergenceResult(electronic, ionic, converged, energy, reason)
