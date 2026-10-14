"""Deterministic project-state scanner: which workflow steps are done, and what comes next.

No LLM here. The agent calls ``project_status`` as a tool, but the answer is computed from
the files in the project folder, so it is the same every time.
"""

import re
from dataclasses import asdict, dataclass, field

from elastemp.base.strain_analysis import MIN_FIT_POINTS
from elastemp.base.symmetry import get_num_constants
from elastemp.base.utils import (
    RESULTS_DIR,
    STRAINS_FILE,
    SYMMETRY_FILE,
    deformation_dir,
    read_strains,
    read_symmetry,
    resolve_workspace,
    strain_dir,
)
from elastemp.dynamic.deformations import HIGHTEMP_DIR, PHONON_DIR
from elastemp.validation.convergence import validate_vasp_run

MAX_FAILED_LISTED = 20  # failed runs listed individually (with reasons) before summarising the rest


@dataclass
class Stage:
    name: str
    kind: str  # "step" (an elastemp command) or "vasp" (runs the user submits)
    status: str  # "done", "partial", "not started", "blocked"
    detail: str = ""
    failed: list = field(default_factory=list)  # "folder: reason" for each failed run
    not_run: list = field(default_factory=list)  # compact summary per deformation folder


def _run_state(folder):
    """'converged', 'failed: <reason>', or 'not run' for one VASP folder."""
    if not (folder / "vasprun.xml").exists() and not (folder / "OSZICAR").exists():
        return "not run"
    result = validate_vasp_run(folder)
    return "converged" if result.converged else f"failed: {result.reason}"


def _vasp_stage(name, folders, ws, *, tolerate_failures=False, per_group=None):
    """Summarise a group of VASP runs.

    :param tolerate_failures: static strain runs may lose a few points (drop-and-refit),
                              so a failure there does not block the next step.
    :param per_group: dict group -> list of folders, used to check the 7-point minimum.
    """
    if not folders:
        return Stage(name, "vasp", "not started", "no run folders created yet")
    states = {f: _run_state(f) for f in folders}
    converged = [f for f, s in states.items() if s == "converged"]
    failed = [f"{f.relative_to(ws)}: {s[len('failed: ') :]}" for f, s in states.items() if s.startswith("failed")]
    not_run = _summarise_by_deformation([f for f, s in states.items() if s == "not run"], ws)
    n_not_run = sum(1 for s in states.values() if s == "not run")
    detail = f"{len(converged)} of {len(folders)} converged, {len(failed)} failed, {n_not_run} not run"

    if len(converged) == len(folders):
        status = "done"
    elif n_not_run:
        status = "partial"
    elif tolerate_failures and per_group:
        short = [g for g, fs in per_group.items() if sum(states[f] == "converged" for f in fs) < MIN_FIT_POINTS]
        status = "blocked" if short else "done"
        if short:
            detail += f"; fewer than {MIN_FIT_POINTS} converged points in: {', '.join(short)}"
        else:
            detail += " (failed points will be dropped from the fit)"
    else:
        status = "blocked"
    if len(failed) > MAX_FAILED_LISTED:
        failed = failed[:MAX_FAILED_LISTED] + [f"... and {len(failed) - MAX_FAILED_LISTED} more failed runs"]
    return Stage(name, "vasp", status, detail, failed, not_run)


def _summarise_by_deformation(folders, ws):
    """Group run folders by deformation: 'deformation-2: all 7 strains' or 'deformation-2: strain-3, strain-5'."""
    groups = {}
    for folder in folders:
        parts = folder.relative_to(ws).parts
        groups.setdefault(parts[0], []).append("/".join(parts[1:]))
    out = []
    for deformation, runs in groups.items():
        total = len(list((ws / deformation).glob("strain-*")))
        if len(runs) == total:
            out.append(f"{deformation}: all {total} strains")
        else:
            out.append(f"{deformation}: {', '.join(sorted(runs, key=_strain_number))}")
    return out


def _strain_number(run):
    match = re.search(r"strain-(\d+)", run)
    return int(match.group(1)) if match else 0


def _step_stage(name, done, detail=""):
    return Stage(name, "step", "done" if done else "not started", detail)


def project_status(workspace=None):
    """Return the state of every workflow stage and the recommended next action, as a dict."""
    ws = resolve_workspace(workspace)
    report = {"workspace": str(ws), "stages": [], "next_action": ""}

    if not (ws / "POSCAR").exists():
        report["next_action"] = "Add POSCAR (and POTCAR) to the project folder."
        return report
    if not (ws / "POTCAR").exists():
        report["warnings"] = ["No POTCAR in the project folder: the make_* steps will need it."]
    if not (ws / SYMMETRY_FILE).exists() or not (ws / STRAINS_FILE).exists():
        report["stages"].append(asdict(_step_stage("make_elastic", False)))
        report["next_action"] = "Run the make_elastic step."
        return report

    symmetry = read_symmetry(ws)
    n_def = get_num_constants(symmetry)
    n_strain = len(read_strains(ws))
    report.update(symmetry=symmetry, independent_constants=n_def, strains=n_strain)

    deformations = [deformation_dir(k) for k in range(1, n_def + 1)]
    all_defs = deformations + [deformation_dir("bulk")]
    strains = [strain_dir(i) for i in range(1, n_strain + 1)]
    results = ws / RESULTS_DIR

    def folders(defs, sub=""):
        out = [ws / d / s / sub if sub else ws / d / s for d in defs for s in strains]
        return [f for f in out if f.exists()]

    static_groups = {d: [ws / d / s for s in strains if (ws / d / s).exists()] for d in all_defs}
    created = all((ws / d).exists() for d in all_defs)
    stages = [
        _step_stage("make_elastic", created),
        _vasp_stage(
            "VASP static runs (deformation-*/strain-*)",
            folders(all_defs),
            ws,
            tolerate_failures=True,
            per_group=static_groups,
        ),
        _step_stage("get_elastic", (results / "elastic_constants_rawdata_0K.txt").exists()),
        _step_stage(
            "make_bulk_dynamic",
            all((ws / deformation_dir("bulk") / s / PHONON_DIR / "POSCAR-unitcell").exists() for s in strains),
        ),
        _vasp_stage("VASP DFPT runs (deformation-bulk)", folders([deformation_dir("bulk")], PHONON_DIR), ws),
        _step_stage("get_bulk_dynamic", (results / "volume_temp.txt").exists()),
        _step_stage(
            "make_temp_static",
            all((ws / d / s / HIGHTEMP_DIR / "POSCAR").exists() for d in deformations for s in strains),
        ),
        _vasp_stage("VASP Hightemp_static runs", folders(deformations, HIGHTEMP_DIR), ws),
        _step_stage(
            "get_temp_static",
            all((ws / d / s / HIGHTEMP_DIR / "ZPE.txt").exists() for d in deformations for s in strains),
        ),
        _step_stage(
            "make_dynamic",
            all((ws / d / s / PHONON_DIR / "POSCAR-unitcell").exists() for d in deformations for s in strains),
        ),
        _vasp_stage("VASP DFPT runs (deformations)", folders(deformations, PHONON_DIR), ws),
        _step_stage(
            "get_dynamic",
            all((ws / d / s / PHONON_DIR / "FORCE_CONSTANTS").exists() for d in deformations for s in strains),
        ),
        _step_stage("get_thermal", (results / "moduli_temp_isothermal.csv").exists()),
    ]
    report["stages"] = [asdict(s) for s in stages]
    report["next_action"] = _next_action(stages)
    return report


def _next_action(stages):
    for stage in stages:
        if stage.status == "done":
            continue
        if stage.kind == "step":
            return f"Run the {stage.name} step."
        if stage.status == "blocked":
            return (
                f"{stage.name}: {stage.detail}. Inspect the failed runs (inspect_run), fix and rerun them "
                "in VASP, then continue."
            )
        if stage.status == "not started":
            return f"{stage.name}: run folders have not been created yet."
        return f"Run VASP (the user submits these jobs): {stage.name}, {stage.detail}."
    return "All steps are complete. Results are in results_dir/."
