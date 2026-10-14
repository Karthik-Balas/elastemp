"""The agent's tools, as plain functions (no LLM or LangGraph imports).

Each returns a JSON string so the model sees structured, verifiable values. ``graph.py`` wraps
them as LangGraph tools bound to one project folder.
"""

import json
import logging
import math
import re

import pandas as pd

from elastemp.agent.status import project_status
from elastemp.base.utils import RESULTS_DIR, deformation_dir, resolve_workspace
from elastemp.cli import OPERATIONS, run_operation
from elastemp.dynamic.deformations import HIGHTEMP_DIR, PHONON_DIR
from elastemp.exceptions import ElastempError
from elastemp.validation.convergence import validate_vasp_run

# Steps that delete or overwrite folders that may hold finished VASP runs:
# step -> glob patterns (relative to the project folder) of the run folders it replaces.
DESTRUCTIVE_STEPS = {
    "make_elastic": ["deformation-*/strain-*"],
    "make_bulk_dynamic": [f"deformation-bulk/strain-*/{PHONON_DIR}"],
    "make_temp_static": [f"deformation-[0-9]*/strain-*/{HIGHTEMP_DIR}"],
    "make_dynamic": [f"deformation-[0-9]*/strain-*/{PHONON_DIR}"],
}

NUMBER = r"-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?"

OUTCAR_WARNING_PATTERNS = re.compile(
    r"WARNING|ERROR|Error|ZBRENT|hermitian|BRMIX|EDDDAV|EDWAV|VERY BAD NEWS|internal error|not converged",
)


def _dump(obj):
    return json.dumps(obj, indent=1, default=str)


def _round(x, sig=6):
    if x is None or not isinstance(x, int | float) or not math.isfinite(x) or x == 0:
        return x
    return round(x, max(0, sig - 1 - int(math.floor(math.log10(abs(x))))))


# --------------------------------------------------------------------------- project_status


def tool_project_status(workspace):
    return _dump(project_status(workspace))


# --------------------------------------------------------------------------- run_step


def _finished_runs_at_risk(ws, step):
    at_risk = []
    for pattern in DESTRUCTIVE_STEPS.get(step, []):
        for folder in sorted(ws.glob(pattern)):
            if (folder / "vasprun.xml").exists() or (folder / "OSZICAR").exists():
                at_risk.append(str(folder.relative_to(ws)))
    return at_risk


class _Capture(logging.Handler):
    def __init__(self):
        super().__init__(logging.INFO)
        self.lines = []

    def emit(self, record):
        self.lines.append(f"{record.levelname}: {record.getMessage()}")


def tool_run_step(workspace, step, confirm=False):
    """Run one workflow step. Refuses steps that would delete finished VASP runs unless confirm=True."""
    ws = resolve_workspace(workspace)
    if step not in OPERATIONS:
        return _dump({"status": "invalid", "message": f"Unknown step '{step}'. Valid steps: {list(OPERATIONS)}"})

    at_risk = _finished_runs_at_risk(ws, step)
    if at_risk and not confirm:
        return _dump(
            {
                "status": "refused",
                "message": (
                    f"{step} would delete or overwrite {len(at_risk)} folders that contain VASP output. "
                    "Do not retry on your own: ask the user to confirm explicitly, and only then call "
                    "run_step again with confirm=true."
                ),
                "examples": at_risk[:5],
            }
        )

    capture = _Capture()
    logger = logging.getLogger("elastemp")
    logger.addHandler(capture)
    previous_level = logger.level
    logger.setLevel(min(previous_level or logging.INFO, logging.INFO))
    try:
        run_operation(step, workspace=ws)
        status, message = "ok", f"{step} finished"
    except ElastempError as exc:
        status, message = "failed", str(exc)
    except Exception as exc:  # noqa: BLE001 - report anything else to the agent instead of crashing the chat
        status, message = "error", f"{type(exc).__name__}: {exc}"
    finally:
        logger.removeHandler(capture)
        logger.setLevel(previous_level)

    log = [line for line in capture.lines if not line.startswith("INFO: Running phonopy")]
    return _dump({"status": status, "message": message, "log": log[-40:]})


# --------------------------------------------------------------------------- inspect_run


def _resolve_inside(ws, folder):
    path = (ws / folder).resolve()
    if path != ws and ws not in path.parents:
        raise ValueError(f"{folder} is outside the project folder")
    return path


def tool_inspect_run(workspace, folder):
    """Convergence verdict plus OSZICAR/OUTCAR evidence for one VASP run folder."""
    ws = resolve_workspace(workspace)
    try:
        path = _resolve_inside(ws, folder)
    except ValueError as exc:
        return _dump({"error": str(exc)})
    if not path.is_dir():
        return _dump({"error": f"{folder} does not exist"})

    result = validate_vasp_run(path)
    report = {
        "folder": str(path.relative_to(ws)),
        "converged": result.converged,
        "electronic_converged": result.electronic,
        "reason": result.reason,
        "final_energy_eV": result.energy,
        "files": sorted(p.name for p in path.iterdir() if p.is_file())[:30],
    }
    oszicar = path / "OSZICAR"
    if oszicar.exists():
        report["oszicar_last_lines"] = oszicar.read_text(errors="replace").splitlines()[-3:]
    outcar = path / "OUTCAR"
    if outcar.exists():
        warnings, seen, finished = [], set(), False
        with open(outcar, errors="replace") as f:
            for line in f:
                if "General timing and accounting" in line:
                    finished = True
                if OUTCAR_WARNING_PATTERNS.search(line):
                    text = line.strip(" |\n")
                    if text and text not in seen:
                        seen.add(text)
                        warnings.append(text)
        report["outcar_finished_normally"] = finished
        report["outcar_warnings"] = warnings[:15]
    return _dump(report)


# --------------------------------------------------------------------------- read_results


def _parse_moduli(text):
    out = {}
    for kind, prop, value in re.findall(rf"(Voigt|Reuss|Hill) (bulk|shear|elastic) modulus =\s*({NUMBER})", text):
        out[f"{prop}_modulus_{kind.lower()}_GPa"] = float(value)
    for key, pattern in (
        ("poisson_ratio", rf"poisson ratio = ({NUMBER})"),
        ("pugh_ratio", rf"pugh's ratio = ({NUMBER})"),
        ("chen_hardness_GPa", rf"Hardness from Chen's model = ({NUMBER})"),
    ):
        m = re.search(pattern, text)
        if m:
            out[key] = float(m.group(1))
    return out


def _parse_stability(text):
    criteria = dict(re.findall(r"condition (\d+) = (.+)", text))
    verdicts = dict(re.findall(r"condition (\d+) constraint:\s*(met|NOT met)", text))
    return [{"criterion": criteria.get(n, ""), "met": verdicts[n] == "met"} for n in sorted(verdicts, key=int)]


def _thermal_section(results, independent):
    iso = pd.read_csv(results / "elastic_constants_temp_isothermal.csv")
    adi = pd.read_csv(results / "elastic_constants_temp_adiabatic.csv")
    mod_iso = pd.read_csv(results / "moduli_temp_isothermal.csv", sep="\t")
    vt = pd.read_csv(results / "volume_temp.txt", sep="\t")
    vt.columns = [c.strip() for c in vt.columns]

    tmax = float(iso["temp"].max())
    temps = [t for t in (300.0, tmax) if t <= tmax]
    base_c = iso[iso["temp"] == 0].iloc[0]
    base_m = mod_iso[mod_iso["T"] == 0].iloc[0]
    section = {"tmax_K": tmax, "at_temperature": {}}
    for T in temps:
        ci, ca = iso[iso["temp"] == T].iloc[0], adi[adi["temp"] == T].iloc[0]
        m, v = mod_iso[mod_iso["T"] == T].iloc[0], vt[vt["temperature"] == T].iloc[0]
        entry = {"isothermal_GPa": {}, "adiabatic_GPa": {}, "percent_change_from_0K_isothermal": {}}
        for name in independent:
            entry["isothermal_GPa"][name] = _round(ci[name])
            entry["adiabatic_GPa"][name] = _round(ca[name])
            if base_c[name]:
                entry["percent_change_from_0K_isothermal"][name] = round(100 * (ci[name] / base_c[name] - 1), 2)
        entry["hill_moduli_isothermal"] = {
            "bulk_GPa": _round(m["K"]),
            "shear_GPa": _round(m["Gh"]),
            "young_GPa": _round(m["Eh"]),
            "poisson_ratio": _round(m["nu"], 4),
            "chen_hardness_GPa": _round(m["hc"]),
        }
        for key, col in (("bulk", "K"), ("shear", "Gh"), ("young", "Eh")):
            entry["percent_change_from_0K_isothermal"][f"{key}_modulus"] = round(100 * (m[col] / base_m[col] - 1), 2)
        entry["linear_thermal_expansion_1e-6_per_K"] = round(float(v["Linear_thermal_expansion"]) * 1e6, 3)
        entry["volume_A3"] = _round(float(v["Volume_T"]))
        section["at_temperature"][f"{T:g}K"] = entry
    v0 = vt.iloc[0]
    section["static_volume_A3"] = _round(float(v0["Volume_0"]))
    section["volume_at_0K_with_zero_point_A3"] = _round(float(v0["Volume_T"]))
    section["zero_point_expansion_percent"] = round(100 * (float(v0["Volume_T"]) / float(v0["Volume_0"]) - 1), 3)
    return section


def tool_read_results(workspace):
    """All computed results and consistency checks, as structured numbers."""
    ws = resolve_workspace(workspace)
    results = ws / RESULTS_DIR
    raw = results / "elastic_constants_rawdata_0K.txt"
    if not raw.exists():
        return _dump({"error": "No results yet: get_elastic has not been run."})

    values = dict(re.findall(rf"c_(\d\d)=\s*({NUMBER})", raw.read_text()))
    independent = [f"c{k}" for k in sorted(values) if int(k[0]) <= int(k[1]) and float(values[k]) != 0]
    out = {"elastic_constants_0K_GPa": {name: float(values[name[1:]]) for name in independent}}

    if (results / "moduli.txt").exists():
        out["moduli_0K"] = _parse_moduli((results / "moduli.txt").read_text())
    if (results / "mechanical_stability_criteria.txt").exists():
        crit = _parse_stability((results / "mechanical_stability_criteria.txt").read_text())
        out["stability_0K"] = {"all_met": all(c["met"] for c in crit), "criteria": crit}

    checks, flags = {}, []
    b_file = ws / deformation_dir("bulk") / "B.txt"
    k_hill = out.get("moduli_0K", {}).get("bulk_modulus_hill_GPa")
    if b_file.exists() and k_hill:
        k_eos = float(re.search(NUMBER, b_file.read_text().split("=")[1]).group(0))
        gap = round(100 * (k_hill - k_eos) / k_eos, 2)
        checks["bulk_modulus_from_Cij_vs_EOS"] = {
            "from_Cij_GPa": k_hill,
            "from_EOS_GPa": k_eos,
            "difference_percent": gap,
        }
        if abs(gap) > 2:
            flags.append(f"Bulk modulus from C_ij and from the EOS fit differ by {gap}%")
    m = out.get("moduli_0K", {})
    if {"shear_modulus_voigt_GPa", "shear_modulus_hill_GPa", "shear_modulus_reuss_GPa"} <= m.keys():
        ok = m["shear_modulus_voigt_GPa"] >= m["shear_modulus_hill_GPa"] >= m["shear_modulus_reuss_GPa"]
        checks["voigt_hill_reuss_order_shear"] = ok
        if not ok:
            flags.append("Shear moduli violate Voigt >= Hill >= Reuss")
    if "stability_0K" in out and not out["stability_0K"]["all_met"]:
        flags.append("At least one Born stability criterion is not met at 0 K")

    dropped = {}
    for energies in sorted(ws.glob("deformation-*/energies.txt")):
        df = pd.read_csv(energies, sep="\t")
        if "converged" in df and not df["converged"].all():
            dropped[energies.parent.name] = [float(s) for s in df.loc[~df["converged"], "strains"]]
    if dropped:
        flags.append(f"Unconverged strain points were dropped from the fits: {dropped}")
    out["checks"] = checks
    out["dropped_strain_points"] = dropped

    if (results / "moduli_temp_isothermal.csv").exists():
        out["temperature_dependence"] = _thermal_section(results, independent)
    out["flags"] = flags
    return _dump(out)
