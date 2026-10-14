"""Command-line entry point: ``elastemp --operations <step> [--workdir PROJECT]``."""

import argparse
import logging
import sys

from elastemp.base.utils import resolve_workspace
from elastemp.exceptions import ElastempError

logger = logging.getLogger("elastemp")

OPERATIONS = (
    "make_elastic",
    "get_elastic",
    "make_bulk_dynamic",
    "get_bulk_dynamic",
    "make_temp_static",
    "get_temp_static",
    "make_dynamic",
    "get_dynamic",
    "get_thermal",
)


def configure_logging(level="INFO", log_file="elastemp.log"):
    """Console logging at ``level``; full DEBUG log (incl. phonopy output) appended to ``log_file``."""
    root = logging.getLogger("elastemp")
    root.setLevel(logging.DEBUG)
    root.handlers.clear()

    console = logging.StreamHandler()
    console.setLevel(getattr(logging, level.upper()))
    console.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
    root.addHandler(console)

    if log_file:
        fh = logging.FileHandler(log_file)
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(name)s: %(message)s"))
        root.addHandler(fh)


def run_operation(operation, symprec=0.05, angle_tolerance=10, workspace=None):
    """Run one workflow step on the project in ``workspace`` (default: current directory)."""
    # Imported lazily so `elastemp --help` stays fast.
    from elastemp.input.make_dynamic_input import get_dynamic_wrapper, get_T_Zpe, make_dynamic_wrapper, make_T_Zpe
    from elastemp.input.make_elastic_input import get_elastic_wrapper, make_elastic_wrapper
    from elastemp.input.make_thermal_input import get_elastic_temp_plots, get_thermal_constants, get_thermal_wrapper

    ws = resolve_workspace(workspace)
    if operation == "make_elastic":
        make_elastic_wrapper(symprec=symprec, angle_tolerance=angle_tolerance, workspace=ws)
    elif operation == "get_elastic":
        get_elastic_wrapper(make_plots=True, workspace=ws)
    elif operation == "make_bulk_dynamic":
        make_dynamic_wrapper(calc_bulk_dynamic=True, workspace=ws)
    elif operation == "get_bulk_dynamic":
        get_dynamic_wrapper(calc_bulk_dynamic=True, workspace=ws)
    elif operation == "make_temp_static":
        make_T_Zpe(ws)
    elif operation == "get_temp_static":
        get_T_Zpe(ws)
    elif operation == "make_dynamic":
        make_dynamic_wrapper(calc_bulk_dynamic=False, workspace=ws)
    elif operation == "get_dynamic":
        get_dynamic_wrapper(calc_bulk_dynamic=False, workspace=ws)
    elif operation == "get_thermal":
        get_thermal_wrapper(ws)
        get_thermal_constants(ws)
        get_elastic_temp_plots(ws)
    else:
        raise ValueError(f"Unknown operation '{operation}'")


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="elastemp",
        description="Quasi-harmonic temperature-dependent elastic constants (VASP + phonopy).",
    )
    parser.add_argument("--operations", required=True, choices=OPERATIONS, help="workflow step to run")
    parser.add_argument(
        "--workdir", default=".", help="project folder containing POSCAR, POTCAR, ... (default: current directory)"
    )
    parser.add_argument("--symprec", type=float, default=0.05, help="symmetry tolerance for lattice parameters")
    parser.add_argument("--angle-tolerance", type=float, default=10, help="symmetry tolerance for angles (deg)")
    parser.add_argument(
        "--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"], help="console log level"
    )
    parser.add_argument(
        "--log-file", default="elastemp.log", help="DEBUG log file, relative to --workdir ('' to disable)"
    )
    args = parser.parse_args(argv)

    workspace = resolve_workspace(args.workdir)
    if not workspace.is_dir():
        parser.error(f"--workdir {workspace} is not a directory")
    configure_logging(args.log_level, workspace / args.log_file if args.log_file else "")
    logger.info("elastemp %s in %s", args.operations, workspace)
    try:
        run_operation(args.operations, args.symprec, args.angle_tolerance, workspace)
    except ElastempError as exc:
        logger.error("%s", exc)
        return 1
    except Exception:
        logger.exception("Unexpected error during %s", args.operations)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
