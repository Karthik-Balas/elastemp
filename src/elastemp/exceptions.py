"""Exception hierarchy for elastemp.

Library code raises these instead of calling ``sys.exit()``, so callers (the CLI,
notebooks, or an orchestration layer) decide how to react to a failure.
"""


class ElastempError(Exception):
    """Base class for all elastemp errors."""


class MissingInputError(ElastempError):
    """A required input file or folder (POSCAR, POTCAR, OSZICAR, vasprun.xml, ...) is missing."""


class ConvergenceError(ElastempError):
    """A VASP calculation did not converge or its output could not be parsed."""


class ExternalToolError(ElastempError):
    """An external program (phonopy) failed or is incompatible."""
