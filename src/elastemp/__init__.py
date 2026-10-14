"""
Elastemp: quasi-harmonic temperature-dependent elastic constants from VASP + Phonopy.

Distributed under the MIT license. Author: Karthik Balasubramanian.

The library only *emits* log records (via ``logging.getLogger(__name__)`` in each
module). It never configures handlers itself; the CLI (``elastemp.cli``) or the
calling application decides where logs go.
"""

import logging

__version__ = "0.2.2"

logging.getLogger(__name__).addHandler(logging.NullHandler())
