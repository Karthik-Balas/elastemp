"""
Backwards-compatible entry point. Equivalent to the ``elastemp`` console script:

    python3 elastic_input.py --operations make_elastic

Elastemp is an open-source python package, distributed under the MIT license, to compute
quasi-harmonic temperature-dependent elastic constants of materials.
Author: Karthik Balasubramanian
"""

import sys

from elastemp.cli import main

if __name__ == "__main__":
    sys.exit(main())
