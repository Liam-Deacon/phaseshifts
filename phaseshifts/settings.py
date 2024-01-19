"""Module for phaseshifts config and settings, typically defined by environment variables.

Notes
-----
Environment variables will be prefixed with :code:`PHASESHIFTS_` to namespace them and
avoid conflicts with other packages.

"""

import os
import json

#: Whether to compile missing libraries on first import
COMPILE_MISSING = bool(json.loads((os.environ.get("PHASESHIFTS_COMPILE_MISSING") or "1").lower()))
