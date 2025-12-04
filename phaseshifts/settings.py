"""Module for phaseshifts config and settings, typically defined by environment variables.

Notes
-----
Environment variables will be prefixed with :code:`PHASESHIFTS_` to namespace them and
avoid conflicts with other packages.

"""

# cspell:ignore TESTRUN

import os
import json


def load_bool_env_var(var_name, default="0"):  # type: (str, str) -> bool
    """Load boolean value from environment variable.

    Parameters
    ----------
    env_var : str
        Environment variable name.
    default : str
        The default value if environment variable is not set.

    Returns
    -------
    bool
        Boolean value of environment variable.

    """
    return bool(json.loads((os.environ.get(var_name) or default).lower()))


#: Internal readthedocs check to help determine execution environment
_READTHEDOCS = load_bool_env_var("READTHEDOCS")

#: Whether to compile missing libraries on first import
# Default to on so legacy CI/users still get a compiled libphsh even when wheels
# omit it; set PHASESHIFTS_COMPILE_MISSING=0 to disable.
COMPILE_MISSING = load_bool_env_var("PHASESHIFTS_COMPILE_MISSING", "1")

#: Whether to show debug messages
DEBUG = load_bool_env_var("PHASESHIFTS_DEBUG")

#: Controls whether to do a test execution run of phsh.py
TESTRUN = load_bool_env_var("PHASESHIFTS_TESTRUN")

#: Controls whether to profile phsh.py using cProfile when executing to understand performance
PROFILE = load_bool_env_var("PHASESHIFTS_PROFILE")
