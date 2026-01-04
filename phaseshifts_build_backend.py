from __future__ import absolute_import

import os
import sys

# Minimal PEP 517 backend shim to keep Python 2.7/3.x builds working.
# - Python <3.8: fall back to setuptools' legacy backend.
# - Python >=3.8: delegate to scikit-build-core for modern CMake builds.
_is_pyodide = os.environ.get("CIBW_PLATFORM") == "pyodide" or sys.platform == "emscripten"

if sys.version_info >= (3, 8) and not _is_pyodide:
    from scikit_build_core.build import *  # type: ignore  # noqa: F401,F403
else:
    from setuptools.build_meta import *  # type: ignore  # noqa: F401,F403
