from __future__ import absolute_import

import os
import sys

# Minimal PEP 517 backend shim to keep Python 2.7/3.x builds working.
# - Python <3.8: fall back to setuptools' legacy backend.
# - Python >=3.8: delegate to scikit-build-core for modern CMake builds.
_is_pyodide = os.environ.get("CIBW_PLATFORM") == "pyodide" or sys.platform == "emscripten"

if sys.version_info >= (3, 8) and not _is_pyodide:
    import scikit_build_core.build as _backend  # type: ignore
else:
    import setuptools.build_meta as _backend  # type: ignore

_BACKEND_EXPORTS = (
    "build_wheel",
    "build_sdist",
    "build_editable",
    "get_requires_for_build_wheel",
    "get_requires_for_build_sdist",
    "get_requires_for_build_editable",
    "prepare_metadata_for_build_wheel",
    "prepare_metadata_for_build_editable",
)

__all__ = [name for name in _BACKEND_EXPORTS if hasattr(_backend, name)]

for _name in __all__:
    globals()[_name] = getattr(_backend, _name)
