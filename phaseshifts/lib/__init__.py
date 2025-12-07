"""The subpackage is for library code, likely compiled from Fortran or C sources."""

import importlib
import importlib.util
import os
import sys

from ._fortran_lib import (
    is_module_importable,
    compile_f2py_shared_library,
    find_shared_library_path,
    FORTRAN_LIBS,
    LIBPHSH_MODULE,
)

from .. import settings

# NOTE: Attempt to compile libphsh library on import if not already available, useful when installing from source dist
if (
    settings.COMPILE_MISSING
    and not os.environ.get("PHASESHIFTS_SKIP_COMPILE_ON_IMPORT")
    and not is_module_importable(LIBPHSH_MODULE)
):
    try:
        compile_f2py_shared_library(**FORTRAN_LIBS[LIBPHSH_MODULE])
    except Exception:  # pragma: no cover - best-effort fallback
        # Avoid failing import just because the optional auto-compile failed.
        # Callers can inspect/log separately if needed.
        pass


def _load_extension_from_location(module_name, location):
    """Load a compiled extension directly from ``location`` if present.

    This supports editable installs where the Python package is imported from
    the source tree but the compiled extension lives in site-packages.
    """

    try:
        spec = importlib.util.spec_from_file_location(module_name, location)
        if spec and spec.loader:
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            sys.modules[module_name] = module
            return module
    except Exception:
        return None

    return None


# Expose libphsh when it exists, even if it was installed outside the source tree.
try:
    libphsh = importlib.import_module(LIBPHSH_MODULE)
except ImportError:
    _lib_location = find_shared_library_path(LIBPHSH_MODULE)
    if _lib_location:
        libphsh = _load_extension_from_location(LIBPHSH_MODULE, _lib_location)
