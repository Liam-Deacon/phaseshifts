"""The subpackage is for library code, likely compiled from Fortran or C sources."""

import importlib
import importlib.util
import os
import sys

try:
    from typing import List
except ImportError:
    pass  # No-op if typing is unavailable (e.g., Python <3.5)

from ._fortran_lib import (
    is_module_importable,
    compile_f2py_shared_library,
    find_shared_library_path,
    FORTRAN_LIBS,
    LIBPHSH_MODULE,
)

from .. import settings


def _add_windows_dll_dirs():
    """Ensure gfortran runtime DLLs are discoverable on Windows runners.

    GitHub-hosted Windows images install gfortran via MSYS/MinGW, leaving the
    runtime DLLs outside Python's default search path. Adding these directories
    prevents ``ImportError: DLL load failed`` when loading the extension.
    """

    added = []
    handles = []
    if sys.platform != "win32":
        return added

    # Python <3.8 lacks add_dll_directory; PATH is the best we can do there.
    add_dir = getattr(os, "add_dll_directory", None)
    if add_dir is None:
        return added

    dll_candidates = [
        "libgfortran-5.dll",
        "libgfortran-3.dll",
        "libgfortran.dll",
        "libgcc_s_seh-1.dll",
        "libquadmath-0.dll",
        "libwinpthread-1.dll",
    ]

    # Common MinGW/MSYS install locations on GitHub Windows runners
    default_dirs = [
        os.environ.get("MINGW_PREFIX"),
        r"C:\\msys64\\mingw64",
        r"C:\\msys64\\ucrt64",
        r"C:\\msys64\\mingw32",
        r"C:\\ProgramData\\chocolatey\\lib\\mingw\\tools\\install\\mingw64",
    ]
    default_dirs = [
        os.path.join(path, "bin") if path else None for path in default_dirs
    ]

    path_dirs = [p for p in os.environ.get("PATH", "").split(os.pathsep) if p]
    search_dirs = path_dirs + [d for d in default_dirs if d]

    seen = set()
    for directory in search_dirs:
        if not os.path.isdir(directory) or directory in seen:
            continue
        for dll in dll_candidates:
            if os.path.isfile(os.path.join(directory, dll)):
                try:
                    handle = add_dir(directory)
                    handles.append(handle)
                    added.append(directory)
                except OSError:
                    pass
                seen.add(directory)
                break

    # Keep handles alive at module scope to prevent DLL dirs from being removed
    if handles:
        globals()["_WIN_DLL_HANDLES"] = handles

    return added


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
    _WIN_DLL_DIRS = _add_windows_dll_dirs()
    libphsh = importlib.import_module(LIBPHSH_MODULE)
except ImportError:
    _lib_location = find_shared_library_path(LIBPHSH_MODULE)
    if _lib_location:
        libphsh = _load_extension_from_location(LIBPHSH_MODULE, _lib_location)
