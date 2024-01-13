"""This module can be used for checking whether Fortran wrapped libraries such as `libphsh` hve compiled successfully.

Notes
-----

The functionality currently targets `libphsh.f` by default and provides utilities to check for the
compiled wrapped `libphsh` fortran shared library and adds the ability to crudely compile
on the fly (useful when installing this package from source distribution).
"""
import importlib
import os
import subprocess
import sys
from typing import Optional

LIBPHSH_MODULE = "phaseshifts.lib.libphsh"

FORTRAN_LIBS = {
    LIBPHSH_MODULE: {"source": "libphsh.f", "module_name": LIBPHSH_MODULE},
}


def is_module_importable(module: str) -> bool:
    """Determine whether `module` is importable."""
    try:
        is_importable = bool(importlib.import_module(module))
    except ModuleNotFoundError:
        is_importable = False
    return is_importable


def compile_f2py_shared_library(
    source: str,
    module_name: Optional[str] = None,
    **f2py_kwargs,
) -> int:
    """Compile `source` into f2py wrapped shared library given by `module_name`.

    Examples
    --------
    >>> compile_f2py_shared_library(**FORTRAN_LIBS[LIBPHSH_MODULE])

    See Also
    --------
    numpy.f2py
    """
    return subprocess.check_call(
        [
            sys.executable,
            "-m",
            "numpy.f2py",
            source,
            "-m",
            module_name or os.path.basename(source),
            "-c",
            *[f"--{key}={val!r}" for key, val in f2py_kwargs.items()],
        ],
        cwd=os.path.dirname(__file__),
    )
