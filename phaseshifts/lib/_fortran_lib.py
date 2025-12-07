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

from ._distutils_compat import ensure_distutils


LIBPHSH_MODULE = "phaseshifts.lib.libphsh"

FORTRAN_LIBS = {
    LIBPHSH_MODULE: {
        "source": "libphsh.f",
        "module_name": LIBPHSH_MODULE.split(".")[-1],
    },
}


def is_module_importable(module):  # (str) -> bool
    """Determine whether `module` is importable."""
    try:
        is_importable = bool(importlib.import_module(module))
    except ImportError:
        is_importable = False
    return is_importable


def compile_f2py_shared_library(source, module_name=None, cwd=None, **f2py_kwargs):
    # type: (str, str|None, str|None, str) -> int
    """Compile `source` into f2py wrapped shared library given by `module_name`.

    Examples
    --------
    >>> compile_f2py_shared_library(**FORTRAN_LIBS[LIBPHSH_MODULE])

    See Also
    --------
    numpy.f2py
    """
    # numpy.f2py relies on distutils; on Python 3.12+ ensure a vendored copy is
    # available so the subprocess does not fail immediately.
    if not ensure_distutils():
        raise ImportError(
            "distutils is unavailable; install setuptools or enable build isolation."
        )

    f2py_args = [
        source,
        "-m",
        module_name or os.path.basename(source),
        "-c",
    ]
    f2py_args += ["--{}={!r}".format(key, val) for key, val in f2py_kwargs.items()]

    # Ensure the subprocess can import the package (for the shim) even when the
    # current working directory is inside phaseshifts/lib.
    env = os.environ.copy()
    pkg_parent = os.path.abspath(
        os.path.join(os.path.dirname(__file__), os.pardir, os.pardir)
    )
    pythonpath_parts = [pkg_parent]
    if env.get("PYTHONPATH"):
        pythonpath_parts.append(env["PYTHONPATH"])
    env["PYTHONPATH"] = os.pathsep.join(pythonpath_parts)

    return subprocess.check_call(
        [sys.executable, "-m", "phaseshifts.lib._f2py_shim"] + f2py_args,
        cwd=cwd or os.path.dirname(__file__),
        env=env,
    )
