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


def _ensure_distutils():
    """Ensure a ``distutils`` implementation is importable.

    Python 3.12+ removes ``distutils`` from the stdlib, but ``numpy.f2py`` still
    expects it. When missing, we fall back to the vendored copy shipped with
    ``setuptools`` so that on-the-fly builds keep working in CI.
    """

    try:  # stdlib distutils (Py <3.12)
        import distutils  # noqa: F401

        return True
    except ImportError:
        pass

    try:  # vendored fallback (Py >=3.12 or slim envs)
        import setuptools._distutils as _du  # type: ignore
    except Exception:
        return False

    sys.modules.setdefault("distutils", _du)
    # Preload the submodules most commonly needed by numpy.f2py
    for _mod in (
        "distutils.ccompiler",
        "distutils.sysconfig",
        "distutils.unixccompiler",
        "distutils.msvccompiler",
    ):
        try:
            __import__(_mod)
        except Exception:
            pass

    return True

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
    if not _ensure_distutils():
        raise ImportError(
            "distutils is unavailable; install setuptools or enable build isolation."
        )

    args = [
        sys.executable,
        "-m",
        "numpy.f2py",
        source,
        "-m",
        module_name or os.path.basename(source),
        "-c",
    ]
    args += ["--{}={!r}".format(key, val) for key, val in f2py_kwargs.items()]
    return subprocess.check_call(
        args,
        cwd=cwd or os.path.dirname(__file__),
    )
