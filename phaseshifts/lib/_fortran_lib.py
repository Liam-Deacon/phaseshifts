"""This module can be used for checking whether Fortran wrapped libraries such as `libphsh` hve compiled successfully.

Notes
-----

The functionality currently targets `libphsh.f` by default and provides utilities to check for the
compiled wrapped `libphsh` fortran shared library and adds the ability to crudely compile
on the fly (useful when installing this package from source distribution).
"""

import importlib
import os
import site
import subprocess
import sys

from ._distutils_compat import ensure_distutils


LIBPHSH_MODULE = "phaseshifts.lib.libphsh"

_F2PY_SOURCE = "libphsh.f"
_F2PY_SIGNATURE = os.path.join(os.path.dirname(__file__), "libphsh.pyf")
if not os.path.exists(_F2PY_SIGNATURE):
    _F2PY_SIGNATURE = None

FORTRAN_LIBS = {
    LIBPHSH_MODULE: {
        "source": _F2PY_SOURCE,
        "signature": _F2PY_SIGNATURE,
        "module_name": LIBPHSH_MODULE.split(".")[-1],
    },
}


def find_shared_library_path(module_name):  # (str) -> str|None
    """Locate a compiled extension for ``module_name`` across common paths.

    In editable installs the pure-Python package is often imported from the
    source tree while the compiled extension lands in site-packages. The
    standard importer then misses the binary. This helper searches sys.path and
    site-packages for any matching shared object.
    """

    try:
        spec = importlib.util.find_spec(module_name)
        if spec and spec.origin and spec.origin != "built-in":
            return spec.origin
    except ImportError:
        pass

    module_basename = module_name.split(".")[-1]
    parent_parts = module_name.split(".")[:-1]
    parent_rel = os.path.join(*parent_parts)

    search_roots = []
    for path_entry in sys.path:
        candidate = os.path.join(path_entry or os.getcwd(), parent_rel)
        if os.path.isdir(candidate):
            search_roots.append(candidate)

    for site_dir in site.getsitepackages() + [site.getusersitepackages()]:
        if not site_dir:
            continue
        candidate = os.path.join(site_dir, parent_rel)
        if os.path.isdir(candidate):
            search_roots.append(candidate)

    for directory in search_roots:
        for filename in os.listdir(directory):
            if filename.startswith(module_basename) and (
                filename.endswith((".so", ".pyd", ".dll")) or ".cpython-" in filename
            ):
                return os.path.join(directory, filename)

    return None


def is_module_importable(module):  # (str) -> bool
    """Determine whether `module` is importable."""
    try:
        is_importable = bool(importlib.import_module(module))
    except ImportError:
        is_importable = False
    return is_importable


def print_package_tree(root_dir=None, prefix=""):
    """Recursively print the directory tree of the top-level phaseshifts package."""
    if root_dir is None:
        root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    entries = sorted(os.listdir(root_dir))
    print(prefix + os.path.basename(root_dir) + "/")
    for i, entry in enumerate(entries):
        path = os.path.join(root_dir, entry)
        is_last = i == len(entries) - 1
        connector = "`-- " if is_last else "|-- "
        print(prefix + connector + entry)
        if os.path.isdir(path) and not entry.startswith("."):
            extension = "    " if is_last else "|   "
            print_package_tree(path + "/", prefix + extension)


def compile_f2py_shared_library(
    source,
    module_name=None,
    cwd=None,
    signature=None,
    **f2py_kwargs,
):
    # type: (str, str|None, str|None, str|None, str) -> int
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
        "-m",
        module_name or os.path.basename(source),
        "-c",
    ]
    f2py_args += ["--{}={!r}".format(key, val) for key, val in f2py_kwargs.items()]
    if signature:
        f2py_args.append(signature)
    f2py_args.append(source)

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
