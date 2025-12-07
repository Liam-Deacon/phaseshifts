"""Compatibility helpers for environments where ``distutils`` is missing."""

import sys
import types


def _create_stub_msvccompiler():
    """Create a stub msvccompiler module for non-Windows platforms.

    numpy.distutils unconditionally imports distutils.msvccompiler even on Linux
    (in mingw32ccompiler.py). On Python 3.12+ where distutils is gone and
    setuptools' vendored version doesn't include msvccompiler on non-Windows,
    we need to provide a stub.
    """
    stub = types.ModuleType("distutils.msvccompiler")
    stub.__doc__ = "Stub msvccompiler for non-Windows platforms"

    # Provide the commonly imported function that numpy.distutils expects
    def get_build_version():
        return None

    stub.get_build_version = get_build_version
    stub.MSVCCompiler = None  # Placeholder class reference
    return stub


def ensure_distutils():
    """Ensure a ``distutils`` implementation is importable."""
    try:
        import setuptools._distutils as _du  # type: ignore
    except Exception:
        try:  # stdlib fallback (Py <3.12)
            import distutils as _du  # type: ignore
        except Exception:
            return False

    # Register the base distutils module first
    sys.modules.setdefault("distutils", _du)

    # Preload critical submodules and force-map them to the vendored copies so
    # imports like ``distutils.msvccompiler`` succeed everywhere (numpy.f2py
    # triggers this even on Linux when probing compilers).
    import importlib

    submodules = ["ccompiler", "sysconfig", "unixccompiler", "msvccompiler"]

    for _name in submodules:
        dist_name = f"distutils.{_name}"

        # Skip if already registered
        if dist_name in sys.modules:
            continue

        vend_name = f"setuptools._distutils.{_name}"
        module = None

        for candidate in (vend_name, dist_name):
            try:
                module = importlib.import_module(candidate)
                break
            except Exception:
                continue

        # If msvccompiler couldn't be found (common on Linux with Python 3.12+),
        # create a stub module to prevent import errors in numpy.distutils
        if module is None:
            if _name == "msvccompiler":
                module = _create_stub_msvccompiler()
            else:
                continue

        sys.modules[dist_name] = module

    return True
