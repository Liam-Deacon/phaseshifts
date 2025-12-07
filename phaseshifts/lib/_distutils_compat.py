"""Compatibility helpers for environments where ``distutils`` is missing."""

import sys


def ensure_distutils():
    """Ensure a ``distutils`` implementation is importable."""
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
