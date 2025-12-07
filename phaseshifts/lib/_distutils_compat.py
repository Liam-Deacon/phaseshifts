"""Compatibility helpers for environments where ``distutils`` is missing."""

import sys


def ensure_distutils():
    """Ensure a ``distutils`` implementation is importable."""
    try:  # stdlib distutils (Py <3.12)
        import distutils  # noqa: F401

        try:
            __import__("distutils.msvccompiler")
        except Exception:
            # Continue to vendored fallback to ensure Windows submodules exist
            raise ImportError
        return True
    except ImportError:
        pass

    try:  # vendored fallback (Py >=3.12 or slim envs)
        import setuptools._distutils as _du  # type: ignore
    except Exception:
        return False

    sys.modules.setdefault("distutils", _du)
    # Preload critical submodules and force-map them to the vendored copies so
    # imports like ``distutils.msvccompiler`` succeed (numpy.f2py triggers
    # this on Windows).
    for _name in ("ccompiler", "sysconfig", "unixccompiler", "msvccompiler"):
        dist_name = f"distutils.{_name}"
        vend_name = f"setuptools._distutils.{_name}"
        m = None
        for candidate in (dist_name, vend_name):
            try:
                m = __import__(candidate)
                break
            except Exception:
                continue
        if m is None:
            continue
        # If the import came from setuptools._distutils, re-register it under
        # the distutils.* namespace so downstream imports resolve.
        sys.modules.setdefault(dist_name, m)

    return True
