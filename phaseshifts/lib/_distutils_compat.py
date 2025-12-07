"""Compatibility helpers for environments where ``distutils`` is missing."""

import sys


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

    for _name in ("ccompiler", "sysconfig", "unixccompiler", "msvccompiler"):
        dist_name = f"distutils.{_name}"
        vend_name = f"setuptools._distutils.{_name}"
        module = None
        for candidate in (vend_name, dist_name):
            try:
                module = importlib.import_module(candidate)
                break
            except Exception:
                continue
        if module is None:
            continue
        sys.modules[dist_name] = module

    return True
