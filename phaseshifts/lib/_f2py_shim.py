"""Run ``numpy.f2py`` with a vendored ``distutils`` fallback."""

import sys

from ._distutils_compat import ensure_distutils

# Ensure distutils (with stubs) is available BEFORE any numpy imports.
# This must happen at module load time because numpy.f2py imports
# distutils.msvccompiler unconditionally even on Linux.
if not ensure_distutils():
    raise ImportError("distutils is unavailable; install setuptools or enable build isolation.")


def main(argv=None):
    from numpy.f2py.f2py2e import main as f2py_main

    args = ["f2py"] + list(argv if argv is not None else sys.argv[1:])
    sys.argv = args
    return f2py_main()


if __name__ == "__main__":
    raise SystemExit(main())
