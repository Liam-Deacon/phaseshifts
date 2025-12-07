"""Run ``numpy.f2py`` with a vendored ``distutils`` fallback."""

import sys

from ._distutils_compat import ensure_distutils


def main(argv=None):
    if not ensure_distutils():
        raise ImportError(
            "distutils is unavailable; install setuptools or enable build isolation."
        )

    from numpy.f2py.f2py2e import main as f2py_main

    args = ["f2py"] + list(argv if argv is not None else sys.argv[1:])
    sys.argv = args
    return f2py_main()


if __name__ == "__main__":
    raise SystemExit(main())
