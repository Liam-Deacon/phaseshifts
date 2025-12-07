"""Locate the Python DLL on Windows runners.

Prints the absolute path to stdout. Exits non-zero if not found.
"""

from __future__ import annotations

import pathlib
import sys
import sysconfig


def find_python_dll() -> pathlib.Path:
    libname = (
        sysconfig.get_config_var("LDLIBRARY")
        or sysconfig.get_config_var("DLLLIBRARY")
        or ""
    )
    if not libname:
        raise RuntimeError("Python DLL name could not be determined from sysconfig")

    bindir = sysconfig.get_config_var("BINDIR") or sys.exec_prefix
    candidates = [
        pathlib.Path(bindir) / libname,
        pathlib.Path(sys.exec_prefix) / libname,
        pathlib.Path(sys.base_prefix) / libname,
    ]

    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()

    raise FileNotFoundError(f"Python DLL {libname!r} not found in {candidates}")


def main() -> int:
    path = find_python_dll()
    print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
