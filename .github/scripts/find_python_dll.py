"""Locate the Python DLL on Windows runners.

Prints the absolute path to stdout. Exits non-zero if not found.
"""

from __future__ import annotations

import pathlib
import sys
import sysconfig
from typing import List, Optional


def find_python_dll() -> pathlib.Path:
    """Find the Python DLL on Windows.

    Returns:
        pathlib.Path: The resolved path to the Python DLL.

    Raises:
        FileNotFoundError: If the Python DLL cannot be located.
    """
    libname: Optional[str] = (
        sysconfig.get_config_var("LDLIBRARY")
        or sysconfig.get_config_var("DLLLIBRARY")
        or sysconfig.get_config_var("PY3DLL")
        or ""
    )
    if not libname:
        libname = "python{}{}.dll".format(
            sys.version_info.major, sys.version_info.minor
        )

    bindir: str = sysconfig.get_config_var("BINDIR") or sys.exec_prefix
    candidates: List[pathlib.Path] = []
    roots = (
        bindir,
        sys.exec_prefix,
        sys.base_prefix,
        sys.prefix,
        pathlib.Path(sys.executable).resolve().parent,
    )
    for base in dict.fromkeys(str(r) for r in roots):
        basepath = pathlib.Path(base)
        candidates.append(basepath / libname)
        candidates.append(basepath / "DLLs" / libname)
        for pattern in (libname, "python*.dll"):
            candidates.extend(basepath.glob(pattern))
            if basepath.parent:
                candidates.extend(basepath.parent.glob(pattern))

    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()

    raise FileNotFoundError(
        "Python DLL {!r} not found in {}".format(libname, candidates)
    )


def main() -> int:
    path = find_python_dll()
    print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
