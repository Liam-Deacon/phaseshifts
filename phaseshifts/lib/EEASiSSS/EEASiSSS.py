from __future__ import absolute_import, division, print_function, unicode_literals

import argparse
import os
import sys
from ctypes import cdll, create_string_buffer
from ctypes.util import find_library

_ext = ".dll" if str(sys.platform).startswith("win") else ".so"
_lib = os.path.join(os.path.dirname(__file__), "lib")
_library_path = find_library("EEASiSSS") or os.path.join(_lib, "libEEASiSSS" + _ext)
_LIB_HANDLE = None


def _load_library():
    global _LIB_HANDLE
    if _LIB_HANDLE is not None:
        return _LIB_HANDLE
    os.environ["PATH"] = _lib + os.pathsep + os.environ.get("PATH", "")
    try:
        _LIB_HANDLE = cdll.LoadLibrary(_library_path)
    except OSError as err:
        raise ImportError(
            "EEASiSSS library not found. Expected at: {} ({})".format(
                _library_path, err
            )
        ) from err
    return _LIB_HANDLE


def eeasisss(input_file="inputX"):
    """Call the EEASiSSS Fortran library using ctypes."""
    lib_eeasisss = _load_library()
    input_bytes = str(input_file).encode("utf-8")
    lib_eeasisss.hartfock_(create_string_buffer(input_bytes, 255))


def main(argv=None):
    """CLI entry point for eeasisss."""
    parser = argparse.ArgumentParser(
        description="Call the EEASiSSS Fortran library using ctypes.",
    )
    parser.add_argument(
        "-i",
        "--input",
        dest="input_file",
        default="inputX",
        help="Input file passed to the EEASiSSS library",
    )
    args = parser.parse_args(argv)
    if not os.path.isfile(args.input_file):
        sys.stderr.write("Error: Input file '{}' not found\n".format(args.input_file))
        return 1
    try:
        eeasisss(args.input_file)
    except Exception as err:  # pylint: disable=broad-except
        sys.stderr.write("Error calling EEASiSSS: {}\n".format(err))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
