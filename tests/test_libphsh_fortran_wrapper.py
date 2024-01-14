import sys
import os
import re

import pytest


def test_libphsh_exists():
    """
    GIVEN a package build including binaries of libphsh.f has been previously executed
    WHEN searching phaseshifts/lib directory
    THEN there should be a compiled shared object for libphsh (.pyd|.dll on Windows, .so otherwise)
    """
    import phaseshifts

    root_dir = os.path.dirname(phaseshifts.__file__)
    lib_dir = os.path.join(root_dir, "lib")
    assert os.path.exists(lib_dir)
    for directory, _, files in os.walk(lib_dir):
        if "__pycache__" in directory:
            continue

        if sys.platform == "win32":
            assert any(
                [
                    re.match(r"^libphsh.*\.(dll|pyd)$", filename, re.IGNORECASE)
                    for filename in files
                ]
            ), "Expected to find a file matching regex 'phaseshifts/lib/libphsh.*\.(dll|pyd)' (case insensitive)"
        else:
            assert any(
                [filename.startswith("libphsh") and filename.endswith(".so") for filename in files]
            ), "Expected to find file matching 'phaseshifts/lib/libphsh*.so' glob pattern"


def test_import_libphsh():
    """
    GIVEN a compiled libphsh.f shared object
    WHEN importing phaseshifts.lib.libphsh compiled wrapper extension
    THEN import should be successful
    """
    ext = ".so" if sys.platform != "win32" else ".pyd"
    try:
        import phaseshifts.lib.libphsh  # type: ignore [import-untyped] # noqa
    except ModuleNotFoundError:
        pytest.fail("libphsh*{} has not been compiled".format(ext))
    except ImportError as err:
        err_message = "{}: ".format(ext).join(str(err).split("{}: ".format(ext))[1:])
        pytest.fail(
            "Unable to import compiled libphsh due to: '{}'".format(err_message)
        )
