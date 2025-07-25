import sys
import os
import re

import pytest

import phaseshifts

PHASESHIFTS_ROOT_DIR = os.path.dirname(phaseshifts.__file__)
PHASESHIFTS_LIB_DIR = os.path.join(PHASESHIFTS_ROOT_DIR, "lib")


def test_libphsh_exists():
    """
    GIVEN a package build including binaries of libphsh.f has been previously executed
    WHEN searching phaseshifts/lib directory
    THEN there should be a compiled shared object for libphsh (.pyd|.dll on Windows, .so otherwise)
    """
    assert os.path.exists(PHASESHIFTS_LIB_DIR), (
        "Directory does not exist: %s\nCurrent working directory: %s\nContents of parent: %s"
        % (
            PHASESHIFTS_LIB_DIR,
            os.getcwd(),
            os.listdir(os.path.dirname(PHASESHIFTS_LIB_DIR)),
        )
    )

    found = False
    details = []
    for directory, _, files in os.walk(PHASESHIFTS_LIB_DIR):
        if "__pycache__" in directory or not directory.endswith("/lib"):
            continue
        details.append("Checked directory: %s\nFiles: %s" % (directory, files))

        if sys.platform == "win32":
            found = any(
                re.match(r"^libphsh.*\.(dll|pyd)$", filename, re.IGNORECASE)
                for filename in files
            )
            expected = "libphsh*.dll or libphsh*.pyd"
        else:
            found = any(
                filename.startswith("libphsh") and filename.endswith(".so")
                for filename in files
            )
            expected = "libphsh*.so"

        if found:
            break

    assert found, (
        "Expected to find file matching '%s' in '%s'\nPlatform: %s\nDirectory walk details:\n%s"
        % (expected, PHASESHIFTS_LIB_DIR, sys.platform, "\n".join(details))
    )


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
        if sys.platform == "win32":
            _ = (
                os.system(
                    "dumpbin /dependents {}{}libphsh.pyd".format(
                        PHASESHIFTS_LIB_DIR, os.path.sep
                    )
                )
                == 0
            )
            if sys.version_info[:2] >= (3, 8):
                os.add_dll_directory(PHASESHIFTS_LIB_DIR)
                import phaseshifts.lib.libphsh  # type: ignore [import-untyped] # noqa

                return  # success after adding DLL directory to path
        err_message = "{}: ".format(ext).join(str(err).split("{}: ".format(ext))[1:])
        pytest.fail(
            "Unable to import compiled libphsh due to: '{}'".format(err_message)
        )
