import sys
import os

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
    shared_obj_candidates = []
    for directory, _, files in os.walk(PHASESHIFTS_LIB_DIR):
        if "__pycache__" in directory or not directory.endswith("/lib"):
            continue
        details.append("Checked directory: %s\nFiles: %s" % (directory, files))

        for filename in files:
            # Accept libphsh, libphsh.so, libphsh.pyd, libphsh.dll
            if (
                filename == "libphsh"
                or filename.startswith("libphsh")
                and (
                    filename.endswith(".so")
                    or filename.endswith(".pyd")
                    or filename.endswith(".dll")
                )
            ):
                shared_obj_candidates.append(os.path.join(directory, filename))

    # Check file type of candidates
    def describe_file_type(filepath):
        if not os.path.exists(filepath):
            return "File does not exist: %s" % filepath
        with open(filepath, "rb") as f:
            header = f.read(4)
        if header.startswith(b"\x7fELF"):
            return "ELF shared object (Linux)"
        elif header[:2] == b"MZ":
            return "Windows DLL"
        elif header[:4] == b"\xcf\xfa\xed\xfe" or header[:4] == b"\xfe\xed\xfa\xcf":
            return "Mach-O (macOS)"
        else:
            return "Unknown file type: %s" % header

    found = False
    file_type_details = []
    for candidate in shared_obj_candidates:
        file_type = describe_file_type(candidate)
        file_type_details.append("%s: %s" % (candidate, file_type))
        if (
            file_type.startswith("ELF shared object")
            or file_type.startswith("Mach-O")
            or file_type.startswith("Windows DLL")
        ):
            found = True
            break

    assert found, (
        "Expected to find a valid libphsh shared object in '%s'\nPlatform: %s\nCandidates: %s\nFile type details:\n%s\nDirectory walk details:\n%s"
        % (
            PHASESHIFTS_LIB_DIR,
            sys.platform,
            shared_obj_candidates,
            "\n".join(file_type_details),
            "\n".join(details),
        )
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
        assert False, "libphsh*{} has not been compiled".format(ext)
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
        assert False, "Unable to import compiled libphsh due to: '{}'".format(
            err_message
        )


def test_libphsh_hb_invocation():
    """
    GIVEN a successfully imported libphsh extension
    WHEN calling a simple Fortran function (hb)
    THEN the function should execute and return a float
    """
    try:
        import phaseshifts.lib.libphsh as libphsh
    except Exception as err:
        assert False, "Could not import libphsh: %s" % str(err)
    try:
        result = libphsh.hb(2.5, 1.0)
    except Exception as err:
        assert False, "Could not call libphsh.hb(2.5, 1.0): %s" % str(err)
    assert isinstance(
        result, float
    ), "Expected float result from libphsh.hb, got %s" % type(result)
