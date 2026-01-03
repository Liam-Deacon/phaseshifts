import sys
import os
import importlib.util

import phaseshifts

PHASESHIFTS_ROOT_DIR = os.path.dirname(phaseshifts.__file__)
PHASESHIFTS_LIB_DIR = os.path.join(PHASESHIFTS_ROOT_DIR, "lib")


def _find_libphsh_location():
    """Find the actual location of the libphsh extension module."""
    # Try to find via importlib (works for both editable and regular installs)
    try:
        spec = importlib.util.find_spec("phaseshifts.lib.libphsh")
        if spec and spec.origin:
            return spec.origin
    except (ImportError, ModuleNotFoundError):
        pass

    # Fallback: search in common locations
    search_dirs = [PHASESHIFTS_LIB_DIR]

    # Also check site-packages for editable installs with redirect mode
    import site

    for site_dir in site.getsitepackages() + [site.getusersitepackages()]:
        if site_dir:
            candidate = os.path.join(site_dir, "phaseshifts", "lib")
            if os.path.isdir(candidate):
                search_dirs.append(candidate)

    for directory in search_dirs:
        if not os.path.isdir(directory):
            continue
        for filename in os.listdir(directory):
            if filename.startswith("libphsh") and (
                filename.endswith(".so")
                or filename.endswith(".pyd")
                or filename.endswith(".dll")
                or ".cpython-" in filename
            ):
                return os.path.join(directory, filename)

    return None


def test_libphsh_exists():
    """
    GIVEN a package build including binaries of libphsh.f has been previously executed
    WHEN searching for the libphsh extension module
    THEN the compiled shared object should be findable (via importlib or filesystem)
    """
    location = _find_libphsh_location()

    if location is None:
        # Gather diagnostic info
        details = []
        details.append("PHASESHIFTS_LIB_DIR: %s" % PHASESHIFTS_LIB_DIR)
        details.append("phaseshifts.__file__: %s" % phaseshifts.__file__)

        if os.path.isdir(PHASESHIFTS_LIB_DIR):
            details.append("Contents of %s: %s" % (PHASESHIFTS_LIB_DIR, os.listdir(PHASESHIFTS_LIB_DIR)))
        else:
            details.append("Directory does not exist: %s" % PHASESHIFTS_LIB_DIR)

        # Check site-packages
        import site

        for site_dir in site.getsitepackages():
            lib_dir = os.path.join(site_dir, "phaseshifts", "lib")
            if os.path.isdir(lib_dir):
                details.append("Contents of %s: %s" % (lib_dir, os.listdir(lib_dir)))

        assert False, (
            "Expected to find a valid libphsh shared object.\n"
            "Platform: %s\n"
            "Diagnostics:\n%s" % (sys.platform, "\n".join(details))
        )

    # Verify it's a valid binary
    with open(location, "rb") as f:
        header = f.read(4)

    valid_binary = (
        header.startswith(b"\x7fELF")  # Linux ELF
        or header[:2] == b"MZ"  # Windows PE
        or header[:4]
        in (
            b"\xcf\xfa\xed\xfe",
            b"\xfe\xed\xfa\xcf",
            b"\xca\xfe\xba\xbe",
        )  # macOS Mach-O
    )

    assert valid_binary, "Found libphsh at %s but it doesn't appear to be a valid binary. Header: %s" % (
        location,
        header,
    )


def test_import_libphsh():
    """
    GIVEN a compiled libphsh.f shared object
    WHEN importing phaseshifts.lib.libphsh compiled wrapper extension
    THEN import should be successful
    """
    ext = ".so" if sys.platform != "win32" else ".pyd"

    # Find where the extension actually is
    location = _find_libphsh_location()

    # On Windows 3.8+, we may need to add DLL search directories
    if sys.platform == "win32" and sys.version_info[:2] >= (3, 8) and location:
        lib_dir = os.path.dirname(location)
        try:
            os.add_dll_directory(lib_dir)
        except OSError:
            pass  # Directory might not exist or already added

    try:
        from phaseshifts.lib import libphsh as _libphsh  # noqa: F401
    except ModuleNotFoundError:
        assert False, "libphsh*{} has not been compiled.\n" "Extension location search result: {}".format(ext, location)
    except ImportError as err:
        # Provide detailed diagnostics
        details = [
            "Platform: %s" % sys.platform,
            "Python version: %s" % sys.version,
            "Extension location: %s" % location,
            "Error: %s" % str(err),
        ]

        if sys.platform == "win32" and location:
            # Try to get DLL dependencies info
            import subprocess

            try:
                result = subprocess.run(
                    ["dumpbin", "/dependents", location],
                    capture_output=True,
                    text=True,
                    timeout=10,
                    check=False,
                )
                details.append("dumpbin /dependents output:\n%s" % result.stdout)
            except (FileNotFoundError, subprocess.TimeoutExpired):
                details.append("dumpbin not available or timed out")

            # Check PATH for common MinGW DLLs
            mingw_dlls = [
                "libgfortran-5.dll",
                "libgcc_s_seh-1.dll",
                "libquadmath-0.dll",
                "libwinpthread-1.dll",
            ]
            path_dirs = os.environ.get("PATH", "").split(os.pathsep)
            dll_found = {}
            for dll in mingw_dlls:
                for path_dir in path_dirs:
                    if os.path.isfile(os.path.join(path_dir, dll)):
                        dll_found[dll] = os.path.join(path_dir, dll)
                        break
                else:
                    dll_found[dll] = "NOT FOUND"
            details.append(
                "MinGW DLL search results:\n%s" % "\n".join("  %s: %s" % (k, v) for k, v in dll_found.items())
            )

        assert False, "Unable to import compiled libphsh due to ImportError.\n" "Diagnostics:\n%s" % "\n".join(details)


def test_libphsh_hb_invocation():
    """
    GIVEN a successfully imported libphsh extension
    WHEN calling a simple Fortran function (hb)
    THEN the function should execute and return a float
    """
    try:
        from phaseshifts.lib import libphsh
    except (ImportError, ModuleNotFoundError) as err:
        assert False, "Could not import libphsh: %s" % str(err)
    try:
        result = libphsh.hb(2.5, 1.0)
    except (TypeError, AttributeError, RuntimeError) as err:
        assert False, "Could not call libphsh.hb(2.5, 1.0): %s" % str(err)
    assert isinstance(result, float), "Expected float result from libphsh.hb, got %s" % type(result)
