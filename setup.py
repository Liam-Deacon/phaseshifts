#!/usr/bin/env python
from __future__ import print_function  # noqa

import os
import sys
import sysconfig
import subprocess  # nosec
import shutil

try:
    from contextlib import suppress
except ImportError:  # TODO: Remove when we drop support for python 2.7
    # For Python <3.4, suppress is not available
    def suppress(*exceptions):  # type: ignore[no-redef]
        """Context manager that suppresses specified exceptions."""

        class SuppressContext:
            def __enter__(self):
                return self

            def __exit__(self, exc_type, exc_value, traceback):
                return isinstance(exc_value, exceptions)

        return SuppressContext()


with suppress(ImportError):
    from typing import Literal  # noqa: F401

sys.path.append(os.path.dirname(__file__))


# Ensure distutils is importable even on slim/modern Python installs where it
# is absent from the stdlib (e.g., Python 3.12+ or Debian/Ubuntu minimal
# builds). numpy.f2py and numpy.distutils expect ``distutils`` to be present.
try:
    from phaseshifts.lib._distutils_compat import ensure_distutils
except Exception:  # pragma: no cover - fallback for very early import failures

    def ensure_distutils():
        try:
            import distutils  # noqa: F401

            return True
        except ImportError:
            pass

        try:
            import setuptools._distutils as _du  # type: ignore
        except Exception:
            return False

        sys.modules.setdefault("distutils", _du)
        for _mod in (
            "distutils.ccompiler",
            "distutils.sysconfig",
            "distutils.unixccompiler",
            "distutils.msvccompiler",
        ):
            try:
                __import__(_mod)
            except Exception:
                pass

        return True


PROJECT_ROOT = os.path.abspath(os.path.dirname(__file__))
F2PY_SOURCE = "libphsh.f"
F2PY_SIGNATURE = "libphsh.pyf"
F2PY_SIGNATURE_PATH = os.path.join(PROJECT_ROOT, "phaseshifts", "lib", F2PY_SIGNATURE)


def _pythonpath_env(extra_path):
    env = os.environ.copy()
    paths = [extra_path]
    existing = env.get("PYTHONPATH")
    if existing:
        paths.append(existing)
    env["PYTHONPATH"] = os.pathsep.join(paths)
    return env


try:
    import phaseshifts
    import phaseshifts.phshift2007
except ModuleNotFoundError:
    phaseshifts = None  # type: ignore [assignment]

try:
    from setuptools import find_packages, setup, Extension  # type: ignore [import-untyped]
except ImportError:
    # distutils is deprecated/removed in Python 3.12+. Use setuptools only.
    try:
        from distutils.core import find_packages
    except ImportError:
        raise ImportError(
            "setuptools is required for building phaseshifts. Please install setuptools."
        )

INCLUDE_DIRS = []

# WARNING: numpy.distutils is completely removed in python 3.12 and deprecated for removal in python 3.11 by Oct 2025
# The project will therefore need to be migrated to use a different build backend, see
# https://numpy.org/doc/stable/reference/distutils_status_migration.html#distutils-status-migration
with suppress(ImportError):
    import skbuild  # noqa: F401

# Select build backend. Modern builds require scikit-build and CMake (default)
BUILD_BACKEND = "skbuild"  # type: Literal["skbuild", "numpy.distutils"]

BUILD_BACKEND = None  # will be set below based on available tooling


# Build logic for legacy and modern Python
if tuple(sys.version_info[:2]) <= (3, 11):
    ensure_distutils()
    # Try modern build first: scikit-build + CMake
    try:
        import skbuild
        from skbuild import setup  # noqa: F811

        BUILD_BACKEND = "skbuild"
    except ImportError:
        print(
            "WARNING: scikit-build not found; falling back to legacy numpy.distutils build.",
            file=sys.stderr,
        )
        BUILD_BACKEND = "numpy.distutils"
        from numpy.distutils.core import setup  # type: ignore  # noqa: F811
        from numpy.distutils.extension import Extension  # type: ignore  # noqa: F811
    except Exception as e:
        print(
            f"WARNING: scikit-build build failed ({e}); falling back to legacy numpy.distutils build.",
            file=sys.stderr,
        )
        BUILD_BACKEND = "numpy.distutils"
        from numpy.distutils.core import setup  # type: ignore  # noqa: F811
        from numpy.distutils.extension import Extension  # type: ignore  # noqa: F811
    if BUILD_BACKEND == "numpy.distutils":
        try:
            import numpy
            import numpy.f2py

            INCLUDE_DIRS += [numpy.get_include(), numpy.f2py.get_include()]
        except ImportError:
            print(
                "WARNING: numpy.f2py not found; Fortran extension will not be built.",
                file=sys.stderr,
            )
        if not any(x in sys.argv for x in ("sdist", "wheel")):
            args = [
                sys.executable,
                "-m",
                "phaseshifts.lib._f2py_shim",
                F2PY_SOURCE,
                "-m",
                "libphsh",
                "-c",
            ]
            try:
                if phaseshifts and hasattr(phaseshifts, "phshift2007"):
                    # Use the platform-filtered flags directly from COMPILER_FLAGS
                    fortran_flags = phaseshifts.phshift2007.COMPILER_FLAGS["gfortran"]
                    # Filter out flags that f2py/numpy.distutils might not handle well
                    f2py_safe_flags = [
                        flag
                        for flag in fortran_flags
                        if flag
                        not in phaseshifts.phshift2007.WINDOWS_INCOMPATIBLE_FLAGS
                    ]
                    if f2py_safe_flags:
                        fortran_flags_str = " ".join(f2py_safe_flags)
                        args.append(f"--f77flags={fortran_flags_str}")
                    else:
                        args.append("--f77flags=-frecursive")
                else:
                    # Default fallback
                    args.append("--f77flags=-frecursive")
            except (NameError, AttributeError, KeyError):
                # Fallback if phaseshifts module or compiler flags not available
                args.append("--f77flags=-frecursive")
            subprocess.check_call(
                args,
                cwd="./phaseshifts/lib",
                env=_pythonpath_env(PROJECT_ROOT),
            )  # nosec
else:
    # Modern build: require scikit-build and CMake for Python 3.12+
    # NOTE: numpy.distutils is completely removed in Python 3.12
    try:
        from skbuild import setup

        BUILD_BACKEND = "skbuild"
    except ImportError:
        # Fallback to setuptools for Python 3.12+ (numpy.distutils is gone)
        print(
            "WARNING: scikit-build not found; falling back to setuptools build.",
            file=sys.stderr,
        )
        print(
            "NOTE: For best results, install scikit-build-core and cmake.",
            file=sys.stderr,
        )
        BUILD_BACKEND = "setuptools"
        from setuptools import setup, Extension  # type: ignore  # noqa: F811

        ensure_distutils()

        # Try to pre-build the Fortran extension using f2py
        try:
            import numpy
            import numpy.f2py

            INCLUDE_DIRS += [numpy.get_include(), numpy.f2py.get_include()]

            # Build the extension if not doing sdist/wheel
            if not any(x in sys.argv for x in ("sdist", "wheel", "--version", "-V")):
                args = [
                    sys.executable,
                    "-m",
                    "phaseshifts.lib._f2py_shim",
                    F2PY_SOURCE,
                    "-m",
                    "libphsh",
                    "-c",
                    "--f77flags=-frecursive -fcheck=bounds -std=legacy",
                ]
                try:
                    print("Building libphsh extension with f2py...")
                    subprocess.check_call(
                        args,
                        cwd="./phaseshifts/lib",
                        env=_pythonpath_env(PROJECT_ROOT),
                    )  # nosec
                except subprocess.CalledProcessError as e:
                    print(
                        "WARNING: f2py build failed ({}); extension may not be available.".format(
                            e
                        ),
                        file=sys.stderr,
                    )
        except ImportError:
            print(
                "WARNING: numpy.f2py not found; Fortran extension will not be built.",
                file=sys.stderr,
            )

if len(sys.argv) == 1:
    sys.argv.append("install")

CMAKE_ARGS = {}

#: True-like command line flags
TRUE_OPTS = {"y", "yes", "on", "true", "1"}

# Read environment variable to control phshift2007 binary build
BUILD_PHSHIFT2007 = (
    os.environ.get("PHASESHIFTS_BUILD_PHSHIFT2007_BINARIES", "OFF").lower() in TRUE_OPTS
)

if BUILD_BACKEND == "skbuild":
    cmake_args = [
        '-DPYTHON_INCLUDE_DIR="{}"'.format(sysconfig.get_path("include")),
        '-DPYTHON_LIBRARY="{}"'.format(sysconfig.get_config_var("LIBDIR")),
        "-DENABLE_PHSHIFT2007_BINARIES={}".format("ON" if BUILD_PHSHIFT2007 else "OFF"),
    ]

    # Windows-specific CMake configuration
    if os.name == "nt":
        # Force MinGW Makefiles generator on Windows
        cmake_args.extend(
            [
                "-G",
                "MinGW Makefiles",
                "-DCMAKE_C_COMPILER=gcc",
                "-DCMAKE_Fortran_COMPILER=gfortran",
            ]
        )

        # Check if MinGW is available in common locations
        mingw_paths = [
            "C:/mingw64/bin/gfortran.exe",
            "C:/tools/mingw64/bin/gfortran.exe",
            "C:/msys64/mingw64/bin/gfortran.exe",
            "C:/ProgramData/chocolatey/lib/mingw/tools/install/mingw64/bin/gfortran.exe",
        ]

        gfortran_found = False
        for path in mingw_paths:
            if os.path.exists(path):
                cmake_args.extend(
                    [
                        f"-DCMAKE_Fortran_COMPILER={path}",
                        f"-DCMAKE_C_COMPILER={path.replace('gfortran.exe', 'gcc.exe')}",
                    ]
                )
                gfortran_found = True
                print(f"Found gfortran at: {path}")
                break

        if not gfortran_found:
            print(
                "WARNING: gfortran not found in common locations. "
                "Please ensure MinGW-w64 is installed and in PATH.",
                file=sys.stderr,
            )

    CMAKE_ARGS = {"cmake_args": cmake_args}

BUILD_EXT_INPLACE_ARGS = ["build_ext", "--inplace"]

# build f2py extensions
libphsh_sources = [
    os.path.join(
        "phaseshifts",
        "lib",
        F2PY_SOURCE,
    ),
]
if os.path.exists(F2PY_SIGNATURE_PATH):
    libphsh_sources.insert(
        0,
        os.path.join(
            "phaseshifts",
            "lib",
            F2PY_SIGNATURE,
        ),
    )
f2py_exts_sources = {"libphsh": libphsh_sources}
try:
    # Detect MSVC toolchain; avoid GNU Fortran flags when it is active
    is_msvc = os.name == "nt" and shutil.which("cl.exe") is not None
    if (not is_msvc) and phaseshifts and hasattr(phaseshifts, "phshift2007"):
        # Use the platform‑filtered flags directly from COMPILER_FLAGS
        gfortran_compiler_args = [
            flag
            for flag in phaseshifts.phshift2007.COMPILER_FLAGS["gfortran"]
            if flag not in phaseshifts.phshift2007.WINDOWS_INCOMPATIBLE_FLAGS
        ]
    else:
        gfortran_compiler_args = []
except (NameError, AttributeError, KeyError):
    gfortran_compiler_args = []

# Platform-specific extra arguments
f2py_platform_extra_args = {
    "darwin": {"extra_link_args": [], "extra_compile_args": []},
    "win32": {
        "extra_link_args": ([] if shutil.which("cl.exe") else gfortran_compiler_args),
        "extra_compile_args": (
            [] if shutil.which("cl.exe") else gfortran_compiler_args
        ),
    },
    "linux": {
        "extra_link_args": gfortran_compiler_args + ["-lgomp"],
        "extra_compile_args": gfortran_compiler_args + ["-fopenmp"],
    },
    "linux2": {
        "extra_link_args": gfortran_compiler_args + ["-lgomp"],
        "extra_compile_args": gfortran_compiler_args + ["-fopenmp"],
    },
}.get(sys.platform, {"extra_link_args": [], "extra_compile_args": []})

build_eeasisss = (
    os.environ.get("PHASESHIFTS_BUILD_EEASISSS", "OFF").lower() in TRUE_OPTS
)
eeasisss_ext = Extension(
    name="phaseshifts.lib.libhartfock",
    extra_compile_args=[],
    language="f90",
    sources=[os.path.join("phaseshifts", "lib", "EEASiSSS", "hf.f90")],
)

f2py_exts = (
    [
        # NOTE: When hacking the build process for Python 3.12, we still want to force wheel to be platform specific
        #       See https://stackoverflow.com/a/53463910/22567758
        Extension(
            name="phaseshifts.lib._native_build",
            sources=[os.path.join("phaseshifts", "lib", "_native_build.c")],
        )
    ]
    if BUILD_BACKEND in ("skbuild", "setuptools")
    else [
        Extension(
            name="phaseshifts.lib.libphsh",
            include_dirs=INCLUDE_DIRS,
            extra_compile_args=f2py_platform_extra_args["extra_compile_args"],
            extra_link_args=f2py_platform_extra_args["extra_link_args"],
            sources=f2py_exts_sources["libphsh"],
        ),
    ]
    + ([eeasisss_ext] if build_eeasisss else [])
)

print("BUILD_BACKEND: {}".format(BUILD_BACKEND))

README = "README.md"
LONG_DESCRIPTION = ""

try:
    with open(README, encoding="utf-8") as readme_file_ptr:
        LONG_DESCRIPTION = readme_file_ptr.read()
except OSError:
    pass

# --- Fallback logic for build ---
setup_args = dict(
    name="phaseshifts",
    packages=find_packages(),
    version=getattr(phaseshifts, "__version__", "0.1.8-dev"),
    author="Liam Deacon",
    author_email="liam.m.deacon@gmail.com",
    license="MIT License",
    url="https://github.com/Liam-Deacon/phaseshifts",
    description=(
        "Python-based version of the Barbieri/Van Hove phase shift calculation package for LEED/XPD modelling"
    ),
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Environment :: X11 Applications :: Qt",  # The end goal is is to have Qt or other GUI frontend
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    extras_require={
        "atorb": ["mendeleev", "elementy", "pydantic; python_version >= '3.9'"],
        "gui": [
            "six",
            "qtpy",
            "pyqt5; python_version < '3.7'",
            "pyside6; python_version >= '3.7'",
            "pyqtgraph",
        ],
        "dev": [
            "black",
            "isort",
            "numpy",
            "pre-commit",
            "ruff",
            "scikit-build-core; python_version > '3.11'",
            "wheel",
        ],
        "test": ["pytest", "pytest-cov"],
        "docs": ["sphinx>=7,<8", "sphinx_rtd_theme", "numpydoc", "ipykernel"],
        "input": ["pyyaml>=6.0", "jsonschema>=4.0"],
        "viperleed": ["viperleed>=0.14.1"],
    },
    keywords="phaseshifts atomic scattering muffin-tin diffraction",
    include_package_data=True,
    package_data=(
        {}
        if BUILD_BACKEND == "skbuild"
        else {
            # If any package contains *.f or *.pyd files, include them:
            "": ["*.f", "*.pyd", "*.so", "*.dll"],
            # If any package contains *.txt or *.rst files, include them:
            ".": ["*.txt", "*.rst", "*.pyw", "ChangeLog"],
            "lib": ["lib/*.f", "lib/*.c", "lib/*.h"] + ["*.so", "*.pyd"],
            "gui": ["gui/*.ui", "gui/*.bat"],
            "gui/res": ["gui/res/*.*"],
        }
    ),
    scripts=["phaseshifts/phsh.py"],
    install_requires=[
        "scipy >= 0.7",
        "numpy >= 1.3",
        "periodictable",
        "typing_extensions",
    ],
    ext_modules=f2py_exts,
    console=[os.path.join("phaseshifts", "phsh.py")],
)

build_failed = False
if BUILD_BACKEND == "skbuild":
    try:
        dist = setup(**setup_args, **CMAKE_ARGS)
    except Exception as e:
        print(f"WARNING: skbuild/CMake build failed ({e})", file=sys.stderr)
        build_failed = True
        # Fallback for Python <3.12
        if tuple(sys.version_info[:2]) <= (3, 11):
            print("Falling back to setuptools/numpy.f2py build.", file=sys.stderr)
            dist = setup(**setup_args)
        else:
            raise
else:
    dist = setup(**setup_args)
