#!/usr/bin/env python
from __future__ import print_function  # noqa

import os
import sys
import sysconfig
from collections import defaultdict

sys.path.append(os.path.dirname(__file__))

try:
    import phaseshifts
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
try:
    import skbuild
except ImportError:
    skbuild = None  # type: ignore [assignment]

BUILD_BACKEND = None
# Build logic for legacy and modern Python
if tuple(sys.version_info[:2]) <= (3, 11):
    # Try modern build first: scikit-build + CMake
    try:
        import skbuild
        from skbuild import setup  # noqa: F811

        BUILD_BACKEND = "skbuild"
    except ImportError:
        print(
            "WARNING: scikit-build not found; falling back to legacy numpy.f2py build.",
            file=sys.stderr,
        )
        BUILD_BACKEND = "numpy.f2py"
    except Exception as e:
        print(
            f"WARNING: scikit-build build failed ({e}); falling back to legacy numpy.f2py build.",
            file=sys.stderr,
        )
        BUILD_BACKEND = "numpy.f2py"
    if BUILD_BACKEND == "numpy.f2py":
        from setuptools import find_packages, setup, Extension  # noqa: F811

        try:
            import numpy.f2py

            INCLUDE_DIRS += [numpy.get_include(), numpy.f2py.get_include()]
        except ImportError:
            print(
                "WARNING: numpy.f2py not found; Fortran extension will not be built.",
                file=sys.stderr,
            )
        # Optionally, run f2py manually if needed
        # subprocess.check_call([...])
else:
    # Modern build: require scikit-build and CMake
    try:
        from skbuild import setup

        BUILD_BACKEND = "skbuild"
    except ImportError:
        raise ImportError(
            "scikit-build is required for building phaseshifts on Python 3.12+. Please install scikit-build and CMake."
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
    CMAKE_ARGS = {
        "cmake_args": [
            '-DPYTHON_INCLUDE_DIR="{}"'.format(sysconfig.get_path("include")),
            '-DPYTHON_LIBRARY="{}"'.format(sysconfig.get_config_var("LIBDIR")),
            "-DENABLE_PHSHIFT2007_BINARIES={}".format(
                "ON" if BUILD_PHSHIFT2007 else "OFF"
            ),
        ]
    }

BUILD_EXT_INPLACE_ARGS = ["build_ext", "--inplace"]

# build f2py extensions
f2py_exts_sources = {
    "libphsh": [
        os.path.join(
            "phaseshifts",
            "lib",
            "libphsh" + (".f" if BUILD_BACKEND == "numpy.distutils" else "module.c"),
        ),
    ]
}
f2py_platform_extra_args = defaultdict(
    dict,
    {
        "darwin": {"extra_link_args": [], "extra_compile_args": []},
        "win32": {"extra_link_args": [], "extra_compile_args": []},
        "linux": {"extra_link_args": ["-lgomp"], "extra_compile_args": ["-fopenmp"]},
        "linux2": {"extra_link_args": ["-lgomp"], "extra_compile_args": ["-fopenmp"]},
    },
)[sys.platform]

f2py_exts = (
    [
        # NOTE: When hacking the build process for Python 3.12, we still want to force wheel to be platform specific
        #       See https://stackoverflow.com/a/53463910/22567758
        Extension(
            name="phaseshifts.lib._native_build",
            sources=[os.path.join("phaseshifts", "lib", "_native_build.c")],
        )
    ]
    if BUILD_BACKEND != "numpy.distutils"
    else [
        Extension(
            name="phaseshifts.lib.libphsh",
            extra_compile_args=f2py_platform_extra_args["extra_compile_args"],
            extra_link_args=f2py_platform_extra_args["extra_link_args"],
            sources=f2py_exts_sources["libphsh"],
        )
    ]
)

print("BUILD_BACKEND: {}".format(BUILD_BACKEND))

README = "README.md"

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
    long_description=(None if not os.path.exists(README) else open(README).read()),
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
        "atorb": ["mendeleev", "elementy"],
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
        "doc": ["sphinx>=7,<8", "sphinx_rtd_theme", "numpydoc", "ipykernel"],
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
