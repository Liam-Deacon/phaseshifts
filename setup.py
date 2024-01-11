#!/usr/bin/env python
import os
import sys

try:
    from setuptools import find_packages
except ImportError:
    from distutils.core import find_packages

# WARNING: numpy.distutils is completely removed in python 3.12 and deprecated for removal in python 3.11 by Oct 2025
# The project will therefore need to be migrated to use a different build backend, see
# https://numpy.org/doc/stable/reference/distutils_status_migration.html#distutils-status-migration
try:
    from numpy.distutils.core import Extension, setup
except ModuleNotFoundError as err:
    if tuple(sys.version_info[:2]) >= (3, 11):
        raise NotImplementedError(
            "numpy.distutils has been removed for python {}".format(
                ".".join(sys.version_info[:2])
            )
        )
    raise

import sys, os

try:
    import py2exe
except ImportError:
    py2exe = None

if len(sys.argv) == 1:
    sys.argv.append("install")

# build f2py extensions
f2py_exts = [
    Extension(
        name="phaseshifts.lib.libphsh",
        extra_compile_args=["-fopenmp"],
        extra_link_args=["-lgomp"],
        sources=[os.path.join("phaseshifts", "lib", "libphsh.f")],
    )
]

dist = setup(
    name="phaseshifts",
    packages=find_packages(),
    version="0.1.6-dev",
    author="Liam Deacon",
    author_email="liam.m.deacon@gmail.com",
    license="MIT License",
    url="https://pypi.python.org/pypi/phaseshifts",
    description="Python-based version of the Barbieri/Van Hove phase "
    "shift calculation package for LEED/XPD modelling",
    long_description=None
    if not os.path.exists("README.rst")
    else open("README.rst").read(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Environment :: X11 Applications :: Qt",  # The end goal is to have Qt or other GUI frontend
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    keywords="phaseshifts atomic scattering muffin-tin diffraction",
    include_package_data=True,
    package_data={
        # If any package contains *.txt or *.rst files, include them:
        "": ["*.txt", "*.rst", "*.pyw", "ChangeLog"],
        "lib": ["lib/*.f", "lib/*.c", "lib/*.h"],
        "gui": ["gui/*.ui", "gui/*.bat"],
        "gui/res": ["gui/res/*.*"],
    },
    scripts=["phaseshifts/phsh.py"],
    install_requires=["scipy >= 0.7", "numpy >= 1.3", "periodictable"],
    ext_modules=f2py_exts,
    console=[os.path.join("phaseshifts", "phsh.py")],
    **(
        {
            "options": {
                "py2exe": {
                    "skip_archive": 1,
                    "compressed": 0,
                    "bundle_files": 2,
                    "dist_dir": os.path.join("dist", "py2exe"),
                    "excludes": ["tcl", "bz2"],
                    "dll_excludes": ["w9xpopen.exe", "tk85.dll", "tcl85.dll"],
                }
            }
        }
        if py2exe
        else {}
    ),
)

if len(sys.argv) < 2:
    # add build_ext and --inplace flags when performing: `python setup.py`
    sys.argv.append("build_ext")
    sys.argv.append("--inplace")
