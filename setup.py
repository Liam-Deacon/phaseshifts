#!/usr/bin/env python
"""Setup script for the phaseshifts package installation."""
from __future__ import print_function  # noqa

import os
import sys

from setuptools import find_packages, setup, Extension  # type: ignore [import-untyped]

sys.path.append(os.path.dirname(__file__))

# pylint: disable=wrong-import-position

import phaseshifts  # noqa
import phaseshifts.build.cmake  # noqa

# WARNING: numpy.distutils is completely removed in python 3.12 and deprecated for removal in python 3.11 by Oct 2025
# The project will therefore need to be migrated to use a different build backend, see
# https://numpy.org/doc/stable/reference/distutils_status_migration.html#distutils-status-migration


README = "README.md"
try:
    with open(README, mode="r") as fh:
        LONG_DESCRIPTION = fh.read()
except IOError:
    LONG_DESCRIPTION = ""


setup(
    name="phaseshifts",
    packages=find_packages(),
    version=getattr(phaseshifts, "__version__", "0.2.0-dev"),
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
        "Development Status :: 3 - Alpha",  # NOTE: Move out of Alpha status once API is stable and well tested
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
        "gui": [],
        "dev": [
            "black",
            "cmake>=3.27",
            "isort",
            "numpy",
            "pre-commit",
            "ruff",
            "wheel",
        ],
        "test": [
            "cmake>=3.27",
            "mock==3.0.5; python_version < '3.3'",
            "pytest",
            "pytest-cov",
        ],
        "doc": ["sphinx>=7,<8", "sphinx_rtd_theme", "numpydoc", "ipykernel"],
    },
    keywords="phaseshifts atomic scattering muffin-tin diffraction",
    include_package_data=True,
    cmdclass={"build_ext": phaseshifts.build.cmake.CMakeBuild},
    data_files=[
        ("bin", ["bin/phsh0", "bin/phsh1", "bin/phsh2cav", "bin/phsh2rel", "bin/phsh2wil", "bin/phsh3"]),
    ],
    package_data={
        "phaseshifts": ["**/*.pyi", "**/*.dll", "**/*.dylib", "**/*.so"],
        "bin": ["bin/*"],
        # If any package contains *.f or *.pyd files, include them:
        "": ["*.f", "*.pyd", "*.so", "*.dll"],
        # If any package contains *.txt or *.rst files, include them:
        ".": ["*.txt", "*.rst", "*.pyw", "ChangeLog"],
        "lib": ["lib/*.f", "lib/*.c", "lib/*.h"] + ["*.so", "*.pyd", "*.dll"],
        "gui": ["gui/*.ui", "gui/*.bat"],
        "gui/res": ["gui/res/*.*"],
    },
    scripts=["phaseshifts/phsh.py"],
    install_requires=[
        "scipy >= 0.7",
        "numpy >= 1.3",
        "periodictable",
        "typing_extensions",
    ],
    ext_modules=[Extension(name="", sources=[])],
    console=[os.path.join("phaseshifts", "phsh.py")],
    zip_safe=False,
)
