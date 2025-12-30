.. _installing_phaseshifts:

**********************************
Installing the phaseshifts package
**********************************

TDLR;
-----

For python 3.11 or older::

  pip install wheel numpy setuptools

  # latest pypi release
  pip install phaseshifts

  # or local version with essential packages for development
  git clone https://github.com/Liam-Deacon/phaseshifts
  cd phaseshifts
  pip install -e .[dev,test]  # !! best do this in a virtualenv

  phsh --help

Details
-------

The `phaseshifts <http://https://pypi.python.org/pypi/phaseshifts/>`_ package
requires CPython 2.7 or later and also uses the `numpy
<http://www.scipy.org/scipylib/download.html>`_, `scipy
<http://www.scipy.org/scipylib/download.html>`_ and `periodictable
<http://https://pypi.python.org/pypi/periodictable>`_ packages.
Currently, it has only been tested extensively with Python 2.7 on Windows,
and Ubuntu Linux so there are no guarantees with other platforms.
To perform a setup follow the steps below.

 1. Install the numpy, scipy and periodictable packages.

    On systems compatible with PyPI this can be done using the command::

      pip install numpy scipy periodictable

    Or if you have the easy_install package::

      easy_install install numpy scipy periodictable

    Older versions of ``numpy`` & ``scipy`` did not allow simultaneous installation -
    if you experience problems then try first installing numpy before
    attempting to install scipy.

    The ``periodictable`` package allows lookup of the most common crystal
    structure for a given element and is instrumental in many of the
    convenience functions contained in the model module.

    Alternatively download and install these packages manually following the
    instructions provided for the respective packages.

 2. To install the phaseshifts package::

      pip install phaseshifts

    Optional extras can be installed as needed::

      pip install "phaseshifts[viperleed]"  # EEASiSSS backend support (alias: viperleed)
      pip install "phaseshifts[input]"     # structured input (YAML/JSON) + validation
      pip install "phaseshifts[atorb]"     # extra element helpers for atorb
      pip install "phaseshifts[gui]"       # Qt GUI dependencies
      pip install "phaseshifts[docs]"      # documentation build deps
      pip install "phaseshifts[test]"      # test deps (pytest, pytest-cov)
      pip install "phaseshifts[dev]"       # dev tooling (black/isort/ruff, etc.)

    .. note:: Until a ``pyproject.toml`` with a working PEP-517 build backend
              is implemented then the user will first need to run
              :code:`pip install numpy setuptools wheel` in order to have the necessary
              python pre-requisites available (along with a fortran compiler) in order
              to compile the FORTRAN source and wrap it to be available via python.

    .. tip:: Running ``make check`` with run a test suite designed to catch issues with
             the installation (however the ``pytest`` package is required).

    With any luck the package has been installed successfully. A set of test scripts
    are provided, however a simple check may suffice using an interactive session of
    the python interpreter:

      >>> import phaseshifts
      >>> from phaseshifts.lib import libphsh  # compiled FORTRAN .pyd or .so using f2py

    If these execute without errors then it is likely that all is well, but in case of
    problems or bugs please use the contact provided below and I will do my best to
    address the problem quickly.

.. tip:: On Windows systems it may be easier to install a scientific python distibution
         rather than install the dependencies from source - `Python(x,y)
         <http://code.google.com/p/pythonxy>`_ or
         `Anaconda <https://www.anaconda.com/download>`_ with mingw (gcc & gfortran)
         installed is highly recommended. Mac OS X users can simply do ``brew install gfortran``
         and Debian/Ubuntu users can do ``sudo apt-get install -y gfortran``.

.. note:: On unix systems, setup the virtualenv on Python 3.10 or lower, activate it and run `make`.

.. note::
   Python 3.12+ is now fully supported! The build system uses `scikit-build-core` and CMake for Fortran extensions.
   - Always use a virtual environment for development and installation.
   - Install CMake and a Fortran compiler (e.g., gfortran) before building.
   - To install locally:

     python3.12 -m venv .venv
     source .venv/bin/activate
     pip install .

   - To install from PyPI:

     pip install phaseshifts

   For Python ≤3.11, the installer will first attempt a modern CMake/scikit-build build (if available), and will fallback to legacy numpy.f2py if CMake/scikit-build is not installed or fails. Just use:

     pip install -e .  # will try CMake/scikit-build first, then fallback to f2py

   A Fortran compiler is required for both build paths.

.. tip::
   If you encounter build errors, ensure you have CMake, scikit-build-core, and a working Fortran compiler installed. If the modern build fails, the installer will automatically fallback to the legacy build if possible.

.. note::
   Optional extras can be installed with:
   - ``pip install "phaseshifts[input]"`` for structured input + validation.
   - ``pip install "phaseshifts[viperleed]"`` to enable the experimental
     EEASiSSS backend (via ViPErLEED, alias: ``viperleed``).

Build distributions
-------------------

To build source distributions and wheels locally, use the standard ``build``
module or the Makefile targets::

  # With build
  python -m pip install build
  python -m build --sdist --wheel --no-isolation

  # Or via Makefile
  make sdist
  make wheel

The EEASiSSS backend is an optional runtime dependency; install the extra
(``phaseshifts[viperleed]``) on systems where you want to run EEASiSSS through
ViPErLEED.

.. tip::
   To optionally disable building the phshift2007 Fortran binaries (e.g., for minimal installs or CI), set the environment variable before installation:

     export PHASESHIFTS_BUILD_PHSHIFT2007_BINARIES=OFF
     pip install .

   By default, the binaries are built and included. Set to "OFF" to skip downloading, building, and packaging the phshift2007 binaries. This is useful for environments where Fortran compilers are unavailable or the binaries are not needed (typically useful for debugging).
