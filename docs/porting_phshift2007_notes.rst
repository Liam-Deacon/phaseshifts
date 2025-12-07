.. _PortingNotes:

==============================
``phshift2007`` Porting  Notes
==============================

This document contains notes on porting the original Barbieri / Van Hove
`phshift2007 <https://www.icts.hkbu.edu.hk/VanHove_files/leed/phshift2007.zip>`_
phase shift package code into the phaseshifts package.

In order to compile the FORTRAN code on a modern system, a number of changes were
made, including:

* The original code was written in FORTRAN 77 with some FORTRAN 66 features
  that are not allowed in common open source f77 compilers such as ``gfortran`` and ``flang``.
  In addition, some language-intrinsic functions were not portable such as ``DFLOAT``
  and were replaced with the equivalent standard functions, e.g. ``DFLOAT``
  (a GNU extension) was replaced with ``DBLE`` [1]_ and ``dexp`` replaced with ``exp``.
* In order to use as a callable library, ``PROGRAM`` units were replaced with ``SUBROUTINE``
  blocks. As part of this process these subroutines gained additional parameters for the input/output
  file names. Without these changes there would be multiple main blocks and the code would not be
  suitable for inclusion in a shared library.
* Data types were updated to be FORTRAN 77 compatible. Notably, ``HOLLERITH`` [2]_ constants
  representing string constants were replaced with the standard ``CHARACTER`` type.
* Leading tabs were replaced with 7 spaces to avoid ``-Wtabs`` warnings on ``gfortran`` due to the
  fact that tabs are not members of the Fortran Character Set [3]_ .

Additional changes were made to improve the readability of the code:

* Translated the comments in ``PHSH3.FOR`` (i.e. the CONPHAS program) from German to English.
  A python implementation can be found under ``phaseshifts.conphas``.
* Trailing whitespace was removed from all lines.
* Continuation lines such as those starting with ``'     1'`` were replaced with ``'     +'`` to
  improve readability by more easily distinguishing continuation lines from labels.
* ``real*8`` was replaced with ``double precision`` (and related casts, i.e. ``dfloat`` to ``dble``).
* ``GOTO`` semantics and DO with labels were refactored to more closely resemble other modern
  programming language constructs.
* ``c$OMP PARALLEL DO`` blocks were added where appropriate to allow for parallel execution of
  loops via OpenMP. This is only enabled when compiling with ``-fopenmp``.

.. [1] https://github.com/Liam-Deacon/phaseshifts/commit/fbd701e20f83d4eca90e5d90ef696d8316717d41
.. [2] https://en.wikipedia.org/wiki/Hollerith_constant#Examples
.. [3] https://gcc.gnu.org/onlinedocs/gfortran/Error-and-Warning-Options.html#index-Wtabs

.. danger::

    Even with the original code, the LEED Calculation Home Page offers no guarantees of correctness
    in the calculated phase shifts. Furthermore, compiling the code with a modern compiler
    against an unknown benchmark means that there are no assurances that the compiled programs
    will produce the same results as the original code authors' intended. There are probably plenty
    of bugs, as such please open an issue if you find any.

.. note::

    A notebook guiding the user through the initial porting process can be found at
    `porting-phshift2007-notes.ipynb <https://github.com/Liam-Deacon/phaseshifts/blob/master/porting-phshift2007-notes.ipynb>`_

    .. image:: https://mybinder.org/badge_logo.svg
     :target: https://mybinder.org/v2/gh/Liam-Deacon/phaseshifts/HEAD?labpath=porting-phshift2007-notes.ipynb

.. warning::

    Significant changes were made to the original code in order to port it for ``f2py``.
    Artifacts of this process are likely present in the code and may cause unexpected
    behaviour that deviates from the original intended purpose. As such, if you are
    looking for a reference implementation of the original code, please run
    ``make phshift2007`` and ``sudo make install``, which will download the original
    phshift2007 package, compile the ``phsh*`` programs and install them to ``$PREFIX/bin``
    (this is ``/usr/local/bin`` by default). They will then be available to run from the
    command line.

.. tip::

    Should you not trust the bundled f2py library, then a future version of ``phsh.py``
    will allow you to run the original phshift2007 programs via wrapped subprocess calls.

Compiler Notes
--------------

When originally porting the code back in 2014, the code was compiled with f2py,
Python 2.7 (32-bit) and the mingw32 toolchain on Windows 7 (installed together
via `Python(x, y) <https://python-xy.github.io/>`_ version <2.7.6.1). This was
a more permissive compiler toolchain than modern GCC-toolchain `gfortran <https://gcc.gnu.org/fortran/>`_
and LLVM-based `(classic) flang <https://github.com/flang-compiler/flang>`_ compilers tested
for the `v0.1.7 release <https://github.com/Liam-Deacon/phaseshifts/releases/tag/v0.1.7>`_.

.. note::

    According to wikipedia [4]_ ``g77`` is no longer included in the GCC project since v4
    as the maintainer decided to no longer support it. Another prominent fortran compiler ``g95``
    was also discontinued in 2012 and has diverged considerably from the original GNU compiler
    collection. As such it is possible that the fortran compiler included in the mingw32 toolchain
    used in the original porting was one of these compilers and this would explain why additional
    changes were required to compile the code with modern compilers.

.. [4] https://en.wikipedia.org/wiki/GNU_Compiler_Collection#Fortran

.. tip::

    Those wishing to perform a Windows build would be advised to use `Anaconda <https://www.anaconda.com/>`_
    and can be installed on Windows 10/11 using :command:`winget install --id Anaconda.Anaconda3`.
    Once installed, the conda environment can be installed with :command:`conda env create -f environment.yml`
    and activated with :command:`conda activate phaseshifts`. The phshift2007 code `could` then be compiled with
    :command:`gfortran -static-libgcc -static-libgfortran ...` (however no modern Windows build has been tried yet)


Compiler Test Matrix
--------------------

The following table compilers provides some summary information on compilers and platforms tested:

+------------------+--------+----------------+--------------+--------+---------------------------------------------------------+--------------+--------------+
| Compiler         | Version| Platform       | Architecture | Status | Notes                                                   | Date Tested  | Commit / Tag |
+==================+========+================+==============+========+=========================================================+==============+==============+
| gfortran         | 11     | Ubuntu 22.04   | x86_64       | ✔      | Built via ``ubuntu-latest`` GitHub Action runner [5]_   | 2024-01-21   | v0.1.8 [6]_  |
+------------------+--------+----------------+--------------+--------+---------------------------------------------------------+--------------+--------------+
| gfortran         | 11     | Mac OS X 12    | x86_64       | ✔      | Built via ``macos-latest`` GitHub Action runner [5]_    | 2024-01-21   | v0.1.8 [6]_  |
+------------------+--------+----------------+--------------+--------+---------------------------------------------------------+--------------+--------------+

.. [5] https://github.com/Liam-Deacon/phaseshifts/actions/workflows/publish-to-pypi.yaml
.. [6] https://github.com/Liam-Deacon/phaseshifts/releases/tag/v0.1.8

Phase Shift Workflow Summary
---------------------------

The following table summarizes the workflow and file dependencies for the Barbieri/Van Hove phase shift package:

+------+-------------------+-------------------------------+-------------------------------+-------------------------------------------------------------+
| Step | Program           | Input Files                   | Output Files                  | Notes                                                       |
+======+===================+===============================+===============================+=============================================================+
| 0    | PhSh0.for         | atorb                         | atelem.i                      | One per element                                             |
+------+-------------------+-------------------------------+-------------------------------+-------------------------------------------------------------+
| 1    | PhSh1.for         | cluster.i, atomic.i           | mufftin.d, check.o            | atomic.i = concat of atelem.i for all inequivalent atoms    |
+------+-------------------+-------------------------------+-------------------------------+-------------------------------------------------------------+
| 2    | PhSh2wil/rel/cav. | mufftin.d                     | phasout, dataph.d, inpdat,    | Use wil for non-rel, rel for relativistic;                 |
|      | for               |                               | leedph.d (wil only)           | phasout is split for next step                              |
+------+-------------------+-------------------------------+-------------------------------+-------------------------------------------------------------+
| 3    | PhSh3.for         | phJ (split from phasout)      | leedph.d, dataph.d            | One phJ per element                                         |
+------+-------------------+-------------------------------+-------------------------------+-------------------------------------------------------------+

Workflow Steps
~~~~~~~~~~~~~~

1. **STEP 0: Atomic Orbitals Calculation (`PhSh0.for`)**
   - Input: `atorb` (e.g., `atorbRh`)
   - Output: `atelem.i` (e.g., `atelemRh.i`)

2. **STEP 1: Muffin-Tin Potential Calculation (`PhSh1.for`)**
   - Input: `cluster.i`, `atomic.i` (concatenation of all `atelemJ.i` files)
   - Output: `mufftin.d`, `check.o`

3. **STEP 2: Phase Shift Calculation (`PhSh2cav.for`, `PhSh2wil.for`, `PhSh2rel.for`)**
   - Input: `mufftin.d`
   - Output: `phasout`, `dataph.d`, `inpdat`, `leedph.d` (wil only)

4. **STEP 3: Phase Shift Postprocessing (`PhSh3.for`)**
   - Input: `phJ` (split from `phasout`)
   - Output: `leedph.d`, `dataph.d`

Known Issues
------------

The following issues are known to exist in the current version of the code:

1. The code is not thread-safe. This is due to the use of global variables
   in the original code as well as large arrays that do not fit into stack memory.
   This is not a major issue if the user is aware of this and the code is not
   used in a multi-threaded context. Should the user need to ensure thread-safety,
   a workaround is to run via ephemeral docker containers, see :ref:`running` section.
2. Many minor compiler warnings have been ignored, such as those related to
   implicit typing of variables. These should be fixed in future releases.
