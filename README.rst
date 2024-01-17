===================
PHASESHIFTS PACKAGE
===================

This package is a Python-based implementation of the Barbieri/Van Hove
phase shift (a.k.a. *phshift*) calculation package needed to produce phase shifts for
various LEED packages (including CLEED), as well as for certain XPD packages.
To quote the original authors' site:

"The phase shift calculation is performed in several steps:

1. Calculation of the radial charge density for a free atom.

2. Calculation of the radial muffin-tin potential for atoms embedded in a
   surface defined by the user (the surface is represented by a slab that
   is periodically repeated in 3 dimensions, within vacuum between the
   repeated slabs); various approximations to the exchange potential
   are available; relativistic effects are taken into account.

3. Calculation of phase shifts from the muffin-tin potential.

4. Elimination of pi-jumps in the energy dependence of the phase shifts."

.. note:: You can get the original Fortran source (& learn more about the *phshift* programs)
   from Michel Van Hove's LEED Calculation Home Page:

   https://www.icts.hkbu.edu.hk/VanHove_files/leed/leedpack.html

   A local copy of the source files can be found under ``phaseshifts/lib/.phsh.orig/phsh[0-2].f``.

Running
=======

The `phsh.py` script (available after installing the package) aims to simplify these
steps with a single command. For more information please read the documentation
at `<http://pythonhosted.org//phaseshifts/>`_ (latest PyPI release) or
`GitHub Pages <https://liam-deacon.github.io/phaseshifts/>`_ (latest master)

The simplest and most reliable cross-platform way to run `phsh.py` is through docker::

  # obtain the image
  docker pull ghcr.io/Liam-Deacon/phaseshifts:latest  # should only need to do this once

  # run phsh.py via the docker image
  docker run ghcr.io/Liam-Deacon/phaseshifts:latest  # will display usage

  # or more generally (adjust as needed)
  docker run ghcr.io/Liam-Deacon//phaseshifts:latest -v /path/to/host/input/data:/data [<phsh-args> ...]


.. tip:: Development docker images can be built locally, e.g.
         :code:`docker build . -t ghcr.io/Liam-Deacon/phaseshifts:dev`

.. warning:: There is a `known possible bug <https://github.com/Liam-Deacon/phaseshifts/issues/6>`_
             where the compiled ``libphsh.f`` is not thread-safe (as ascertained by the fortran compiler),
             as such if you anticipate using this library in concurrent environments then it is advised to
             run ``phsh.py`` via :code:`docker run ghcr.io/Liam-Deacon/phaseshifts:latest` as this works around
             this limitation due to the emphereal nature of container instances created using ``docker run``.

Install
=======

TDLR;
-----

For python 3.11 or older::

  pip install wheel numpy setuptools
  pip install -e .
  phsh --help

Details
-------

The `phaseshifts <http://https://pypi.python.org/pypi/phaseshifts/>`_ package
requires CPython 2.7 or later and also uses the `numpy
<http://www.scipy.org/scipylib/download.html>`_, `scipy
<http://www.scipy.org/scipylib/download.html>`_ and `periodictable
<http://https://pypi.python.org/pypi/periodictable>`_ packages.
Currently, it has only been tested extensively with Python 2.7 on Windows, so
there are no guarantees with other platforms. To perform a setup follow the
steps below.

 1. Install the numpy, scipy and periodictable packages.

    On systems compatible with PyPI this can be done using the command::

      pip install numpy scipy periodictable

    Or if you have the easy_install package::

      easy_install install numpy scipy periodictable

    Older versions of numpy & scipy did not allow simultaneous installation -
    if you experience problems then try first installing numpy before
    attempting to install scipy.

    The periodictable package allows lookup of the most common crystal
    structure for a given element and is instrumental in many of the
    convenience functions contained in the model module.

    Alternatively download and install these packages manually following the
    instructions provided for the respective packages.

 2. To install the phaseshifts package::

      python setup.py install

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

.. warning:: Python 3.12 compatibility is a work in progress due to the removal of ``numpy.distuils``
             build backend for ``f2py`` preventing simple installation via ``pip install``,
             `this github issue <https://github.com/Liam-Deacon/phaseshifts/issues/8>`_
             tracks progress on fixing this known issue.


About the code
==============

The example source codes provided in this package are intended to be
instructional in calculating phase shifts. While it is not recommended to
use the example code in production, the code
should be sufficient to explain the general use of the library.

If you aren't familiar with the phase shift calculation process, you can
read further information in ``doc/`` folder:

+ ``phshift2007.rst`` - a brief user guide/documentation concerning the input files
  (& details of the original fortran `phshift` package).
+ ``phaseshifts.pdf`` - a more detailed overview of the library functions and how to
  calculate phase shifts using the convenience functions in this package. This is not
  yet finished and so the reader is referred to the above document for the time being.

For those wanting a crash course of the Van Hove / Tong programs, I advise reading the
phsh2007.txt document.
See the ``examples/`` directory to get an idea of the structure of the input files
(for a random selection of models & elements). In particular see the ``cluster_Ni.i``
file for helpful comments regarding each line of input.

Those of you who are eager to generate phase shifts - first look at the example
cluster files for a bulk and slab calculation, noting that the atoms in the model
are in fractional units of the *a* basis vector for the unitcell (SPA units). Next,
after creating a bulk and slab model in the ``cluster.i`` format, simply use
the following python code:

   >>> from phaseshifts.phsh import Wrapper as phsh
   >>> phsh.autogen_from_inputs(bulk_file, slab_file)

This will hopefully produce the desired phase shift output files (at least for
simple models) and works by assessing the two models to determine what output to
produce. For more detailed documentation and function use refer to the pdf manual.

.. tip:: A standalone command line utility **phsh.py** is provided as a way of
         automating the generation of phase shifts as part of a script. For more
         information use:

         .. code:: bash

            phsh.py --help

.. note:: The `phaseshifts.leed` module provides a conversion class for CLEED ``.inp`` and
          ``.bul`` files. This is included as part of the `phsh.py` module,
          however the file extension is important (needs ``.inp``, ``.pmin``, ``.bul``,
          or ``.bmin``) and error checking is limited. There are also plans to include a
          validator to check the files for malformatted input at some point in the
          future.

Alternatives
------------

A number of alternatives are available, notably the following:

1. `AQuaLEED <https://physics.mff.cuni.cz/kfpp/povrchy/software>`_ (with a useful
   `poster overview of phaseshifts calculations <https://physics.mff.cuni.cz/kfpp/povrchy/files/1179-Poster.pdf>`_).
   This is an officially mentioned piece of software on Michel Van Hove's
   `LEEDPACK webpage <https://www.icts.hkbu.edu.hk/VanHove_files/leed/leedpack.html>`_,
   however when tested as of January 2024 the link appears to be dead (with a ``500 INTERNAL_SERVER_ERROR``).
   Furthermore, although the poster mentions that the software is written in python,
   this software is not (currently) distributed on https://PyPI.org and therefore harder to
   intergrate with other python LEED-related projects such as `CLEED <https://github.com/Liam-Deacon/CLEED>`_
   and `cleedpy <https://github.com/empa-scientific-it/cleedpy>`_.

.. note:: Should you know of alternatives, please either
          `open an issue <https://Liam-Deacon/phaseshifts/issues>`_ or
          (better yet) create a PR with changes to this documentation
          to keep this list up to date.


Acknowledgements
================

As with all scientific progress, we stand on the shoulders of giants. If this
package is of use to you in publishing papers then please acknowledge the
following people who have made this package a reality:

 - **A. Barbieri** and **M.A. Van Hove** - who developed most of the original
   fortran code. Use *A. Barbieri and M.A. Van Hove, private communication.*
   (see ``doc/phsh2007.txt`` for further details).

 - **E.L. Shirley** - who developed part of the fortran code during work towards his
   PhD thesis (refer to the thesis: *E.L. Shirley, "Quasiparticle calculations in
   atoms and many-body core-valence partitioning", University of Illinois, Urbana, 1991*).

 - **Christoph Gohlke** - who developed the elements.py module used extensively throughout
   for the modelling convenience functions (see 'elements.py' for license details).

 I would also be grateful if you acknowledge this python package (*phaseshifts*) as:
 *L.M. Deacon, private communication.*


Thanks
------

I wish to personally add a heart-felt thanks to both Eric Shirley and Michel Van Hove
who have kindly allowed the use of their code in the ``libphsh.f`` file needed for the
underlying low-level functions in this package.

Contact
=======

This package is developed/maintained in my spare time so any bug reports, patches,
or other feedback are very welcome and should be sent to: liam.deacon@diamond.ac.uk

The project is in the early developmental stages and so anyone who wishes to get
involved are most welcome (simply contact me using the email above).

To Do
=====

 1. Documentation - the manual has been started, but is not complete and thus is a
    high priority. The current aim is to use sphinx to generate html and latex documents
    for semi-automated generation of both the tutorial and supporting website. If
    you have the phaseshifts source and the `sphinx <https://pypi.python.org/pypi/Sphinx>`_
    and the `numpydoc <https://pypi.python.org/pypi/numpydoc>`_ PyPi packages then you
    can try making html or latex manuals using ``make html`` or ``make latexpdf`` commands
    from the ``doc/`` directory.

 2. Test suit to verify the package is working as expected.

 3. GUI frontend (Qt ui files are provided in the ``gui/`` directory for anyone
    wishing to undertake this challenge). Other frontends are welcome (I use Qt
    due to familiarity/experience). For those wishing a sneak preview, try executing
    ``main.pyw``

See ``TODO.rst`` for more information.

Contacts
========

  - `Liam Deacon <mailto://liam.m.deacon@gmail.com>`_ - *current maintainer*
  - `Michel Van Hove <mailto://vanhove@cityu.edu.hk>`_ - Contact for original LEEDPACK ``phsh[0-3].f`` programs
