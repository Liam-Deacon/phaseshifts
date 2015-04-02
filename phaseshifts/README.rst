===================
PHASESHIFTS PACKAGE
===================

This package is a Python package which produces atomic phase shifts for 
various LEED packages (including CLEED), as well as for certain XPD packages. 

Currently, it uses the Barbieri/Van Hove phase shift calculation package and 
John Rundgren's EEASiSSS package as backends that are wrapped using f2py and a 
few Python modules to provide a user-friendly object-orientated interface for 
the end user.

The `phsh.py` script aims to simplify/unify the steps needed for the phase 
shift calculations into a single command intended for user. For more 
information please read the documentation at 
`<http://pythonhosted.org//phaseshifts/>`_

-------------------------
Barbieri/Van Hove backend
-------------------------

The original phase shift package developed by A. Barbieri & M. A. Van Hove 
during the 1970's & 1980's. To quote the authors' site: 

"The phase shift calculation is performed in several steps:

1. Calculation of the radial charge density for a free atom.

2. Calculation of the radial muffin-tin potential for atoms embedded in a 
   surface defined by the user (the surface is represented by a slab that 
   is periodically repeated in 3 dimensions, within vacuum between the 
   repeated slabs); various approximations to the exchange potential 
   are available; relativistic effects are taken into account.

3. Calculation of phase shifts from the muffin-tin potential.

4. Elimination of pi-jumps in the energy dependence of the phase shifts."

.. note:: You can get the original Fortran source (& learn more about the 
   *phshift* programs) from `here
   <http://www.icts.hkbu.edu.hk/surfstructinfo/SurfStrucInfo_files/leed/leedpack.html>`_

Please contact Michel Van Hove <vanhove@hkbu.edu.hk> regarding this package.

----------------
EEASiSSS backend
----------------

The Elastic Electron-Atom Scattering in Solids and Surface Slabs (EEASiSSS) 
package [#]_ [#]_ [#]_ can also calculate phase shifts and is used 
in a number of recent works on oxides [#]_ [#]_. In the words of the package's 
author, John Rundgren, the main qualifications of the program are:

+ The program accepts three sources of atomic potentials: 
    
    1. E. L. Shirley’s atomic program [#]_ applied together with Mattheiss’s 
    superposition method.
           
    2. The DFT-SCF program of V. Eyert using the full-potential 
    Augmented Spherical Wave method [#]_ .
          
    3. The DFT-SCF program `WIEN2k <http://www.wien2k.at/>`_ using the 
    full-potential Augmented Plane Wave method.
           
+ The exchange-correlation interaction between the scattered electron and the 
  crystal’s electron gas generates an energy-dependent inner potential. 
  The phase shifts are referred to the in-crystal kinetic energy, and it is 
  supposed that an associated LEED code uses the same standard.
+ The crystal potential is everywhere continuous so as to exclude fortuitous 
  standing-wave electron resonances in the muffin-tin spheres and pertaining 
  fortuitous wobblings in the phase shift versus energy curves.
+ The optimization of the muffin-tin radii is made using the method of 
  `Differential Evolution <http://www.physik.uni-augsburg.de/~eyert/ASWhome.shtml/>`_ 
  , an extremely efficient minimizer.

.. Note:: A short EEASiSSS users guide is appended to the input template files
          `inputA` and `inputX` distributed with the program package.

.. [#] J Rundgren, Phys. Rev. B 68, 125405 (2003).
.. [#] J Rundgren, Phys. Rev. B 76, 195441 (2007).
.. [#] E. A. Soares, C. M. C. De Castillho, and V. E. Carvalho, J. Phys.: Condens. Matter
   23,303001 (2011).
.. [#] R. Pentcheva, W. Moritz, J. Rundgren, S. Frank, D. Schrupp, and M. Scheffler, Surf. Sci
   602, 1299 (2008).
.. [#] V.B. Nascimento, R.G. Moore, J. Rundgren, J. Zhang, L. Cai, R. Jin, D.G. Mandrus,
   and E.W. Plummer, Phys. Rev. B 75, 035408 (2007).
.. [#] S. Kotochigova, Z. H. Levine, E. L. Shirley, M. D. Stiles, and C. W. Clark, Phys. Rev. B
   55, 191 (1997).
.. [#] R Storn and K. Price, J. Global Optimization 11, 341 (1997).

Please contact John Rundgren <jru@kth.se> for queries, comments or suggestions 
related to this package.


Install
=======

The `phaseshifts <http://https://pypi.python.org/pypi/phaseshifts/>`_ package 
requires CPython 2.6 or later and also uses the `numpy 
<http://www.scipy.org/scipylib/download.html>`_, `scipy 
<http://www.scipy.org/scipylib/download.html>`_ and `periodictable 
<http://https://pypi.python.org/pypi/periodictable>`_ packages. 
Currently, it has only been tested extensively with Python 2.7 on Windows
and Ubuntu Linux, so there are no guarantees with other platforms. 
To perform a setup follow the steps below.

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

    With any luck the package has been installed successfully. A set of test scripts
    are provided, however a few simple checks may suffice using the command and an interactive session of the python interpreter::

      phsh.py --help
      python
      >>> import phaseshifts
      >>> from phaseshifts.lib import libphsh  # compiled FORTRAN .pyd or .so using f2py
      >>> libphsh.hartfock
      <fortran object>
      >>> exit(0) # okay if no errors found above

    If these execute without errors then it is likely that all is well, but in case of 
    problems or bugs please use the contact provided below and I will do my best to 
    address the problem quickly.

.. tip:: On Windows systems it may be easier to install a scientific python distibution 
         rather than install the dependencies from source - `Python(x,y) 
         <http://code.google.com/p/pythonxy>`_ with mingw (gcc & gfortran) 
         installed is highly recommended.

Linux Users
-----------

Installation on Linux hosts is not as rigorously tested compared with Windows, 
therefore the following is meant as a guide only. If you experience any 
difficulties then please don't hesitate to contact me using 
liam.deacon@diamond.ac.uk

.. note:: You may need to specify an explicit version of ``python`` when trying 
          the commands below (e.g. ``python2.7``)

Ubuntu
++++++

On Ubuntu (14.04.01 LTS) some dependencies were needed, which can be installed 
using the following bash commands::

   $ sudo apt-get install python-pip python-numpy python-scipy python-dev \
   libblas-dev gfortran # should install f2py as a dependency of numpy/scipy 
   $ sudo pip install phaseshifts==0.1.5-dev | tee phaseshifts.log # explicit version needed
   $ sudo ln -s /usr/bin/python2.7 /usr/bin/python # can skip if only one version of python is installed 

You may need to verify whether the ``phsh.py`` script is installed on 
the system $PATH and therefore please check by typing::

   $ python /usr/local/bin/phsh.py --help # should work okay
   $ phsh.py --help # may not work: phsh.py should be on the $PATH and executable


OpenSUSE
++++++++

On OpenSUSE (13.2 x86-64) the required dependencies can be installed using::

    $ sudo zypper install python-pip python-devel python-numpy-devel \ 
    python-scipy-devel libblas-devel gcc-fortran
    $ sudo pup install phaseshifts==0.1.5-dev | tee phaseshifts.log
    
Again, you may need to verify whether the ``phsh.py`` script is installed on the 
system $PATH and therefore please check by typing::

   $ python /usr/local/bin/phsh.py --help # should work okay
   $ phsh.py --help # may not work: phsh.py should be on the $PATH and executable
    

Development
===========

For those wishing to see the latest code, please visit: 
`<https://bitbucket.org/Liam_Deacon/phaseshifts/overview>`_ 

The latest code is actively developed and may potentially be bug ridden... 
It is therefore recommended that you use it with caution and a great deal of 
scepticism! Please raise the 
`alarm here <https://bitbucket.org/Liam_Deacon/phaseshifts/issues?status=new&status=open>`_ 
if you find any bugs, or better yet submit a patch into the development branch!
Your contributions will be greatly appreciated.

About the code
==============

The example source codes provided in this package are intended to be 
instructional in calculating phase shifts. While it is not recommended to 
use the example code in production (unless you are lucky enough to have a 
matching model), the code should be sufficient to explain the general use of 
the library.

If you aren't familiar with the phase shift calculation process, you can 
read further information in ``doc/`` folder:

+ ``phshift2007.rst`` - a brief user guide/documentation for performig 
  calculations using the Barbieri/Van Hove phase shift package and details 
  the input files (for the original fortran `phshift` package).
+ ``phaseshifts.pdf`` - a more detailed overview of the library functions and 
  how to calculate phase shifts using the convenience functions in this 
  package. This is not yet finished and so there may be gaps in the supplied 
  information.

For those wanting a crash course of the Van Hove / Tong programs, I advise 
reading the phsh2007.txt document. See the ``examples/`` directory to get an 
idea of the structure of the input files (for a random selection of models 
& elements). In particular see the ``cluster_Ni.i`` file for helpful comments 
regarding each line of input.

Likewise it is advised to look at ``inputA`` and ``inputX`` files to get an 
idea of how the EEASiSSS input is structured. Notice that whilst there are 
similarities between both packages, the two types of input are not compatible. 
This may be something that is worked upon for a future release.

Those of you who are eager to generate phase shifts - first look at the example
cluster files for a bulk and slab calculation, noting that the atoms in the model
are in fractional units of the *a* basis vector for the unitcell (SPA units). Next, 
after creating a bulk and slab model in the ``cluster.i`` format, simply use 
the following code within the python interpreter:
 
   >>> from phaseshifts.wrappers import VHTWrapper as phsh
   >>> phsh.autogen_from_inputs(bulk_file, slab_file)

This will hopefully produce the desired phase shift output files (at least for 
simple models) and works by assessing the two models to determine what output to
produce. For more detailed documentation and function use refer to the pdf manual.  

.. tip:: A standalone command line utility **phsh.py** is provided as a way of 
         automating the generation of phase shifts as part of a script. For more 
         information use:
         
         .. code:: bash
            
            phsh.py --help
         
.. note:: The `leed.py` module provides a conversion class for CLEED .inp and 
          .bul files. This is included as part of the `phsh.py` module, 
          however the file extension is important (needs .inp, .pmin, .bul, or .bmin) 
          and error checking is limited. There are also plans to include a 
          validator to check the files for malformatted input at some point in the 
          future.
         
Acknowledgements
================

As with all scientific progress, we stand on the shoulders of giants. If this 
package is of use to you in publishing papers then please acknowledge the 
following people who have made this package a reality:

 - **A. Barbieri** and **M.A. Van Hove** - who developed most of the original 
   fortran code. Use *A. Barbieri and M.A. Van Hove, private communication.* 
   (see ``doc/phsh2007.txt`` for further details).
 
 - **E.L. Shirley** - who developed part of the fortran code during work 
   towards his PhD thesis (refer to the thesis: *E.L. Shirley, "Quasiparticle 
   calculations in atoms and many-body core-valence partitioning", 
   University of Illinois, Urbana, 1991*).

 - **J. Rundgren** - who developed the EEASiSSS package and collaborated on 
   numerous improvements to the underlying FORTRAN code base. Please cite
   *J. Rundgren, Phys. Rev. B 68 125405 (2003).*

 - **Christoph Gohlke** - who developed the elements.py module used extensively 
   throughout for the modelling convenience functions (see 'elements.py' 
   for license details). 

 I would also be grateful if you acknowledge use of this python package 
 (*phaseshifts*) as: *L.M. Deacon, private communication.*


Thanks
------

I wish to personally add a heart-felt thanks to Eric Shirley, John Rundgren and 
Michel Van Hove who have kindly allowed the use of their code needed for the 
underlying low-level functions in this package. 

Contact
=======

This package is developed/maintained in my spare time so any bug reports, 
patches, or other feedback are very desirable and should be sent to: 
liam.deacon@diamond.ac.uk

The project is in the early developmental stages and so anyone who wishes to get 
involved are most welcome (simply contact me using the email above).

To do
=====

 1. Documentation - the manual has been started, but is not complete and thus 
    is a high priority. The current aim is to use sphinx to generate html and 
    latex documents for semi-automated generation of both the tutorial and 
    supporting website. If you have the phaseshifts source and the 
    `sphinx <https://pypi.python.org/pypi/Sphinx>`_ and the 
    `numpydoc <https://pypi.python.org/pypi/numpydoc>`_ PyPi packages then you 
    can try making html or latex manuals using ``make html`` or 
    ``make latexpdf`` commands from the ``doc/`` directory.

 2. Test suit to verify the package is working as expected.

 3. GUI frontend (Qt ui files are provided in the ``gui/`` directory for anyone 
    wishing to undertake this challenge). Other frontends are welcome (I use Qt 
    due to familiarity/experience). For those wishing a sneak preview, try executing
    ``main.pyw``

See ``TODO.rst`` for more information.

Author list
===========

  - `Liam Deacon <liam.deacon@diamond.ac.uk>`_ - *current maintainer*