.. _installing_phaseshifts:

**********************************
Installing the phaseshifts package
**********************************

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