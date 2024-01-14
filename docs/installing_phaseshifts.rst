.. _installing_phaseshifts:

**********************************
Installing the phaseshifts package
**********************************

The `phaseshifts <http://https://pypi.python.org/pypi/phaseshifts/>`_ package 
requires CPython 2.6 or later and also uses the `numpy 
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
         <http://code.google.com/p/pythonxy>`_ with mingw (gcc & gfortran) 
         installed is highly recommended.

