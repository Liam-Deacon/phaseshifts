.. _api:

===============
phaseshifts API
===============

Package Contents
----------------

This chapter covers the main modules of the phaseshifts and provides some 
API documentation for those wishing to incorporate this package into 
their own projects. 

.. automodule:: phaseshifts
    :members:
    :undoc-members:
    :show-inheritance:
    
Sub-packages
------------

The main sub-packages are listed below, however some notable sub-packages that 
deserve a brief mention are: 

* :py:mod:`phaseshifts.gui` - includes all the necessary files for the graphical user interface.
* :py:mod:`phaseshifts.lib` - contains the Fortran libphsh library and the python wrappings.
* :py:mod:`phaseshifts.doc` - source documentation for the phaseshifts package.
* :py:mod:`phaseshifts.test` - modules for testing the phaseshift package.

Submodules
----------

phaseshifts.atorb
+++++++++++++++++

.. automodule:: phaseshifts.atorb
    :members:
    :private-members:
    :undoc-members:
    :show-inheritance:
    
    .. autoclass:: Atorb 
      :members:
      :private-members:
      :undoc-members:
      :show-inheritance:

phaseshifts.conphas
+++++++++++++++++++

.. automodule:: phaseshifts.conphas
    :members:
    :show-inheritance:
    
    .. autoclass:: Conphas
      :members:
      :show-inheritance:
      :private-members:

phaseshifts.elements
++++++++++++++++++++

.. automodule:: phaseshifts.elements
    :members:
    :private-members:
    :undoc-members:
    :show-inheritance:
    
    .. autoclass:: Element
      :members:
      :private-members:
      :undoc-members:
      :show-inheritance:
      
    .. autodata:: ELEMENTS
      :annotation: = list of known elements

phaseshifts.factories
+++++++++++++++++++++

.. automodule:: phaseshifts.factories
    :members:
    :private-members:
    :undoc-members:
    :show-inheritance:

phaseshifts.leed
++++++++++++++++

.. automodule:: phaseshifts.leed
    :members:
    :private-members:
    :undoc-members:
    :show-inheritance:
    
phaseshifts.model
+++++++++++++++++

.. automodule:: phaseshifts.model
    :members:
    :private-members:
    :undoc-members:
    :show-inheritance:
            
phaseshifts.phsh
++++++++++++++++

.. automodule:: phaseshifts.phsh
    :members:
    :private-members:
    :undoc-members:
    :show-inheritance:

phaseshifts.utils
+++++++++++++++++

.. automodule:: phaseshifts.utils
    :members:
    :private-members:
    :undoc-members:
    :show-inheritance:
    
phaseshifts.wrappers
++++++++++++++++++++

.. automodule:: phaseshifts.wrappers
    :members:
    :private-members:
    :undoc-members:
    :show-inheritance: