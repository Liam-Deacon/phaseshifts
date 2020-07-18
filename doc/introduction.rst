.. _introduction:

************
Introduction
************

This package is a Python-based implementation of the Barbieri/Van Hove 
phase shift (*phsh*) calculation package needed to produce phase shifts for 
various LEED packages (including CLEED), as well as for certain XPD packages. 

To quote the original authors site: 

"The phase shift calculation is performed in several steps:

1. Calculation of the radial charge density for a free atom.

2. Calculation of the radial muffin-tin potential for atoms embedded in a 
   surface defined by the user (the surface is represented by a slab that 
   is periodically repeated in 3 dimensions, within vacuum between the 
   repeated slabs); various approximations to the exchange potential 
   are available; relativistic effects are taken into account.

3. Calculation of phase shifts from the muffin-tin potential.

4. Elimination of pi-jumps in the energy dependence of the phase shifts."

.. note:: You can get the original Fortran source 
 (& learn more about the *phsh* programs) from:

   http://www.icts.hkbu.edu.hk/surfstructinfo/SurfStrucInfo_files/leed/leedpack.html

The aim of this package is to both automate and simplify the generation of 
phase shift files in a manner that is easy for the computational hitch-hiker, 
but powerful for those that wish to extend the package for particular needs.