#!/usr/bin/env python
#encoding: utf-8

##############################################################################
# Author: Liam Deacon                                                        #
#                                                                            #
# Contact: liam.deacon@diamond.ac.uk                                         #
#                                                                            #
# Copyright: Copyright (C) 2013-2014 Liam Deacon                             #
#                                                                            #
# License: MIT License                                                       #
#                                                                            #
# Permission is hereby granted, free of charge, to any person obtaining a    #
# copy of this software and associated documentation files (the "Software"), #
# to deal in the Software without restriction, including without limitation  #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,   #
# and/or sell copies of the Software, and to permit persons to whom the      #
# Software is furnished to do so, subject to the following conditions:       #
#                                                                            #
# The above copyright notice and this permission notice shall be included in #
# all copies or substantial portions of the Software.                        #
#                                                                            #
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    #
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING    #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER        #
# DEALINGS IN THE SOFTWARE.                                                  #
#                                                                            #
##############################################################################
"""
**atorb.py**

Provides convenience functions for generating input and calculating 
atomic charge densities for use with the Barbieri/Van Hove phase 
shift calculation package.

:See: http://www.icts.hkbu.edu.hk/surfstructinfo/SurfStrucInfo_files/leed/

:Requires: f2py (for libphsh fortran wrapper generation) 

.. note::
   To generate libphsh fortran wrappers (libphsh.pyd) for your platform
   then use 'python setup.py' in the lib directory of this package to 
   install into your python distribution. Alternatively, use::
           
     f2py -c -m libphsh libphsh.f
             
   Windows users may have to add appropriate compiler switches, e.g. ::
           
    f2py -c -m libphsh --fcompiler=gfortran --compiler=mingw-32 libphsh.f

"""

import os
import elements
from lib.libphsh import hartfock
from collections import OrderedDict

elements_dict = OrderedDict([
('H', 'Hydrogen'), 
('He', 'Helium'), 
('Li', 'Lithium'), 
('Be', 'Beryllium'), 
('B', 'Boron'), 
('C', 'Carbon'), 
('N', 'Nitrogen'), 
('O', 'Oxygen'), 
('F', 'Fluorine'), 
('Ne', 'Neon'), 
('Na', 'Sodium'), 
('Mg', 'Magnesium'), 
('Al', 'Aluminium'), 
('Si', 'Silicon'), 
('P', 'Phosphorus'), 
('S', 'Sulfur'), 
('Cl', 'Chlorine'), 
('Ar', 'Argon'), 
('K', 'Potassium'), 
('Ca', 'Calcium'), 
('Sc', 'Scandium'), 
('Ti', 'Titanium'), 
('V', 'Vanadium'), 
('Cr', 'Chromium'), 
('Mn', 'Manganese'), 
('Fe', 'Iron'), 
('Co', 'Cobalt'), 
('Ni', 'Nickel'), 
('Cu', 'Copper'), 
('Zn', 'Zinc'), 
('Ga', 'Gallium'), 
('Ge', 'Germanium'), 
('As', 'Arsenic'), 
('Se', 'Selenium'), 
('Br', 'Bromine'), 
('Kr', 'Krypton'), 
('Rb', 'Rubidium'), 
('Sr', 'Strontium'), 
('Y', 'Yttrium'), 
('Zr', 'Zirconium'), 
('Nb', 'Niobium'), 
('Mo', 'Molybdenum'), 
('Tc', 'Technetium'), 
('Ru', 'Ruthenium'), 
('Rh', 'Rhodium'), 
('Pd', 'Palladium'), 
('Ag', 'Silver'), 
('Cd', 'Cadmium'), 
('In', 'Indium'), 
('Sn', 'Tin'), 
('Sb', 'Antimony'), 
('Te', 'Tellurium'), 
('I', 'Iodine'), 
('Xe', 'Xenon'), 
('Cs', 'Cesium'), 
('Ba', 'Barium'), 
('La', 'Lanthanum'), 
('Ce', 'Cerium'), 
('Pr', 'Praseodymium'), 
('Nd', 'Neodymium'), 
('Pm', 'Promethium'), 
('Sm', 'Samarium'), 
('Eu', 'Europium'), 
('Gd', 'Gadolinium'), 
('Tb', 'Terbium'), 
('Dy', 'Dysprosium'), 
('Ho', 'Holmium'), 
('Er', 'Erbium'), 
('Tm', 'Thulium'), 
('Yb', 'Ytterbium'), 
('Lu', 'Lutetium'), 
('Hf', 'Hafnium'), 
('Ta', 'Tantalum'), 
('W', 'Tungsten'), 
('Re', 'Rhenium'), 
('Os', 'Osmium'), 
('Ir', 'Iridium'), 
('Pt', 'Platinum'), 
('Au', 'Gold'), 
('Hg', 'Mercury'), 
('Tl', 'Thallium'), 
('Pb', 'Lead'), 
('Bi', 'Bismuth'), 
('Po', 'Polonium'), 
('At', 'Astatine'), 
('Rn', 'Radon'), 
('Fr', 'Francium'), 
('Ra', 'Radium'), 
('Ac', 'Actinium'), 
('Th', 'Thorium'), 
('Pa', 'Protactinium'), 
('U', 'Uranium'), 
('Np', 'Neptunium'), 
('Pu', 'Plutonium'), 
('Am', 'Americium'), 
('Cm', 'Curium'), 
('Bk', 'Berkelium'), 
('Cf', 'Californium'), 
('Es', 'Einsteinium'), 
('Fm', 'Fermium'), 
('Md', 'Mendelevium'), 
('No', 'Nobelium'), 
('Lr', 'Lawrencium'), 
('Rf', 'Rutherfordium'), 
('Db', 'Dubnium'), 
('Sg', 'Seaborgium'), 
('Bh', 'Bohrium'), 
('Hs', 'Hassium'), 
('Mt', 'Meitnerium'), 
('Ds', 'Darmstadtium'), 
('Rg', 'Roentgenium'), 
('Cn', 'Copernicium'), 
('Uut', 'Ununtrium'), 
('Fl', 'Flerovium'), 
('Uup', 'Ununpentium'), 
('Lv', 'Livermorium'), 
('Uus', 'Ununseptium'), 
('Uuo', 'Ununoctium'), 
])


class Atorb(object):
    '''
    Description
    -----------
     A python wrapper for the atorb program by Eric Shirley for use in
     calculating atomic scattering for different elements
    
    Notes
    -----
    Original author: Eric Shirley  

    There are nr grid points, and distances are in Bohr radii
    :math:`a_0 \simeq 0.539 \mathrm{\AA}`
    
    :math:`r(i) = r_{min} \cdot (r_{max} / r_{min})^{(i/n_r)}`, 
    :math:`i=1,2,3,...n_r-1,n_r`
    
    The orbitals are stored in phe(), first index goes :math:`1...n_r`, the
    second index is the orbital index (:math:`i...n_{el}`)
    
    Look at the atomic files after printing this out to see everything...
    Suffice it to say, that the charge density at radius :math:`r(i)`
    in units of electrons per cubic Bohr radius is given by:
    
    :math:`\sum_{j-1}^{n_el}{occ(j) \cdot phe(i,j)^2 / (4.0\,\pi\,{r(i)^2)}}` 
    
    Think of the phe functions as plotting the radial wave-functions
    as a function of radius on a logarithmic mesh...  
    
    The Dirac equation is solved for the orbitals, whereas their density
    is treated by setting :math:`phe(i,j)` to Dirac's 
    :math:`\sqrt{F(i,j)^2 + G(i,j)^2}` times the sign of :math:`G(i,j)`...
    
    So we are doing Dirac-Fock, except that we are not treating exchange 
    exactly, in terms of working with major and minor components of the 
    orbitals, and the phe's give the CORRECT CHARGE DENSITY...
    
    The above approximation ought to be very small for valence states,
    so you need not worry about it...
    
    The Breit interaction has been neglected altogether...it should not 
    have a huge effect on the charge density you are concerned with...
    
    '''

    def __init__(self, **kwargs):
        '''
        Constructor
        '''
        self.__dict__.update(kwargs)
        
    @staticmethod
    def get_quantum_info(shell):
        r"""
        Description
        -----------
        Get a tuple of quantum information for a given orbital 's', 'p', 'd' 
        or 'f' from a given subshell string.
        
        Returns
        =======
        tuple : (int, int, list[float, float], list[float, float])
            (n, l, j=[l-s, l+s], occ=[:math:`n^-_r`, :math:`n^+_r`])
            
        Notes
        -----
        - *n* is the principle quantum number (:math:`n > 0`).
        - *l* is the azimuthal quantum number (:math:`0 \leq l \leq n-1`).
        - *s* is the spin quantum number (:math:`s \pm \frac{1}{2}`).
        - *j* is the total angular momentum quantum numbers for both 
          :math:`l-s` or :math:`l+s`, respectively.
        - :math:`n_r` is the occupancy of the spin-split :math:`l-s` 
          and :math:`l+s` levels, respectively. 
        
        Example
        -------
        >>> Atorb.get_quantum_info('3d6')
         (3, 2, [1.5, 2.5], [2.4, 3.6])
        
        """
        
        subshell = "".join([s for s in shell if s.isalpha()])
        try:
            (n, nelectrons) = [t(s) for t, s in zip((int, int),
                                shell.replace(subshell, ' ').split())]
        except ValueError:  # assume 1 electron in shell
            n = int(shell.replace(subshell, ' ').split()[0])
            nelectrons = 1
            
        s = 0.5
        if subshell == 's':
            l = 0
            occ = [nelectrons / 1.0]
            j = [l + s]
            return (n, l, j, occ)
        elif subshell == 'p':
            # 3 subshells
            l = 1
            max_occ = 6 
            occ = []
            for j in [l - s, l + s]:
                occ.append(((2.0 * j) + 1) * nelectrons / max_occ)
            return(n, l, [l - s, l + s], occ)
        elif subshell == 'd':
            # 5 subshells
            l = 2
            max_occ = 10 
            occ = []
            for j in [l - s, l + s]:
                occ.append(((2.0 * j) + 1) * nelectrons / max_occ)
            return(n, l, [l - s, l + s], occ)
        elif subshell == 'f':
            # 7 subshells!
            l = 3
            max_occ = 14
            occ = []
            for j in [l - s, l + s]:
                occ.append(((2.0 * j) + 1) * nelectrons / max_occ)
            return(n, l, [l - s, l + s], occ)

    @staticmethod
    def replace_core_config(electron_config):
        """
        Description
        -----------
        Replace nobel gas core with equivalent electronic shell configuration
        
        Parameters
        ----------
        electron_config : str
            String containing the electronic configuration of the given 
            element.
        
        Returns
        -------
        str :
            A substituted string where the nobel gas core has been replaced.
        
        Examples
        --------
        >>> Atorb.replace_core_config('[Ar] 4s2')
         '1s2 2s2 2p6 3s2 3p6 4s2'

        >>> Atorb.replace_core_config('[Xe] 6s2 5d1')
         '1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 5s2 4d10 5p6 6s2 5d1'
            
        """
        cores = {'[He]': '1s2', '[Ne]': '1s2 2s2 2p6', 
                 '[Ar]': '1s2 2s2 2p6 3s2 3p6', 
                 '[Kr]': '1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6',
                 '[Xe]': '1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 5s2 4d10 5p6', 
                 '[Rn]': '1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 5s2 4d10 5p6'
                          '4f14 5d10 6s2 6p6'}
        core = electron_config.split()[0]

        if core in cores:
            return electron_config.replace(core, cores.get(core))
        else:
            return electron_config

    @staticmethod
    def gen_input(element, **kwargs):
        """
        Description
        -----------
        Generate atorb input file from <element> and optional **kwargs 
        arguments. Returns filename of input file once generated.
        
        Parameters
        ----------
        element : int or str
            Either the atomic number, symbol or name for a given element
        output : str, optional
            File string for atomic orbital output (default: 'at_<symbol>.i')
        ngrid : int, optional 
            Number of points in radial grid (default: 1000)
        rel : bool, optional
            Specify whether to consider relativistic effects
        filename : str, optional
            Name for generated input file (default: 'atorb')
        header : str, optional
            Comment at beginning of input file
        method : str, optional
            Exchange correlation method using either 0.0=Hartree-Fock,
            1.0=LDA, -alpha = xalpha (default: 0.0)
        relic : float, optional
            Relic value for calculation (default: 0)
        mixing_SCF : float, optional
            Self consisting field value (default: 0.5)
        tolerance : float, optional
            Eigenvalue tolerance (default: 0.0005)
        ech : float, optional 
            (default: 100)
        
        Example
        -------
        >>> Atorb.gen_input('H',rel=False,filename="atorb.txt",output='at_H.i')
        >>> with open('atorb_H.txt', 'r') as f: print("".join(f.readlines())
         C*********************************************************************
         C  atorb input file: atorb_H.txt.
         C*********************************************************************
         i
         1 1000                   ! Z NR (number of points in radial grid)
         d
         0                        ! 1=rel, 0=n.r.
         x
         0.d0                     ! 0.d0=HF, 1.d0=LDA, -alfa = xalfa...
         a
         0 1 0.5 0.0005 100       ! relic,levels,mixing SCF, eigen. tol,for ech
         1 0 0 -0.5 1 1.0         ! n, l, l, -j, <1>, occupation
         w
         at_H.i
         q
         

        """
        
        ele = elements.ELEMENTS[element]
        Z = ele.protons
        
        # get full electronic configuration
        config = Atorb.replace_core_config(ele.eleconfig)

        # get quantum numbers & occupancy for each electronic obrital in atom
        electrons = []
        nlevels = 0
        for shell in config.split():
            (n, l, J, occ) = Atorb.get_quantum_info(shell)
            for i, j in enumerate(J):
                electrons.append((n, l, l, -j, 1, occ[i]))
                nlevels += 1
        
        # test kwargs and generate output arguments
        if 'output' in kwargs:
            output = kwargs.get('output')
        else:
            output = "at_{0}.i".format(ele.symbol)
        
        if 'ngrid' in kwargs:
            NR = kwargs.get('ngrid')
        else:
            NR = 1000  # default grid resolution
            
        if 'rel' in kwargs:
            rel = int(kwargs.get('rel'))
        else:
            rel = 1  # default is relativistic 
        
        if 'filename' in kwargs:
            filename = kwargs.get('filename')
        else:
            filename = "atorb_{0}.txt".format(ele.symbol)
            
        if 'header' in kwargs:
            header = kwargs.value('header')
        else:
            header = "  atorb input file: {0}.".format(
                                            os.path.basename(filename))
            
        if 'method' in kwargs:
            method = str(kwargs.value('method'))
            methods_dict = {'HF': '0.d0', 'LDA': '1.d0', 'xalpha': '-alpha'}
            if method in methods_dict:
                method = methods_dict.get('method')
            else:
                method = '0.d0'
        else:
            method = '0.d0'
            
        if 'relic' in kwargs:
            relic = float(kwargs.get('relic'))
        else:
            relic = 0
        
        if 'mixing_SCF' in kwargs:
            mixing_SCF = float(kwargs.get('mixing_SCF'))
        else:
            mixing_SCF = 0.5
            
        if 'tolerance' in kwargs:
            eigen_tol = float(kwargs.get('tolerance'))
        else:
            eigen_tol = 0.0005
            
        if 'ech' in kwargs:
            ech = int(kwargs.get('ech'))
        else:
            ech = 100
            
        # produce output file
        with open(filename, 'w') as f:
            f.write("C".ljust(70, '*') + "\n")
            f.write("C" + str(header) + "\n")
            f.write("C".ljust(70, '*') + "\n")
            f.write('i\n')
            f.write('{0} {1}'.format(Z, NR).ljust(30, ' ')
                    + '! Z NR (number of points in radial grid)\n')
            f.write('d\n')
            f.write('{0}'.format(rel).ljust(30) + '! 1=rel, 0=n.r.\n')
            f.write('x\n')
            f.write('{0}'.format(method).ljust(30) +
                        '! 0.d0=HF, 1.d0=LDA, -alfa = xalfa...\n')
            f.write('a\n')
            f.write('{0} {1} {2} {3} {4}'.format(
                        relic, nlevels, mixing_SCF, eigen_tol, ech).ljust(30)
                        + '! relic,levels,mixing SCF, eigen. tol,for ech.\n')
            for i in range(0, nlevels):
                f.write('{0} {1} {2} {3} {4} {5}'.format(*electrons[i]
                            ).ljust(30) + '! n, l, l, -j, <1>, occupation\n')
            f.write('w\n')
            f.write('{0}\n'.format(output))
            f.write('q\n')
        
        return filename  # return output filename for further use
        
    @staticmethod
    def calculate_Q_density(**kwargs):
        """
        Description
        -----------
        Calculate the radial charge density of a given element or atorb input 
        file.
        
        Usage
        -----
        Atorb.calculate_Q_density(**kwargs)
        
        Parameters
        ----------
        kwargs may be any of the following.
        
        element : int or str, optional
            Generate element atorb input file on the fly. Additional
            kwargs may be used to govern the structure of the input
            file - please use ``help(phaseshifts.Atorb.gen_input)`` 
            for more information. 
        input : str, optional
            Specify atorb input file otherwise will use the class
            instance value.
        output_dir : str, optional
            Specify the output directory for the `at_*.i` file
            generated, otherwise the default current working directory 
            is used.
        
        Returns
        -------
        str : filename
        
        Examples
        --------
        >>> Atorb.calculate_Q_density(input='atorb_C.txt')
              18.008635    -33.678535
               4.451786    -36.654271
               1.569616    -37.283660
               0.424129    -37.355634
               0.116221    -37.359816
               0.047172    -37.360317
               0.021939    -37.360435
               0.010555    -37.360464
               0.005112    -37.360471
               0.002486    -37.360473
               0.001213    -37.360473
               0.000593    -37.360473
               0.000290    -37.360474
            N L M J S OCC.
            1   0 0  -1/2   1    2.0000        -11.493862
            2   0 0  -1/2   1    2.0000         -0.788618
            2   1 1  -1/2   1    0.6667         -0.133536
            2   1 1  -3/2   1    1.3333         -0.133311
         TOTAL ENERGY =      -37.360474  -1016.638262

        >>> Atorb.calculate_Q_density(element='H')
               0.500007     -0.343752
               0.152392     -0.354939
               0.065889     -0.357254
               0.028751     -0.357644
               0.012732     -0.357703
               0.005743     -0.357711
               0.002641     -0.357712
               0.001236     -0.357713
               0.000587     -0.357713
               0.000282     -0.357713
         N L M J S OCC.
            1   0 0  -1/2   1    1.0000         -0.229756
         TOTAL ENERGY =       -0.357713     -9.733932


        """
        if 'input' in kwargs:
            inp = os.path.abspath(kwargs.pop('input'))
            
        if 'element' in kwargs:
            inp = os.path.abspath(
                    Atorb.gen_input(element=kwargs.pop('element'), **kwargs))
        
        current_dir = os.path.curdir
        if 'output_dir' in kwargs:
            output_dir = kwargs.pop('output_dir')    
            if os.path.isdir(output_dir):
                os.chdir(output_dir)
            else:
                try:
                    os.makedirs(output_dir)
                    os.chdir(output_dir)
                except:
                    raise IOError
                
        hartfock(inp)  # calculates atomic orbital charge densities for atom
        
        # get output filename
        lines = []
        output_filename = 'atorb'
        
        try:
            with open(inp, 'r') as f:
                lines = [line for line in f]
        except IOError:
            raise IOError
        
        lines = [line.replace('\n', '').replace('\r', '').lstrip(' ').split(
                          '!')[0].split('#')[0].rstrip(' ') for line in lines]
        
        for i in range(len(lines)):
            if lines[i].startswith('w'):
                output_filename = lines[i + 1]
                break
                     
        os.chdir(current_dir)  # return to original directory
        
        return (os.path.join(output_dir, output_filename) 
                    if output_dir != None else output_filename)
