#!/usr/bin/env python
# encoding: utf-8

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
**conphas.py**

Provides a native python version of the conphas (phsh3) FORTRAN program 
by W. Moritz, which is distributed as part of the SATLEED code 
(see "Barbieri/Van Hove phase shift calculation package" section) and can
be found at: http://www.icts.hkbu.edu.hk/surfstructinfo/SurfStrucInfo_files/
leed/leedpack.html

The Conphas() class also provides a number of convenience functions (see 
docstrings below). 

Examples
--------
>>> from os.path import join
>>> from phaseshifts.conphas import Conphas
>>> con = Conphas(output_file=join('testing', 'leedph_py.d'),
                 lmax=10)
>>> con.set_input_files([join('testing', 'ph1')])
>>> con.set_format('cleed')
>>> con.calculate()

"""

import ntpath
import os
import sys
import platform
import re
from numpy import loadtxt
from math import pi
from getpass import getuser
from time import gmtime, strftime
from copy import copy, deepcopy

try:
    import StringIO
except ImportError:
    import io.StringIO as StringIO


# Constants
HARTREE = 27.21  # 139 eV in Van Hove LEED program
VERBOSE = 1


class Conphas():
    """Class Conphas
    
    Description
    -----------
    Generates continuous phase shifts (as a function of energy) 
    from phases with discontinuities (jumps UM +/- PI). It reformats 
    scattered phases and energies to use as input datasets for LEED programs. 
    
    Notes
    -----
    This work is based on the original conphas (phsh3) FORTRAN program by 
    W. Moritz, which is distributed as part of the SATLEED code (see 
    "Barbieri/Van Hove phase shift calculation package" section) and can be 
    found at: http://www.icts.hkbu.edu.hk/surfstructinfo/SurfStrucInfo_files/
    leed/leedpack.html

    """

    def __init__(self, input_files=[], output_file=[], formatting=None, 
                 lmax=10, **kwargs):
        """
        Parameters
        ----------
        input_files : list(str)
            a list of the input files (default None)
        output_file : str
            output filepath - this file will have all pi/2 jumps 
            removed to form continuous phase shifts
        lmax (optional) : int
            Maximum angular momentum quantum number to calculate and 
            must be in the range 0 <= lmax <= 18 .
        formatting (optional) : str 
            output style - use either 'CLEED' or None (default None)

        """
        self.input_files = [filename for filename in input_files 
                            if os.path.isfile(filename)]
        self.output_file = ntpath.abspath(str(output_file))
        if int(lmax) >= 0 and int(lmax) <= 18:
            self.lmax = int(lmax)
        self.set_format(formatting)
        self.__dict__.update(kwargs)

    # Fix for escape characters in Windows-based paths
    # also removes invalid characters
    def __fix_path(self, file_path):
        """Fix escaped characters in filepath"""
        if platform.system() is 'Windows':
      
            file_path = os.path.abspath(file_path)
            fix_list = {'\a': '\\a', '\b': '\\b',
                        '\f': '\\f', '\n': '\\n',
                        '\r': '\\r', '\t': '\\t',
                        '\v': '\\v', '\\\\': '\\'}
            for fix in fix_list:
                file_path = file_path.replace(fix, fix_list[fix])
        
            for fix in fix_list:
                file_path = file_path.replace(fix, fix_list[fix])
        
        return "".join(x for x in file_path if x.isalnum() or x in ':\\/-_.')
    
    def read_datafile(self, filename):
        """
        Read in discontinuous phase shift file
        
        Parameters
        ----------
        filename : str
            The path to the discontinuous phase shift file

        """
        if not ntpath.isfile(filename):
            return
        try:
            with open(filename) as f:
                data = []
                data = [data.append(line.replace('-', ' -').replace('\n', 
                                                    '').split()) for line in f]
                data = "".join(line.replace('-', ' -').rstrip() for line in f)
                filename = "C:\\Users\\Liam\\Desktop\\leedph.d"
                data = "".join(line.rstrip() for line in f)
            self.data = loadtxt("C:\\Users\\Liam\\Desktop\\leedph.d", 
                                dtype=float, comments="#")
        except IOError:
            assert IOError
            
    # Set internal data for conphas
    def __set_data(self, data=None):
        if data != None:
            self.data = data
    
    # Load phase shift data from file
    def load_data(self, filename):
        """
        Load (discontinuous) phase shift data from file
        
        Parameters
        ----------
        file : str
           Path to phase shift file.
            
        Returns
        -------
        tuple: (double, double, int, int, ndarray)
           (initial_energy, energy_step, n_phases, lmf, data)
            
        Notes
        -----
        + `initial_energy` is the starting energy of the phase shifts.
        + `energy_step` is the change in energy between consecutive values.
        + `n_phases` is the number of phase shifts contained in the file.
        + `lmf` is the maximum azimuthal quantum number considered.
        + `data` is a (2 x n_phases) array containing the phase shift data.

        """
        with open(filename, 'r') as f:
            title = f.readline()  # skip first line
            (initial_energy, energy_step, n_phases, lmf) = [t(s) 
                for (t, s) in zip((float, float, int, int),
                    f.readline().replace('-', ' -').split())]
            # get parameters
            data_lines = [line.replace('-', ' -').replace('\n', '')
                            for line in f.readlines()]
            data = [float(number) for number in "".join(data_lines).split()]
        return (initial_energy, energy_step, n_phases, lmf, data)
    
    @staticmethod
    def split_phasout(filename, output_filenames=[]):
        """split phasout input file into separate files"""
        try:
            with open(filename, 'r') as phasout:
                lines = phasout.readlines()
        
        except IOError:
            assert IOError("Cannot open file '%s'" % filename)
        
        # get list of phase shifts in phasout
        phsh_list = []
        [phsh_list.append(i) for (i, line) in enumerate(lines)
         if re.match("^[A-Za-z]", line.replace(' ', ''))]
        
        # try to guess filenames from header lines in file
        guessed_filenames = [lines[i].split('#')[0].split(' ')[-1].replace(
                                '\n', '').replace('\r', '') + '.ph' 
                             for i in phsh_list]
        
        # determine list of output filenames
        phsh_filenames = []
        if isinstance(output_filenames, list):
            # generate list of filenames from list
            phsh_filenames = [name for name in output_filenames]
            if len(phsh_filenames) < len(phsh_list):  # not enough names
                for i in range(len(phsh_filenames), len(phsh_list)):
                    phsh_filenames.append(guessed_filenames[i])
        elif isinstance(output_filenames, str):
            # generate list of filenames from trunk filename
            output_filenames = os.path.splitext(output_filenames)[0]
            phsh_filenames = [output_filenames + '_%i.ph' % i 
                              for i, name in enumerate(phsh_filenames)]
        else:
            # try to guess from header lines in file
            phsh_filenames = guessed_filenames
            
        # write separate files for each phase shift
        phsh_list.append(len(lines))
        for i_phsh in range(1, len(phsh_list)):
            try:
                with open(phsh_filenames[i_phsh - 1], 'w') as phsh_file:
                    [phsh_file.write(lines[i]) for i in range(
                        phsh_list[i_phsh - 1], phsh_list[i_phsh])]
            except IOError:
                raise IOError
        
        return phsh_filenames[:len(phsh_list) - 1]  # return written files
         
    def set_input_files(self, input_files=[]):
        """set list of input filenames"""
        if input_files:
            input_files = [self.__fix_path(filename) 
                           for filename in input_files]
            temp_input_files = [filename for filename 
                                in input_files if ntpath.isfile(file)]
            if temp_input_files != None and temp_input_files != []:
                self.input_files = temp_input_files
        
    def set_output_file(self, output_file):
        """set output filename"""
        self.output_file = output_file
                    
    # Set max orbital angular momentum
    def set_lmax(self, lmax):
        """
        Set max orbital angular momentum (azimuthal quantum number)
        
        Parameters
        ----------
        lmax : int
            Maximum azimuthal quantum number to be considered in calculations.
            
        """
        try:
            self.lmax = int(lmax)
        except:
            pass
    
    # set appropriate format from available options 
    def set_format(self, formatting=None):
        """
        Set appropriate format from available options
        
        Parameters
        ----------
        format : str, optional
            The format identifier for different packages; can be 'cleed'
            or None. 

        """
        if str(formatting).lower() in ['cleed', 'curve']:
            self.format = str(formatting).lower()
        else:
            self.format = None
    
    # process to create continuous phase shifts
    def calculate(self):
        """
        Calculates continuous phase shifts from input file(s).
            
        Examples
        --------
        >>> con = Conphas(output_file=r'testing\leedph_py.d', lmax=10)
        >>> con.set_input_files([r'testing\\ph1'])
        >>> con.set_format('cleed')
        >>> con.calculate()
         L = 0
         jump between 25.0 eV and 30.0 eV; IFAK = -1
         L = 1
         jump between 65.0 eV and 70.0 eV; IFAK = -1
         L = 2
         jump between 20.0 eV and 25.0 eV; IFAK = 1
         jump between 80.0 eV and 85.0 eV; IFAK = 0
         L = 3
         L = 4
         jump between 275.0 eV and 280.0 eV; IFAK = 1
         L = 5
         L = 6
         L = 7
         L = 8
         L = 9
         L = 10

        """
        neng = n_phases = 0
        
        energy = []
        phas = []
        conpha = []
        
        # read phase scattering
        for i, input_file in enumerate(self.input_files):
            
            (initial_energy, energy_step, 
                n_phases, lmf, data) = self.load_data(input_file)
            
            if n_phases > 250:
                n_phases = 250
            if i == 0:
                initial_energy0 = copy(initial_energy)
                energy_step0 = copy(energy_step)
                n_phases0 = copy(n_phases)
                lmf0 = copy(lmf)
            
            # increase dimensions of array for input
            phas.append([])
    
            # populate data arrays
            index = 0
            try:
                while index < len(data):
                    for ie in range(n_phases):
                        energy.append(data[index])
                        index += 1
                        phas[i].append([])  # now 3-dimensional
                        for l in range(lmf + 1):
                            phas[i][ie].append(data[index])
                            index += 1
                            
            except IndexError:
                raise IndexError("Malformatted or incomplete input")
            
            neng += n_phases
            conpha = deepcopy(phas)
            
            # produce continuous scattering phases
            for l in range(0, self.lmax + 1):
                if l > 18:
                    break
                ifak = 0
    
                conpha[i][0][l] = phas[i][0][l]
                
                if VERBOSE: 
                    print("L = {0}".format(l))
                for ie in range(1, n_phases):
                    dif = phas[i][ie][l] - phas[i][ie - 1][l]
                    if dif >= pi / 2.:
                        ifak -= 1
                        if VERBOSE: 
                            print("jump between {} eV and {} eV; IFAK = {}"
                                  .format(energy[ie - 1], energy[ie], ifak))
                    elif dif <= -pi / 2.:
                        ifak += 1
                        if VERBOSE: 
                            print("jump between {} eV and {} eV; IFAK = {}"
                                  .format(energy[ie - 1], energy[ie], ifak))
                    conpha[i][ie][l] = phas[i][ie][l] + (pi * ifak)
            
            # get root name of output
            if ntpath.exists(self.input_files[i]):
                root = ntpath.dirname(self.input_files[i])
                name = ntpath.splitext(ntpath.basename(self.input_files[i]))[0]
                dataph = ntpath.join(root, str('dataph_' + name + '.d'))
                leedph = ntpath.join(root, str('leedph_' + name + '.d'))
            else:
                root = ntpath.dirname(self.output_file)
                name = ntpath.splitext(ntpath.basename(self.output_file))[0]
                dataph = ntpath.join(root, 'dataph_{0}_{1}.d'.format(name, i))
                leedph = ntpath.join(root, 'leedph_{0}_{1}.d'.format(name, i))
                
            # write datafile 'dataph.d'
            try:
                with open(dataph, 'w') as f:
                    for kk in range(0, self.lmax + 1):
                        f.write("\"L = {0}\n".format(kk))
                        for ii in range(0, n_phases):
                            f.write("%9.7f\t%9.7f\n" % (
                                    energy[ii] / HARTREE, conpha[i][ii][kk]))
                        f.write("\n")
            except IOError:
                from tempfile import gettempdir
                base = ntpath.basename(dataph)
                with open(os.path.join(gettempdir(), 'phsh', base), 'w') as f:
                    for kk in range(self.lmax + 1):
                        f.write("\"L = {0}\n".format(kk))
                        for ii in range(0, n_phases):
                            f.write("%9.7f\t%9.7f\n" % (
                                    energy[ii] / HARTREE, conpha[i][ii][kk]))
                        f.write("\n")
    
        # write final output
        with open(self.output_file, 'w') as f:
            if str(self.format).lower() == 'cleed': 
                # add formatted header for Held CLEED package
                f.write("{0} {1} neng lmax (calculated by {2} on {3})\n"
                            .format(
                                    neng, self.lmax, getuser(), 
                                    strftime("%Y-%m-%d at %H:%M:%S", gmtime())
                                    )
                        )
            if str(self.format).lower() != 'curve':
                # data format is kept the same for compatibility reasons
                for ie in range(n_phases):
                    f.write("%7.4f\n" % (energy[ie] / HARTREE))
                    # append phase shifts from all files
                    for ii in range(len(self.input_files)):
                        for l in range(self.lmax + 1):
                            f.write("%7.4f" % (conpha[ii][ie][l]))
                    f.write("\n")
            else:  # write output in "x y y ..." format
                # write header
                f.write('# phase shift curve: %s\n' % 
                        os.path.basename(self.output_file))
                f.write('#energy ')
                for l in range(self.lmax + 1):
                    f.write('l=%i ' % l)
                f.write('\n')
                
                # write data
                for ie in range(n_phases):
                    f.write("%7.4f " % (energy[ie] / HARTREE))
                    # append phase shifts from all files
                    for ii in range(len(self.input_files)):
                        for l in range(self.lmax + 1):
                            f.write("%7.4f " % (conpha[ii][ie][l]))
                    f.write('\n')
