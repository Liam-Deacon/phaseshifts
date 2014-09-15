#!/usr/bin/env python
# encoding: utf-8

##############################################################################
# Author: Liam Deacon                                                        #
#                                                                            #
# Contact: liam.deacon@diamond.ac.uk                                         #
#                                                                            #
# Copyright: Copyright (C) 2014 Liam Deacon                                  #
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
'''
Provides CLEED validator and Converter classes.

The CLEED_validator() class provides a method for checking 
the input files for errors, whereas the Converter.import_CLEED()
method allows importing CLEED input files as a MTZ_model class 

'''

import os
import sys

import model

try:
    from io import IOBase
except ImportError:
    from io import iobase as IOBase
    
from math import sqrt, pow
from glob import glob


class CLEED_validator(object):
    '''
    Class for validation of CLEED input files
    '''
    def __init__(self):
        '''
        Constructor
        '''
        pass
    
    @staticmethod
    def is_CLEED_file(filename):
        """
        Determine if file is a CLEED input file
        
        Returns
        -------
        True
            if a valid filename ending in any of .bul, .inp, .bsr, .bmin
            
        False
            otherwise 
        """
        
        # expand filename string and match to first file found
        filename = glob(os.path.expanduser(os.path.expandvars(filename)))[0]
        
        # check filename exists
        if not os.path.isfile(filename):
            return False
        
        # check extension is valid
        ext = os.path.splitext(filename)[1]
        if ext not in ['.bmin', '.bul', '.inp', '.pmin']:
            if ext != '.i':
                sys.stderr.write('Warning: %s is not a valid CLEED extension\n' 
                                 % filename)
                sys.stderr.flush()
            return False
        else:
            return True
    
    def validate(self, filename, aoi=False):
        '''
        Checks CLEED input file for errors
        
        Parameters
        ----------
        filename : str
            Path to input file. Should be \*.bul , \*.ctr, \*.inp or \*.bmin
        aoi : bool
            Check for angle of incidence parameters

        '''
        filename = glob(os.path.expanduser(os.path.expandvars(filename)))[0]
        
        basename, ext = os.path.splitext(filename)
        if isinstance(ext, int):
            aoi = True
            basename, ext = os.path.splitext(basename)
        
        if ext in ['.bul', '.bmin', '.bsr']:  # bulk input file
            filetype = 'bulk'
        elif ext in ['.inp', '.pmin']:  # surface input file
            filetype = 'surface'
        elif ext == '.ctr':  # control file
            filetype = 'control'
        elif ext == '.par':  # parameter file
            filetype = 'parameter'
        elif ext == '.mkiv':  # mkiv file
            filetype = 'make_iv'
        else:  # try to determine filetype from file data
            with open(filename, 'r') as f:
                # get input lines, stripping left whitespace and comments
                lines = "".join([line.lstrip() for line in f 
                                 if not line.lstrip().startswith('#')])
                # check for control parameters
                if ('ef=' in lines and 'ti=' in lines and 
                    'id=' in lines and 'wt=' in lines):  
                    # probably control file
                    filetype = 'control'
                elif ('pb:' in lines and ('ei:' in lines or 'ef:' in lines
                      or 'es:' in lines or 'ep:' in lines or 'vi:' in lines
                      or 'vr:' in lines or 'lm:' in lines or 'm1:' in lines
                      or 'm2:' in lines)):  # probably bulk file
                    filetype = 'bulk'
                elif ('po:' in lines and ('a1:' in lines or 'a2:' in lines
                      or 'm1:' in lines or 'm2:' in lines)):  # surface file?
                    filetype = 'surface'
                elif ('po:' in lines and not ('a1:' in lines or 'a2:' in lines
                      or 'm1:' in lines or 'm2:' in lines)):
                    filetype = 'parameter'
                else:
                    raise IOError("unknown CLEED format for file '%s'"
                                  % filename)
              
        with open(filename, 'r') as f:
            # get input lines, stripping left whitespace and comments
            lines = [line for line in f]
            
        with open(filename, 'w') as f:
            for (i, line) in enumerate(lines):
                line = line.lstrip().rstrip()
                if not line.startswith('#') or line.startswith(''):
                    # process line based on file type
                    if line.startswith('po:'):
                        if filetype not in ['surface', 'parameter']:
                            # comment out line
                            sys.stderr.write('warning: commented out line '
                                             '%i\n' % i)
                            sys.stderr.flush()
                            line = '#' + line
                            
                    else:  # line unknown - comment out
                            sys.stderr.write('warning: commented out line '
                                             '%i\n' % i)
                            sys.stderr.flush()
                            line = '#' + line
                    if filetype == 'bulk' or filetype == 'surface':
                        pass
                f.write(line + '\r\n')  # playing it safe with line endings
                sys.stderr.flush()
        # unit cell
        try:
            basis = [line for line in lines if line.startswith('a') and 
                     int(line[1]) <= 3 and int(line[1]) > 0]

        except ValueError:
            raise ValueError("'%s' is not a valid basis vector" 
                             % line[:2])
        
        # try to get minimum radii for mufftin input
        try:
            minimum_radii = ["".join(line.split(':')[1].split()[:4]) 
                            for line in lines if line.startswith('rm:')]
        except ValueError:
            raise ValueError("'%s' is not a valid minimum radius input" 
                             % line[:2]) 
        radii_dict = {}
        for radius in minimum_radii:
            tag, value = zip((str, float), radius.split()[:2])
            radii_dict[tag] = value 
        
        try:
            overlayer_atoms = ["".join(line.split(':')[1].split()[:4]) 
                               for line in lines 
                               if line.startswith('po:')]
        except ValueError:
            raise ValueError("'%s' is not a valid overlayer atom input" 
                             % line[:2])
            
        try:
            bulk_atoms = ["".join(line.split(':')[1].split()[:4]) 
                           for line in lines if line.startswith('pb:')]
        except ValueError:
            raise ValueError("'%s' is not a valid overlayer atom input" 
                             % line[:2])


class Converter(object):
    '''
    Convert different input into phaseshift compatible input
    
    '''
    def __init__(self):
        '''
        Constructor
        '''
        pass
    
    @staticmethod
    def import_CLEED(filename):
        '''
        Imports CLEED input file and converts model to muffin-tin input.
         
        It assumes the following:
          * the basis vectors a1, a2, & a3 are x,y,z cartezian coordinates
          * if no a3 is found, the maximum z distance between atoms multiplied
            by four is given
          * the unitcell is converted from cartezian to fractional coordinates
          * atom sites are converted from Angstrom to Bohr units
          * additional info from the phase shift filename is provided by 
            splitting the '_' chars:
              1. First string segment is element or symbol, e.g. Ni
              2. Second string segment is the oxidation, e.g. +2
          * lines with '**rm:**' provide the radii dictionary of the atomic 
            species
          * if no '**rm:**' found for that species, the atomic radius is used 
            for zero valence, otherwise the covalent radius is used. 
           
        Additional information can, however, be provided using 'phs:' at the 
        start of a line within the input file and may have the following 
        formats:
          1. "**phs:** c *<float>* nh *<int>* nform *<int>* exchange *<float>*"
          2. "**phs:** *<phase_shift>* valence *<float>* radius *<float>*"
          
        The identifiers ``exchange``, ``nform``, ``valence`` and ``radius`` may 
        be abbreviated to ``exc``, ``nf``, ``val`` and ``rad``, respectively.
         
        
        Parameters
        ----------
        filename : str
            Path to input file.
            
        Returns
        -------
        phaseshifts.model.MTZ_model
        
        Raises
        ------
        IOError : filename invalid
        ValueError : bad formatting of input

        '''
        
        if not isinstance(filename, str):
            # careful: passing file object not string
            raise IOError("'%s' is not a string" % filename)  
        
        filename = glob(os.path.expanduser(os.path.expandvars(filename)))[0]
        
        if os.path.isdir(filename):
            # check directory for unique .bul or .bmin files
            raise IOError("'%s' is a directory, not a file" % filename)
        
        elif os.path.isfile(filename):
            # determine type of file
            ext = os.path.splitext(filename)[1]
            
            if ext in ['.bul', '.bmin']:  # bulk input file
                pass
            elif ext in ['.inp', '.pmin']:  # surface input file
                pass
            
        else:
            raise IOError("cannot open '%s'" % filename)
            
        with open(filename, 'r') as f:
            # get input lines, stripping left whitespace and comments
            lines = [line.lstrip() for line in f 
                     if not line.lstrip().startswith('#')]
            
        # initialise model
        atoms = []
        unitcell = model.Unitcell(1, 2, [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        z_dist = 4.0  # calculate 'c' for slab from atom coordinates
        z_min = sys.float_info.max
        z_max = sys.float_info.min

        # unit cell
        basis = [line for line in lines if line.startswith('a') and 
                 int(line[1]) <= 3 and int(line[1]) > 0]
                
        a1 = a2 = a3 = [0.0, 0.0, 0.0]  # intialise basis vectors
        for vector_line in basis:
            line = " ".join(vector_line.replace(':', '').split()[:4])
            
            try:
                (vector_str, x, y, z) = [t(s) for (t, s) in zip(
                                        (str, float, float, float), 
                                          line.split())]
            except ValueError:
                raise ValueError("'%s' is not a valid basis vector" 
                             % vector_line.replace('\n', ''))
       
            exec("%s = [%f, %f, %f]" % (vector_str, x, y, z))
        
        zmin = sys.float_info.max
        a = zmax = sys.float_info.min
        for vector in [a1, a2, a3]:
                    
            # compare magnitude of xy components
            len_a = sqrt(pow(vector[0], 2) + pow(vector[1], 2))
            if len_a > a:
                a = len_a
            
            # find min and max z
            if vector[2] < zmin:
                zmin = vector[2]
                
            if vector[2] > zmax:
                zmax = vector[2]     

        # try to get minimum radii for mufftin input
        try:
            minimum_radii = [" ".join(line.split(':')[1].split()[:4]) 
                            for line in lines if line.startswith('rm:')]
            
            radii_dict = {}
            for radius in minimum_radii:
                tag, value = zip((str, float), radius.split()[:2])
                radii_dict[tag] = value 
            
        except ValueError:
            raise ValueError("'%s' is not a valid minimum radius input" 
                             % line[:2]) 
                    
        # try to get the overlayer atoms
        overlayer_atoms = [" ".join(line.split(':')[1].split()[:4]) 
                           for line in lines if line.startswith('po:')]
            
        # extract details for each overlayer atom
        for overlayer_atom in overlayer_atoms:
            try:
                id_str, x, y, z = [t(s) for t, s in zip((str, float, 
                                   float, float), overlayer_atom.split())]

            except ValueError:
                raise ValueError("'%s' line input is invalid")
            
            # extract further information from id string
            id_str = os.path.basename(os.path.expanduser(
                                        os.path.expandvars(id_str))
                                      )  # ensure path not included
            element = id_str.split('_')[0]  # assume element name / symbol
            oxidation = str(id_str.split('_')[1])  # assume valence is 2nd
            if (oxidation.startswith('-') or oxidation.startswith('+') 
                or oxidation.isdigit()):
                try:
                    oxidation = float(oxidation)
                except ValueError:
                    oxidation = 0.0  # invalid format
            else:
                oxidation = 0.0  # oxidation not specified
            try:
                rm = radii_dict[id_str]
                atom = model.Atom(element, [x, y, z], tag=id_str, 
                              valence=oxidation, radius=rm) 
            except KeyError:
                atom = model.Atom(element, [x, y, z], tag=id_str, 
                              valence=oxidation)

            # update z limits                    
            if z < z_min:
                z_min = z
            if z > z_max:
                z_max = z
            
            # add atom to list
            atoms.append(atom)
        
        # try to get the bulk atoms    
        bulk_atoms = [" ".join(line.split(':')[1].split()[:4]) 
                       for line in lines if line.startswith('pb:')]
        
        # extract details for each bulk atom
        for bulk_atom in bulk_atoms:
            # split line data
            try:
                (id_str, x, y, z) = [t(s) for t, s in zip(
                                      (str, float, float, float), 
                                      bulk_atom.split())]
            except ValueError:
                raise ValueError("'%s' line input is invalid")
            
            # extract further information from id string
            id_str = os.path.basename(os.path.expanduser(
                                        os.path.expandvars(id_str))
                                      )  # ensure path not included
            element = id_str.split('_')[0]  # assume element name / symbol
            oxidation = str(id_str.split('_')[1])  # assume valence is 2nd
            if oxidation.startswith('-') or oxidation.startswith('+'):
                try:
                    oxidation = float(oxidation)
                except ValueError:
                    oxidation = 0.0  # invalid format
            else:
                oxidation = 0.0  # oxidation not specified
            try:
                rm = radii_dict[id_str]
                atom = model.Atom(element, [x, y, z], tag=id_str, 
                              valence=oxidation, radius=rm) 
            except KeyError:
                atom = model.Atom(element, [x, y, z], tag=id_str, 
                              valence=oxidation)

            # update z limits                    
            if z < z_min:
                z_min = z
            if z > z_max:
                z_max = z
            
            # add atom to list
            atoms.append(atom)
            
        # change basis to fractional coordinates (SPA units)
        # i.e. in terms of magnitude 'a'
        for vector in [a1, a2, a3]:
            vector = [x / a for x in vector]
        
        if a3 == [0.0, 0.0, 0.0]:  # no a3 in input file => slab
            z_dist = (z_max - z_min) / float(a) / 0.529
            c = z_dist * 5  # slab calculation
            a3[2] = c
        else:
            c = a3[2] / float(a) / 0.529
        
        # create unitcell
        unitcell = model.Unitcell(a, c, [a1, a2, a3])

        mtz_model = model.MTZ_model(unitcell, atoms)
        
        # read additional phaseshift statement lines 'phs:'
        phs_lines = [line.lower().lstrip("phs:").rstrip("#") 
                     for line in lines if line.lower().startswith("phs:")]
        
        for line in phs_lines:
            strings = line.split()
            phase_shift = strings[0]
            phs_atoms = [atom.tag for atom in mtz_model.atoms
                               if atom.tag == phase_shift] 
            if phase_shift not in phs_atoms:
                for i in range(len(strings)):
                    try:
                        if strings[i] == 'c':
                            c = strings[i] / float(a) / 0.529
                            mtz_model.unitcell.set_c(c)
                        
                        elif strings[i] == 'nh':
                            mtz_model.set_nh(strings[i + 1])
                        
                        elif (strings[i] == 'exc' 
                              or strings[i] == 'exchange'):
                            mtz_model.set_exchange(strings[i + 1])
                        
                        elif (strings[i] == 'nf' or 
                              strings[i] == 'nform'):
                            mtz_model.set_nform(strings[i + 1])
                            
                    except IndexError:
                        break
            else:

                # apply commands to atom 
                for i in range(1, len(strings)):
                    try:
                        if (strings[i] == 'val' or
                            strings[i] == 'valence'):
                            for atom in phs_atoms:
                                atom.valence = float(strings[i + 1])
                        
                        elif (strings[i] == 'rad' or
                              strings[i] == 'radius'):
                            for atom in phs_atoms:
                                atom.radius = float(strings[i + 1])
                        
                    except IndexError:
                        break
                        
        return mtz_model  # converted muffin-tin potential model
