#!/usr/bin/env python
# encoding: utf-8
##############################################################################
# Author: Liam Deacon                                                        #
#                                                                            #
# Contact: liam.m.deacon@gmail.com                                           #
#                                                                            #
# Copyright: Copyright (C) 2014-2016 Liam Deacon                             #
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
**wrappers.py**

Provides wrapper classes for phaseshift calculation backends.
"""
from __future__ import print_function, unicode_literals
from __future__ import absolute_import, division, with_statement

import sys
import os
import tempfile

from glob import glob
from abc import ABCMeta, abstractmethod
from getpass import getuser
from time import gmtime, strftime
from re import compile
from copy import deepcopy
from tempfile import gettempdir

from phaseshifts import model
from phaseshifts import atorb
from phaseshifts.leed import Converter, CLEED_validator
from phaseshifts.lib.libphsh import phsh_rel, phsh_wil, phsh_cav
from phaseshifts.conphas import Conphas
from phaseshifts.elements import Element, ELEMENTS
from phaseshifts.model import Atom, Model, MTZModel
from phaseshifts.utils import FileUtils, expand_filepath, stringify


class PhaseShiftWrapper(object):
    """
    Abstract base wrapper class for generating phase shifts
    
    Attributes
    ----------
    lmax : int
        Maximum angular momentum quantum number where 
        :math:`1 \le l_{max} \le 18`.
    format : str
        Formatting of output files.
    range : list of float or tuple of float
        Energy range for phase shift calculations. Takes three values:
        :math:`E_{min}`, :math:`E_{max}` and :math:`\deltaE`.
    min_energy : float
        Minimum energy :math:`E_{min}` for phase shift calculation range.
    max_energy : float
        Maximum energy :math:`E_{max}` for phase shift calculation range.
    energy_step : float
        Energy step :math:`\deltaE` for phase shift calculations.
    tmp_dir : str
        File directory to store intermediate/temporary files.
    slab : str or Model
        Either filepath to slab input file or else a :py:class:`Model` instance 
        describing the model.
    bulk : str or Model
        Either filepath to bulk input file or else a :py:class:`Model` instance 
        describing the model.    
    model_name : str
        Name of the model. Used for naming any generated files.
    store_all_files : bool
        Flag to state whether to keep intermediate files. Used primarily for 
        debugging purposes.
        
    Methods
    -------
    calculate_phaseshifts(model, tmp_dir=None, lmax=None)
        Calculates phase shifts for `model`.
    autogen_from_input(kwargs)
        Generates phase shifts from an input file (either `slab`, `bulk` 
        or both).
    """
    __metaclass__ = ABCMeta
    
    FORMATS = [None, 'bvh', 'eeasisss', 'cleed', 'curve']
    DEFAULTS = {'format': None, 
                'lmax': 10, 
                'range': (20., 600., 5.),
                'slab_file': None,
                'bulk_file': None,
                'tmp_dir': None,
                'store_all_files': False, 
                'model_name': None
                }
    
    @abstractmethod
    def calculate_phaseshifts(self, model, tmp_dir=None, lmax=None):
        """
        Abstract base method for calculating phase shifts for model.
        """
        pass
    
    @abstractmethod
    def autogen_from_input(self, **kwargs):
        """
        Abstract base method for generating phase shifts from an input file.
        """
        pass
    
    @abstractmethod
    def autogen_atorbs(self, elements=[], output_dir='.'):
        """ 
        Abstract base method for generating atomic orbital input for 
        Eric Shirley's hartfock program.
        
        Parameters
        ----------
        elements : set of :py:class:`elements.Element`
            Set of elements.
        output_dir : str
            Destination directory for generated files.
        """
        pass
    
    def __init_defaults__(self):
        for key in self.DEFAULTS:
            try:
                val = self.DEFAULTS[key]
                if not isinstance(val, str):
                    exec('self.{property} = {value}'.format(property=key, 
                                                            value=val))
                else:
                    exec('self.{property} = "{value}"'.format(property=key, 
                                                              value=val))
            except:
                pass
            
    def __init__(self, slab, 
                 bulk=None, 
                 lmax=10, 
                 fmt=None,  
                 energy_range=(20., 600., 5.), 
                 tmp_dir=None, 
                 store=False, 
                 **kwargs):
        self.__init_defaults__()
        
        self.lmax = lmax
        self.fmt = fmt
        self.range = energy_range
        self.tmp_dir = tmp_dir
        self.store_all_files = store
        self.bulk = bulk
        self.slab = slab
        
        self.__dict__.update(kwargs)
    
    @property
    def lmax(self):
        """
        Returns the default maximum angular momentum quantum number used in the 
        phase shift calculations.
        """
        return self._lmax
    
    @lmax.setter
    def lmax(self, lmax):
        """
        Sets the default maximum angular momentum quantum number used in the 
        phase shift calculations.
        """
        try:
            self._lmax = lmax if lmax > 0 and lmax < 19 else self._lmax
        except ValueError:
            pass
    
    @property
    def format(self):
        """
        Returns the format of the generated output. 
        """
        return self._format or None
    
    @format.setter
    def format(self, fmt):
        """
        Sets the format of the generated output. (default: None)
        """
        try:
            self._format = fmt if fmt in self.FORMATS else None
        except:
            pass
    
    @property
    def range(self):
        """
        Returns a tuple of the minimum energy, maximum energy and energy step.
        """
        return self._range
    
    @range.setter
    def range(self, energy_range):
        """
        Sets a tuple of the minimum energy, maximum energy and energy step.
        """
        try:
            self._range = tuple([t(s) for t, s in zip((float, float, float), 
                                                      energy_range)])
        except:
            pass
    
    @property
    def min_energy(self):
        """ 
        Returns the minimum energy used in the phase shift calculations. 
        """
        return self.range[0]
    
    @property
    def max_energy(self):
        """
        Returns the maximum energy used in the phase shift calculations.
        """
        return self.range[1]
    
    @property
    def energy_step(self):
        """
        Returns the energy step used in the phase shift calculations.
        """
        return self.range[2]
    
    @property
    def tmp_dir(self):
        """ 
        Returns the directory that is used to store temporary files. 
        """
        return self._tmp_dir
    
    @tmp_dir.setter
    def tmp_dir(self, tmp_dir):
        """ 
        Sets the directory that is used to store temporary files. 
        (default: :py:meth:`tempfile.gettempdir()`) 
        """
        try:
            self._tmp_dir = (tmp_dir if isinstance(tmp_dir, str) 
                             else gettempdir())
        except:
            pass
    
    @property
    def slab(self):
        """ 
        Returns the slab model or filename
        """
        return self._slab
    
    @slab.setter
    def slab(self, model):
        """ 
        Sets the slab model or filename
        """
        self._set_model(model, "slab")
        
    @property
    def bulk(self):
        """ 
        Returns the bulk model or filename
        """
        return self._bulk
    
    @bulk.setter
    def bulk(self, model):
        """ 
        Returns the bulk model or filename
        """
        self._set_model(model, "bulk")
        
    def _set_model(self, model, model_type='slab'):
        """ 
        Helper method for setting the specified model data member.
        
        Raises
        ------
        Exception
            If `model` is neither a valid input filename nor a Model
            instance.
        ValueError
            If `model_type` is invalid.
        """
        model_types = ['slab', 'bulk']
        if model_type in model_types:
            try:
                if isinstance(model, Model):
                    exec("self._" + model_type + " = model")
                    return
                elif isinstance(model, str):
                    if os.path.isfile(expand_filepath(model)): 
                        exec("self._" + model_type + " = model")
                        return 
                raise Exception
            except:
                raise Exception("%s is neither a valid input filename "
                                "nor Model" % model) 
        else:
            raise ValueError("model_type must be one of: '" + 
                             "' '".join(model_types) + "'")        
        
    @property
    def model_name(self):
        """
        Returns the name of the model.
        """
        return self._model_name
    
    @model_name.setter
    def model_name(self, model):
        """
        Sets the name of the model.
        """
        try:
            if isinstance(model, str):
                self._model_name = model
            elif isinstance(model, Model):
                self._model_name = model.name 
        except:
            pass
        
    @property
    def store_all_files(self):
        """
        Returns boolean value stating whether to store all 
        intermediate generated files during the phase shifts 
        calculation(s). Used primarily for debugging purposes.
        """
        return self._store_all_files
    
    @store_all_files.setter
    def store_all_files(self, store):
        """
        Sets boolean value stating whether to store all 
        intermediate generated files during the phase shifts 
        calculation(s). Used primarily for debugging purposes.
        """
        try:
            self._store_all_files = bool(store)
        except:
            pass
    
    @staticmethod
    def _copy_phsh_files(phsh_files, store='.', out_format=None):
        """
        Copies a list of phase shift files
        """ 
        out_format = None if out_format is None else out_format.lower()
        
        # change file extension for CLEED phase shifts
        if out_format == 'cleed':
            phsh_files = [os.path.splitext(f)[0] + '.phs' for f in phsh_files]
        
        if out_format != 'cleed':
            dst = expand_filepath(str(store)) if store != '.' else '.'
            dst = dst if os.path.isdir(dst) else '.'
            dst = os.path.abspath(dst)
            FileUtils.copy_files(phsh_files, dst, verbose=True)
            
        elif 'CLEED_PHASE' in os.environ and out_format == 'cleed':
            dst = os.path.abspath(expand_filepath('$CLEED_PHASE'))
            FileUtils.copy_files(phsh_files, dst, verbose=True)

        else:
            FileUtils.copy_files(phsh_files, 
                                 os.path.abspath('.'), verbose=True)
            
    def _add_header(self, phsh_file=None, fmt=None):
        """
        Prepends a header line to the beginning of the phase shift file. 
        It assumes that the number of lines in the file are the number of phase 
        shifts (excluding any trailing blank lines).
        
        Parameters
        ----------
        phsh_file : str
            Path to raw phase shift file to prepend header to. If no file is 
            given then the function simply returns a string containing the 
            header line.
        fmt : str
            Format of the added header line. Currently only :py:obj:`'cleed'` 
            and :py:obj:`None` are supported. 
            
        Returns
        -------
        str
            Header line string.
        
        Raises
        ------
        IOError
            If `phsh_file` cannot be read or is empty.
            
        Notes
        -----
        If `neng` and `lmax` cannot be found in the object dictionary then 
        a best estimate will be made from the input file given.
        
        """
        if str(self.format).lower() == 'cleed': 
            # add formatted header for Held CLEED package
            header = ''
            neng = self.neng if 'neng' in self.__dict__ else 0
            lmax = self.lmax if 'lmax' in self.__dict__ else 0
            if phsh_file is not None:
                lines = []
                
                # get old content from file
                with open(phsh_file, 'r') as f:
                    lines = f.readlines()
                    
                if lines.count('\n') == 0:
                    raise IOError("phase shift file '%s' is empty" % phsh_file)
                    
                # remove trailing lines
                while lines[-1] == '\n':
                    lines.pop(-1)  # removes last element
                
                # remove comments and empty lines for input
                stripped_lines = [s for s in lines 
                                  if s != '\n' and s not in 
                                  compile('[a-zA-Z]').findall(lines)]
                
                if lmax == 0:
                    # very crude count of l elements in input line
                    # uses multiple input lines to get best estimate
                    n = 10 if len(lines) > 10 else len(lines)
                    for i in range(n):
                        num = self._get_phaseshift_input(stripped_lines[i])
                        num = len(num.split())
                        if num > lmax:
                            lmax = num
                
                if neng == 0:
                    # guess that neng is half number of non-empty input lines    
                    neng = "".join(stripped_lines).count('\n') / 2
                    
                header = ("%i %i neng lmax (calculated by %s on %s)\n" 
                          % (neng, lmax, getuser(), 
                             strftime("%Y-%m-%d at %H:%M:%S", gmtime()))
                          )
                lines = lines.insert(0, header)
                
                # write new contents to file
                with open(phsh_file, 'w') as f:
                    f.writelines(lines)
                    
            return header
        
        return ""  # zero length string
    
    def _get_phaseshift_input(self, line):
        """
        Parses a Fortran formatted string containing (multiple) decimal numbers
        and returns a Python-friendly line string to be processed further.
        """
        regex = '[-+0-9]{1,5}\.\d{1,6}'  # basic decimal number
        decimal = compile('({})'.format(regex))
        decimal_and_exponent = compile('({}[eED][+-\d]\d{0,6})'.format(regex))
        
        output = decimal.findall(line)
        output = output if output != [] else decimal_and_exponent.findall(line)
        output = [s.replace('D', 'E') for s in output]
        return " ".join(output)
    

class EEASiSSSWrapper(PhaseShiftWrapper):
    """
    Wrapper class to easily generate phase shifts using EEASiSSS backend
    """
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    def autogen_atorbs(self, elements=[], output_dir='.'):
        """ 
        Generates chgden input files for a set of elements and calculates their 
        atomic charge densities using the EEASiSS variant of Eric Shirley's 
        hartfock program.
        
        Parameters
        ----------
        elements : set of :py:class:`elements.Element`
            Set of elements.
        output_dir : str 
            Destination directory for generated files.
        
        Returns
        -------
        dict   
            Dictionary of atorb filepaths for all elements
        """
        EEASiSSSWrapper.calculate_Q_density(elements, output_dir=output_dir)
        atomic_dict = {}
        for symbol in [element.symbol for element in elements]:
            atomic_dict[symbol] = os.path.join(output_dir, "chgden%s" % symbol)
        return atomic_dict

    @staticmethod
    def autogen_from_input(bulk_file, slab_file, tmp_dir=None,
                           model_name=None, lmax=10, verbose=False, **kwargs):
        """
        Generate phase shifts from a slab/cluster input file.

        Parameters
        ----------
        slab_file : str
            Path to the cluster slab input file.
        bulk_file : str
            Path to the cluster bulk input file.
        tmp_dir : str
            Temporary directory for intermediate files.
        store : bool or int
            Specify whether to keep generated files.
        format : str
            Specify formatting of generated phase shift files
        range : tuple(float, float, float)
            Specify the energy of the start, stop and step in eV.
        model_name : str
            Name of model.
        verbose : bool
            Determines whether to print output to screen.

        Returns
        -------
        output_files : list of str
           A list of phase shift output filenames
        """
        dummycell = model.Unitcell(1, 2, [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        
        # check formatting of arguments
        out_format = kwargs['format'].lower() if 'format' in kwargs else None
        model_name = model_name or 'atomic'
        lmax = lmax if isinstance(lmax, int) else 10

        # check for intermediate storage directory, temp folder otherwise
        tmp_dir = tmp_dir or tempfile.gettempdir() 

        # Load bulk model and calculate MTZ
        bulk_mtz = model.MTZModel(dummycell, atoms=[])
        if CLEED_validator.is_CLEED_file(bulk_file):
            bulk_mtz = Converter.import_CLEED(bulk_file, verbose=verbose)
            full_fname = glob(expand_filepath(bulk_file))[0]
            bulk_root = os.path.splitext(os.path.basename(full_fname))[0]
            bulk_file = os.path.join(tmp_dir, bulk_root + '_bulk.i')
            bulk_mtz.gen_input(filename=bulk_file)
        else:
            bulk_mtz.load_from_file(bulk_file)

        # Load slab model and calculate MTZ
        slab_mtz = model.MTZModel(dummycell, atoms=[])
        if CLEED_validator.is_CLEED_file(slab_file):
            slab_mtz = Converter.import_CLEED(slab_file)
            full_fname = glob(expand_filepath(slab_file))[0]
            slab_root = os.path.splitext(os.path.basename(full_fname))[0]
            slab_file = os.path.join(tmp_dir, slab_root + '_slab.i')
            slab_mtz.gen_input(filename=slab_file)
        else:
            slab_mtz.load_from_file(slab_file)

        # generate atomic charge densities for each element in bulk model
        if not isinstance(bulk_mtz, model.MTZModel):
            raise AttributeError("bulk_mtz is not an MTZModel() instance")

        # get unique elements in bulk and slab
        bulk_elements = [atom.element.symbol for atom in bulk_mtz.atoms]
        slab_elements = [atom.element.symbol for atom in slab_mtz.atoms]
        all_elements = set(bulk_elements + slab_elements)
        atomic_dict = EEASiSSSWrapper.autogen_atorbs(elements=all_elements,
                                                     output_dir=tmp_dir)  

        if verbose:
            print("\nModel")
            print("bulk atoms: %s" % [s for s in bulk_mtz.atoms])
            print("slab atoms: %s" % [s for s in slab_mtz.atoms])

        phsh_files = EEASiSSSWrapper.autogen_from_input(inputX, 
                                                        tmp_dir=tmp_dir, 
                                                        **kwargs)

        # copy files to storage location
        if 'store' in kwargs and out_format != 'cleed':
            if kwargs['store'] != '.':
                dst = os.path.abspath(expand_filepath(kwargs['store']))
            else:
                dst = os.path.abspath('.')
            EEASiSSSWrapper._copy_files(phsh_files, dst, verbose=True)

        elif 'CLEED_PHASE' in os.environ and out_format == 'cleed':
            dst = os.path.abspath(expand_filepath('$CLEED_PHASE'))
            EEASiSSSWrapper._copy_files(phsh_files, dst, verbose=True)

        else:
            EEASiSSSWrapper._copy_files(phsh_files, os.path.abspath('.'), 
                                        verbose=True)

        return phsh_files


class BVHWrapper(PhaseShiftWrapper):
    """
    Wrapper class to easily generate phase shifts 
    using the Barbieri - Van Hove backend.
    """
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    def _adjust_range(self, energy_range, mufftin_file):
        """
        Adjusts the minimum energy, maximum energy and energy step of a 
        mufftin.d file generated using the Barbieri/Van Hove package.
        
        Parameters
        ----------
        energy_range : array_like
            Array containing the min, max and energy step.
        mufftin_file : str
            Filepath to 'mufftin.d' output file to edit.
        
        Returns
        -------
        tuple 
            Tuple containing :math:`E_i`, :math:`E_f` and :math:`\delta E` 
            used.
        """
        try:
            # get values
            with open(mufftin_file, 'r') as f:
                lines = [line for line in f]

            ei, de, ef, lsm, vc = [t(s) for t, s in 
                                   zip((float, float, float, int, float),
                                       lines[1].replace('D', 'E').split()[:5])]

            # assign new values
            (ei, ef, de) = [t(s) for t, s in zip((float, float, float), 
                                                 energy_range)]

            # edit energy range
            lines[1] = str('%12.4f%12.4f%12.4f    %3i    %12.4f\n'
                           % (ei, de, ef, lsm, vc)).replace('e', 'D')

            with open(mufftin_file, 'w') as f:
                f.write("".join([str(line) for line in lines]))

        except any:
            sys.stderr.write('Unable to change phase shift energy '
                             'range - using Barbieri/Van Hove '
                             'default of 20-300eV in 5eV steps\n')
            sys.stderr.flush()
        finally:
            return (ei, ef, de)

    def _get_lmax_dict(self, atoms, default_lmax=10):
        lmax_dict = {}
        for atom in set([atom for atom in atoms if isinstance(atom, Atom)]):
            try:
                lmax_dict[atom.tag] = atom.lmax
            except AttributeError:
                print(atom.tag, 'default lmax used:', default_lmax)
                lmax_dict[atom.tag] = default_lmax
        return lmax_dict
    
    def _phsh_cav(self, mufftin_filepath, phasout_filepath, 
                  dataph_filepath, phasout_files, phaseshifts, lmax_dict,
                  zph_filepath='zph.o', out_format=None, **kwargs):
        """
        Perform Cavendish phase shift calculation 
        """
        # calculate phase shifts
        phsh_cav(mufftin_filepath, phasout_filepath,
                 dataph_filepath, zph_filepath)

        # split phasout
        phasout_files = Conphas.split_phasout(filename=phasout_filepath,
                                              output_filenames=phasout_files)

        phsh_files = []

        # eliminate pi-jumps
        for i, phaseshift in enumerate(phaseshifts):
            filename = os.path.splitext(phasout_files[i])[0]
            if out_format == 'curve':
                filename += '.cur'
            else:
                filename += '.phs'
            phsh_files.append(filename)
            print("\nRemoving pi/2 jumps in '%s':\n"
                  % os.path.basename(filename))
            phsh = Conphas(input_files=[phasout_files[i]],
                           output_file=filename, formatting=out_format,
                           lmax=lmax_dict[phaseshift])
            phsh.calculate()
            
        return phsh_files
    
    def _phsh_wil(self, filepath, mufftin_filepath, phasout_filepath, 
                  dataph_filepath, phasout_files, energy_range, 
                  fmt=None, phaseshifts=[],
                  **kwargs):
        """
        Perform William's phase shift calculation 
        """
        # calculate phase shifts
        phsh_wil(mufftin_filepath, phasout_filepath,
                 dataph_filepath, filepath + '_zph.o')

        # split phasout
        phasout_files = Conphas.split_phasout(filename=phasout_filepath,
                                              output_filenames=phasout_files)
        
    def _phsh_rel(self, 
                  filepath, 
                  mufftin_filepath, 
                  phasout_filepath, 
                  dataph_filepath, 
                  phasout_files, 
                  energy_range, 
                  fmt=None, 
                  phaseshifts=[], 
                  lmax_dict={}, 
                  tmp_dir=gettempdir(), 
                  **kwargs):
        """
        Perform relativistic phase shift calculation 
        """
        # calculate phase shifts
        # print("Current time " + time.strftime("%X"))
        phsh_rel(mufftin_filepath, phasout_filepath,
                 dataph_filepath, filepath + '_inpdat.txt')
        # print("Current time " + time.strftime("%X"))

        # split phasout
        phasout_files = Conphas.split_phasout(filename=phasout_filepath,
                                              output_filenames=phasout_files)

        phsh_files = []

        # eliminate pi-jumps
        for i, phaseshift in enumerate(phaseshifts):
            filename = os.path.splitext(phasout_files[i])[0]
            if fmt == 'curve':
                filename += '.cur'
            else:
                filename += '.phs'
            phsh_files.append(filename)
            print("\nRemoving pi/2 jumps in '%s':\n"
                  % os.path.basename(filename))
            phsh = Conphas(input_files=[phasout_files[i]],
                           output_file=filename, formatting=fmt,
                           lmax=lmax_dict[phaseshift])
            phsh.calculate()
        
        return phsh_files or []

    def autogen_atorbs(self, elements=[], output_dir='.'):
        """ 
        Generates atomic orbital input files for a set of elements and 
        calculates their atomic charge densities according to the Barbieri /
        Van Hove variant of Eric Shirley's hartfock program.
        
        Parameters
        ----------
        elements : set of elements.Element objects.
        output_dir : destination directory for generated files.
        
        Returns
        -------
        dict
            Dictionary of atorb filepaths for all elements.
        """
        atomic_dict = {}
        for elem in [element.symbol if isinstance(element, Element) 
                     else ELEMENTS[element] for element in elements]:
            at_file = os.path.join(output_dir, "at_%s.i" % elem)
            if not os.path.isfile(at_file):
                print('\nCalculating atomic charge density for %s...' % elem)
                chgden = atorb.Atorb.calculate_Q_density
                atomic_dict[elem] = chgden(element=elem, output_dir=output_dir)
            else:
                atomic_dict[elem] = at_file
        return atomic_dict

    @staticmethod
    def autogen_from_input(bulk_file, slab_file, tmp_dir=None,
                           model_name=None, verbose=False, **kwargs):
        """
        Generate phase shifts from a slab/cluster input file.

        Parameters
        ----------
        slab_file : str
            Path to the cluster slab MTZ input file.
        bulk_file : str
            Path to the cluster bulk MTZ input file.
        tmp_dir : str
            Temporary directory for intermediate files.
        store : bool or int
            Specify whether to keep generated files.
        format : str
            Specify formatting of generated phase shift files
        range : tuple(float, float, float)
            Specify the energy of the start, stop and step in eV.
        model_name : str
            Name of model.

        Returns
        -------
        output_files : list(str)
           A list of phase shift output filenames
        """
        dummycell = model.Unitcell(1, 2, [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        
        model_name = model_name or 'atomic'

        # check formatting
        store = kwargs['store'] if 'store' in kwargs else '.'
        out_format = kwargs['format'].lower() if 'format' in kwargs else None
        lmax = kwargs['lmax'] if 'lmax' in kwargs else 10

        # check for intermediate storage directory, temp folder otherwise
        tmp_dir = str(tmp_dir)  # ensure string does not have escape chars
        if os.path.isdir(str(tmp_dir)):
            tmp_dir = tempfile.gettempdir()

        # Load bulk model and calculate MTZ
        bulk_mtz = model.MTZModel(dummycell, atoms=[])
        if CLEED_validator.is_CLEED_file(bulk_file):
            bulk_mtz = Converter.import_CLEED(bulk_file, verbose=verbose)
            full_fname = glob(expand_filepath(bulk_file))[0]
            bulk_root = os.path.splitext(os.path.basename(full_fname))[0]
            bulk_file = os.path.join(tmp_dir, bulk_root + '_bulk.i')
            bulk_mtz.gen_input(filename=bulk_file)
        else:
            bulk_mtz.load_from_file(bulk_file)

        # Load slab model and calculate MTZ
        slab_mtz = model.MTZModel(dummycell, atoms=[])
        if CLEED_validator.is_CLEED_file(slab_file):
            slab_mtz = Converter.import_CLEED(slab_file)
            full_fname = glob(expand_filepath(slab_file))[0]
            slab_root = os.path.splitext(os.path.basename(full_fname))[0]
            slab_file = os.path.join(tmp_dir, slab_root + '_slab.i')
            slab_mtz.gen_input(filename=slab_file)
        else:
            slab_mtz.load_from_file(slab_file)

        # generate atomic charge densities for each element in bulk model
        if not isinstance(bulk_mtz, model.MTZModel):
            raise AttributeError("bulk_mtz is not an MTZModel")

        # get unique elements in bulk and slab
        bulk_elements = [atom.element.symbol for atom in bulk_mtz.atoms]
        slab_elements = [atom.element.symbol for atom in slab_mtz.atoms]
        atomic_dict = BVHWrapper.autogen_atorbs(elements=set(bulk_elements + 
                                                             slab_elements),
                                                output_dir=tmp_dir)

        # prepare at files for appending into atomic file
        bulk_at_files = [atomic_dict[atom.element.symbol]
                         for atom in set(bulk_mtz.atoms)]

        # create atomic.i input file from mtz model
        bulk_model_name = os.path.basename(os.path.splitext(bulk_file)[0])
        model_root = os.path.join(tmp_dir, bulk_model_name) 
        bulk_atomic_file = bulk_mtz.gen_atomic(input_files=bulk_at_files,
                                               output_file=(model_root + 
                                                            '_bulk.i')
                                               )

        if verbose:
            print("\nModel")
            print("bulk atoms: %s" % [s for s in bulk_mtz.atoms])
            print("slab atoms: %s" % [s for s in slab_mtz.atoms])

        # calculate muffin-tin potential for bulk model
        print('\nCalculating bulk muffin-tin potential...')
        if verbose:
            print("\tcluster file: '%s'" % bulk_file)
            print("\tatomic file: '%s'" % bulk_atomic_file)
            print("\tslab calculation: '%s'" % str(False))
            print("\toutput file: '%s'" 
                  % os.path.join(tmp_dir, bulk_model_name + '.bmtz'))
            print("\tmufftin file: '%s'" 
                  % os.path.join(tmp_dir, bulk_model_name + '_mufftin.d'))

        bulk_mtz_file = bulk_mtz.calculate_MTZ(cluster_file=bulk_file,
                                               atomic_file=bulk_atomic_file,
                                               slab=False,
                                               output_file=(model_root + 
                                                            '.bmtz'),
                                               mufftin_file=(model_root + 
                                                             '_mufftin.d')
                                               )
        print('Bulk MTZ = %f' % bulk_mtz.mtz)

        # prepare at files for appending into atomic file
        slab_at_files = [atomic_dict[atom.element.symbol]
                         for atom in set(slab_mtz.atoms)]

        # create atomic.i input file from mtz model
        slab_model_name = os.path.basename(os.path.splitext(slab_file)[0])
        model_root = os.path.join(tmp_dir, slab_model_name) 
        slab_atomic_file = slab_mtz.gen_atomic(input_files=slab_at_files,
                                               output_file=(model_root + 
                                                            '_slab.i')
                                               )

        # calculate muffin-tin potential for slab model
        mufftin_filepath = os.path.join(tmp_dir,
                                        slab_model_name + '_mufftin.d')
        print('\nCalculating slab muff-tin potential...')
        if verbose:
            print("\tcluster file: '%s'" % slab_file)
            print("\tatomic file: '%s'" % slab_atomic_file)
            print("\tslab calculation: %s" % str(True))
            print("\toutput file: '%s'" 
                  % os.path.join(tmp_dir, slab_model_name + '.bmtz'))
            print("\tmufftin file: '%s'" 
                  % os.path.join(tmp_dir, mufftin_filepath))
            print("\tmtz value: %s" % str(bulk_mtz.mtz))

        slab_mtz_file = slab_mtz.calculate_MTZ(cluster_file=slab_file,
                                               output_file=model_root + '.mtz',
                                               atomic_file=slab_atomic_file,
                                               mufftin_file=mufftin_filepath,
                                               mtz_string=str(bulk_mtz.mtz), 
                                               slab=True)

        # create raw phase shift files
        print("\nGenerating phase shifts from '%s'..."
              % os.path.basename(mufftin_filepath))

        # create overall model
        model = deepcopy(slab_mtz) + deepcopy(bulk_mtz)
        model.name = slab_model_name
        
        # calculate phase shifts
        phsh_files = list(BVHWrapper._calculate_phaseshifts(self, model, 
                                                            tmp_dir, 
                                                            lmax, 
                                                            mufftin_filepath)) 

        # copy files to storage location
        BVHWrapper.__bases__[-1]._copy_phsh_files(phsh_files, 
                                                  store, 
                                                  out_format)
        return phsh_files
         
    def _calculate_phaseshifts(self, model, 
                               tmp_dir=None, 
                               lmax=None,
                               energy_range=None,
                               mufftin_filepath=None):
        """
        Calculates phase shifts for model.
        
        Parameters
        ----------
        lmax : int
            Maximum angular momentum quantum number (1 <= `lmax` <= 18).
        tmp_dir : str
            Temporary directory storage path.
            
        Notes
        -----
        If `lmax` or `tmp_dir` is :py:obj:`None` then the class instance 
        data member value is used instead. 
        """
        filepath = os.path.join(tmp_dir, model.name)
        
        # check arguments
        tmp_dir = lmax or self.tmp_dir
        lmax = lmax or self.lmax
        mufftin_filepath = mufftin_filepath or filepath + '_mufftin.d'
        
        # set other paths
        phasout_filepath = filepath + '_phasout.i'
        dataph_filepath = filepath + '_dataph.d'

        atoms = set(model.atoms)
        phaseshifts = [atom.tag for atom in atoms]
        phasout_files = [os.path.join(tmp_dir, phaseshift + '.ph') 
                         for phaseshift in phaseshifts]

        # assign phase shift specific lmax values
        lmax_dict = self._get_lmax_dict(atoms, default_lmax=lmax)
    
        # calculate phase shifts
        nform = str(model.nform)
        
        # select type of calculation
        phsh_method = self._select_method(nform)
        
        return phsh_method(filepath, 
                           mufftin_filepath=mufftin_filepath, 
                           phasout_filepath=phasout_filepath, 
                           dataph_filepath=dataph_filepath, 
                           energy_range=energy_range, 
                           fmt=None, 
                           phaseshifts=phaseshifts,
                           phasout_files=phasout_files,
                           lmax_dict=lmax_dict 
                           )
            
    def _select_method(self, nform):
        """
        Returns the method to invoke for the phase shift calculations
        depending on `nform`.
        
        Raises
        ------
        ValueError
            If `nform` is not one of: {nform_dict}
        """.format(stringify(nform_dict=MTZModel.nforms))
        nform = str(nform)[:3]
        if MTZModel.nforms[nform] == MTZModel.cav:
            return self._phsh_cav
        elif MTZModel.nforms[nform] == MTZModel.wil:
            return self._phsh_wil
        elif MTZModel.nforms[nform] == MTZModel.rel:
            return self._phsh_rel
        else:
            raise ValueError('Invalid MTZModel nform (0-2). Value should be '
                             'one of: {}'.format(stringify(MTZModel.nforms)))
