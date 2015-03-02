#!/usr/bin/env python
# encoding: utf-8

##############################################################################
# Author: Liam Deacon                                                        #
#                                                                            #
# Contact: liam.deacon@diamond.ac.uk                                         #
#                                                                            #
# Copyright: Copyright (C) 2014-2015 Liam Deacon                                  #
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
Provides wrapper classes for phaseshift calculation backends
'''

from __future__ import print_function, unicode_literals
from __future__ import absolute_import, division, with_statement

import sys
import os
import tempfile

from glob import glob
from shutil import copy

from phaseshifts import model, atorb
from phaseshifts.leed import Converter, CLEED_validator, CSearch
from phaseshifts.lib.libphsh import phsh_rel, phsh_wil, phsh_cav
from phaseshifts.conphas import Conphas

class VHTWrapper(object):
    '''Wrapper class to easily generate phase shifts'''
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    @staticmethod
    def autogen_from_input(bulk_file, slab_file, tmp_dir=None,
                           model_name=None, **kwargs):
        '''
        Description
        -----------
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
        '''
        dummycell = model.Unitcell(1, 2, [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        if model_name == None:
            model_name = 'atomic'

        # check formatting
        if 'format' in kwargs:
            out_format = kwargs['format'].lower()
        else:
            out_format = None

        # check lmax
        if 'lmax' in kwargs:
            lmax = kwargs['lmax']
        else:
            lmax = 10

        # check for intermediate storage directory, temp folder otherwise
        tmp_dir = str(tmp_dir)  # ensure string does not have escape chars
        if not os.path.isdir(tmp_dir):
            tmp_dir = tempfile.gettempdir()

        # Load bulk model and calculate MTZ
        bulk_mtz = model.MTZ_model(dummycell, atoms=[])
        if CLEED_validator.is_CLEED_file(bulk_file):
            bulk_mtz = Converter.import_CLEED(bulk_file, verbose=VERBOSE)
            full_fname = glob(os.path.expanduser(os.path.expandvars(
                                                bulk_file)))[0]
            bulk_file = os.path.join(tmp_dir,
                            os.path.splitext(os.path.basename(full_fname))[0]
                            + '_bulk.i')
            bulk_mtz.gen_input(filename=bulk_file)
        else:
            bulk_mtz.load_from_file(bulk_file)

        # Load slab model and calculate MTZ
        slab_mtz = model.MTZ_model(dummycell, atoms=[])
        if CLEED_validator.is_CLEED_file(slab_file):
            slab_mtz = Converter.import_CLEED(slab_file)
            full_fname = glob(os.path.expanduser(os.path.expandvars(
                                                slab_file)))[0]
            slab_file = os.path.join(tmp_dir,
                            os.path.splitext(os.path.basename(full_fname))[0]
                            + '_slab.i')
            slab_mtz.gen_input(filename=slab_file)
        else:
            slab_mtz.load_from_file(slab_file)

        # generate atomic charge densities for each element in bulk model
        if not isinstance(bulk_mtz, model.MTZ_model):
            raise AttributeError("bulk_mtz is not an MTZ_model")

        # get unique elements in bulk and slab
        atomic_dict = {}
        bulk_elements = [atom.element.symbol for atom in bulk_mtz.atoms]
        slab_elements = [atom.element.symbol for atom in slab_mtz.atoms]
        for elem in set(bulk_elements + slab_elements):
            at_file = os.path.join(tmp_dir, "at_%s.i" % elem)
            if not os.path.isfile(at_file):
                print('\nCalculating atomic charge density for %s...' % elem)
                atomic_dict[elem] = atorb.Atorb.calculate_Q_density(
                                        element=elem, output_dir=tmp_dir)
            else:
                atomic_dict[elem] = at_file

        # prepare at files for appending into atomic file
        bulk_at_files = [atomic_dict[atom.element.symbol]
                         for atom in set(bulk_mtz.atoms)]

        # create atomic.i input file from mtz model
        bulk_model_name = os.path.basename(os.path.splitext(bulk_file)[0])
        bulk_atomic_file = bulk_mtz.gen_atomic(input_files=bulk_at_files,
                                            output_file=os.path.join(tmp_dir,
                                                bulk_model_name + '_bulk.i'))

        if VERBOSE:
            print("\nModel")
            print("bulk atoms: %s" % [s for s in bulk_mtz.atoms])
            print("slab atoms: %s" % [s for s in slab_mtz.atoms])

        # calculate muffin-tin potential for bulk model
        print('\nCalculating bulk muffin-tin potential...')
        if VERBOSE:
            print("\tcluster file: '%s'" % bulk_file)
            print("\tatomic file: '%s'" % bulk_atomic_file)
            print("\tslab calculation: '%s'" % str(False))
            print("\toutput file: '%s'" % os.path.join(tmp_dir,
                                                bulk_model_name + '.bmtz'))
            print("\tmufftin file: '%s'" % os.path.join(tmp_dir,
                                            bulk_model_name + '_mufftin.d'))

        bulk_mtz_file = bulk_mtz.calculate_MTZ(
                                    cluster_file=bulk_file,
                                    atomic_file=bulk_atomic_file,
                                    slab=False,
                                    output_file=os.path.join(tmp_dir,
                                            bulk_model_name + '.bmtz'),
                                    mufftin_file=os.path.join(tmp_dir,
                                            bulk_model_name + '_mufftin.d'))
        print('Bulk MTZ = %f' % bulk_mtz.mtz)

        # prepare at files for appending into atomic file
        slab_at_files = [atomic_dict[atom.element.symbol]
                         for atom in set(slab_mtz.atoms)]

        # create atomic.i input file from mtz model
        slab_model_name = os.path.basename(os.path.splitext(slab_file)[0])
        slab_atomic_file = slab_mtz.gen_atomic(input_files=slab_at_files,
                                            output_file=os.path.join(tmp_dir,
                                            slab_model_name + '_slab.i'))

        # calculate muffin-tin potential for slab model
        mufftin_filepath = os.path.join(tmp_dir,
                                            slab_model_name + '_mufftin.d')
        print('\nCalculating slab muff-tin potential...')
        if VERBOSE:
            print("\tcluster file: '%s'" % slab_file)
            print("\tatomic file: '%s'" % slab_atomic_file)
            print("\tslab calculation: %s" % str(True))
            print("\toutput file: '%s'" % os.path.join(tmp_dir,
                                                slab_model_name + '.bmtz'))
            print("\tmufftin file: '%s'" % os.path.join(tmp_dir,
                                            mufftin_filepath))
            print("\tmtz value: %s" % str(bulk_mtz.mtz))

        slab_mtz_file = slab_mtz.calculate_MTZ(cluster_file=slab_file,
                                output_file=os.path.join(tmp_dir,
                                            slab_model_name + '.mtz'),
                                atomic_file=slab_atomic_file,
                                mufftin_file=mufftin_filepath,
                                mtz_string=str(bulk_mtz.mtz), slab=True)

        # create raw phase shift files
        print("\nGenerating phase shifts from '%s'..."
              % os.path.basename(mufftin_filepath))
        filepath = os.path.join(tmp_dir, slab_model_name)
        phasout_filepath = filepath + '_phasout.i'
        dataph_filepath = filepath + '_dataph.d'

        phaseshifts = [atom.tag for atom in set(slab_mtz.atoms+bulk_mtz.atoms)]
        phasout_files = [os.path.join(tmp_dir, atom.tag + '.ph')
                         for atom in set(slab_mtz.atoms+bulk_mtz.atoms)]
        phsh_files = []

        # assign phase shift specific lmax values
        lmax_dict = {}
        for atom in set(slab_mtz.atoms + bulk_mtz.atoms):
            try:
                lmax_dict[atom.tag] = atom.lmax
            except AttributeError:
                print(atom.tag, 'default lmax used:', lmax)
                lmax_dict[atom.tag] = lmax

        try:
            if (slab_mtz.nform == 0
                or str(slab_mtz.nform).lower().startswith('cav')):
                # calculate phase shifts
                phsh_cav(mufftin_filepath, phasout_filepath,
                         dataph_filepath, filepath + '_zph.o')

                # split phasout
                phasout_files = Conphas.split_phasout(
                                    filename=phasout_filepath,
                                    output_filenames=phasout_files)

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

            if (slab_mtz.nform == 1
                or str(slab_mtz.nform).lower().startswith('wil')):
                # calculate phase shifts
                phsh_wil(mufftin_filepath, phasout_filepath,
                         dataph_filepath, filepath + '_zph.o')

                # split phasout
                phasout_files = Conphas.split_phasout(
                                    filename=phasout_filepath,
                                    output_filenames=phasout_files)

            if (slab_mtz.nform == 2
                or str(slab_mtz.nform).lower().startswith('rel')):
                # check energy range
                if 'range' in kwargs:
                    try:

                        # get values
                        with open(mufftin_filepath, 'r') as f:
                            lines = [line for line in f]

                        ei, de, ef, lsm, vc = [t(s) for t, s in zip(
                                    (float, float, float, int, float),
                                    lines[1].replace('D', 'E').split()[:5])]

                        # assign new values
                        (ei, ef, de) = [t(s) for t, s in zip(
                                       (float, float, float),
                                       kwargs['range'])]

                        # edit energy range
                        lines[1] = str('%12.4f%12.4f%12.4f    %3i    %12.4f\n'
                                    % (ei, de, ef, lsm, vc)).replace('e', 'D')
#                         lines[1] = ff.FortranRecordWriter(
#                                         '(3D12.4,4X,I3,4X,D12.4)'
#                                         ).write([ei, de, ef, lsm, vc]) + '\n'

                        with open(mufftin_filepath, 'w') as f:
                            f.write("".join([str(line) for line in lines]))

                    except any as e:
                        sys.stderr.write('Unable to change phase shift energy '
                                         'range - using Barbieri/Van Hove '
                                         'default of 20-300eV in 5eV steps\n')
                        sys.stderr.flush()

                # calculate phase shifts
                #print("Current time " + time.strftime("%X"))
                phsh_rel(mufftin_filepath, phasout_filepath,
                         dataph_filepath, filepath + '_inpdat.txt')
                #print("Current time " + time.strftime("%X"))

                # split phasout
                phasout_files = Conphas.split_phasout(
                                    filename=phasout_filepath,
                                    output_filenames=phasout_files)

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

        except AttributeError:
            raise AttributeError("MTZ_model has no NFORM (0-2) specified!")

        # copy files to storage location
        if 'store' in kwargs and out_format != 'cleed':
            if kwargs['store'] != '.':
                dst = os.path.abspath(os.path.expanduser(os.path.expandvars(
                        kwargs['store'])))
            else:
                dst = os.path.abspath('.')
            Wrapper._copy_files(phsh_files, dst, verbose=True)

        elif 'CLEED_PHASE' in os.environ and out_format == 'cleed':
            dst = os.path.abspath(os.path.expanduser(os.path.expandvars(
                        '$CLEED_PHASE')))
            Wrapper._copy_files(phsh_files, dst, verbose=True)

        else:
            Wrapper._copy_files(phsh_files, os.path.abspath('.'), verbose=True)

        return phsh_files

    @staticmethod
    def _copy_files(files, dst, verbose=False):
        '''copy list of files into destination directory'''
        # check if using native Windows Python with cygwin
        env = ''
        if platform.system() == 'Windows' and dst.startswith('/cygdrive'):
            if os.environ['CLEED_PHASE'] == dst:
                env = 'CLEED_PHASE='
            dst = '"%s"' % (dst.split('/')[2] + ':' +
                             os.path.sep.join(dst.split('/')[3:]))

        # do check and create directory if needed
        if os.path.isfile(dst):
            dst = os.path.dirname(dst)
        if not os.path.exists(dst):
            try:
                os.makedirs(dst)
            except WindowsError:
                pass

        # copy each phase shift file to directory
        if verbose:
            print("\nCopying files to %s'%s'" % (env, dst))
        for filename in files:
            try:
                copy(filename, dst)
                if verbose:
                    print(os.path.basename(filename))
            except IOError:
                sys.stderr.write("Cannot copy file '%s'\n" % filename)
                sys.stderr.flush()


