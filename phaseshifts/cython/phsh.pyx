#!/usr/local/bin/python2.7
# encoding: utf-8

##############################################################################
# Author: Liam Deacon                                                        #
#                                                                            #
# Contact: liam.m.deacon@gmail.com                                           #
#                                                                            #
# Copyright: Copyright (C) 2013-2016 Liam Deacon                             #
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
**phsh.py** - quickly generate phase shifts

phsh provides convenience functions to create phase shifts files
suitable for input into LEED-IV programs such as SATLEED and CLEED.

Examples
--------
.. code:: bash
   
   phsh.py -i *.inp -b *.bul -f CLEED -S phase_dir


'''
from __future__ import print_function, unicode_literals
from __future__ import absolute_import, division, with_statement

import sys
import os
import tempfile

from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
from glob import glob
from shutil import copy
from subprocess import Popen
from sys import platform
from ctypes import c_char_p, c_int, cdll
import time

# local imports
import atorb
import model
from conphas import Conphas
from leed import Converter, CLEED_validator

# try:  # loading Fortran shared object
#     if platform == 'win32':
#         dll_path = os.path.abspath('./libphsh.dll')
#     else:
#         dll_path = os.path.abspath('./libphsh.so')
#          
#     libphsh = cdll.LoadLibrary(dll_path)
#      
#     phsh_cav = libphsh.phsh_cav_
#     phsh_wil = libphsh.phsh_wil_
#     phsh_rel = libphsh.phsh_rel_
#      
#     phsh_cav.argtypes = [c_char_p, c_char_p, c_char_p, c_char_p]
#     phsh_wil.argtypes = [c_char_p, c_char_p, c_char_p, c_char_p]
#     phsh_rel.argtypes = [c_char_p, c_char_p, c_char_p, c_char_p]
#      
#      
# except:  # fallback to pyd object
#     sys.stderr.write('Cannot load %s - falling back to pyd version\n' % dll_path)
#     sys.stderr.flush()
#     try:
#         from libphsh import phsh_cav, phsh_wil, phsh_rel
#     except ImportError:
#         sys.stderr.write('Cannot import libphsh - Exiting...\n')
#         sys.exit(1)

from libphsh import phsh_cav, phsh_wil, phsh_rel

__all__ = []
__version__ = '0.1.1'
__date__ = '2013-11-15'
__updated__ = '2014-02-23'
__contact__ = 'liam.m.deacon@gmail.com'

DEBUG = 0
TESTRUN = 0
PROFILE = 0


class Wrapper(object):
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
            out_format = kwargs['format']
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
            bulk_mtz = Converter.import_CLEED(bulk_file)
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

        # calculate muffin-tin potential for bulk model
        print('\nCalculating bulk muffin-tin potential...')
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
        phasout_files = [os.path.join(tmp_dir, atom.tag + '.ph') 
                         for atom in set(slab_mtz.atoms)]
        phsh_files = []
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
                for phaseshift in phasout_files:
                    filename = os.path.splitext(phaseshift)[0] + '.phs'
                    phsh_files.append(filename)
                    print("\nRemoving pi/2 jumps in '%s':\n" 
                          % os.path.basename(filename))
                    phsh = Conphas(input_files=[phaseshift], 
                                   output_file=filename, 
                                   formatting=out_format, lmax=lmax)
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
                for phaseshift in phasout_files:
                    filename = os.path.splitext(phaseshift)[0] + '.phs'
                    phsh_files.append(filename)
                    print("\nRemoving pi/2 jumps in '%s':" 
                          % os.path.basename(filename))
                    phsh = Conphas(input_files=[phaseshift], 
                                   output_file=filename, 
                                   formatting=out_format, lmax=lmax)
                    phsh.calculate()
                    
        except AttributeError:
            raise AttributeError("MTZ_model has no NFORM (0-2) specified!")

        # copy files to storage location
        if 'store' in kwargs and kwargs['store'] != '.':
            dst = os.path.abspath(os.path.expanduser(os.path.expandvars(
                        kwargs['store'])))
            
            # do check and create directory if needed
            if os.path.isfile(dst):
                dst = os.path.dirname(dst)
            if not os.path.exists(dst):
                try:
                    os.makedirs(dst)
                except WindowsError:
                    pass
            
            # copy each phase shift file to directory
            print("\nCopying files to '%s'" % dst)
            for phsh in phsh_files:
                copy(phsh, dst)
                print(os.path.basename(phsh))
                
        elif 'CLEED_PHASE' in os.environ:
            if (kwargs['format'].lower() == 'cleed' and 
                os.path.isdir(os.environ['CLEED_PHASE'])):
                dst = os.path.abspath(os.path.expanduser(os.path.expandvars(
                        '$CLEED_PHASE')))
                
                # do check and create directory if needed
                if os.path.isfile(dst):
                    dst = os.path.dirname(dst)
                if not os.path.exists(dst):
                    try:
                        os.makedirs(dst)
                    except WindowsError:
                        pass
                
                # copy each phase shift file to directory    
                print("\nCopying files to $CLEED_PHASE:")
                for phsh in phsh_files:
                    copy(phsh, dst)
                    print(os.path.basename(phsh))
        
        elif 'CLEED_PHASE' not in os.environ and kwargs['store'] == '.':
            dst = os.path.abspath('.')
            
            # do check and create directory if needed
            if os.path.isfile(dst):
                dst = os.path.dirname(dst)
            if not os.path.exists(dst):
                try:
                    os.makedirs(dst)
                except WindowsError:
                    pass
            
            # copy each phase shift file to directory
            print("\nCopying files to '%s'" % dst)
            for phsh in phsh_files:
                copy(phsh, dst)
                print(os.path.basename(phsh))
        

class CLIError(Exception):
    '''Generic exception to raise and log different fatal errors.'''
    def __init__(self, msg):
        super(CLIError).__init__(type(self))
        self.msg = "E: %s" % msg
        
    def __str__(self):
        return self.msg
    
    def __unicode__(self):
        return self.msg


def main(argv=None):
    '''Command line options.'''

    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    # display help if no arguments
    if len(argv) == 1:
        argv.append('--help')
    
    program_name = os.path.basename(sys.argv[0])
    program_version = "v%s" % __version__
    program_build_date = str(__updated__)
    program_version_message = '%%(prog)s %s (%s)' % (program_version, 
                                                     program_build_date)
    program_shortdesc = '''phsh.py - quickly generate phase shifts'''
    program_license = '''%s

      Created by Liam Deacon on %s.
      Copyright 2013-2016 Liam Deacon. All rights reserved.

      Licensed under the MIT license (see LICENSE file for details)

      Please send your feedback, including bugs notifications
      and fixes, to: %s

    usage:-
    ''' % (program_shortdesc, str(__date__), __contact__)

    try:
        # Setup argument parser
        parser = ArgumentParser(description=program_license, 
                                formatter_class=RawDescriptionHelpFormatter)
        parser.add_argument('-b', '--bulk', dest='bulk', metavar='<bulk_file>', 
                            help="path to MTZ bulk or CLEED *.bul input file")
        parser.add_argument('-i', '--slab', dest='slab', metavar='<slab_file>', 
                            help="path to MTZ slab or CLEED *.inp input file",
                            required=True)
        parser.add_argument('-t', '--tmpdir', dest='tmpdir', 
                            metavar='<temp_dir>', 
                            help="temporary directory for intermediate "
                            "file generation")
        parser.add_argument('-l', '--lmax', dest='lmax', metavar='<lmax>', 
                            default=10, help="Maximum angular momentum "
                            "quantum number [default: %(default)s]")
        parser.add_argument('-f', '--format', dest='format', 
                            metavar='<format>', default="CLEED",
                            help="Use specific phase shift format "
                            "i.e. 'CLEED' [default: %(default)s]")
        parser.add_argument('-g', '--generate-only', dest='generate', metavar='', 
                            default=False, help="Exit after generating "
                            "phaseshifts; do not launch subprocess using "
                            "PHASESHIFTS_LEED environment variable. "
                            "[default: %(default)s]")
        parser.add_argument('-S', '--store', dest='store', metavar='<subdir>', 
                            default=False,
                            help="Keep intermediate files in subdir when done")
        parser.add_argument("-v", "--verbose", dest="verbose", action="count",
                            help="set verbosity level [default: %(default)s]")
        parser.add_argument('-V', '--version', action='version', 
                            version=program_version_message)

        # Process arguments
        args, unknown = parser.parse_known_args()

        verbose = False
        try:
            verbose = args.verbose
        except:
            pass

        if verbose > 0 and len(unknown) > 0:
            for arg in unknown:
                sys.stderr.write("phsh - warning: Unknown option '%s'\n" 
                                 % arg)
            sys.stderr.flush()
        
        if args.bulk == None:
            args.bulk = str(os.path.splitext(args.slab)[0] + '.bul')
        
        if args.store == False:
            args.store = '.'
        
    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0
    except Exception, e:
        if DEBUG or TESTRUN:
            raise(e)
        indent = len(program_name) * " "
        sys.stderr.write(program_name + ": " + repr(e) + "\n")
        sys.stderr.write(indent + "  for help use --help")
        return 2

    # create phase shifts (warning: black magic within)
    Wrapper.autogen_from_input(args.bulk, args.slab, tmp_dir=args.tmpdir, 
                               lmax=int(args.lmax),
                               format=args.format, store=args.store)
    
    # chain loop commands to next program
    if 'PHASESHIFTS_LEED' in os.environ and not args.generate:
        if os.path.isfile(os.environ['PHASESHIFTS_LEED']):
            leed_cmd = [os.environ['PHASESHIFTS_LEED']]
            for arg in argv:
                leed_cmd.append(arg) 
            
            if verbose:
                print("phsh - starting subprocess: '%s'..." % " ".join(leed_cmd))
            
            # execute subprocess
            Popen(leed_cmd)
    
    return 0

if __name__ == "__main__":
    if DEBUG:
        sys.argv.append("-h")
        sys.argv.append("-v")
        sys.argv.append("-r")
    if TESTRUN:
        import doctest
        doctest.testmod()
    if PROFILE:
        import cProfile
        import pstats
        profile_filename = 'wrapper_profile.txt'
        cProfile.run('main()', profile_filename)
        statsfile = open("profile_stats.txt", "wb")
        p = pstats.Stats(profile_filename, stream=statsfile)
        stats = p.strip_dirs().sort_stats('cumulative')
        stats.print_stats()
        statsfile.close()
        sys.exit(0)
    sys.exit(main())
