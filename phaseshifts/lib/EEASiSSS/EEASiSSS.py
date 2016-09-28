#!/usr/bin/env python
# encoding: utf-8
#
##############################################################################
# Author: Liam Deacon                                                        #
#                                                                            #
# Contact: liam.m.deacon@gmail.com                                           #
#                                                                            #
# Created on 18 Mar 2015                                                     #
#                                                                            #
# Copyright: Copyright (C) 2015 Liam Deacon                                  #
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
**eeasisss.py** - calculate EEASiSSS phase shifts.

Description
-----------
eeasisss.py is a command line script which calculates phase shifts
using John Rundgren's ``EEASiSSS()`` subroutine.

Warning
-------
The EEASiSSS implementation currently relies on making the subroutine call to
a compiled shared library of the :ref:`EEASiSSS` code (minus the main program),
using ctypes to achieve this. As such the function call may be outdated and is
subject to change.

To do
-----
Implement an f2py wrapped version of the EEASiSSS code that will be safer and
more tightly integrated into the overall :ref:`phaseshifts` package.

Examples
--------
The following calculates the phase shifts from the modelling parameters within
'inputX' and from the atomic charge density file(s) found in '~/atlib'
and places the generated files into 'results/'. Note that if the necessary
chgden* files are missing then they will also be generated and placed into
the location given by the '-a' option.

.. code:: bash

    eeasisss.py -i inputX -a ~/atlib/ -o results/

The above command will print the results to ``stdout`` i.e. to the terminal,
however the log can be saved to a file (e.g. 'ilogX') by using:

.. code:: bash

    eeasisss.py -i inputX -a ~/atlib/ -o results/ -l ilogX

The following does the same, but uses file redirection for the calculation log.
Note the subtle difference in command line usage.

.. code:: bash

    eeasisss.py -i inputX -a ~/atlib/ -o results/ > ilogX

On Linux there is also the possibility to use the 'tee' command pipe to print
the log to both the terminal screen and to a file.

.. code:: bash

    eeasisss.py -i inputX -a ~/atlib/ -o 'results/' | tee ilogX

'''
from __future__ import absolute_import, division, with_statement
from __future__ import print_function, unicode_literals

from ctypes import cdll, create_string_buffer
from ctypes.util import find_library
import os
import sys

from argparse import RawDescriptionHelpFormatter, ArgumentParser


_ext = '.dll' if str(sys.platform).startswith('win') else '.so'
_lib = os.path.join(os.path.dirname(__file__), 'lib')

os.environ['PATH'] = _lib + ';' + os.environ['PATH']
_library_path = (find_library('EEASiSSS') or
                 os.path.join(_lib, 'libEEASiSSS' + _ext))
libEEASiSSS = cdll.LoadLibrary(_library_path)

__date__ = '2015-03-18'
__updated__ = '2016-09-18'
__author__ = "Liam Deacon"
__contact__ = 'liam.m.deacon@gmail.com'

DEBUG = 0
TESTRUN = 0
PROFILE = 0


def eeasisss(input_file='inputX',
             output_dir=os.path.abspath('.'),
             chgden_dir=(os.path.expandvars('ATLIB')
                         if os.path.expandvars('ATLIB') != ''
                         else os.path.expanduser('~/atlib')),
             log_file=None, ini_file=None, verbose=False):
    '''
    Description
    -----------
    Wrapper function to call EEASiSSS Fortran library using ctypes

    Parameters
    ----------
    input_file : str
        Path to the input file containing the model geometry and associated
        calculation parameters. Default is 'inputX'.

    output_dir : str
        Path to place all the generated phase shifts and their plots into.
        If the directory does not exist then the function will attempt to
        create it. Default is './results/'.

    chgden_dir : str
        Path to the stored chgden* atomic charge density files calculated
        using Eric Shirley's :code:`hartfock()` subroutine. Default is
        either :envvar:`ATLIB` environment variable if it exists,
        or :code:`~/atlib/`. If neither directory exists then the necessary
        chgden files will be created for each element that is missing.

    log_file : str
        Path to logfile. If None then the calculation log will be printed
        to stdout. Default is to print to stdout.

    ini_file : str
        Path to INI file containing the user-specified default parameters for
        :code:`hartfock()` when automatically generating the ``chgden*`` files.
        If None (the default), then the files will be generated with the
        inbuilt default parameter values.

    Returns
    -------

    '''
    log_file = '' if log_file is None else log_file

    if not os.path.isfile(input_file):
        sys.stderr.write("eeasisss - error: input file '%s' not found\n"
                         % input_file)
        sys.stderr.flush()
        sys.exit(1)

    chgden_dir = os.path.expanduser(os.path.expandvars(chgden_dir))
    if not os.path.isdir(chgden_dir):
        if verbose:
            sys.stderr.write("eeasisss - warning: '%s' does not exist. "
                             "Attempting to read from './'" % chgden_dir)
            sys.stderr.flush()
        chgden_dir = os.path.abspath('.')

    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    libEEASiSSS.eeasisss_(create_string_buffer(str(input_file), size=255),
                          create_string_buffer(str(output_dir), size=255),
                          create_string_buffer(str(chgden_dir), size=255),
                          create_string_buffer(str(log_file), size=255))


def main(argv=None):
    '''Command line options to hf'''
    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    # display help if no arguments
    if len(argv) == 1:
        argv.append('--help')

    program_name = os.path.basename(sys.argv[0])
    program_version = "v%s" % __date__
    program_build_date = str(__updated__)
    program_version_message = '%%(prog)s %s (%s)' % (program_version,
                                                     program_build_date)
    program_shortdesc = __import__('__main__').__doc__.split("\n")[1]
    program_license = '''%s

      Please contact John Rundgren <jru@kth.se> for queries and comments
      relating to the native EEASiSSS or %s <%s> regarding this python package.

    usage:-
    ''' % (program_shortdesc, __author__, __contact__)

    try:
        # Setup argument parser
        parser = ArgumentParser(description=program_license,
                                formatter_class=RawDescriptionHelpFormatter)
        parser.add_argument('-i', '--input', dest='input', metavar='<inputA>',
                            help="path to the EEASiSSS input file which "
                            "contains the model geometry and calculation "
                            "parameters", required=True)
        parser.add_argument('-o', '--output', dest='output_dir',
                            metavar='<output_dir>', default='./result/',
                            help="storage path for the calculated phase "
                            "shifts and postscript plots. "
                            "[default: %(default)s")
        parser.add_argument('-l', '--log', dest='log', metavar='<log_file>',
                            default=None, help="Optionally redirects the "
                            "EEASiSSS calculation log to the file specified "
                            "by <log_file>. The default is to print to stdout "
                            " if the argument is None. [default: %(default)s]")
        parser.add_argument('-a', '--atom_dir', dest='chgden_dir',
                            metavar='<chgden_dir>',
                            default=(os.path.expandvars('ATLIB')
                                     if os.path.expandvars('ATLIB') != ''
                                     else os.path.expanduser('~/atlib')),
                            help="Specifies the directory filepath to where "
                            "the calculated atomic charge density files "
                            "reside. If no argument is given, the 'ATLIB' "
                            "environment variable is used as the output "
                            "directory, else '~/atlib' is used. Use '.' "
                            "to use the current working directory instead. "
                            "[default: %(default)s]")
        parser.add_argument("-v", "--verbose", dest="verbose",
                            action="store_true", default=False,
                            help="Verbosity flag for log output. Note that "
                            "setting this will also always generate "
                            "postscript plots of the calculated phase shifts. "
                            "[default: %(default)s]")
        parser.add_argument('-V', '--version', action='version',
                            version=program_version_message)

        # Process arguments
        args, unknown = parser.parse_known_args()

        if args.verbose > 0 and len(unknown) > 0:
            for arg in unknown:
                sys.stderr.write("eeasisss - warning: Unknown option '%s'\n"
                                 % arg)
            sys.stderr.flush()

    except KeyboardInterrupt:
        # handle keyboard interrupt #
        return 0
    except Exception as e:
        if DEBUG or TESTRUN:
            raise e
        indent = len(program_name) * " "
        sys.stderr.write(program_name + ": " + repr(e) + "\n")
        sys.stderr.write(indent + "  for help use --help")
        return 2

    # call underlying hartfock() subroutine to create charge density files
    eeasisss(args.input, args.output_dir, args.chgden_dir,
             args.log_file, args.verbose)

if __name__ == '__main__':
    # use this file as a script
    main()
