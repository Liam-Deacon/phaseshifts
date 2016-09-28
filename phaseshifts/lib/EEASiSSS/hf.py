#!/usr/bin/env python
# encoding: utf-8
#
##############################################################################
# Author: Liam Deacon                                                        #
#                                                                            #
# Contact: liam.m.deacon@gmail.com                                           #
#                                                                            #
# Created on 26 Apr 2015                                                     #
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
**hf.py** - calculate atomic charge densities

Description
-----------
hf.py is a command line script which calculates atomic charge densities
for a given element or elements using a modified version of Eric Shirley's
``hartfock()`` subroutine.

Examples
--------
The following calculates the the atomic charge density file(s) specified in
the input file (e.g. 'inputA') and places them into the directory given by
the '-a' option, which is  '~/atlib/' in this example.

.. code:: bash

    hf.py -i inputA -a ~/atlib/

The above command will print the results to ``stdout`` i.e. to the terminal,
however the log can be saved to a file (e.g. 'ilogA') by using:

.. code:: bash

    hf.py -i inputA -a ~/atlib/ -l ilogA

The following does the same, but uses file redirection for the calculation log.
Note the subtle difference in command line usage.

.. code:: bash

    hf.py -i inputA -a ~/atlib/ > ilogA

On Linux there is also the possibility to use the 'tee' command pipe to print
the log to both the terminal screen and to a file.

.. code:: bash

    hf.py -i inputA -a ~/atlib/ | tee ilogA

'''
from __future__ import absolute_import, division, with_statement
from __future__ import print_function, unicode_literals

import os
import sys

from argparse import ArgumentParser, RawDescriptionHelpFormatter

from phaseshifts.lib.libhartfock import hartfock
from phaseshifts.utils import expand_filepath


__date__ = '2015-04-26'
__updated__ = '2016-09-18'
__author__ = 'Liam Deacon'
__contact__ = 'liam.m.deacon@gmail.com'

DEBUG = 0
TESTRUN = 0
PROFILE = 0


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

      Please contact Eric Shirley <eric.shirley@nist.gov> for queries
      and comments relating to the FORTRAN calculations or %s <%s> relating
      to the Python package.

    usage:-
    ''' % (program_shortdesc, __author__, __contact__)

    try:
        # Setup argument parser
        parser = ArgumentParser(description=program_license,
                                formatter_class=RawDescriptionHelpFormatter)
        parser.add_argument('-i', '--input', dest='input', metavar='<inputA>',
                            help="path to the atomic orbital input file which "
                            "is an appended list of Barbieri/Van Hove atorb "
                            "files (with some modifications).",
                            required=True)
        parser.add_argument('-l', '--log', dest='log', metavar='<log_file>',
                            default=None, help="Optionally redirects the "
                            "hartfock() calculation log to the file specified "
                            "by <log_file>. The default is to print to stdout "
                            " if the argument is None. [default: %(default)s]")
        parser.add_argument('-a', '--atom_dir', dest='chgden_dir',
                            metavar='<chgden_dir>',
                            default=(os.path.expandvars('ATLIB')
                                     if os.path.expandvars('ATLIB') != ''
                                     else os.path.expanduser('~/atlib')),
                            help="Specifies the output directory to place "
                            "the calculated atomic charge density files into. "
                            "If no argument is given, the 'ATLIB' environment "
                            "variable is used as the output directory, else "
                            "'~/atlib' is used. Use '.' to place the output "
                            "files into the current working directory. "
                            "[default: %(default)s]")
        parser.add_argument("-v", "--verbose", dest="verbose", action="count",
                            help="Set verbosity level. [default: %(default)s]")
        parser.add_argument('-V', '--version', action='version',
                            version=program_version_message)

        # Process arguments
        args, unknown = parser.parse_known_args()

        verbose = args.verbose

        if verbose > 0 and len(unknown) > 0:
            for arg in unknown:
                sys.stderr.write("hf - warning: Unknown option '%s'\n" % arg)
            sys.stderr.flush()

        args.log = '' if args.log is None else args.log

        if not os.path.isfile(args.input):
            sys.stderr.write("hf - error: input file '%s' not found\n"
                             % args.input)
            sys.stderr.flush()

        args.chgden_dir = expand_filepath(args.chgden_dir)
        if not os.path.isdir(args.chgden_dir):
            if verbose:
                sys.stderr.write("hf - warning: '%s' does not exist. "
                                 "Creating directory..." % args.chgden_dir)
                sys.stderr.flush()
            os.makedirs(args.chgden_dir)

    except KeyboardInterrupt:
        # handle keyboard interrupt #
        return 0
    except Exception, e:
        if DEBUG or TESTRUN:
            raise(e)
        indent = len(program_name) * " "
        sys.stderr.write(program_name + ": " + repr(e) + "\n")
        sys.stderr.write(indent + "  for help use --help")
        return 2

    # call underlying hartfock() subroutine to create charge density files
    hartfock(args.input, args.log, args.chgden_dir)


if __name__ == '__main__':
    main()
