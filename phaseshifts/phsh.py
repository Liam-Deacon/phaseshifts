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
"""
phsh.py - quickly generate phase shifts

phsh provides convenience functions to create phase shifts files
suitable for input into LEED-IV programs such as SATLEED and CLEED.

Examples
--------
.. code:: bash

   phsh.py -i *.inp -b *.bul -f CLEED -S phase_dir

   phsh.py --gui


"""
from __future__ import absolute_import, division, with_statement
from __future__ import print_function, unicode_literals

from subprocess import Popen
import datetime
import os
import platform
import sys

from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
import argparse


try:
    from . import __version__, __author_email__
    from .factories import PhaseshiftFactory
    from .utils import FileUtils, stringify
    from .leed import CSearch
except ValueError:
    from phaseshifts import __version__, __author_email__
    from phaseshifts.factories import PhaseshiftFactory
    from phaseshifts.utils import FileUtils, stringify
    from phaseshifts.leed import CSearch

try:
    from .gui import MainWindow
except ValueError:
    try:
        from phaseshifts.gui import MainWindow
    except ImportError:
        MainWindow = lambda x=None: sys.stderr.write("GUI not supported\n")
except ImportError:
    MainWindow = lambda x=None: sys.stderr.write("GUI not supported\n")


__all__ = []
__date__ = '2013-11-15'
__updated__ = '2016-09-18'

DEBUG = 0
TESTRUN = 0
PROFILE = 0
VERBOSE = 0

PHASESHIFT_FORMATS = ['cleed', 'curve', 'none']


def required_length(nmin, nmax):
    """custom action to check range"""
    class RequiredLength(argparse.Action):

        def __call__(self, parser, args, values, option_string=None):
            if not nmin <= len(values) <= nmax:
                msg = 'argument "{f}" requires between '
                '{nmin} and {nmax} arguments'.format(
                    f=self.dest, nmin=nmin, nmax=nmax)
                raise argparse.ArgumentTypeError(msg)
            setattr(args, self.dest, values)
    return RequiredLength


class CLIError(Exception):
    """Generic exception to raise and log different fatal errors."""

    def __init__(self, msg):
        super(CLIError).__init__(type(self))
        self.msg = "E: %s" % msg

    def __str__(self):
        return self.msg

    def __unicode__(self):
        return self.msg


def main(argv=None):
    """Command line options."""

    global VERBOSE

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
    program_shortdesc = 'phsh.py - quickly generate phase shifts'
    try:
        program_shortdesc = __import__('__main__').__doc__.split("\n")[1]
    except (ImportError, AttributeError):
        pass
    program_license = """{short_description}

      Created by Liam Deacon on {build_date}.
      Copyright 2013-{year} Liam Deacon. All rights reserved.

      Licensed under the MIT license (see LICENSE file for details)

      Please send your feedback, including bug notifications
      and fixes, to: {contact}

    usage:-
    """.format(short_description=program_shortdesc,
               build_date=str(__date__),
               year=datetime.datetime.now().year,
               contact=__author_email__)

    if '--gui' in argv and '-i' not in argv:
        argv.append('-i')
        argv.append('.')

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
                            default=10, type=int,
                            help="Maximum angular momentum "
                            "quantum number [default: %(default)s]")
        parser.add_argument('-f', '--format', dest='format',
                            metavar='<format>', default="cleed",
                            help="Use specific phase shift format "
                            "options are: {list}. [default: %(default)s]"
                            "".format(list=stringify(PHASESHIFT_FORMATS)))
        parser.add_argument('-r', '--range', dest='range', nargs='+',
                            action=required_length(2, 3), type=float,
                            metavar='<energy>', default=(20., 600., 5.),
                            help="Energy range in eV with the format: "
                            "'<start> <stop> [<step>]'. The <step> "
                            "value is optional. Valid for relativistic "
                            "calculations only. [default: %(default)s]")
        parser.add_argument('-g', '--generate-only', dest='generate',
                            action='store_true', default=False,
                            help="Exit after generating "
                            "phaseshifts; do not launch subprocess using "
                            "PHASESHIFTS_LEED environment variable. "
                            "[default: %(default)s]")
        parser.add_argument('-a', '--atorbs-only', dest='atorbs_only',
                            action='store_true', default=False,
                            help="Only generate atomic orbitals of elements "
                            "found in the input files using Eric Shirley's "
                            "hartfock routine, then exit. "
                            "[default: %(default)s]")
        parser.add_argument('-p', '--package', dest='package',
                            metavar='<package>', default='bvh',
                            help="Selects package to use for phase shift "
                            "calculations. Choices are: {}. "
                            "[default: %(default)s]"
                            "".format(stringify(PhaseshiftFactory.BACKENDS)))
        parser.add_argument('-S', '--store', dest='store', metavar='<subdir>',
                            default=False,
                            help="Keep intermediate files in subdir when done")
        parser.add_argument("-v", "--verbose", dest="verbose", action="count",
                            help="Set verbosity level. Note this will also "
                            "produce postscript graphs when using the EEASiSSS"
                            "backend. [default: %(default)s]")
        parser.add_argument('-V', '--version', action='version',
                            version=program_version_message)
        parser.add_argument("--gui", dest="gui", action="store_true",
                            help="Starts GUI frontend (experimental)")

        # Process arguments
        args, unknown = parser.parse_known_args()

        verbose = False
        verbose = args.verbose
        VERBOSE = verbose

        if verbose > 0 and len(unknown) > 0:
            for arg in unknown:
                sys.stderr.write("phsh - warning: Unknown option '%s'\n"
                                 % arg)
            sys.stderr.flush()

        if args.bulk is None:
            args.bulk = str(os.path.splitext(args.slab)[0] + '.bul')

        if args.store is False:
            args.store = '.'

        if args.lmax < 1 or args.lmax > 18:
            raise argparse.ArgumentError("lmax is not between 1 and 18")

        if len(args.range) < 3:  # add default step to list
            args.range = list(args.range).append(5)

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

    # create phase shifts (warning: black magic within - needs testing)
    if verbose:
        print("Phase shift auto-generation parameters")
        print("\tbulk input file: %s" % args.bulk)
        print("\tslab input file: %s" % args.slab)
        print("\tformat: %s" % args.format)
        print("\tlmax: %s" % args.lmax)
        print("\trange: %s eV" % [s for s in args.range])

    try:
        package = args.package
    except:
        package = 'vht'

    if args.atorbs_only is True:
        # only produce atomic orbital input files for Eric Shirley's hartfock
        sys.stderr.write("option '-a' or '--atorb-only' is not implemented\n")
        sys.exit(0)

    if args.gui:
        sys.exit(MainWindow.main(sys.argv))

    phaseshifts = PhaseshiftFactory(package,
                                    bulk_file=args.bulk,
                                    slab_file=args.slab,
                                    tmp_dir=args.tmpdir,
                                    lmax=int(args.lmax),
                                    format=args.format,
                                    store=args.store,
                                    range=args.range
                                    )

    phsh_files = phaseshifts.getPhaseShiftFiles()

    # chain loop commands to next program
    if 'PHASESHIFTS_LEED' in os.environ and not args.generate:

        # copy files into sub-directory
        csearch = CSearch(os.path.splitext(args.slab)[0])
        last_iteration = csearch.getIteration(-1)
        if last_iteration is not None:
            it = str(last_iteration).split('par:')[0].replace(' ', '')
            it = it.replace('#', '').rjust(3, '0')
            model = os.path.splitext(os.path.basename(args.slab))[0]
            parent = os.path.dirname(args.slab)
            name, ext = os.path.splitext(os.path.basename(args.slab))
            dest = os.path.join(parent, 'phsh_' + model, name + '_' + it + ext)
            FileUtils.copy_files(phsh_files, dest, verbose)

        # create subprocess
        leed_cmd = [os.environ['PHASESHIFTS_LEED']]
        # check if using native Windows Python with cygwin
        if (platform.system() == 'Windows' and
                leed_cmd.startwith('/cygdrive')):
            leed_cmd = '"%s"' % (leed_cmd.split('/')[2] + ':' +
                                 os.path.sep.join(leed_cmd.split('/')[3:]))

        for arg in argv:
            leed_cmd.append(arg)

        if verbose:
            print("phsh - starting subprocess: '%s'..." % " ".join(leed_cmd))

        # execute subprocess
        Popen(leed_cmd)


def gui_main(argv=None):
    """ Simple GUI wrapper function """
    argv = argv or sys.argv
    if '--gui' not in argv:
        argv.insert(0, '--gui')
    import subprocess
    subprocess.Popen(['phsh'] + argv, shell=True)


if __name__ == "__main__":
    if DEBUG:
        #  sys.argv.append("-h")
        sys.argv.append("-v")
        #  sys.argv.append("-r")

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
