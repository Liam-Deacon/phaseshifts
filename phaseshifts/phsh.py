#! /usr/bin/env python3
# -*- coding: utf-8 -*-

# pylint: disable=consider-using-f-string

from __future__ import (
    absolute_import,
    division,
    print_function,
    unicode_literals,
    with_statement,
)

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

__doc__ = """
phsh.py - quickly generate phase shifts

Orchestrates the full phase shift calculation workflow for LEED (Low-Energy Electron Diffraction)
and related surface science applications, following the Barbieri/Van Hove methodology.

This module provides a command-line interface and Python API for generating phase shift files
suitable for input into LEED-IV programs such as SATLEED and CLEED. It automates the process of:
- Atomic charge density calculation (Dirac-Fock/Hartree-Fock)
- Muffin-tin potential generation
- Phase shift calculation (relativistic/non-relativistic)
- Removal of phase discontinuities (pi-jumps)
- Output formatting for LEED analysis

Scientific Context
------------------
Phase shifts are central to electron scattering theory and LEED surface structure analysis. The workflow
implemented here follows the Barbieri/Van Hove package, ensuring that all physical and computational steps
are performed in the correct sequence, with appropriate file formats and conventions.

References
----------
- `Barbieri, G., & Van Hove, M. A. (1979). "Phase shift calculation package for LEED." Surf. Sci. 90, 1-25.
     <https://doi.org/10.1016/0039-6028(79)90470-7>`_
- Moritz, W. (SATLEED code): http://www.icts.hkbu.edu.hk/surfstructinfo/SurfStrucInfo_files/leed/leedpack.html
- See also: https://phaseshifts.readthedocs.io/en/latest/phshift2007.html

Examples
--------
.. code:: bash

   phsh.py -i *.inp -b *.bul -f CLEED -S phase_dir

"""

# pylint: disable=wrong-import-position
import argparse
import datetime
import os
import platform
import subprocess
import sys
import tempfile
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from glob import glob

import phaseshifts
import phaseshifts.backends as _backend_registry
import phaseshifts.settings
from phaseshifts import atorb, model
from phaseshifts.conphas import Conphas
from phaseshifts.leed import CLEEDInputValidator, Converter, CSearch
from phaseshifts.utils import FileUtils

BMTZ_SUFFIX = ".bmtz"
MUFFTIN_SUFFIX = "_mufftin.d"

try:
    from phaseshifts.lib.libphsh import (  # type: ignore [import-untyped]
        phsh_cav,
        phsh_rel,
        phsh_wil,
    )
except ImportError:
    # TODO: Create pure-python implementations of the fortran libphsh.f code
    def function_not_implemented(*args, **kwargs):
        """Simple placeholder function, which raises a NotImplementedError."""
        raise NotImplementedError

    phsh_rel = phsh_wil = phsh_cav = function_not_implemented  # type: ignore

# pylint: enable=wrong-import-position


def required_length(nmin, nmax):
    """custom action to check range"""

    class RequiredLength(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            if not nmin <= len(values) <= nmax:
                msg = 'argument "{arg}" requires between {nmin} and {nmax} arguments'.format(
                    arg=self.dest, nmin=nmin, nmax=nmax
                )
                raise argparse.ArgumentTypeError(msg)
            setattr(args, self.dest, values)

    return RequiredLength


class Wrapper(object):
    """Class for automating the Barbieri/Van Hove phase shift calculation workflow.

    This class provides methods to generate phase shift files from atomic and cluster input data,
    handling all intermediate steps: atomic charge density, muffin-tin potential, phase shift calculation,
    phase-jump removal, and output formatting. It supports both Python and Fortran backends, and is
    compatible with SATLEED and CLEED file formats.

    Scientific Rationale
    -------------------
    Automating the workflow ensures reproducibility and correct sequencing of physical and computational steps,
    as required for LEED surface structure analysis.

    References
    ----------
    - `Barbieri, G., & Van Hove, M. A. (1979). "Phase shift calculation package for LEED." Surf. Sci. 90, 1-25.
     <https://doi.org/10.1016/0039-6028(79)90470-7>`_
    - Moritz, W. (SATLEED code): http://www.icts.hkbu.edu.hk/surfstructinfo/SurfStrucInfo_files/leed/leedpack.html
    - See also: https://phaseshifts.readthedocs.io/en/latest/phshift2007.html

    """

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    @staticmethod
    def autogen_from_input(bulk_file, slab_file, tmp_dir=None, model_name=None, **kwargs):
        """
           Generate phase shifts from slab/cluster and bulk input files, following the Barbieri/Van Hove workflow.

           This method automates the full sequence described in:
           - `Barbieri, G., & Van Hove, M. A. (1979). "Phase shift calculation package for LEED." Surf. Sci. 90, 1-25.
        <https://doi.org/10.1016/0039-6028(79)90470-7>`_
           - Moritz, W. (SATLEED code): http://www.icts.hkbu.edu.hk/surfstructinfo/SurfStrucInfo_files/leed/leedpack.html
           - See also: https://phaseshifts.readthedocs.io/en/latest/phshift2007.html

           Steps:
           1. Atomic charge density calculation (Dirac-Fock/Hartree-Fock)
           2. Muffin-tin potential generation
           3. Phase shift calculation (relativistic/non-relativistic)
           4. Removal of phase discontinuities (pi-jumps)
           5. Output formatting for LEED analysis

           Parameters
           ----------
           bulk_file : str
               Path to the bulk MTZ or CLEED input file.
           slab_file : str
               Path to the slab/cluster MTZ or CLEED input file.
           tmp_dir : str, optional
               Temporary directory for intermediate files.
           model_name : str, optional
               Name of the model.
           store : bool or str, optional
               Whether to keep generated files and where to store them.
           format : str, optional
               Output format for phase shift files ('cleed', 'curve', etc.).
           range : tuple(float, float, float), optional
               Energy range for phase shift calculation (start, stop, step in eV).

           Returns
           -------
           output_files : list of str
               List of generated phase shift output filenames.

           Scientific Context
           ------------------
           This workflow implements the full Barbieri/Van Hove methodology for LEED phase shift generation,
           ensuring all physical and computational steps are performed in the correct order. For further details,
           see the references above and the user guide at https://phaseshifts.readthedocs.io/en/latest/phshift2007.html

        """
        verbose = kwargs.get("verbose", False)
        dummycell = model.Unitcell(1, 2, [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        if model_name is None:
            model_name = "atomic"

        # check formatting
        if "format" in kwargs:
            out_format = kwargs["format"].lower()
        else:
            out_format = None

        # check lmax
        if "lmax" in kwargs:
            lmax = kwargs["lmax"]
        else:
            lmax = 10

        # check for intermediate storage directory, temp folder otherwise
        tmp_dir = str(tmp_dir)  # ensure string does not have escape chars
        if not os.path.isdir(tmp_dir):
            tmp_dir = tempfile.gettempdir()

        # Load bulk model and calculate MTZ
        bulk_mtz = model.MTZ_model(dummycell, atoms=[])
        if CLEEDInputValidator.is_cleed_file(bulk_file):
            bulk_mtz = Converter.import_CLEED(bulk_file, verbose=verbose)
            full_fname = glob(os.path.expanduser(os.path.expandvars(bulk_file)))[0]
            bulk_file = os.path.join(tmp_dir, os.path.splitext(os.path.basename(full_fname))[0] + "_bulk.i")
            bulk_mtz.gen_input(filename=bulk_file)
        else:
            bulk_mtz.load_from_file(bulk_file)

        # Load slab model and calculate MTZ
        slab_mtz = model.MTZ_model(dummycell, atoms=[])
        if CLEEDInputValidator.is_cleed_file(slab_file):
            slab_mtz = Converter.import_CLEED(slab_file)
            full_fname = glob(os.path.expanduser(os.path.expandvars(slab_file)))[0]
            slab_file = os.path.join(tmp_dir, os.path.splitext(os.path.basename(full_fname))[0] + "_slab.i")
            slab_mtz.gen_input(filename=slab_file)
        else:
            slab_mtz.load_from_file(slab_file)

        # generate atomic charge densities for each element in bulk model
        if not isinstance(bulk_mtz, model.MTZ_model):
            raise AttributeError("bulk_mtz is not an MTZ_model")

        # get unique elements in bulk and slab
        atomic_dict = {}
        bulk_elements = [getattr(atom.element, "symbol", str(atom.element)) for atom in bulk_mtz.atoms]
        slab_elements = [getattr(atom.element, "symbol", str(atom.element)) for atom in slab_mtz.atoms]
        for elem in set(bulk_elements + slab_elements):
            at_file = os.path.join(tmp_dir, "at_%s.i" % elem)
            if not os.path.isfile(at_file):
                print("\nCalculating atomic charge density for %s..." % elem)
                atomic_dict[elem] = atorb.Atorb.calculate_Q_density(element=elem, output_dir=tmp_dir)
            else:
                atomic_dict[elem] = at_file

        # prepare at files for appending into atomic file
        bulk_at_files = [
            atomic_dict[getattr(atom.element, "symbol", str(atom.element))] for atom in set(bulk_mtz.atoms)
        ]

        # create atomic.i input file from mtz model
        bulk_model_name = os.path.basename(os.path.splitext(bulk_file)[0])
        bulk_atomic_file = bulk_mtz.gen_atomic(
            input_files=bulk_at_files,
            output_file=os.path.join(tmp_dir, bulk_model_name + "_bulk.i"),
        )

        if verbose:
            print("\nModel")
            print("bulk atoms: %s" % [s for s in bulk_mtz.atoms])
            print("slab atoms: %s" % [s for s in slab_mtz.atoms])

        # calculate muffin-tin potential for bulk model
        print("\nCalculating bulk muffin-tin potential...")
        if verbose:
            print("\tcluster file: '%s'" % bulk_file)
            print("\tatomic file: '%s'" % bulk_atomic_file)
            print("\tslab calculation: '%s'" % str(False))
            print("\toutput file: '%s'" % os.path.join(tmp_dir, bulk_model_name + BMTZ_SUFFIX))
            print("\tmufftin file: '%s'" % os.path.join(tmp_dir, bulk_model_name + MUFFTIN_SUFFIX))

        # bulk_mtz_file = bulk_mtz.calculate_MTZ(
        bulk_mtz.calculate_MTZ(
            cluster_file=bulk_file,
            atomic_file=bulk_atomic_file,
            slab=False,
            output_file=os.path.join(tmp_dir, bulk_model_name + BMTZ_SUFFIX),
            mufftin_file=os.path.join(tmp_dir, bulk_model_name + MUFFTIN_SUFFIX),
        )
        print("Bulk MTZ = %f" % bulk_mtz.mtz)

        # prepare at files for appending into atomic file
        slab_at_files = [atomic_dict[atom.element.symbol] for atom in set(slab_mtz.atoms)]

        # create atomic.i input file from mtz model
        slab_model_name = os.path.basename(os.path.splitext(slab_file)[0])
        slab_atomic_file = slab_mtz.gen_atomic(
            input_files=slab_at_files,
            output_file=os.path.join(tmp_dir, slab_model_name + "_slab.i"),
        )

        # calculate muffin-tin potential for slab model
        mufftin_filepath = os.path.join(tmp_dir, slab_model_name + MUFFTIN_SUFFIX)
        print("\nCalculating slab muff-tin potential...")
        if verbose:
            print("\tcluster file: '%s'" % slab_file)
            print("\tatomic file: '%s'" % slab_atomic_file)
            print("\tslab calculation: %s" % str(True))
            print("\toutput file: '%s'" % os.path.join(tmp_dir, slab_model_name + BMTZ_SUFFIX))
            print("\tmufftin file: '%s'" % os.path.join(tmp_dir, mufftin_filepath))
            print("\tmtz value: %s" % str(bulk_mtz.mtz))

        # slab_mtz_file = slab_mtz.calculate_MTZ(
        slab_mtz.calculate_MTZ(
            cluster_file=slab_file,
            output_file=os.path.join(tmp_dir, slab_model_name + ".mtz"),
            atomic_file=slab_atomic_file,
            mufftin_file=mufftin_filepath,
            mtz_string=str(bulk_mtz.mtz),
            slab=True,
        )

        # create raw phase shift files
        print("\nGenerating phase shifts from '%s'..." % os.path.basename(mufftin_filepath))
        filepath = os.path.join(tmp_dir, slab_model_name)
        phasout_filepath = filepath + "_phasout.i"
        dataph_filepath = filepath + "_dataph.d"

        phaseshifts = [atom.tag for atom in set(slab_mtz.atoms)]
        phasout_files = [os.path.join(tmp_dir, atom.tag + ".ph") for atom in set(slab_mtz.atoms)]
        phsh_files = []

        # assign phase shift specific lmax values
        lmax_dict = {}
        for atom in set(slab_mtz.atoms + bulk_mtz.atoms):
            try:
                lmax_dict[atom.tag] = atom.lmax
            except AttributeError:
                print(atom.tag, "default lmax used:", lmax)
                lmax_dict[atom.tag] = lmax

        try:
            if slab_mtz.nform == 0 or str(slab_mtz.nform).lower().startswith("cav"):
                # calculate phase shifts
                phsh_cav(
                    mufftin_filepath,
                    phasout_filepath,
                    dataph_filepath,
                    filepath + "_zph.o",
                )

                # split phasout
                phasout_files = Conphas.split_phasout(filename=phasout_filepath, output_filenames=phasout_files)

                # eliminate pi-jumps
                for i, phaseshift in enumerate(phaseshifts):
                    filename = os.path.splitext(phasout_files[i])[0]
                    if out_format == "curve":
                        filename += ".cur"
                    else:
                        filename += ".phs"
                    phsh_files.append(filename)
                    print("\nRemoving pi/2 jumps in '%s':\n" % os.path.basename(filename))
                    phsh = Conphas(
                        input_files=[phasout_files[i]],
                        output_file=filename,
                        formatting=out_format,
                        lmax=lmax_dict[phaseshift],
                    )
                    phsh.calculate()

            if slab_mtz.nform == 1 or str(slab_mtz.nform).lower().startswith("wil"):
                # calculate phase shifts
                phsh_wil(
                    mufftin_filepath,
                    phasout_filepath,
                    dataph_filepath,
                    filepath + "_zph.o",
                )

                # split phasout
                phasout_files = Conphas.split_phasout(filename=phasout_filepath, output_filenames=phasout_files)

            if slab_mtz.nform == 2 or str(slab_mtz.nform).lower().startswith("rel"):
                # check energy range
                if "range" in kwargs:
                    try:
                        # get values
                        with open(mufftin_filepath, mode="r", encoding="utf-8") as f:
                            lines = [line for line in f]

                        ei, de, ef, lsm, vc = [
                            t(s)
                            for t, s in zip(
                                (float, float, float, int, float),
                                lines[1].replace("D", "E").split()[:5],
                            )
                        ]

                        # assign new values
                        (ei, ef, de) = [t(s) for t, s in zip((float, float, float), kwargs["range"])]

                        # edit energy range
                        lines[1] = str("%12.4f%12.4f%12.4f    %3i    %12.4f\n" % (ei, de, ef, lsm, vc)).replace(
                            "e", "D"
                        )
                        #                         lines[1] = ff.FortranRecordWriter(
                        #                                         '(3D12.4,4X,I3,4X,D12.4)'
                        #                                         ).write([ei, de, ef, lsm, vc]) + '\n'

                        with open(mufftin_filepath, mode="w", encoding="utf-8") as f:
                            f.write("".join([str(line) for line in lines]))

                    except Exception:
                        sys.stderr.write(
                            "Unable to change phase shift energy "
                            "range - using Barbieri/Van Hove "
                            "default of 20-300eV in 5eV steps\n"
                        )
                        sys.stderr.flush()

                # calculate phase shifts
                # print("Current time " + time.strftime("%X"))
                phsh_rel(
                    mufftin_filepath,
                    phasout_filepath,
                    dataph_filepath,
                    filepath + "_inpdat.txt",
                )
                # print("Current time " + time.strftime("%X"))

                # split phasout
                phasout_files = Conphas.split_phasout(filename=phasout_filepath, output_filenames=phasout_files)

                # eliminate pi-jumps
                for i, phaseshift in enumerate(phaseshifts):
                    filename = os.path.splitext(phasout_files[i])[0]
                    if out_format == "curve":
                        filename += ".cur"
                    else:
                        filename += ".phs"
                    phsh_files.append(filename)
                    print("\nRemoving pi/2 jumps in '%s':\n" % os.path.basename(filename))
                    phsh = Conphas(
                        input_files=[phasout_files[i]],
                        output_file=filename,
                        formatting=out_format,
                        lmax=lmax_dict[phaseshift],
                    )
                    phsh.calculate()

        except AttributeError:
            raise AttributeError("MTZ_model has no NFORM (0-2) specified!")

        # copy files to storage location
        if "store" in kwargs and out_format != "cleed":
            if kwargs["store"] != ".":
                dst = os.path.abspath(os.path.expanduser(os.path.expandvars(kwargs["store"])))
            else:
                dst = os.path.abspath(".")
            Wrapper._copy_files(phsh_files, dst, verbose=True)

        elif "CLEED_PHASE" in os.environ and out_format == "cleed":
            dst = os.path.abspath(os.path.expanduser(os.path.expandvars("$CLEED_PHASE")))
            Wrapper._copy_files(phsh_files, dst, verbose=True)

        else:
            Wrapper._copy_files(phsh_files, os.path.abspath("."), verbose=True)

        return phsh_files

    @staticmethod
    def _copy_files(files, dst, verbose=False):
        """Copy list of files into destination directory."""
        FileUtils.copy_files(files, dst, verbose=verbose)


class CLIError(Exception):
    """Generic exception to raise and log different fatal errors."""

    def __init__(self, msg):
        super(CLIError).__init__(type(self))
        self.msg = "E: %s" % msg

    def __str__(self):
        return self.msg

    def __unicode__(self):
        return self.msg


def _generate_atomic_orbitals(bulk_file, slab_file, tmp_dir=None, verbose=False):
    """Generate atomic charge density files for elements in bulk and slab inputs."""
    if not bulk_file or not slab_file:
        raise CLIError("--atorbs-only requires both bulk and slab inputs.")

    tmp_dir = tmp_dir or tempfile.gettempdir()
    if not os.path.isdir(tmp_dir):
        os.makedirs(tmp_dir)

    dummycell = model.Unitcell(1, 2, [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    bulk_mtz = model.MTZ_model(dummycell, atoms=[])
    if CLEEDInputValidator.is_cleed_file(bulk_file):
        bulk_mtz = Converter.import_CLEED(bulk_file, verbose=verbose)
    else:
        bulk_mtz.load_from_file(bulk_file)

    slab_mtz = model.MTZ_model(dummycell, atoms=[])
    if CLEEDInputValidator.is_cleed_file(slab_file):
        slab_mtz = Converter.import_CLEED(slab_file, verbose=verbose)
    else:
        slab_mtz.load_from_file(slab_file)

    elements = [getattr(atom.element, "symbol", str(atom.element)) for atom in bulk_mtz.atoms + slab_mtz.atoms]

    atomic_dict = {}
    for elem in sorted(set(elements)):
        at_file = os.path.join(tmp_dir, "at_%s.i" % elem)
        if not os.path.isfile(at_file):
            if verbose:
                sys.stdout.write("\nCalculating atomic charge density for %s...\n" % elem)
            atomic_dict[elem] = atorb.Atorb.calculate_Q_density(element=elem, output_dir=tmp_dir)
        else:
            atomic_dict[elem] = at_file

    return atomic_dict


def main(argv=None):
    """Run the command-line interface for phase shift generation.

    Parses user arguments, sets up the workflow, and executes the full Barbieri/Van Hove phase shift
    calculation sequence. Supports options for input files, output formatting, energy range, and
    intermediate file management.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments (default: sys.argv).

    Returns
    -------
    int
        Exit code (0 for success).

    Scientific Context
    ------------------
    This CLI enables reproducible, automated phase shift generation for LEED analysis, following
    best practices in surface science computation.

    References
    ----------
    - `Barbieri, G., & Van Hove, M. A. (1979). "Phase shift calculation package for LEED." Surf. Sci. 90, 1-25.
     <https://doi.org/10.1016/0039-6028(79)90470-7>`_
    - Moritz, W. (SATLEED code): http://www.icts.hkbu.edu.hk/surfstructinfo/SurfStrucInfo_files/leed/leedpack.html
    - See also: https://phaseshifts.readthedocs.io/en/latest/phshift2007.html

    """
    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    # display help if no arguments
    if len(argv) == 1:
        argv.append("--help")

    program_name = os.path.basename(sys.argv[0])
    program_version = "v%s" % phaseshifts.__version__
    program_version_message = "%%(prog)s %s" % program_version
    indent = len(program_name) * " "
    # Use this module's docstring or a static string for shortdesc
    program_shortdesc = (
        (__doc__ or "phsh.py - quickly generate phase shifts").split("\n")[1]
        if (__doc__ and len((__doc__ or "").split("\n")) > 1)
        else "phsh.py - quickly generate phase shifts"
    )
    program_license = """%s

      Created by Liam Deacon using LEEDPACK code (kindly permitted by Micheal Van Hove).
      Copyright 2013-%s Liam Deacon. All rights reserved.

      Licensed under the MIT license (see LICENSE file for details)

      Please send your feedback, including bug notifications
      and fixes, to: %s

    usage:-
    """ % (
        program_shortdesc,
        datetime.date.today().year,
        phaseshifts.__contact__,
    )

    try:
        # Setup argument parser
        parser = ArgumentParser(description=program_license, formatter_class=RawDescriptionHelpFormatter)
        parser.add_argument(
            "-b",
            "--bulk",
            dest="bulk",
            metavar="<bulk_file>",
            help="path to MTZ bulk or CLEED *.bul input file",
        )
        parser.add_argument(
            "-i",
            "--input",
            dest="input",
            metavar="<input.json|input.yml>",
            help="structured input (cleedpy-style) describing bulk and slab in JSON or YAML",
        )
        parser.add_argument(
            "-s",
            "--slab",
            dest="slab",
            metavar="<slab_file>",
            help="path to MTZ slab or CLEED *.inp input file (required unless --input is used)",
            required=False,
        )
        parser.add_argument(
            "-t",
            "--tmpdir",
            dest="tmpdir",
            metavar="<temp_dir>",
            help="temporary directory for intermediate " "file generation",
        )
        parser.add_argument(
            "-l",
            "--lmax",
            dest="lmax",
            metavar="<lmax>",
            default=10,
            type=int,
            help="Maximum angular momentum " "quantum number [default: %(default)s]",
        )
        parser.add_argument(
            "-f",
            "--format",
            dest="format",
            metavar="<format>",
            default="CLEED",
            help="Use specific phase shift format "
            "i.e. 'cleed', 'curve', 'viper', or 'viperleed' "
            "[default: %(default)s]",
        )
        parser.add_argument(
            "--backend",
            dest="backend",
            metavar="<backend>",
            default="bvh",
            help="Phase shift backend to use: 'bvh' (default), " "'eeasisss' (native or ViPErLEED; alias: viperleed).",
        )
        parser.add_argument(
            "--backend-params",
            dest="backend_params",
            metavar="<parameters>",
            help="Backend-specific parameters file. For viperleed mode, pass "
            "the ViPErLEED PARAMETERS file. For native eeasisss, pass an inputX "
            "file.",
        )
        parser.add_argument(
            "--backend-workdir",
            dest="backend_workdir",
            metavar="<dir>",
            help="Backend working directory (eeasisss/viperleed uses it for " "EEASiSSS files).",
        )
        parser.add_argument(
            "-r",
            "--range",
            dest="range",
            nargs="+",
            action=required_length(2, 3),
            type=float,
            metavar="<energy>",
            default=(20.0, 600.0, 5.0),
            help="Energy range in eV with the format: "
            "'<start> <stop> [<step>]'. The <step> "
            "value is optional. Valid for relativistic "
            "calculations only. [default: %(default)s]",
        )
        parser.add_argument(
            "-g",
            "--generate-only",
            dest="generate",
            action="store_true",
            default=False,
            help="Exit after generating "
            "phaseshifts; do not launch subprocess using "
            "PHASESHIFTS_LEED environment variable. "
            "[default: %(default)s]",
        )
        parser.add_argument(
            "-a",
            "--atorbs-only",
            dest="atorbs_only",
            action="store_true",
            default=False,
            help="Only generate atomic orbitals of elements "
            "found in the input files using Eric Shirley's "
            "hartfock routine, then exit. "
            "[default: %(default)s]",
        )
        parser.add_argument(
            "-p",
            "--package",
            dest="package",
            metavar="<package>",
            default="VHT",
            help="Selects package to use for phase shift "
            "calculations. Choices are 'VHT' (van Hove-Tong) "
            "or 'Rundgren' (EEASiSSS). [default: %(default)s]",
        )
        parser.add_argument(
            "-S",
            "--store",
            dest="store",
            metavar="<subdir>",
            default=False,
            help="Keep intermediate files in subdir when done",
        )
        parser.add_argument(
            "-v",
            "--verbose",
            dest="verbose",
            action="count",
            help="Set verbosity level. Note this will also "
            "produce postscript graphs when using the EEASiSSS "
            "backend. [default: %(default)s]",
        )
        parser.add_argument("-V", "--version", action="version", version=program_version_message)

        # Process arguments
        args, unknown = parser.parse_known_args()

        verbose = 0
        try:
            verbose = args.verbose
        except AttributeError:
            pass

        if verbose and len(unknown) > 0:
            for arg in unknown:
                sys.stderr.write("phsh - warning: Unknown option '%s'\n" % arg)
            sys.stderr.flush()

        default_backend = parser.get_default("backend")
        backend_name = str(getattr(args, "backend", None) or default_backend).lower()
        package_name = str(getattr(args, "package", "") or "").lower()
        if package_name and backend_name == str(default_backend).lower():
            if package_name in ("vht", "bvh", "van hove", "barbieri"):
                backend_name = "bvh"
            elif package_name in ("rundgren", "eeasisss"):
                backend_name = "eeasisss"

        # Structured input support: auto-generate bulk and slab inputs from geometry
        if getattr(args, "input", None):
            if args.bulk or args.slab:
                sys.stderr.write("phsh: --input provided; ignoring --bulk/--slab arguments\n")
                sys.stderr.flush()
            try:
                bulk_file, slab_file, metadata = Converter.cleedpy_to_inputs(args.input, tmp_dir=args.tmpdir)
            except ImportError as exc:
                raise CLIError(str(exc))
            args.bulk = bulk_file
            args.slab = slab_file
            if metadata.get("maximum_angular_momentum") and args.lmax == parser.get_default("lmax"):
                args.lmax = int(metadata["maximum_angular_momentum"])

        if not args.input and args.slab is None:
            raise CLIError("slab input is required unless using --input.")

        if args.store is False:
            args.store = "."

        if args.lmax < 1 or args.lmax > 18:
            raise ValueError("lmax is not between 1 and 18")

        if len(args.range) < 3:  # add default step to list
            args.range = list(args.range) + [5]

    except KeyboardInterrupt:
        # handle keyboard interrupt
        return 0
    except Exception as err:
        sys.stderr.write("{}: '{}'\n".format(program_name, err))
        sys.stderr.write("{} for help use --help".format(indent))
        return 2

    def _fatal(err):
        sys.stderr.write("{}: '{}'\n".format(program_name, err))
        sys.stderr.write("{} for help use --help".format(indent))
        return 2

    try:
        backend = _backend_registry.get_backend(backend_name)
    except _backend_registry.BackendError as err:
        return _fatal(err)

    backend_name = getattr(backend, "name", backend_name)

    if backend_name in ("eeasisss", "viperleed") and getattr(args, "input", None):
        return _fatal(CLIError("--input is not supported with the eeasisss backend."))

    if (
        args.bulk is None
        and not args.input
        and backend_name
        not in (
            "eeasisss",
            "viperleed",
        )
    ):
        args.bulk = str(os.path.splitext(args.slab)[0] + ".bul")

    if args.atorbs_only:
        try:
            atomic_dict = _generate_atomic_orbitals(args.bulk, args.slab, tmp_dir=args.tmpdir, verbose=verbose)
        except Exception as err:  # pylint: disable=broad-exception-caught
            return _fatal(err)
        if verbose:
            sys.stdout.write("\nGenerated atomic orbitals:\n")
            for at_file in sorted(atomic_dict.values()):
                sys.stdout.write("\t{}\n".format(at_file))
        return 0

    # create phase shifts (warning: black magic within - needs testing)
    if verbose:
        sys.stdout.write("Phase shift auto-generation parameters\n")
        sys.stdout.write("\tbulk input file: {}\n".format(args.bulk))
        sys.stdout.write("\tslab input file: {}\n".format(args.slab))
        sys.stdout.write("\tformat: {}\n".format(args.format))
        sys.stdout.write("\tbackend: {}\n".format(args.backend))
        sys.stdout.write("\tlmax: {}\n".format(args.lmax))
        sys.stdout.write("\trange: {} eV\n".format(args.range))

    backend_workdir = args.backend_workdir or args.tmpdir or args.store
    output_file = os.path.join(args.store, "PHASESHIFTS") if backend_name in ("eeasisss", "viperleed") else None

    try:
        phsh_files = backend.autogen_from_input(
            args.bulk,
            args.slab,
            tmp_dir=args.tmpdir,
            lmax=int(args.lmax),
            format=args.format,
            store=args.store,
            range=args.range,
            verbose=verbose,
            backend_params=args.backend_params,
            backend_workdir=backend_workdir,
            output_file=output_file,
        )
    except _backend_registry.BackendError as err:
        return _fatal(err)

    # chain loop commands to next program
    if "PHASESHIFTS_LEED" in os.environ and not args.generate:
        # copy files into subdirectory
        csearch = CSearch(os.path.splitext(args.slab)[0])
        last_iteration = csearch.getIteration(-1)
        if last_iteration is not None:
            it = str(last_iteration).split("par:", 1)[0].replace(" ", "").replace("#", "").rjust(3, "0")
            model = os.path.splitext(os.path.basename(args.slab))[0]
            parent = os.path.dirname(args.slab)
            name, ext = os.path.splitext(os.path.basename(args.slab))
            dest = os.path.join(parent, "phsh_" + model, name + "_" + it + ext)
            Wrapper._copy_files(phsh_files, dest, verbose)

        # create subprocess
        leed_cmd = [os.environ.get("PHASESHIFTS_LEED") or "cleed"]
        # check if using native Windows Python with cygwin
        if platform.system() == "Windows" and leed_cmd[0].startswith("/cygdrive"):
            leed_cmd[0] = '"%s"' % (leed_cmd[0].split("/")[2] + ":" + os.path.sep.join(leed_cmd[0].split("/")[3:]))

        leed_cmd.extend(argv)

        if verbose:
            print("phsh - starting subprocess: '{}'...".format(" ".join(leed_cmd)))

        # execute subprocess
        subprocess.check_call(leed_cmd)


if __name__ == "__main__":
    if phaseshifts.settings.DEBUG:
        sys.argv.append("-v")

    if phaseshifts.settings.TESTRUN:
        import doctest

        doctest.testmod()

    if phaseshifts.settings.PROFILE:
        import cProfile
        import pstats

        profile_filename = "wrapper_profile.txt"
        cProfile.run("main()", profile_filename)
        statsfile = open("profile_stats.txt", "wb")
        p = pstats.Stats(profile_filename, stream=statsfile)
        stats = p.strip_dirs().sort_stats("cumulative")
        stats.print_stats()
        statsfile.close()
        sys.exit(0)

    sys.exit(main())
