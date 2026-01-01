#!/usr/bin/env python
# encoding: utf-8

##############################################################################
# Author: Liam Deacon                                                        #
#                                                                            #
# Contact: liam.deacon@diamond.ac.uk                                         #
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
"""
eeasisss.py - calculate EEASiSSS phase shifts.

eeasisss calculates phase shifts using John Rundgren's EEASiSSS() subroutine.

Examples
--------
.. code:: bash

    eeasisss.py -i inputX -a ~/atlib/ -l ilogA


"""
from __future__ import absolute_import, division, print_function, unicode_literals

import argparse
import os
import sys
from ctypes import cdll, create_string_buffer
from ctypes.util import find_library


_ext = ".dll" if str(sys.platform).startswith("win") else ".so"
_lib = os.path.join(os.path.dirname(__file__), "lib")

os.environ["PATH"] = _lib + os.pathsep + os.environ.get("PATH", "")
_library_path = find_library("EEASiSSS") or os.path.join(_lib, "libEEASiSSS" + _ext)
try:
    lib_eeasisss = cdll.LoadLibrary(_library_path)
except OSError as err:
    raise ImportError(
        "EEASiSSS library not found. Expected at: {} ({})".format(_library_path, err)
    ) from err


def eeasisss(input_file="inputX"):
    """Call the EEASiSSS Fortran library using ctypes."""
    if not os.path.isfile(input_file):
        raise FileNotFoundError("Input file '{}' not found".format(input_file))
    input_bytes = str(input_file).encode("utf-8")
    lib_eeasisss.hartfock_(create_string_buffer(input_bytes, 255))


def main(argv=None):
    """CLI entry point for eeasisss."""
    parser = argparse.ArgumentParser(
        description="Call the EEASiSSS Fortran library using ctypes.",
    )
    parser.add_argument(
        "-i",
        "--input",
        dest="input_file",
        default="inputX",
        help="Input file passed to the EEASiSSS library",
    )
    args = parser.parse_args(argv)
    if not os.path.isfile(args.input_file):
        sys.stderr.write("Error: Input file '{}' not found\n".format(args.input_file))
        return 1
    try:
        eeasisss(args.input_file)
    except Exception as err:  # pylint: disable=broad-except
        sys.stderr.write("Error calling EEASiSSS: {}\n".format(err))
        return 1
    return 0


if __name__ == "__main__":
    # use this file as a script
    main()
