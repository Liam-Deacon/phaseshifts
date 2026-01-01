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

"""Provides an abstract factory class for phase shift calculations"""

import sys

try:
    from typing import List  # noqa: F401
except ImportError:
    pass

from .wrappers import BVHWrapper, EEASiSSSWrapper


class PhaseshiftFactory(object):
    """Class for backend selection."""

    backend = None  # type: object

    def __init__(self, backend, **kwargs):
        """Initialize the factory with a backend identifier."""
        package = str(backend).lower()
        self.__dict__.update(kwargs)
        self.phsh_files = []  # type: List[str]
        try:
            valid_bvh = ["vht", "bvh", "van hove", "barbieri", "bvt"]
            valid_eeasisss = ["eeasisss", "rundgren", "viperleed", "viper"]
            if package not in valid_bvh + valid_eeasisss:
                sys.stderr.write("Invalid package selected - using default (BVH)\n")
                sys.stderr.flush()
                self.backend = BVHWrapper()
            else:
                if package in valid_bvh:
                    self.backend = BVHWrapper()
                elif package in valid_eeasisss:
                    self.backend = EEASiSSSWrapper()
        except KeyError:
            sys.stderr.write("Invalid phaseshifts backend\n")
            sys.stderr.flush()
            sys.exit(-2)

    def _require_attrs(self, *names):
        missing = [name for name in names if not hasattr(self, name)]
        if missing:
            raise AttributeError(
                "Missing required phaseshift inputs: {}".format(", ".join(missing))
            )

    def createAtorbFiles(self):  # noqa: N802
        """Generate atomic orbital input files for the configured backend."""
        raise NotImplementedError("createAtorbFiles is not implemented.")

    def getPhaseShiftFiles(self):  # noqa: N802
        """Return a list of generated phase shift files."""
        self._require_attrs(
            "bulk_file",
            "slab_file",
            "tmp_dir",
            "lmax",
            "format",
            "store",
            "range",
        )
        return self.backend.autogen_from_input(
            self.bulk_file,
            self.slab_file,
            tmp_dir=self.tmp_dir,
            lmax=int(self.lmax),
            format=self.format,
            store=self.store,
            range=self.range,
        )
