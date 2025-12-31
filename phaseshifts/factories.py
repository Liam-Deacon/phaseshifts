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
    from typing import Any, List
except ImportError:
    pass

from .wrappers import BVHWrapper, EEASiSSSWrapper


class PhaseshiftFactory(object):
    """Class for backend selection"""

    backend = object
    phsh_files = []  # type: List[str]

    def __init__(self, backend, **kwargs):
        # type: (str, **Any) -> None
        self.__dict__.update(kwargs)
        self.phsh_files = []
        package = str(backend).lower()
        bvh_aliases = ("bvh", "vht", "van hove", "barbieri")
        eeasisss_aliases = ("eeasisss", "rundgren")
        if package in eeasisss_aliases:
            self.backend = EEASiSSSWrapper
        elif package in bvh_aliases:
            self.backend = BVHWrapper
        else:
            sys.stderr.write("Invalid package selected - using default (BVH)\n")
            sys.stderr.flush()
            self.backend = BVHWrapper

    def _require_attrs(self, *names):
        missing = [name for name in names if not hasattr(self, name)]
        if missing:
            raise AttributeError(
                "Missing required phaseshift inputs: {}".format(", ".join(missing))
            )

    def createAtorbFiles(self):
        pass

    def getPhaseShiftFiles(self):
        """Returns a list of generated phase shift files"""
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
