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
""" Provides an abstract factory class for phase shift calculations """

from phaseshifts.wrappers import BVHWrapper, EEASiSSSWrapper, PhaseShiftWrapper


class PhaseshiftFactory(object):
    """
    Abstract factory class for backend selection of 
    phase shift package wrapper.
    """
    
    BACKENDS = {None: PhaseShiftWrapper,
                'bvh': BVHWrapper, 
                'eeasisss': EEASiSSSWrapper}
    """:py:obj:`dict` of available backends and their corresponding wrappers"""
    
    def __init__(self, backend='bvh', **kwargs):
        """ Barbieri/Van Hove package is currently the default """
        self.backend = backend or 'bvh'  
        self.__dict__.update(kwargs)
    
    @property
    def backend(self):
        """
        Returns the backend package wrapper class used for calculating the 
        phase shifts.
        """
        return self._backend
    
    @backend.setter
    def backend(self, package):
        """
        Sets the backend package wrapper to use for the phase shift 
        calculations. The current backends supported are: {} 
        
        Raises
        ------
        ValueError if ``package`` is not a known and supported backend.
        
        """.format(str("'" + "' '".join([k for k in self.BACKENDS]) + "'")) 
        
        if str(package).lower() not in self.BACKENDS:
            raise ValueError("Invalid package selected - please use one of: '"
                             "' '".join([key for key in self.BACKENDS] + "'"))
        else:
            self._backend = self.BACKENDS[package]
    
    def getPhaseShiftFiles(self):
        """Returns a list of generated phase shift files"""
        return self.backend.autogen_from_input(self.bulk_file, 
                                               self.slab_file, 
                                               tmp_dir=self.tmp_dir, 
                                               lmax=int(self.lmax),
                                               format=self.format, 
                                               store=self.store,
                                               range=self.range
                                               )     