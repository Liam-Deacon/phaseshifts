#!/usr/bin/env python
# encoding: utf-8
##############################################################################
# Author: Liam Deacon                                                        #
#                                                                            #
# Contact: liam.deacon@diamond.ac.uk                                         #
#                                                                            #
# Created on 12 May 2015                                                       #
#                                                                            #
# Copyright: Copyright (C) 2015 Liam Deacon                                #
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
Unit tests for :py:mod:`phaseshifts.model` module
'''
from __future__ import print_function, unicode_literals
from __future__ import absolute_import, division, with_statement

import unittest

from phaseshifts.tests.tests import TestCase


class TestAtomClass(TestCase):
    
    def test_instance(self):
        ''' Test Atom constructor '''
        from phaseshifts.model import Atom
        atom = Atom('H')
        self.assertTrue(isinstance(atom, Atom))
    
    def test_coordinates(self):
        from phaseshifts.model import Atom, CoordinatesError
        atom = Atom('H')
        
        # test coordinate assignment and retrieval for different types
        self.assertTrue(atom.coordinates, [0., 0., 0.])
        atom.coordinates = [1., 2., 3.]
        self.assertTrue(atom.coordinates, [1., 2., 3.])
        atom.coordinates = (3., 2., 1.)
        self.assertTrue(atom.coordinates, [3., 2., 1.])
        try:
            import numpy as np
            atom.coordinates = np.array([1., 2., 3.])
            self.assertTrue(atom.coordinates, [1., 2., 3.])
        except ImportError:
            pass
        atom.coordinates = (0, '1.', '2')
        self.assertTrue(atom.coordinates, [0., 1., 2.])
        
        # test raising of errors
        with self.assertRaises(CoordinatesError):
            atom.coordinates = [0.]
        with self.assertRaises(CoordinatesError):
            atom.coordinates = [0., 1., 2., 3.]
        with self.assertRaises(TypeError):
            atom.coordinates = '1 2 3'
        with self.assertRaises(ValueError):
            atom.coordinates = (3., 2., 'one')
            
    def test_bohr_coordinates(self):
        from phaseshifts.model import Atom
        atom = Atom('H')
    
    def test_valence(self):
        from phaseshifts.model import Atom
        atom = Atom('H')
    
    def test_tag(self):
        from phaseshifts.model import Atom
        atom = Atom('H')
        self.assertTrue(atom.tag == 'H')
    
    def test_radius(self):
        from phaseshifts.model import Atom
        atom = Atom('H', radius=10.)
        self.assertTrue(atom.radius == 10.)
        atom.radius = 20.
        self.assertTrue(atom.radius == 20.)
    
    def test_occupancy(self):
        from phaseshifts.model import Atom
        atom = Atom('H')
    
    def test_element(self):
        from phaseshifts.model import Atom
        from phaseshifts.elements import ELEMENTS
        self.assertEqual(Atom('H'), ELEMENTS[1])
        self.assertEqual(Atom('1'), ELEMENTS[1])
        self.assertEqual(Atom('Hydrogen'), ELEMENTS[1])
    
    def test_Z(self):
        from phaseshifts.model import Atom
        atom = Atom('H')
        self.assertTrue(atom.Z == 1)
        self.assertTrue(atom.Z == atom.element.protons)
        atom.element = Atom('Ar')
        self.assertFalse(atom.Z == 1)
        self.assertTrue(atom.Z == atom.element.protons)


class TestModelModule(TestCase):


    def testName(self):
        pass


if __name__ == "__main__":
    unittest.main()
