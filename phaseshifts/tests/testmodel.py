#!/usr/bin/env python
# encoding: utf-8
##############################################################################
# Author: Liam Deacon                                                        #
#                                                                            #
# Contact: liam.m.deacon@gmail.com                                           #
#                                                                            #
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
Unit tests for :py:mod:`phaseshifts.model` module
'''
from __future__ import print_function, unicode_literals
from __future__ import absolute_import, division, with_statement

import unittest

from unittest import TestSuite
from phaseshifts.tests.tests import TestCase


class TestAtomClass(TestCase):
    
    def test_instance(self):
        ''' Test Atom constructor '''
        from phaseshifts.model import Atom
        atom = Atom('H')
        self.assertIsInstance(atom, Atom)
        with self.assertRaises(TypeError):
            atom = Atom()  # cannot determine element for atom
    
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
        atom = Atom('H', coordinates=[0., 1., 2.])
        self.assertTrue(atom.bohr_coordinates == 
                        [0., 1. / Atom.bohr, 2. / Atom.bohr])
    
    def test_valence(self):
        from phaseshifts.model import Atom
        atom = Atom('H', valence=2.)
        self.assertTrue(atom.valence == 2.)
        atom.valence = -2
        self.assertTrue(atom.valence == -2)
        with self.assertRaises(ValueError):
            atom.valence = 'NaN'
        
    def test_tag(self):
        from phaseshifts.model import Atom
        atom = Atom('H')
        self.assertTrue(atom.tag == 'H')
    
    def test_radius(self):
        from phaseshifts.model import Atom
        atom = Atom('H', radius=(10. * Atom.bohr))
        self.assertTrue(atom.bohr_radius == 10.)
            
    def test_bohr_radius(self):
        from phaseshifts.model import Atom
        atom = Atom('H', radius=10.)
        self.assertTrue(atom.radius == 10.)
        atom.radius = 20.
        self.assertTrue(atom.radius == 20.)
        with self.assertRaises(TypeError):
            atom.radius = 'tiny'
    
    def test_occupancy(self):
        from phaseshifts.model import Atom
        atom = Atom('H', occupancy=1.)
        self.assertTrue(atom.occupancy == 1.)
        atom.occupancy = 0.
        self.assertTrue(atom.occupancy == 0.)
        with self.assertRaises(ValueError):
            atom.occupancy = 2.
        with self.assertRaises(ValueError):
            atom.occupancy = -0.5
    
    def test_element(self):
        from phaseshifts.model import Atom
        from phaseshifts.elements import ELEMENTS
        self.assertEqual(Atom('H'), ELEMENTS[1])
        self.assertEqual(Atom('1'), ELEMENTS[1])
        self.assertEqual(Atom('Hydrogen'), ELEMENTS[1])
        with self.assertRaises(KeyError):
            Atom('Krytonite')
        with self.assertRaises(KeyError):
            Atom(-1)
        with self.assertRaises(KeyError):
            Atom(9001)  # its over nine thousand!  
    
    def test_Z(self):
        from phaseshifts.model import Atom
        atom = Atom('H')
        self.assertTrue(atom.Z == 1)
        self.assertTrue(atom.Z == atom.element.protons)
        atom.element = Atom('Ar')
        self.assertFalse(atom.Z == 1)
        self.assertTrue(atom.Z == atom.element.protons)
        
    def test_uniqueness(self):
        from phaseshifts.model import Atom
        atoms = [Atom('H'), Atom('H', valence=1), Atom('H'),
                 Atom('C'), Atom('H', tag='?'), Atom('H', occupancy=0.5)]
        self.assertNotEqual(atoms, set(atoms))


class TestUnitcellClass(TestCase):
    
    def test_instance(self):
        from phaseshifts.model import Unitcell
        uc = Unitcell()
        self.assertIsInstance(uc, Unitcell)
        self.assertTrue(uc.a == 1.)
        uc = Unitcell(a=1., c=1.6)
        self.assertTrue(uc.a == 1. and uc.b == uc.a and uc.c == 1.6)

    def test_a(self):
        from phaseshifts.model import Unitcell
        uc = Unitcell(a=10)
        self.assertTrue(uc.a == 10.)
        uc.a = 20.
        self.assertTrue(uc.a == 20.)

    def test_b(self):
        from phaseshifts.model import Unitcell
        uc = Unitcell(a=10)
        self.assertTrue(uc.b == uc.a and uc.b == 10.)
        uc.b = 20.
        self.assertTrue(uc.b == 20.)
        
    def test_c(self):
        from phaseshifts.model import Unitcell
        uc = Unitcell(a=10, c=5.)
        self.assertTrue(uc.b == uc.a and uc.a == 10. and uc.c == 5.)
        uc.c = 20.
        self.assertTrue(uc.b == uc.a and uc.a == 10. and uc.c == 20.)
        
    def test_alpha(self):
        from phaseshifts.model import Unitcell
        uc = Unitcell()
        self.assertTrue(uc.alpha == 90.)
        uc.alpha = 45.
        self.assertTrue(uc.alpha == 45.)
        uc.alpha = -45.
        self.assertTrue(uc.alpha == 315.)
        
    def test_beta(self):
        from phaseshifts.model import Unitcell
        uc = Unitcell()
        self.assertTrue(uc.beta == 90.)
        uc.beta = 45.
        self.assertTrue(uc.beta == 45.)
        
    def test_gamma(self):
        from phaseshifts.model import Unitcell
        uc = Unitcell()
        self.assertTrue(uc.gamma == 90.)
        uc.gamma = 45.
        self.assertTrue(uc.gamma == 45.)
        
    def test_basis(self):
        from phaseshifts.model import Unitcell
        uc = Unitcell()
        self.assertTrue(uc.basis == [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
        new_basis = [[1., 1., 0.], [0., 1., 1.], [1., 0., 1.]]
        uc.basis = new_basis
        self.assertTrue(uc.basis == new_basis)
        with self.assertRaises(IndexError):
            uc.basis = [[0.], [1., 2.]]
        with self.assertRaises(TypeError):
            uc.basis = [[0.]]
        with self.assertRaises(ValueError):
            uc.basis = [[1., 0., 'zero'], [0., 1., 0.], [0., 0., 1.]]


class TestModelClass(TestCase):
    
    def test_instance(self):
        from phaseshifts.model import Model, Unitcell
        model = Model()
        self.assertIsInstance(model, Model)
        
        # check defaults
        self.assertTrue(model.unitcell == Unitcell)
        self.assertTrue(model.atoms == [])
    
    def test_add_atom(self):
        from phaseshifts.model import Atom, Model
        model = Model(atoms=[Atom('H')])
        model.add_atom('C', [0., 0., 0.])
        
    def test_atoms(self):
        from phaseshifts.model import Atom, Model
        model = Model(atoms=(Atom('C'), Atom('H'), Atom('O')))
        self.assertTrue(model.atoms == list(Atom('C'), Atom('H'), Atom('O')))
        with self.assertRaises(TypeError):
            Model(atoms=Atom('H'))
    
    def test_name(self):
        from phaseshifts.model import Model
        model = Model(name='The Standard Model')  # Pun intended
        self.assertTrue(model.name == 'The Standard Model')
        model = Model(atoms=['C', 'O', 'O'])
        self.assertTrue(model.name == model.formula)
        model.name = 'The Greenhouse Effect'
        self.assertTrue(model.name == 'The Greenhouse Effect')
        
    def test_formula(self):
        from phaseshifts.model import Model
        model = Model(atoms=['C', 'O', 'O'])
        self.assertTrue(model.formula == 'CO2')
    
    def test_unitcell(self):
        from phaseshifts.model import Model, Unitcell
        model = Model()
        self.assertEqual(model.unitcell, Unitcell())
        model.unitcell = Unitcell(a=1., b=2., c=3., 
                                  alpha=120., beta=60., gamma=90.)
        self.assertEqual(model.unitcell, 
                         Unitcell(a=1., b=2., c=3., 
                                  alpha=120., beta=60., gamma=90.))
        with self.assertRaises(TypeError):
            model.unitcell = None
      
    def test_uniqueness(self):
        pass
        
    def test_check_coordinates(self):
        from phaseshifts.model import Atom, Model, CoordinatesError
        CoZnO = [Atom('Zn', coordinates=[0., 0., 0.], occupancy=0.5), 
                 Atom('Co', coordinates=[0., 0., 0.], occupancy=0.5), 
                 Atom('O', coordinates=[1., 0., 0.], occupancy=1.0)]
        model = Model(atoms=CoZnO)
        self.assertIsInstance(model.check_coordinates(), dict)
        with self.assertRaises(CoordinatesError):
            CoZnO[1].occupancy = 0.999  
            # occupancy for site [0., 0., 0.] is now > 1. 
            model.check_coordinates()


class TestMTZModelClass(TestCase):
    
    @classmethod
    def setUpClass(cls):
        from phaseshifts.model import MTZModel
        super(TestMTZModelClass, cls).setUpClass()
        cls._model = MTZModel()
    
    def test_nh(self):
        self.assertTrue(self._model.nh == 10)
        self._model.nh = 12
        self.assertTrue(self._model.nh == 12)    
    
    def test_nform(self):
        from phaseshifts.model import MTZModel
        self.assertTrue(self._model.nform == MTZModel.rel)
        self._model.nform = MTZModel.cav
        self.assertTrue(self._model.nform == MTZModel.cav)
        self._model.nform = MTZModel.wil
        self.assertTrue(self._model.nform == MTZModel.wil)
        self._model.nform = 'rel'
        self.assertTrue(self._model.nform == MTZModel.rel)
        self._model.nform = 'wil'
        self.assertTrue(self._model.nform == MTZModel.wil)
        self._model.nform = 'cav'
        self.assertTrue(self._model.nform == MTZModel.cav)        
        with self.assertRaises(ValueError):
            self._model.nform = None
        
    def test_exchange(self):
        pass
        

class TestModelModule(TestSuite):

    def testName(self):
        pass


if __name__ == "__main__":
    unittest.main()
