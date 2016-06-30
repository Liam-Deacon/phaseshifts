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
Master testing module for phaseshifts package.
'''
from __future__ import print_function, unicode_literals
from __future__ import absolute_import, division, with_statement

import unittest
import sys
import glob


class TestCase(unittest.TestCase):

    def msg(self, text, newline=True):
        ch = '\n' if newline else ''
        sys.stderr.write(ch + str(text) + '...')
        sys.stderr.flush()


class TestSuite(unittest.TestSuite):
    
    test_dict = {'atorb': ['test.testatorb.TestAtorbModule'],
                 'conphas': ['test.testconphas.TestConphasModule'],
                 'elements': ['test.testelements.TestElementsModule'],
                 'utils': ['test.testutils.TestUtilsModule']
                 }
    
    def test_all(self, modules=[], test_imports=True):
        for module in modules:
            try:
                self.msg('Testing {} module'.format(module), newline=True)
                exec('import {}'.format(self.test_dict[module]))
            except any as e:
                raise e

    def __init__(self):
        import phaseshifts.tests as tests
        self.addTests([tests.testimports.TestImports('Test '
                                                     'phaseshifts imports'),
                       tests.testatorb.TestAtorbModule('Test atorb module'),
                       tests.testmodel.TestAtomClass('Test Atom class'),
                       tests.testmodel.TestModelModule('Test model module')
                       ])


if __name__ == "__main__":
    testSuite = unittest.TestSuite()
    module_strings = [test_unit[0:len(test_unit) - 3] 
                      for test_unit in glob.glob('test_*.py')]
    unittest.main()
