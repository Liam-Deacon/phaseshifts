#!/usr/bin/env python
# encoding: utf-8
##############################################################################
# Author: Liam Deacon                                                        #
#                                                                            #
# Contact: liam.deacon@diamond.ac.uk                                         #
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
Tests for phaseshifts.utils module
'''
from __future__ import print_function, unicode_literals
from __future__ import absolute_import, division, with_statement

import unittest
import os
import sys

from phaseshifts.tests.tests import TestCase


class TestUtilsModule(TestCase):
    def test_fix_path(self):
        ''' Test fix_path() '''
        self.msg(self.shortDescription())
        from phaseshifts.utils import fix_path
        self.assertEqual(fix_path("C:\test\newfile.txt"), 
                         "C:\\test\\newfile.txt")
        self.assertNotEqual(fix_path("C:\test\newfile.txt"), 
                            "C:\test\newfile.txt")
        sys.stderr.write("SUCCESS")
        
    def test_expand_filepath(self):
        ''' Test expand_filepath() '''
        self.msg(self.shortDescription())
        from phaseshifts.utils import expand_filepath
        self.assertNotEqual(expand_filepath("~/"), "~/")
        self.assertEqual(expand_filepath("C:\\Users\\Liam"), "C:\\Users\\Liam")
        self.assertEqual(expand_filepath("$OS"), os.environ['OS'])
        sys.stderr.write("SUCCESS")

    def test_stringify(self):
        ''' Test stringify() '''
        self.msg(self.shortDescription())
        from phaseshifts.utils import stringify
        self.assertEqual(stringify(None), 'None')
        self.assertEqual(stringify([1, 'two', None]), 
                         u"'1', 'two', or 'None'")
        # note dictionaries seems to order keys in alphabetical-numeric 
        # order for list comprehension
        self.assertEqual(stringify({3: 'three', 'four': 4, False: None}),
                         u"'four', 'False', or '3'")
        self.assertRaises(TypeError, stringify)
        sys.stderr.write("SUCCESS")

if __name__ == "__main__":
    unittest.main()
