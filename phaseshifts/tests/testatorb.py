'''
testatorb.py

Copyright (C) 2013 Liam Deacon
MIT License (see LICENSE file for details)

Test Atorb - checks whether the Atorb class functions as expected
'''
import os
import sys
import unittest

import tests.TestCase


class TestAtorbModule(tests.TestCase):

    def test_atorb(self):
        """ Test Atorb class instance """
        self.msg(self.shortDescription())
        from phaseshifts.atorb import Atorb
        at = Atorb()
        self.failUnless(isinstance(at, Atorb))

    def test_info(self):
        """ Test static method Atorb.get_quantum_info """
        sys.stderr.write('\n')
        self.msg(self.shortDescription())
        from phaseshifts.atorb import Atorb
        self.failUnlessEqual(Atorb.get_quantum_info('3d6'),
                             (3, 2, [1.5, 2.5], [2.4, 3.6]))

    def test_core_replace(self):
        """ Test static method Atorb.replace_core_config """
        sys.stderr.write('\n')
        self.msg(self.shortDescription())
        from phaseshifts.atorb import Atorb
        self.failUnless(Atorb.replace_core_config('[Ar]'),
                        '1s2 2s2 2p6 3s2 3p6')

    def test_calc(self):
        """ Test static method Atorb.gen_input """
        sys.stderr.write('\n')
        self.msg(self.shortDescription())
        from phaseshifts.atorb import Atorb
        filename = Atorb.gen_input('C', rel=True, output=os.tmpfile())
        self.failIf(not os.path.isfile(filename))

if __name__ == "__main__":
    sys.stderr.write('=' * 72 + '\n')
    sys.stderr.write('TESTING: %s\n' % os.path.basename(__file__))
    sys.stderr.write('=' * 72 + '\n')
    unittest.main()
