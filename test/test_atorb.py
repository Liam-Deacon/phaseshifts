"""
testatorb.py

Copyright (C) 2013 Liam Deacon
MIT License (see LICENSE file for details)

Test Atorb - checks whether the Atorb class functions as expected
"""

import unittest, sys, os


def print_test_message(text):
    sys.stderr.write(str(text) + "...")
    sys.stderr.flush()


class Test(unittest.TestCase):
    def test_atorb(self):
        """Test Atorb class instance"""
        print_test_message(self.shortDescription())
        from phaseshifts.atorb import Atorb

        at = Atorb()
        self.assertTrue(isinstance(at, Atorb))

    def test_info(self):
        """Test static method Atorb.get_quantum_info"""
        sys.stderr.write("\n")
        print_test_message(self.shortDescription())
        from phaseshifts.atorb import Atorb

        self.assertEqual(Atorb.get_quantum_info("3d6"), (3, 2, [1.5, 2.5], [2.4, 3.6]))

    def test_core_replace(self):
        """Test static method Atorb.replace_core_config"""
        sys.stderr.write("\n")
        print_test_message(self.shortDescription())
        from phaseshifts.atorb import Atorb

        self.assertTrue(Atorb.replace_core_config("[Ar]"), "1s2 2s2 2p6 3s2 3p6")

    def test_calc(self):
        """Test static method Atorb.gen_input"""
        sys.stderr.write("\n")
        print_test_message(self.shortDescription())
        from phaseshifts.atorb import Atorb

        import tempfile

        file = Atorb.gen_input("C", rel=True, output=tempfile.NamedTemporaryFile().name)
        self.assertTrue(os.path.isfile(file))


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    sys.stderr.write(
        "======================================================================\n"
    )
    sys.stderr.write("TESTING: %s\n" % os.path.basename(__file__))
    sys.stderr.write(
        "======================================================================\n"
    )
    unittest.main()
