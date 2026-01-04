import unittest
import os
import shutil
import pytest

phsh2 = pytest.importorskip("phaseshifts.phsh2")


class TestPhsh2(unittest.TestCase):

    def setUp(self):
        self.test_dir = os.path.join(os.path.dirname(__file__), "temp_phsh2_test")
        os.makedirs(self.test_dir, exist_ok=True)
        os.chdir(self.test_dir)

    def tearDown(self):
        os.chdir(os.path.join(os.path.dirname(__file__), ".."))
        shutil.rmtree(self.test_dir)

    def test_cav(self):
        # Create mock input files for cav program
        with open("mufftin.d", "w") as f:
            f.write("1\n")
            f.write("Test Atom\n")
            f.write("1.0 1.0 1.0\n")
            f.write("2\n")
            f.write("0.1 0.2\n")
            f.write("0.3 0.4\n")

        # Run the cav program
        phsh2.run_cav()

        # Check if output files are created
        self.assertTrue(os.path.exists("zph.o"))
        self.assertTrue(os.path.exists("phasout"))
        self.assertTrue(os.path.exists("dataph.d"))

    def test_rel(self):
        # Create mock input files for rel program
        with open("mufftin.d", "w") as f:
            f.write("Test Atom\n")
            f.write("0.0 1.0 10.0 1 0.0\n")
            f.write("1 1.0 100 1.0 1.0 1.0 1.0 1.0 1.0\n")
            f.write("0.1 0.2 0.3 0.4 0.5\n")

        # Run the rel program
        phsh2.run_rel()

        # Check if output files are created
        self.assertTrue(os.path.exists("inpdat"))
        self.assertTrue(os.path.exists("phasout"))
        self.assertTrue(os.path.exists("dataph.d"))

    def test_wil(self):
        # Create mock input files for wil program
        with open("mufftin.d", "w") as f:
            f.write(
                "&NL16 CS=0.0, Z=1.0, E1=4.0, E2=24.0, NE=30, NL=9, NR=101, IX=1, RT=1.0, NEUI=1, NEUO=2, POTYP=2 &END\n"
            )
            for _ in range(200):
                f.write("0.1 0.2\n")

        # Run the wil program
        phsh2.run_wil()

        # Check if output files are created
        self.assertTrue(os.path.exists("zph.o"))
        self.assertTrue(os.path.exists("leedph.d"))
        self.assertTrue(os.path.exists("dataph.d"))


if __name__ == "__main__":
    unittest.main()
