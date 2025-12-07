import os
import tempfile
import unittest

from phaseshifts.conphas import Conphas


class TestConphasViperFormat(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmpdir.cleanup()

    def _write_phasout(self):
        phasout_path = os.path.join(self.tmpdir.name, "sample.ph")
        with open(phasout_path, "w") as f:
            f.write("HEADER\n")
            f.write("0.0 5.0 2 1\n")
            # energy, l=0, l=1 for two energies
            f.write("0.0 0.1 0.2 5.0 0.3 0.4\n")
        return phasout_path

    def test_writes_viper_phaseshifts_format(self):
        phasout_path = self._write_phasout()
        output_path = os.path.join(self.tmpdir.name, "PHASESHIFTS")

        conphas = Conphas(
            input_files=[phasout_path],
            output_file=output_path,
            formatting="viper",
            lmax=1,
        )
        conphas.calculate()

        with open(output_path) as f:
            lines = [line.strip() for line in f.readlines() if line.strip()]

        header_tokens = lines[0].split()
        self.assertEqual(header_tokens[0], "1")
        self.assertEqual(header_tokens[1:5], ["-10.17", "-0.08", "-74.19", "19.18"])

        self.assertAlmostEqual(float(lines[1]), 0.0, places=4)
        block0 = [float(value) for value in lines[2].split()]
        self.assertEqual(block0, [0.1, 0.2])

        self.assertAlmostEqual(float(lines[3]), 5.0 / 27.21, places=4)
        block1 = [float(value) for value in lines[4].split()]
        self.assertEqual(block1, [0.3, 0.4])


if __name__ == "__main__":
    unittest.main()
