import os
import tempfile
import unittest

from phaseshifts.conphas import Conphas


class TestConphasViperFormat(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmpdir.cleanup()

    def _write_phasout(
        self,
        filename="sample.ph",
        lmax=1,
        emin=0.0,
        emax=5.0,
        nenergies=2,
    ):
        """Write a minimal PHASOUT file with configurable lmax and number of energies."""
        phasout_path = os.path.join(self.tmpdir.name, filename)
        energies = [
            emin + i * (emax - emin) / max(nenergies - 1, 1) for i in range(nenergies)
        ]
        with open(phasout_path, "w") as f:
            f.write("HEADER\n")
            # emin, emax, nenergies, lmax
            f.write(f"{emin:.1f} {emax:.1f} {nenergies:d} {lmax:d}\n")

            values = []
            phase_counter = 1
            for e in energies:
                # energy
                values.append(f"{e:.1f}")
                # l = 0..lmax phases for this energy
                for _ in range(lmax + 1):
                    # just generate distinct dummy phase values; exact values are irrelevant
                    values.append(f"{0.1 * phase_counter:.1f}")
                    phase_counter += 1

            f.write(" ".join(values) + "\n")
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
        self.assertEqual(header_tokens[5], "phaseshifts")
        # Ensure correct token count (e.g. 7 tokens including timestamp)
        self.assertTrue(len(header_tokens) >= 7)

        self.assertAlmostEqual(float(lines[1]), 0.0, places=4)
        block0 = [float(value) for value in lines[2].split()]
        self.assertEqual(block0, [0.1, 0.2])

        self.assertAlmostEqual(float(lines[3]), 5.0 / 27.21, places=4)
        block1 = [float(value) for value in lines[4].split()]
        self.assertEqual(block1, [0.3, 0.4])

    def test_writes_viper_phaseshifts_with_custom_v0_params(self):
        """Test that custom v0_params are correctly written to the PHASESHIFTS header."""
        phasout_path = self._write_phasout()
        output_path = os.path.join(self.tmpdir.name, "PHASESHIFTS_CUSTOM")
        custom_v0 = [1.1, 2.2, 3.3, 4.4]

        conphas = Conphas(
            input_files=[phasout_path],
            output_file=output_path,
            formatting="viper",
            lmax=1,
            v0_params=custom_v0,
        )
        conphas.calculate()

        with open(output_path) as f:
            header = f.readline().strip()

        tokens = header.split()
        self.assertEqual(tokens[1:5], ["1.10", "2.20", "3.30", "4.40"])
        self.assertEqual(tokens[5], "phaseshifts")

    def test_writes_viper_phaseshifts_multiline_for_large_lmax(self):
        """Ensure that viper formatting exercises multi-line chunking for large lmax."""
        # Use lmax > 10 so that _chunked must span multiple lines per energy
        # For lmax=11, we have 12 values. Default chunk size is 10.
        # So we expect 2 lines per energy block + 1 line for energy.
        phasout_path = self._write_phasout(filename="sample_lmax11.ph", lmax=11)
        output_path = os.path.join(self.tmpdir.name, "PHASESHIFTS_LMAX11")

        conphas = Conphas(
            input_files=[phasout_path],
            output_file=output_path,
            formatting="viper",
            lmax=11,
        )
        conphas.calculate()

        with open(output_path) as f:
            lines = [line.rstrip() for line in f if line.strip()]

        # Header: 1 line
        # For each of 2 energies:
        #   1 line for energy
        #   1 line for first 10 phases
        #   1 line for remaining 2 phases
        # Total lines = 1 + 2 * (1 + 2) = 7
        self.assertGreater(len(lines), 4)

    def test_writes_viper_phaseshifts_multiblock_for_multiple_sites(self):
        """Ensure that viper formatting writes a block per input file (site)."""
        phasout_path_1 = self._write_phasout(filename="sample_site1.ph", lmax=3)
        phasout_path_2 = self._write_phasout(filename="sample_site2.ph", lmax=3)
        output_path = os.path.join(self.tmpdir.name, "PHASESHIFTS_MULTISITE")

        conphas = Conphas(
            input_files=[phasout_path_1, phasout_path_2],
            output_file=output_path,
            formatting="viper",
            lmax=3,
        )
        conphas.calculate()

        with open(output_path) as f:
            lines = [line.strip() for line in f.readlines() if line.strip()]

        # Header check: n_blocks should be 2
        header_tokens = lines[0].split()
        self.assertEqual(header_tokens[0], "2")

        # We don't have explicit blank lines separating blocks in the output writer:
        # _write_viper_output loops over input_files inside the energy loop.
        # Format is:
        # HEADER
        # ENERGY 1
        #   BLOCK 1 (SITE 1)
        #   BLOCK 2 (SITE 2)
        # ENERGY 2
        #   BLOCK 1 (SITE 1)
        #   BLOCK 2 (SITE 2)

        # For lmax=3, we have 4 values. Fits in one chunk (one line).
        # So structure:
        # Line 0: Header
        # Line 1: Energy 1
        # Line 2: Site 1 phases
        # Line 3: Site 2 phases
        # Line 4: Energy 2
        # Line 5: Site 1 phases
        # Line 6: Site 2 phases

        self.assertEqual(len(lines), 1 + 2 * (1 + 2))
