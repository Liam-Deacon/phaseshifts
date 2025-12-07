import os
from pathlib import Path

import pytest

from phaseshifts.atorb import Atorb, validate_atorb_file


def test_gen_input_cl_has_safe_comment_spacing(tmp_path):
    dest = tmp_path / "atorb_Cl.txt"
    Atorb.gen_input("Cl", filename=str(dest))

    lines = dest.read_text().splitlines()
    orbital_lines = [line for line in lines if "occupation" in line]
    assert orbital_lines, "Expected generated atorb file to contain orbital entries"

    for line in orbital_lines:
        assert "!" in line
        bang_index = line.index("!")
        assert line[
            bang_index - 1
        ].isspace(), "Inline comments must be separated from data by whitespace"

    validated = Atorb.validate_input_file(str(dest))
    assert validated.nlevels == len(validated.orbitals) == 7
    occ_sum = sum(o.occ for o in validated.orbitals if o.n == 3 and o.l == 1)
    assert pytest.approx(occ_sum, rel=1e-9) == 5.0


def test_validate_input_file_warns_on_missing_comment_spacing(tmp_path):
    bad_file = tmp_path / "atorb_bad.txt"
    bad_file.write_text(
        "\n".join(
            [
                "i",
                "17 1000",
                "d",
                "1",
                "x",
                "0.d0",
                "a",
                "0 1 0.5 0.0005 100",
                "3 1 1 -0.5 1 1.6667!bad comment",
                "w",
                "at_Cl.i",
                "q",
                "",
            ]
        )
    )

    with pytest.warns(UserWarning, match="whitespace"):
        validate_atorb_file(str(bad_file))
