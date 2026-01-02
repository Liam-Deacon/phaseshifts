from pathlib import Path

import pytest

import phaseshifts.wrappers as wrappers


@pytest.mark.parametrize(
    "out_format, expected_ext",
    [
        ("curve", ".cur"),
        ("cleed", ".phs"),
    ],
)
def test_remove_pi_jumps_generates_output(
    monkeypatch, tmp_path, out_format, expected_ext
):
    calls = []

    class DummyConphas(object):
        def __init__(self, input_files, output_file, formatting, lmax):
            calls.append(
                {
                    "input_files": input_files,
                    "output_file": output_file,
                    "formatting": formatting,
                    "lmax": lmax,
                }
            )
            self.output_file = output_file

        def calculate(self):
            Path(self.output_file).write_text("data")

    monkeypatch.setattr(wrappers, "Conphas", DummyConphas)

    phasout = tmp_path / "phasout.dat"
    phasout.write_text("data")

    result = wrappers.Wrapper._remove_pi_jumps(
        ["Ni"], [str(phasout)], {"Ni": 4}, out_format
    )

    expected = str(phasout.with_suffix(expected_ext))

    assert result == [expected]
    assert len(calls) == 1
    assert Path(result[0]).is_file()
    assert calls[0]["formatting"] == out_format
    assert calls[0]["lmax"] == 4
