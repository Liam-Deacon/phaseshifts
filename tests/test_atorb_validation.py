import io
from configparser import ConfigParser

import pytest

from phaseshifts.atorb import Atorb, get_substr_positions, validate_atorb_file
from phaseshifts.validation.atorb import (
    AtorbElectron,
    AtorbInputModel,
    ValidationError,
    render_atorb_file,
)


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


def test_atorb_electron_validation():
    # Valid
    e = AtorbElectron(n=1, l=0, m=0, j=0.5, s=1, occ=1.0)
    e.ensure_valid()

    # Invalid n
    with pytest.raises(ValidationError, match="n must be positive"):
        AtorbElectron(n=0, l=0, m=0, j=0.5, s=1, occ=1.0).ensure_valid()

    # Invalid l
    with pytest.raises(ValidationError, match="l must satisfy"):
        AtorbElectron(n=1, l=1, m=0, j=0.5, s=1, occ=1.0).ensure_valid()

    # Invalid s
    with pytest.raises(ValidationError, match="spin multiplicity"):
        AtorbElectron(n=1, l=0, m=0, j=0.5, s=2, occ=1.0).ensure_valid()

    # Invalid occ
    with pytest.raises(ValidationError, match="occ must be non-negative"):
        AtorbElectron(n=1, l=0, m=0, j=0.5, s=1, occ=-1.0).ensure_valid()

    # Invalid j
    with pytest.raises(ValidationError, match="j must be a multiple of 0.5"):
        AtorbElectron(n=1, l=0, m=0, j=0.3, s=1, occ=1.0).ensure_valid()


def test_atorb_input_model_validation():
    # Setup valid orbitals
    orb = AtorbElectron(n=1, l=0, m=0, j=0.5, s=1, occ=2.0)

    # Valid model
    m = AtorbInputModel(
        z=1,
        nr=100,
        rel=1,
        method="HF",
        relic=0.0,
        nlevels=1,
        mixing_scf=0.5,
        eigen_tol=1e-5,
        ech=100,
        orbitals=[orb],
        output="out.i",
    )
    m.ensure_valid()

    # Invalid Z
    with pytest.raises(ValidationError, match="Atomic number"):
        m_bad = AtorbInputModel(
            z=0,
            nr=100,
            rel=1,
            method="HF",
            relic=0.0,
            nlevels=1,
            mixing_scf=0.5,
            eigen_tol=1e-5,
            ech=100,
            orbitals=[orb],
            output="out.i",
        )
        m_bad.ensure_valid()

    # Invalid NR
    with pytest.raises(ValidationError, match="Radial grid points"):
        m_bad = AtorbInputModel(
            z=1,
            nr=0,
            rel=1,
            method="HF",
            relic=0.0,
            nlevels=1,
            mixing_scf=0.5,
            eigen_tol=1e-5,
            ech=100,
            orbitals=[orb],
            output="out.i",
        )
        m_bad.ensure_valid()

    # Invalid rel
    with pytest.raises(ValidationError, match="rel must be 0"):
        m_bad = AtorbInputModel(
            z=1,
            nr=100,
            rel=2,
            method="HF",
            relic=0.0,
            nlevels=1,
            mixing_scf=0.5,
            eigen_tol=1e-5,
            ech=100,
            orbitals=[orb],
            output="out.i",
        )
        m_bad.ensure_valid()

    # Mismatch nlevels
    with pytest.raises(ValidationError, match="nlevels.*does not match"):
        m_bad = AtorbInputModel(
            z=1,
            nr=100,
            rel=1,
            method="HF",
            relic=0.0,
            nlevels=2,
            mixing_scf=0.5,
            eigen_tol=1e-5,
            ech=100,
            orbitals=[orb],
            output="out.i",
        )
        m_bad.ensure_valid()

    # Empty output
    with pytest.raises(ValidationError, match="output filename"):
        m_bad = AtorbInputModel(
            z=1,
            nr=100,
            rel=1,
            method="HF",
            relic=0.0,
            nlevels=1,
            mixing_scf=0.5,
            eigen_tol=1e-5,
            ech=100,
            orbitals=[orb],
            output=" ",
        )
        m_bad.ensure_valid()


def test_parser_failure_modes(tmp_path):
    # Helper to write and attempt validation
    def check_failure(content, match):
        f = tmp_path / "bad.txt"
        f.write_text(content)
        with pytest.raises(ValueError, match=match):
            validate_atorb_file(str(f))

    # Empty
    check_failure("", "empty")

    # Missing i
    check_failure("z\n", "start with an 'i'")

    # Bad Z/NR
    check_failure("i\nbad 100", "Unable to parse Z and NR")

    # Missing d
    check_failure("i\n1 100\nz", "Expected 'd' line")

    # Unexpected EOF before d
    check_failure("i\n1 100", "Unexpected end of file while expecting 'd'")

    # Bad rel
    check_failure("i\n1 100\nd\nbad", "Unable to parse relativistic flag")

    # Missing x
    check_failure("i\n1 100\nd\n1\nz", "Expected 'x' line")

    # Unexpected EOF before x
    check_failure("i\n1 100\nd\n1", "Unexpected end of file while expecting 'x'")

    # Missing a
    check_failure("i\n1 100\nd\n1\nx\n0.d0\nz", "Expected 'a' line")

    # Unexpected EOF before a
    check_failure(
        "i\n1 100\nd\n1\nx\n0.d0", "Unexpected end of file while expecting 'a'"
    )

    # Bad SCF
    check_failure("i\n1 100\nd\n1\nx\n0.d0\na\nbad", "Unable to parse SCF")

    # Bad orbital (too few tokens)
    check_failure(
        "i\n1 100\nd\n1\nx\n0.d0\na\n0 1 0.5 1e-5 100\n1 0",
        "Unable to parse orbital entry",
    )

    # Missing w
    check_failure(
        "i\n1 100\nd\n1\nx\n0.d0\na\n0 1 0.5 1e-5 100\n1 0 0 0.5 1 1.0\nz",
        "Expected 'w' line",
    )

    # Unexpected EOF before w
    check_failure(
        "i\n1 100\nd\n1\nx\n0.d0\na\n0 1 0.5 1e-5 100\n1 0 0 0.5 1 1.0",
        "Unexpected end of file while expecting 'w'",
    )

    # Missing output filename
    check_failure(
        "i\n1 100\nd\n1\nx\n0.d0\na\n0 1 0.5 1e-5 100\n1 0 0 0.5 1 1.0\nw",
        "Output filename missing",
    )


def test_render_with_header(tmp_path):
    orb = AtorbElectron(n=1, l=0, m=0, j=0.5, s=1, occ=2.0)
    m = AtorbInputModel(
        z=1,
        nr=100,
        rel=1,
        method="HF",
        relic=0.0,
        nlevels=1,
        mixing_scf=0.5,
        eigen_tol=1e-5,
        ech=100,
        orbitals=[orb],
        output="out.i",
        header="Custom Header",
    )
    dest = tmp_path / "rendered.txt"
    render_atorb_file(m, str(dest))

    content = dest.read_text()
    assert "Custom Header" in content
    assert (
        "atorb input file" not in content
    )  # Should override default if present, or at least appear


def test_validate_comment_spacing_at_start(tmp_path):
    # ! at start of line should be fine (no preceding char to check)
    f = tmp_path / "comment_start.txt"
    f.write_text(
        "i\n!comment\n1 100\nd\n1\nx\n0.d0\na\n0 1 0.5 1e-5 100\n1 0 0 0.5 1 1.0\nw\nout.i\nq"
    )

    # We just want to ensure it doesn't warn.
    import warnings

    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        validate_atorb_file(str(f))
        # Filter for UserWarning related to whitespace
        relevant = [w for w in record if "whitespace" in str(w.message)]
        assert len(relevant) == 0


def test_clean_atorb_lines_ignores_comments(tmp_path):
    f = tmp_path / "comments.txt"
    f.write_text(
        "C comment\ni\nC another comment\n1 100\nd\n1\nx\n0.d0\na\n0 1 0.5 1e-5 100\n1 0 0 0.5 1 1.0\nw\nout.i\nq"
    )
    model = validate_atorb_file(str(f))
    assert model.z == 1


def test_render_without_header(tmp_path):
    orb = AtorbElectron(n=1, l=0, m=0, j=0.5, s=1, occ=2.0)
    m = AtorbInputModel(
        z=1,
        nr=100,
        rel=1,
        method="HF",
        relic=0.0,
        nlevels=1,
        mixing_scf=0.5,
        eigen_tol=1e-5,
        ech=100,
        orbitals=[orb],
        output="out.i",
        header="",
    )
    dest = tmp_path / "rendered_default.txt"
    render_atorb_file(m, str(dest))

    content = dest.read_text()
    assert "atorb input file" in content


def test_get_substr_positions():
    assert get_substr_positions("a\nb\nc", "\n") == [1, 3]
    assert get_substr_positions("abc", "\n") == []
    assert get_substr_positions("", "\n") == []
    assert get_substr_positions("a,b,c", ",") == [1, 3]


def test_gen_conf_file_uses_rel_value(tmp_path):
    conf_path = tmp_path / "hf.conf"
    at = Atorb(rel=False, ngrid=123)
    at.gen_conf_file(str(conf_path))

    config = ConfigParser(allow_no_value=True)
    config.read(str(conf_path))

    assert config["DEFAULT"]["rel"] == "False"

    conf_true = tmp_path / "hf_true.conf"
    at_true = Atorb(rel=True, ngrid=123)
    at_true.gen_conf_file(str(conf_true))

    config_true = ConfigParser(allow_no_value=True)
    config_true.read(str(conf_true))

    assert config_true["DEFAULT"]["rel"] == "True"


def test_render_atorb_file_returns_handle_name():
    orb = AtorbElectron(n=1, l=0, m=0, j=0.5, s=1, occ=2.0)
    model = AtorbInputModel(
        header="test",
        z=1,
        nr=10,
        rel=0,
        method="HF",
        relic=0,
        nlevels=1,
        mixing_scf=0.05,
        eigen_tol=0.0005,
        ech=100,
        orbitals=[orb],
        output="at_Cl.i",
    )

    handle = io.StringIO()
    result = render_atorb_file(model, "unused.txt", file_handle=handle)

    assert result == ""
    assert "C test" in handle.getvalue()
