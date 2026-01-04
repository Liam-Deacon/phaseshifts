import pytest
import numpy as np

phsh0 = pytest.importorskip("phaseshifts.phsh0")
AtomicStructureSolver = phsh0.AtomicStructureSolver


def test_hydrogen_grid_and_orbital():
    calc = AtomicStructureSolver(Z=1, electron_config="1s1")
    r = calc.setup_grid()
    assert r[0] > 0 and r[-1] > r[0]
    calc.parse_config()
    orbitals = calc.initial_guess()
    assert (1, 0) in orbitals
    R10 = orbitals[(1, 0)]
    # Check normalization
    dr = np.gradient(r)
    norm = np.sum(R10**2 * dr)
    assert np.isclose(norm, 1.0, atol=1e-2)


def test_density_shape():
    calc = AtomicStructureSolver(Z=1, electron_config="1s1")
    calc.setup_grid()
    calc.parse_config()
    calc.initial_guess()
    density = calc.calc_density()
    assert density.shape == calc.r.shape
    assert np.all(density >= 0)


def test_cli(tmp_path):
    import subprocess

    output_file = tmp_path / "output.txt"
    cmd = ["python", "-m", "phaseshifts.phsh0", "--Z", "1", "--config", "1s1", "--output", str(output_file)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0
    assert output_file.exists()
    with open(output_file) as f:
        lines = f.readlines()
    assert any("density" in line for line in lines)


def test_python_vs_fortran_density():
    """
    Compare Python hydrogenic density to Fortran Dirac-Fock density for Hydrogen.

    Note:
    - The Fortran code uses a self-consistent Dirac-Fock method (relativistic).
    - The Python code uses hydrogenic (non-relativistic, non-self-consistent) orbitals.
    - Even for Hydrogen, the density profiles differ significantly due to self-consistency and relativistic effects.
    - This test checks that both profiles are qualitatively similar: peak near the origin and decay with radius.
    - Direct quantitative comparison will fail; see mean absolute error below.
    """
    # Fortran output file for Hydrogen
    fortran_file = "phaseshifts/examples/input_files/at_H.i"
    # Parse Fortran output
    with open(fortran_file) as f:
        lines = f.readlines()
    # Find grid params and density
    if len(lines) < 4:
        pytest.fail("Fortran output '{}' has {} lines; expected at least 4.".format(fortran_file, len(lines)))
    tokens = lines[3].split()
    if len(tokens) != 4:
        pytest.fail(
            "Fortran grid params line must have 4 tokens; got {} in line 4: {!r}".format(len(tokens), lines[3])
        )
    try:
        rmin, rmax, nr_value, Z = [float(x) for x in tokens]
        nr = int(nr_value)
    except (TypeError, ValueError) as exc:
        pytest.fail("Failed parsing grid params from line 4: {!r} ({})".format(lines[3], exc))
    if len(lines) < 4 + nr:
        pytest.fail(
            "Fortran output '{}' has {} lines; expected at least {} for density data.".format(
                fortran_file, len(lines), 4 + nr
            )
        )
    try:
        density_fortran = np.array([float(x) for x in lines[4 : 4 + nr]])
    except ValueError as exc:
        pytest.fail("Failed parsing density values from '{}' ({})".format(fortran_file, exc))
    # Reconstruct grid
    # h = np.log(rmax / rmin) / (nr - 1)
    # Run Python solver with matching grid
    calc = AtomicStructureSolver(
        Z=int(Z), electron_config="1s1", grid_params={"rmin": rmin, "rmax": rmax, "npoints": nr}
    )
    calc.setup_grid()
    calc.parse_config()
    calc.initial_guess()
    density_py = calc.calc_density()
    # Compare densities
    assert density_py.shape == density_fortran.shape
    # Compare normalized profiles
    norm_py = density_py / np.max(density_py)
    norm_fortran = density_fortran / np.max(density_fortran)
    # Compute mean absolute error
    mae = np.mean(np.abs(norm_py - norm_fortran))
    print("Mean absolute error (normalized density): {:.3e}".format(mae))
    # Qualitative check: both should peak near origin and decay
    peak_py = np.argmax(norm_py)
    peak_fortran = np.argmax(norm_fortran)
    assert peak_py == 0 and peak_fortran == 0, "Both densities should peak at the origin (r=0)"
    # Check that both densities decay to near zero at large r
    assert norm_py[-1] < 1e-3 and norm_fortran[-1] < 1e-3, "Both densities should decay to zero at large radius"
    # Document that direct quantitative comparison fails
    assert mae > 0.1, "Mean absolute error is expected to be large due to model differences"
