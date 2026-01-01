import pytest


yaml = pytest.importorskip("yaml")

from phaseshifts.io import (  # noqa: E402
    AtomParameters,
    CleedpyYamlAdapter,
    EnergyRange,
    InputParameters,
    Position,
    SuperstructureMatrix,
    UnitCell,
)


def test_cleedpy_yaml_roundtrip(tmp_path):
    adapter = CleedpyYamlAdapter()
    unit_cell = UnitCell((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
    superstructure = SuperstructureMatrix((1.0, 0.0), (0.0, 1.0))
    energy = EnergyRange(20.0, 200.0, 5.0)
    atom = AtomParameters("Ni", Position(0.0, 0.0, 0.0), ["dr3", 0.01, 0.01, 0.01])

    model = InputParameters(
        system_name="Ni(111)",
        unit_cell=unit_cell,
        superstructure_matrix=superstructure,
        overlayers=[atom],
        bulk_layers=[atom],
        minimum_radius={"Ni": 1.2},
        energy_range=energy,
    )

    path = tmp_path / "input.yml"
    adapter.write(model, str(path))
    loaded = adapter.read(str(path))

    assert loaded.system_name == "Ni(111)"
    assert loaded.unit_cell.a1 == unit_cell.a1
