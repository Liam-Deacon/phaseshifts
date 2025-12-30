from phaseshifts.io import (
    AtomParameters,
    EnergyRange,
    InputParameters,
    Position,
    SuperstructureMatrix,
    UnitCell,
)


def test_io_model_to_dict():
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

    data = model.to_dict()
    assert data["system_name"] == "Ni(111)"
    assert data["unit_cell"]["a1"] == unit_cell.a1
    assert data["superstructure_matrix"]["m1"] == superstructure.m1
