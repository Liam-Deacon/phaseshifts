from __future__ import absolute_import, division, print_function, unicode_literals

from . import IOAdapter
from .model import (
    AtomParameters,
    EnergyRange,
    InputParameters,
    Position,
    SuperstructureMatrix,
    UnitCell,
)

try:
    import yaml  # type: ignore
except ImportError:  # pragma: no cover - optional dependency
    yaml = None


class CleedpyYamlAdapter(IOAdapter):
    def _require_yaml(self):
        if yaml is None:
            raise ImportError("PyYAML is required for cleedpy YAML support")

    def read(self, path):
        self._require_yaml()
        with open(path, "r") as handle:
            data = yaml.safe_load(handle)
        return self.to_model(data)

    def write(self, model, path):
        self._require_yaml()
        data = self.from_model(model)
        with open(path, "w") as handle:
            yaml.safe_dump(data, handle, sort_keys=False)
        return path

    def to_model(self, data):
        if not isinstance(data, dict):
            raise ValueError("Expected mapping data for cleedpy YAML input")
        required_keys = [
            "unit_cell",
            "superstructure_matrix",
            "overlayers",
            "bulk_layers",
            "energy_range",
            "system_name",
            "minimum_radius",
        ]
        missing = [key for key in required_keys if key not in data]
        if missing:
            raise ValueError(
                "Missing required keys in YAML: {}".format(", ".join(missing))
            )
        unit_cell = UnitCell(
            data["unit_cell"]["a1"],
            data["unit_cell"]["a2"],
            data["unit_cell"]["a3"],
        )
        superstructure = SuperstructureMatrix(
            data["superstructure_matrix"]["m1"],
            data["superstructure_matrix"]["m2"],
        )
        overlayers = [self._parse_atom(item) for item in data["overlayers"]]
        bulk_layers = [self._parse_atom(item) for item in data["bulk_layers"]]
        energy_range = EnergyRange(
            data["energy_range"]["initial"],
            data["energy_range"]["final"],
            data["energy_range"]["step"],
        )

        return InputParameters(
            system_name=data["system_name"],
            unit_cell=unit_cell,
            superstructure_matrix=superstructure,
            overlayers=overlayers,
            bulk_layers=bulk_layers,
            minimum_radius=data["minimum_radius"],
            optical_potential=data.get("optical_potential", (8.0, 4.0)),
            energy_range=energy_range,
            polar_incidence_angle=data.get("polar_incidence_angle", 0.0),
            azimuthal_incidence_angle=data.get("azimuthal_incidence_angle", 0.0),
            epsilon=data.get("epsilon", 0.01),
            maximum_angular_momentum=data.get("maximum_angular_momentum", 8),
            sample_temperature=data.get("sample_temperature", 300.0),
        )

    def from_model(self, model):
        return model.to_dict()

    def _parse_atom(self, data):
        position = self._parse_position(data["position"])
        return AtomParameters(
            phase_file=data["phase_file"],
            position=position,
            vibrational_displacement=data["vibrational_displacement"],
        )

    def _parse_position(self, data):
        if isinstance(data, (list, tuple)):
            return Position(data[0], data[1], data[2])
        return Position(data["x"], data["y"], data["z"])
