from __future__ import absolute_import, division, print_function, unicode_literals

try:
    from typing import Sequence, Tuple, Union
    VectorLike = Union[Tuple[float, float, float], Sequence[float]]
except ImportError:
    VectorLike = object


class Position(object):
    def __init__(self, x, y, z):
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)

    def to_dict(self):
        return {"x": self.x, "y": self.y, "z": self.z}


class UnitCell(object):
    def __init__(self, a1, a2, a3):
        self.a1 = tuple(a1)
        self.a2 = tuple(a2)
        self.a3 = tuple(a3)

    def to_dict(self):
        return {"a1": self.a1, "a2": self.a2, "a3": self.a3}


class SuperstructureMatrix(object):
    def __init__(self, m1, m2):
        self.m1 = tuple(m1)
        self.m2 = tuple(m2)

    def to_dict(self):
        return {"m1": self.m1, "m2": self.m2}


class AtomParameters(object):
    def __init__(self, phase_file, position, vibrational_displacement):
        self.phase_file = phase_file
        self.position = position
        self.vibrational_displacement = vibrational_displacement

    def to_dict(self):
        return {
            "phase_file": self.phase_file,
            "position": (
                self.position.to_dict()
                if hasattr(self.position, "to_dict")
                else self.position
            ),
            "vibrational_displacement": self.vibrational_displacement,
        }


class EnergyRange(object):
    def __init__(self, initial, final, step):
        self.initial = float(initial)
        self.final = float(final)
        self.step = float(step)

    def to_dict(self):
        return {"initial": self.initial, "final": self.final, "step": self.step}


class InputParameters(object):
    def __init__(
        self,
        system_name,
        unit_cell,
        superstructure_matrix,
        overlayers,
        bulk_layers,
        minimum_radius,
        energy_range,
        optical_potential=(8.0, 4.0),
        polar_incidence_angle=0.0,
        azimuthal_incidence_angle=0.0,
        epsilon=0.01,
        maximum_angular_momentum=8,
        sample_temperature=300.0,
    ):
        self.system_name = system_name
        self.unit_cell = unit_cell
        self.superstructure_matrix = superstructure_matrix
        self.overlayers = list(overlayers)
        self.bulk_layers = list(bulk_layers)
        self.minimum_radius = dict(minimum_radius)
        self.optical_potential = tuple(optical_potential)
        self.energy_range = energy_range
        self.polar_incidence_angle = float(polar_incidence_angle)
        self.azimuthal_incidence_angle = float(azimuthal_incidence_angle)
        self.epsilon = float(epsilon)
        self.maximum_angular_momentum = int(maximum_angular_momentum)
        self.sample_temperature = float(sample_temperature)

    def to_dict(self):
        return {
            "system_name": self.system_name,
            "unit_cell": self.unit_cell.to_dict(),
            "superstructure_matrix": self.superstructure_matrix.to_dict(),
            "overlayers": [item.to_dict() for item in self.overlayers],
            "bulk_layers": [item.to_dict() for item in self.bulk_layers],
            "minimum_radius": self.minimum_radius,
            "optical_potential": self.optical_potential,
            "energy_range": self.energy_range.to_dict(),
            "polar_incidence_angle": self.polar_incidence_angle,
            "azimuthal_incidence_angle": self.azimuthal_incidence_angle,
            "epsilon": self.epsilon,
            "maximum_angular_momentum": self.maximum_angular_momentum,
            "sample_temperature": self.sample_temperature,
        }
