"""
Type stubs for phaseshifts.io.model.
"""

from typing import Any, Dict, List, Tuple

class Position(object):
    x: float
    y: float
    z: float

    def __init__(self, x: float, y: float, z: float) -> None: ...
    def to_dict(self) -> Dict[str, float]: ...

class UnitCell(object):
    a1: Tuple[float, float, float]
    a2: Tuple[float, float, float]
    a3: Tuple[float, float, float]

    def __init__(self, a1: Any, a2: Any, a3: Any) -> None: ...
    def to_dict(self) -> Dict[str, Any]: ...

class SuperstructureMatrix(object):
    m1: Tuple[float, float, float]
    m2: Tuple[float, float, float]

    def __init__(self, m1: Any, m2: Any) -> None: ...
    def to_dict(self) -> Dict[str, Any]: ...

class AtomParameters(object):
    phase_file: str
    position: Position
    vibrational_displacement: float

    def __init__(
        self, phase_file: str, position: Position, vibrational_displacement: float
    ) -> None: ...
    def to_dict(self) -> Dict[str, Any]: ...

class EnergyRange(object):
    initial: float
    final: float
    step: float

    def __init__(self, initial: float, final: float, step: float) -> None: ...
    def to_dict(self) -> Dict[str, float]: ...

class InputParameters(object):
    system_name: str
    unit_cell: UnitCell
    superstructure_matrix: SuperstructureMatrix
    overlayers: List[AtomParameters]
    bulk_layers: List[AtomParameters]
    minimum_radius: Dict[str, float]
    optical_potential: Tuple[float, float]
    energy_range: EnergyRange
    polar_incidence_angle: float
    azimuthal_incidence_angle: float
    epsilon: float
    maximum_angular_momentum: int
    sample_temperature: float

    def __init__(
        self,
        system_name: str,
        unit_cell: UnitCell,
        superstructure_matrix: SuperstructureMatrix,
        overlayers: List[AtomParameters],
        bulk_layers: List[AtomParameters],
        minimum_radius: Dict[str, float],
        energy_range: EnergyRange,
        optical_potential: Tuple[float, float] = ...,
        polar_incidence_angle: float = ...,
        azimuthal_incidence_angle: float = ...,
        epsilon: float = ...,
        maximum_angular_momentum: int = ...,
        sample_temperature: float = ...,
    ) -> None: ...
    def to_dict(self) -> Dict[str, Any]: ...
