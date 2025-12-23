"""
Type stubs for phaseshifts.io.
"""

from typing import Any, Optional, Type

from .model import (  # noqa: F401
    AtomParameters,
    EnergyRange,
    InputParameters,
    Position,
    SuperstructureMatrix,
    UnitCell,
)
from .cleedpy_yaml import CleedpyYamlAdapter as _CleedpyYamlAdapter

class IOAdapter(object):
    def read(self, path: str) -> Any: ...
    def write(self, model: Any, path: str) -> str: ...
    def to_model(self, data: Any) -> Any: ...
    def from_model(self, model: Any) -> Any: ...

CleedpyYamlAdapter: Optional[Type[_CleedpyYamlAdapter]]
