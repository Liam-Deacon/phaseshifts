"""
Type stubs for phaseshifts.io.cleedpy_yaml.
"""

from typing import Any

from . import IOAdapter
from .model import AtomParameters, InputParameters, Position


class CleedpyYamlAdapter(IOAdapter):
    def _require_yaml(self) -> None:
        ...

    def read(self, path: str) -> InputParameters:
        ...

    def write(self, model: InputParameters, path: str) -> str:
        ...

    def to_model(self, data: Any) -> InputParameters:
        ...

    def from_model(self, model: InputParameters) -> Any:
        ...

    def _parse_atom(self, data: Any) -> AtomParameters:
        ...

    def _parse_position(self, data: Any) -> Position:
        ...
