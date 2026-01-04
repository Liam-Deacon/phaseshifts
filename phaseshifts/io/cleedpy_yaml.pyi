"""
Type stubs for phaseshifts.io.cleedpy_yaml.
"""

from typing import Any

from . import IOAdapter
from .model import InputParameters

class CleedpyYamlAdapter(IOAdapter):
    def read(self, path: str) -> InputParameters: ...
    def write(self, model: InputParameters, path: str) -> str: ...
    def to_model(self, data: Any) -> InputParameters: ...
    def from_model(self, model: InputParameters) -> Any: ...
