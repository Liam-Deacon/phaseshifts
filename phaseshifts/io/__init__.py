from __future__ import absolute_import, division, print_function, unicode_literals


class IOAdapter(object):
    """Base adapter interface for LEED input/output formats."""

    def read(self, path):
        raise NotImplementedError

    def write(self, model, path):
        raise NotImplementedError

    def to_model(self, data):
        raise NotImplementedError

    def from_model(self, model):
        raise NotImplementedError


from .model import (  # isort:skip
    AtomParameters,
    EnergyRange,
    InputParameters,
    Position,
    SuperstructureMatrix,
    UnitCell,
)

try:
    from .cleedpy_yaml import CleedpyYamlAdapter
except ImportError:  # pragma: no cover - optional dependency
    CleedpyYamlAdapter = None

__all__ = [
    "AtomParameters",
    "CleedpyYamlAdapter",
    "EnergyRange",
    "InputParameters",
    "IOAdapter",
    "Position",
    "SuperstructureMatrix",
    "UnitCell",
]
