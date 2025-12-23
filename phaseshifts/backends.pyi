"""
Type stubs for phaseshifts.backends.
"""

from typing import Any, Dict, Optional, Type

class PhaseShiftBackend(object):
    name: Optional[str]

    def autogen_from_input(self, *args: Any, **kwargs: Any) -> Any: ...

class VHTBackend(PhaseShiftBackend):
    name: str

    def autogen_from_input(self, *args: Any, **kwargs: Any) -> Any: ...

class EEASiSSSBackend(PhaseShiftBackend):
    name: str

    def autogen_from_input(self, *args: Any, **kwargs: Any) -> Any: ...

DEFAULT_BACKENDS: Dict[str, Type[PhaseShiftBackend]]

def get_backend(
    name: Optional[str], registry: Optional[Dict[str, Type[PhaseShiftBackend]]] = ...
) -> PhaseShiftBackend: ...
