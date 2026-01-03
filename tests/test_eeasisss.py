import importlib
import sys
from types import SimpleNamespace

import ctypes
import pytest


def test_eeasisss_missing_input_raises(monkeypatch, tmp_path):
    module_name = "phaseshifts.lib.EEASiSSS.EEASiSSS"
    if module_name in sys.modules:
        del sys.modules[module_name]

    dummy_lib = SimpleNamespace(hartfock_=lambda *args, **kwargs: 0)
    monkeypatch.setattr(ctypes.cdll, "LoadLibrary", lambda *args, **kwargs: dummy_lib)

    module = importlib.import_module(module_name)

    with pytest.raises(FileNotFoundError):
        module.eeasisss(str(tmp_path / "missing_input"))
