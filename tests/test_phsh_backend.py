import pytest

import phaseshifts.backends as backends
from phaseshifts.backends import (
    BVHBackend,
    BackendError,
    EEASiSSSBackend,
    ViperLeedBackend,
    get_backend,
)


def test_get_backend_default():
    backend = get_backend(None)
    assert isinstance(backend, BVHBackend)


def test_get_backend_case_insensitive():
    backend = get_backend("EEASiSSS")
    assert isinstance(backend, EEASiSSSBackend)


def test_get_backend_invalid():
    with pytest.raises(BackendError):
        get_backend("unknown")


def test_eeasisss_backend_requires_params():
    backend = get_backend("eeasisss")
    with pytest.raises(BackendError):
        backend.autogen_from_input("bulk", "POSCAR")


def test_eeasisss_backend_requires_native_or_viperleed(monkeypatch):
    backend = get_backend("eeasisss")
    monkeypatch.setattr(backends, "_load_native_eeasisss", lambda: None)
    monkeypatch.setattr(backends, "_load_viperleed_modules", lambda: None)
    with pytest.raises(BackendError):
        backend.autogen_from_input("bulk", None, backend_params="inputX")


def test_viperleed_backend():
    backend = get_backend("viperleed")
    assert isinstance(backend, ViperLeedBackend)
