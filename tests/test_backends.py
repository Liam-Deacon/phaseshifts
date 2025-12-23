import pytest

from phaseshifts.backends import EEASiSSSBackend, VHTBackend, get_backend


def test_get_backend_default():
    backend = get_backend(None)
    assert isinstance(backend, VHTBackend)


def test_get_backend_unknown():
    with pytest.raises(ValueError):
        get_backend("unknown-backend")


def test_get_backend_eeasisss_constructs():
    backend = get_backend("eeasisss")
    assert isinstance(backend, EEASiSSSBackend)
