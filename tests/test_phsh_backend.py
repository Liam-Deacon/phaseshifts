import pytest

from phaseshifts.backends import (
    BVHBackend,
    BackendError,
    ViperLeedBackend,
    get_backend,
)


def test_get_backend_default():
    backend = get_backend(None)
    assert isinstance(backend, BVHBackend)


def test_get_backend_case_insensitive():
    backend = get_backend("EEASiSSS")
    assert isinstance(backend, ViperLeedBackend)


def test_get_backend_invalid():
    with pytest.raises(BackendError):
        get_backend("unknown")


def test_eeasisss_backend_requires_params():
    backend = get_backend("eeasisss")
    with pytest.raises(BackendError):
        backend.autogen_from_input("bulk", "POSCAR")


def test_eeasisss_backend_requires_slab():
    backend = get_backend("eeasisss")
    with pytest.raises(BackendError):
        backend.autogen_from_input("bulk", None, backend_params="PARAMETERS")


def test_viperleed_alias_maps_to_eeasisss():
    backend = get_backend("viperleed")
    assert isinstance(backend, ViperLeedBackend)
