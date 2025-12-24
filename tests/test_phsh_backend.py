import pytest

from phaseshifts.phsh import BVHBackend, CLIError, ViperLeedBackend, get_backend


def test_get_backend_default():
    backend = get_backend(None)
    assert isinstance(backend, BVHBackend)


def test_get_backend_case_insensitive():
    backend = get_backend("ViPeRLeed")
    assert isinstance(backend, ViperLeedBackend)


def test_get_backend_invalid():
    with pytest.raises(CLIError):
        get_backend("unknown")


def test_viperleed_backend_requires_params():
    backend = get_backend("viperleed")
    with pytest.raises(CLIError):
        backend.autogen_from_input("bulk", "POSCAR")


def test_viperleed_backend_requires_slab():
    backend = get_backend("viperleed")
    with pytest.raises(CLIError):
        backend.autogen_from_input("bulk", None, backend_params="PARAMETERS")
