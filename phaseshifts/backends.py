from __future__ import absolute_import, division, print_function, unicode_literals


class PhaseShiftBackend(object):
    """Base class for phase shift backends."""

    name = None

    def autogen_from_input(self, *args, **kwargs):
        raise NotImplementedError


class VHTBackend(PhaseShiftBackend):
    """Barbieri/Van Hove backend (default)."""

    name = "vht"

    def autogen_from_input(self, *args, **kwargs):
        from phaseshifts.wrappers import BVHWrapper

        return BVHWrapper.autogen_from_input(*args, **kwargs)


class EEASiSSSBackend(PhaseShiftBackend):
    """Rundgren EEASiSSS backend (optional)."""

    name = "eeasisss"

    def autogen_from_input(self, *args, **kwargs):
        try:
            from phaseshifts.wrappers import EEASiSSSWrapper
        except ImportError as err:
            raise RuntimeError(
                "EEASiSSS backend is not available. Install phaseshifts[eeasisss] "
                "and ensure the Fortran build is present."
            ) from err

        return EEASiSSSWrapper.autogen_from_input(*args, **kwargs)


DEFAULT_BACKENDS = {
    VHTBackend.name: VHTBackend,
    EEASiSSSBackend.name: EEASiSSSBackend,
}


def get_backend(name, registry=None):
    backend_name = (name or VHTBackend.name).lower()
    registry = registry or DEFAULT_BACKENDS
    try:
        backend_cls = registry[backend_name]
    except KeyError as err:
        available = ", ".join(sorted(registry.keys()))
        raise ValueError(
            "Unknown backend '{}'. Available: {}".format(backend_name, available)
        ) from err
    return backend_cls()
