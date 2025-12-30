from __future__ import absolute_import, division, print_function, unicode_literals

import os


class BackendError(Exception):
    """Error raised when backend selection or execution fails."""


class PhaseShiftBackend(object):
    """Base interface for phase shift backends."""

    name = None

    def autogen_from_input(self, bulk_file, slab_file, tmp_dir=None, **kwargs):
        raise NotImplementedError


class BVHBackend(PhaseShiftBackend):
    """Barbieri/Van Hove backend using the built-in Wrapper."""

    name = "bvh"

    def autogen_from_input(self, bulk_file, slab_file, tmp_dir=None, **kwargs):
        from phaseshifts.phsh import Wrapper

        return Wrapper.autogen_from_input(
            bulk_file, slab_file, tmp_dir=tmp_dir, **kwargs
        )


class ViperLeedBackend(PhaseShiftBackend):
    """EEASiSSS backend using the ViPErLEED toolchain."""

    name = "eeasisss"

    def autogen_from_input(self, bulk_file, slab_file, tmp_dir=None, **kwargs):
        # ViPErLEED only needs the slab POSCAR input.
        _ = bulk_file
        parameters_file = kwargs.get("backend_params")
        workdir = kwargs.get("backend_workdir") or tmp_dir or os.getcwd()
        output_file = kwargs.get("output_file") or "PHASESHIFTS"

        if not parameters_file:
            raise BackendError(
                "eeasisss backend requires --backend-params <PARAMETERS>."
            )
        if not slab_file:
            raise BackendError("eeasisss backend requires --slab <POSCAR>.")

        try:
            from viperleed.calc.files import parameters as viper_params
            from viperleed.calc.files import phaseshifts as viper_phaseshifts
            from viperleed.calc.files import poscar as viper_poscar
            from viperleed.calc.psgen import runPhaseshiftGen
        except ImportError:
            raise BackendError(
                "eeasisss backend requires 'phaseshifts[viperleed]' to be installed."
            )

        rparams = viper_params.read(parameters_file)
        viper_params.interpret(rparams)
        slab = viper_poscar.read(slab_file)

        if not os.path.isdir(workdir):
            try:
                os.makedirs(workdir)
            except OSError:
                pass

        output_path = os.path.join(workdir, output_file)
        cwd = os.getcwd()
        try:
            os.chdir(workdir)
            result = runPhaseshiftGen(slab, rparams)
            firstline = result[0]
            phaseshifts = result[1]
            viper_phaseshifts.writePHASESHIFTS(
                firstline, phaseshifts, file_path=output_path
            )
        finally:
            os.chdir(cwd)

        return [output_path]


DEFAULT_BACKENDS = {
    BVHBackend.name: BVHBackend,
    ViperLeedBackend.name: ViperLeedBackend,
}

ALIASES = {
    "vht": "bvh",
    "vanhove": "bvh",
    "barbieri": "bvh",
    "bvt": "bvh",
    "easisss": "eeasisss",
    "easiss": "eeasisss",
    "viperleed": "eeasisss",
    "viper": "eeasisss",
}


def register_backend(name, backend_cls, registry=None):
    """Register a backend class in the provided registry."""
    registry = registry or DEFAULT_BACKENDS
    registry[str(name).lower()] = backend_cls


def get_backend(name, registry=None, aliases=None):
    """Return a backend instance by name."""
    registry = registry or DEFAULT_BACKENDS
    aliases = aliases or ALIASES
    key = str(name or BVHBackend.name).lower()
    key = aliases.get(key, key)
    backend_cls = registry.get(key)
    if backend_cls is None:
        available = ", ".join(sorted(registry.keys()))
        raise BackendError(
            "Unknown backend '{}'. Available: {}.".format(name, available)
        )
    return backend_cls()
