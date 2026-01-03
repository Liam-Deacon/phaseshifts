from __future__ import absolute_import, division, print_function, unicode_literals

import os


class BackendError(Exception):
    """Error raised when backend selection or execution fails."""


class PhaseShiftBackend(object):
    """Base interface for phase shift backends."""

    name = None

    def __init__(self, **kwargs):
        _ = kwargs

    def autogen_from_input(self, bulk_file, slab_file, tmp_dir=None, **kwargs):
        raise NotImplementedError


class BVHBackend(PhaseShiftBackend):
    """Barbieri/Van Hove backend using the built-in Wrapper."""

    name = "bvh"

    def autogen_from_input(self, bulk_file, slab_file, tmp_dir=None, **kwargs):
        from phaseshifts.phsh import Wrapper

        return Wrapper.autogen_from_input(bulk_file, slab_file, tmp_dir=tmp_dir, **kwargs)


def _load_viperleed_modules():
    try:
        from viperleed.calc.files import parameters as viper_params
        from viperleed.calc.files import phaseshifts as viper_phaseshifts
        from viperleed.calc.files import poscar as viper_poscar
        from viperleed.calc.psgen import runPhaseshiftGen
    except ImportError:
        return None
    return viper_params, viper_phaseshifts, viper_poscar, runPhaseshiftGen


def _load_native_eeasisss():
    try:
        from phaseshifts.lib.EEASiSSS import EEASiSSS as native_eeasisss
    except ImportError:
        return None
    return native_eeasisss


def _run_viperleed(parameters_file, slab_file, workdir, output_file):
    modules = _load_viperleed_modules()
    if modules is None:
        raise BackendError("eeasisss backend requires 'phaseshifts[viperleed]' to be installed.")

    viper_params, viper_phaseshifts, viper_poscar, runPhaseshiftGen = modules

    rparams = viper_params.read(parameters_file)
    viper_params.interpret(rparams)
    slab = viper_poscar.read(slab_file)

    os.makedirs(workdir, exist_ok=True)

    output_path = os.path.join(workdir, output_file)
    cwd = os.getcwd()
    try:
        os.chdir(workdir)
        result = runPhaseshiftGen(slab, rparams)
        firstline = result[0]
        phaseshifts = result[1]
        viper_phaseshifts.writePHASESHIFTS(firstline, phaseshifts, file_path=output_path)
    finally:
        os.chdir(cwd)

    return [output_path]


def _run_native_eeasisss(parameters_file, workdir, output_file):
    native = _load_native_eeasisss()
    if native is None:
        raise BackendError(
            "Native EEASiSSS library not available; install it or use "
            "'phaseshifts[viperleed]' with --backend viperleed."
        )

    if not os.path.isfile(parameters_file):
        raise BackendError("eeasisss backend input file not found: {}".format(parameters_file))

    os.makedirs(workdir, exist_ok=True)

    output_path = os.path.join(workdir, output_file)
    cwd = os.getcwd()
    try:
        os.chdir(workdir)
        native.eeasisss(parameters_file)
    finally:
        os.chdir(cwd)

    if not os.path.exists(output_path):
        raise BackendError("EEASiSSS did not produce '{}' in '{}'.".format(output_file, workdir))
    return [output_path]


class EEASiSSSBackend(PhaseShiftBackend):
    """EEASiSSS backend using native library or ViPErLEED."""

    name = "eeasisss"
    _valid_modes = ("auto", "native", "viperleed")

    def __init__(self, mode="auto", **kwargs):
        _ = kwargs
        self.mode = mode or "auto"

    def autogen_from_input(self, bulk_file, slab_file, tmp_dir=None, **kwargs):
        _ = bulk_file
        parameters_file = kwargs.get("backend_params")
        workdir = kwargs.get("backend_workdir") or tmp_dir or os.getcwd()
        output_file = kwargs.get("output_file") or "PHASESHIFTS"
        verbose = kwargs.get("verbose")
        mode = str(self.mode or "auto").lower()

        if mode not in self._valid_modes:
            raise BackendError("eeasisss backend mode must be one of: {}.".format(", ".join(sorted(self._valid_modes))))

        if not parameters_file:
            raise BackendError("eeasisss backend requires --backend-params <PARAMETERS|inputX>.")

        if mode == "viperleed":
            if not slab_file:
                raise BackendError("eeasisss backend in viperleed mode requires --slab <POSCAR>.")
            if _load_viperleed_modules() is None:
                raise BackendError("eeasisss backend requires 'phaseshifts[viperleed]' to be installed.")
            return _run_viperleed(parameters_file, slab_file, workdir, output_file)

        if mode == "native":
            return _run_native_eeasisss(parameters_file, workdir, output_file)

        if slab_file:
            modules = _load_viperleed_modules()
            if modules is not None:
                return _run_viperleed(parameters_file, slab_file, workdir, output_file)
            if _load_native_eeasisss() is not None:
                if verbose:
                    print("eeasisss: ViPErLEED not available; falling back to native.")
                return _run_native_eeasisss(parameters_file, workdir, output_file)
            raise BackendError("eeasisss backend requires 'phaseshifts[viperleed]' to be installed.")

        if _load_native_eeasisss() is not None:
            return _run_native_eeasisss(parameters_file, workdir, output_file)

        if _load_viperleed_modules() is not None:
            raise BackendError("eeasisss backend requires --slab <POSCAR> when using ViPErLEED.")

        raise BackendError("eeasisss backend requires a native EEASiSSS library or " "'phaseshifts[viperleed]'.")


DEFAULT_BACKENDS = {
    BVHBackend.name: BVHBackend,
    EEASiSSSBackend.name: EEASiSSSBackend,
}

ALIASES = {
    "vht": "bvh",
    "vanhove": "bvh",
    "barbieri": "bvh",
    "bvt": "bvh",
    "easisss": "eeasisss",
    "easiss": "eeasisss",
    "viperleed": ("eeasisss", {"mode": "viperleed"}),
    "viper": ("eeasisss", {"mode": "viperleed"}),
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
    init_kwargs = {}
    alias = aliases.get(key)
    if alias is not None:
        if isinstance(alias, tuple):
            key = alias[0]
            if len(alias) > 1 and alias[1]:
                init_kwargs = alias[1]
        else:
            key = alias
    backend_cls = registry.get(key)
    if backend_cls is None:
        available = ", ".join(sorted(registry.keys()))
        raise BackendError("Unknown backend '{}'. Available: {}.".format(name, available))
    if init_kwargs:
        return backend_cls(**init_kwargs)
    return backend_cls()
