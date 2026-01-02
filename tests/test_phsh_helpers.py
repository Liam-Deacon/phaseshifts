import os
from pathlib import Path
from types import SimpleNamespace

import pytest

from phaseshifts import model
import phaseshifts.phsh as phsh


def test_wrapper_copy_files(tmp_path):
    src = tmp_path / "source.phs"
    src.write_text("data")
    dest = tmp_path / "dest"

    phsh.Wrapper._copy_files([str(src)], str(dest), verbose=True)

    assert (dest / "source.phs").read_text() == "data"


def test_generate_atomic_orbitals_requires_inputs():
    with pytest.raises(phsh.CLIError):
        phsh._generate_atomic_orbitals(None, "slab")


def test_generate_atomic_orbitals_collects_elements(monkeypatch, tmp_path):
    bulk_file = tmp_path / "bulk.bul"
    slab_file = tmp_path / "slab.inp"
    bulk_file.write_text("bulk")
    slab_file.write_text("slab")

    bulk_model = SimpleNamespace(atoms=[model.Atom("Ni")])
    slab_model = SimpleNamespace(atoms=[model.Atom("Ni"), model.Atom("Cu")])

    def fake_import(path, verbose=False):
        return bulk_model if path == str(bulk_file) else slab_model

    monkeypatch.setattr(phsh.Converter, "import_CLEED", fake_import)

    calls = []

    def fake_calc(element, output_dir, **kwargs):
        calls.append(element)
        output_path = Path(output_dir) / ("at_%s.i" % element)
        output_path.write_text("data")
        return str(output_path)

    monkeypatch.setattr(
        phsh.atorb.Atorb, "calculate_Q_density", staticmethod(fake_calc)
    )

    result = phsh._generate_atomic_orbitals(
        str(bulk_file), str(slab_file), tmp_dir=str(tmp_path)
    )

    assert set(result.keys()) == {"Ni", "Cu"}
    assert sorted(calls) == ["Cu", "Ni"]
    assert all(os.path.isfile(path) for path in result.values())


def test_main_package_maps_backend(monkeypatch, tmp_path):
    bulk_file = tmp_path / "bulk.bul"
    slab_file = tmp_path / "slab.inp"
    bulk_file.write_text("bulk")
    slab_file.write_text("slab")

    seen = {}

    class DummyBackend(object):
        name = "eeasisss"

    def fake_get_backend(name):
        # Backend registry maps package names (e.g., "rundgren" -> "eeasisss").
        seen["name"] = name
        return DummyBackend()

    monkeypatch.setattr(phsh._backend_registry, "get_backend", fake_get_backend)
    monkeypatch.setattr(
        phsh, "_generate_atomic_orbitals", lambda *args, **kwargs: {"Ni": "at_Ni.i"}
    )
    monkeypatch.setattr(
        phsh.sys,
        "argv",
        [
            "phsh.py",
            "--atorbs-only",
            "--package",
            "rundgren",
            "--bulk",
            str(bulk_file),
            "--slab",
            str(slab_file),
        ],
    )

    result = phsh.main()

    assert result == 0
    assert seen["name"] == "eeasisss"
