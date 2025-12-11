import os

import pytest

from phaseshifts.leed import Converter


def test_cleedpy_yaml_to_inputs(tmp_path):
    pytest.importorskip("yaml")

    yaml_content = """
system_name: test_structure
unit_cell:
  a1: [1.0, 0.0, 0.0]
  a2: [0.0, 1.0, 0.0]
  a3: [0.0, 0.0, 5.0]
overlayers:
  - phase_file: "Cu_extra.dat"
    position: [0.5, 0.5, 2.5]
bulk_layers:
  - phase_file: "Ni_pot.inp"
    position: [0.0, 0.0, 0.0]
minimum_radius:
  Cu: 1.0
  Ni: 0.9
maximum_angular_momentum: 4
"""
    yaml_path = tmp_path / "input.yml"
    yaml_path.write_text(yaml_content)

    bulk_file, slab_file, metadata = Converter.cleedpy_to_inputs(
        str(yaml_path), tmp_dir=tmp_path
    )

    assert os.path.exists(bulk_file)
    assert os.path.exists(slab_file)
    assert metadata["maximum_angular_momentum"] == 4

    # use explicit loader (dependency injection)
    import yaml

    bulk_model, slab_model, parsed = Converter.import_cleedpy_input(
        str(yaml_path), yaml_loader=yaml.safe_load
    )
    assert len(bulk_model.atoms) == 1
    assert len(slab_model.atoms) == 2
    assert getattr(slab_model.atoms[0], "lmax") == 4
    assert parsed["system_name"] == "test_structure"
