.. _input_output_reference:

**************************
Input and Output Reference
**************************

Overview
--------

This page summarizes the input and output files used by the CLI, the
Barbieri/Van Hove backend, and the EEASiSSS backend. It is a quick
reference; the detailed formats live in the backend guides and examples.

Inputs
------

Command-line inputs (phsh.py)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``--bulk`` and ``--slab`` accept either:
  - CLEED input files (``*.bul``, ``*.inp``)
  - MTZ cluster input files (``*_bulk.i``, ``*_slab.i``, or ``cluster.i``)
- ``--input`` accepts a cleedpy-style JSON/YAML document. The schema lives in
  :py:data:`phaseshifts.leed._CLEEDPY_SCHEMA` and is validated when
  ``jsonschema`` is installed.
  Required top-level keys: ``unit_cell``, ``bulk_layers``, ``overlayers``.
- ``--backend`` selects the calculation backend (``bvh``, ``eeasisss``, ``viperleed``).
  Use ``--backend-params`` for backend-specific input files:
  - EEASiSSS native: ``inputX`` (model file) and ``inputA`` (atomic orbital list).
  - ViPErLEED mode: ``PARAMETERS`` file.

Input file inventory (Barbieri/Van Hove)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table:: Simplified file flow
   :header-rows: 1

   * - Stage
     - Inputs
     - Outputs
     - Notes
   * - Atorb (PhSh0)
     - ``atorb`` per element
     - ``at_*.i`` / ``atelem*.i``
     - Atomic orbital charge densities
   * - Muffin-tin (PhSh1 / cavpot)
     - ``cluster.i``, ``atomic.i``
     - ``mufftin.d``, ``*.bmtz``
     - ``atomic.i`` is concatenated ``at_*.i``
   * - Phase shifts (PhSh2*)
     - ``mufftin.d``
     - ``phasout``, ``dataph.d``, ``inpdat``
     - Algorithm depends on ``nform``
   * - Postprocess (PhSh3 / conphas)
     - ``phasout`` or ``*.ph``
     - ``*.phs`` / ``*.cur``
     - Removes pi-jumps and formats output

Example inputs live under ``phaseshifts/examples/input_files/`` and
``phaseshifts/examples/CLEED/``.

Outputs
-------

Primary outputs
~~~~~~~~~~~~~~~

- ``phsh.py`` writes one phase shift file per element. The extension and
  format are controlled by ``--format``:
  - ``cleed`` writes ``*.phs`` and copies to ``$CLEED_PHASE`` when set.
  - ``curve`` writes ``*.cur`` for plotting.
  - ``viperleed`` follows ViPErLEED conventions.
  - ``none`` keeps the backend raw output filenames.

Intermediate outputs
~~~~~~~~~~~~~~~~~~~~

- The default workflow uses a temporary directory. Use ``--tmpdir`` to
  override and ``--store`` to keep intermediate files.
- Intermediate files commonly include:
  - ``*_bulk.i``, ``*_slab.i`` (generated MTZ inputs)
  - ``*_mufftin.d`` (muffin-tin potential)
  - ``*_phasout.i``, ``*_dataph.d``, ``*_inpdat.txt``
  - ``*.ph`` per element before formatting

See Also
--------

- :doc:`scripts` for CLI usage
- :doc:`bvh_package` for the detailed Barbieri/Van Hove input format
- :doc:`eeasisss_package` for EEASiSSS input files
- :doc:`porting_phshift2007_notes` for step-by-step file lists
