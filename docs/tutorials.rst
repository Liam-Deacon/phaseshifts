.. _tutorials:

=========
Tutorials
=========

Welcome to the phaseshifts tutorials! These interactive Jupyter notebooks will guide you through
calculating electron scattering phase shifts for Low-Energy Electron Diffraction (LEED) analysis.

.. note::

   All tutorials are available as interactive Jupyter notebooks. You can:

   - Run them locally after installing phaseshifts
   - Launch them on Binder for cloud-based execution (click the Binder badge)
   - Read them as static documentation below

Quick Start
-----------

If you're new to phase shift calculations, start with Tutorial 1 which covers the basics
using the ``phsh.py`` high-level interface.

For a deeper understanding of the underlying physics and the four-step Barbieri/Van Hove
methodology, work through Tutorial 2.

Available Tutorials
-------------------

.. toctree::
   :maxdepth: 2
   :caption: Tutorial Notebooks

   tutorials/01_ni111_quickstart
   tutorials/02_zno_wurtzite_phshift2007

Tutorial 1: Ni(111) Quick Start
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Difficulty:** Beginner

This tutorial demonstrates the end-to-end phase shift calculation workflow for a
simple Ni(111) surface using the high-level ``phsh.py`` wrapper.

**Topics covered:**

- Understanding LEED phase shifts
- Creating bulk (.bul) and slab (.inp) input files
- Running calculations with ``phaseshifts.phsh.Wrapper``
- Visualizing phase shift curves
- Key parameters: energy range, l_max, muffin-tin radii

**Launch on Binder:** |binder-ni111|

.. |binder-ni111| image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/Liam-Deacon/phaseshifts/HEAD?labpath=docs/tutorials/01_ni111_quickstart.ipynb


Tutorial 2: ZnO Wurtzite 4-Step Workflow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Difficulty:** Intermediate

This tutorial provides a detailed walkthrough of the complete Barbieri/Van Hove
phase shift calculation methodology using Wurtzite ZnO (0001) as the example system.

**Topics covered:**

- Wurtzite crystal structure visualization
- **Step 0 (phsh0):** Atomic orbital charge density calculation
- **Step 1 (phsh1):** Muffin-tin potential generation
- **Step 2 (phsh2):** Phase shift calculation (relativistic vs non-relativistic)
- **Step 3 (conphas):** Removing pi-discontinuities
- Understanding output formats (CLEED, curve, ViPErLEED)
- Integration with cleedpy

**Launch on Binder:** |binder-zno|

.. |binder-zno| image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/Liam-Deacon/phaseshifts/HEAD?labpath=docs/tutorials/02_zno_wurtzite_phshift2007.ipynb


Running Tutorials Locally
-------------------------

To run the tutorials on your local machine:

1. Install phaseshifts with tutorial dependencies:

   .. code-block:: bash

      pip install phaseshifts[docs]
      # or
      pip install phaseshifts jupyter matplotlib ase

2. Clone the repository and navigate to the tutorials:

   .. code-block:: bash

      git clone https://github.com/Liam-Deacon/phaseshifts.git
      cd phaseshifts/docs/tutorials

3. Launch Jupyter:

   .. code-block:: bash

      jupyter notebook

4. Open the desired tutorial notebook.

.. note::

   Some tutorials require the compiled Fortran library (``libphsh``). If you encounter
   import errors, ensure you have built the library:

   .. code-block:: bash

      make libphsh
      # or
      python setup.py build_ext --inplace


Using Binder
------------

Binder provides a free, cloud-based environment for running the tutorials without
any local installation. Simply click the "launch binder" badge on any tutorial.

**Advantages:**

- No installation required
- Works on any device with a web browser
- Pre-configured environment

**Limitations:**

- Sessions are temporary (save your work!)
- Limited computational resources
- May take 1-2 minutes to start


Contributing New Tutorials
--------------------------

We welcome community contributions! To add a new tutorial:

1. Create a Jupyter notebook in ``docs/tutorials/``
2. Follow the naming convention: ``NN_descriptive_name.ipynb``
3. Include extensive markdown explanations
4. Add the tutorial to this index page
5. Submit a pull request

See our :doc:`contributing` guide for more details.


Further Reading
---------------

- :ref:`Van_Hove_Phase_Shift_Package_Guide` - Detailed documentation of the original phshift2007 package
- :doc:`input_files` - Input file format reference
- :doc:`running` - Command-line usage guide
- `cleedpy project <https://github.com/empa-scientific-it/cleedpy>`_ - Python LEED analysis package
