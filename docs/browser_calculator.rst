Browser Calculator
==================

Interactive Phase Shift Calculator
-----------------------------------

This browser-based calculator uses WebAssembly to run the Fortran
phase shift calculation code directly in your browser — no installation required.

.. raw:: html

   <div style="margin-bottom: 1em;">
     <a href="../calculator/" target="_blank" rel="noopener noreferrer"
        style="display: inline-block; padding: 10px 20px; background: #2980b9; color: white;
               text-decoration: none; border-radius: 4px; font-weight: bold;">
       ↗ Open Calculator in New Tab
     </a>
     <span style="color: #666; margin-left: 10px; font-size: 0.9em;">
       (Recommended for mobile users)
     </span>
   </div>

   <iframe
       src="../calculator/"
       style="width: 100%; height: 800px; border: 1px solid #ccc; border-radius: 4px;"
       title="Phase Shifts Calculator"
       loading="lazy">
       <p>Your browser does not support iframes.
          <a href="../calculator/">Open the calculator directly</a>.</p>
   </iframe>

.. note::

   If the calculator doesn't load in the embedded view above, please use the
   "Open Calculator in New Tab" button or access it directly at
   `calculator/index.html <../calculator/>`_.


Features
--------

The browser calculator supports three phase shift calculation methods:

**Relativistic (phsh_rel)**
   Full Dirac equation treatment with spin-orbit coupling. Best for heavy
   elements (Z > 30) where relativistic effects are significant.

**Cavity LEED (phsh_cav)**
   Traditional cavity method using Loucks grid. Suitable for standard
   LEED-IV analysis.

**Williams Method (phsh_wil)**
   A.R. Williams' phase shift calculation approach. Alternative method
   for comparison.


Supported Elements
------------------

All elements from Hydrogen (Z=1) to Uranium (Z=92) are supported, with
presets available for common surface science elements:

- **Metals**: Cu, Ag, Au, Ni, Fe, Pt, Pd, Al
- **Semiconductors**: Si, Ge, GaAs
- **Adsorbates**: C, N, O, S


Output Formats
--------------

Calculated phase shifts can be exported in several formats:

CLEED Format
   Standard format for the CLEED package (Georg Held's implementation).

ViPErLEED Format
   Compatible with the Vienna Package for LEED analysis.

CSV Format
   Comma-separated values for use in spreadsheets or custom analysis.


Command-Line Alternative
------------------------

For batch processing or integration into scripts, see the :doc:`running`
page for command-line usage of the ``phsh`` programs.


Technical Details
-----------------

The browser calculator works by:

1. Compiling the original Fortran code (``libphsh.f``) to WebAssembly
   using Emscripten
2. Using MEMFS (in-memory filesystem) to handle Fortran file I/O
3. Providing a JavaScript API wrapper for browser interaction
4. Rendering results with Chart.js for interactive visualization

When WebAssembly is not available (unsupported browser or build issues),
the calculator falls back to a demo mode with synthetic phase shift data.

Virtual Filesystem Contract
---------------------------

The WebAssembly build uses MEMFS (an in-memory filesystem) with fixed file
paths for inputs and outputs. These paths are part of the public API and
are especially useful when driving calculations from notebooks (e.g.,
Jupyter + Pyodide) or other custom tooling:

- Inputs
  - ``/input/atorb.i`` — atorb input (charge density setup)
  - ``/input/mufftin.o`` — optional potential override
  - ``/input/phsh.i`` — phase shift input
- Outputs
  - ``/output/atorb.o`` — charge density output
  - ``/output/phasout.o`` — phase shift output
  - ``/output/dataph.o`` — phase shift diagnostics
  - ``/output/inpdat.o`` — input echo data

Roadmap: the Fortran routines currently use fixed filenames. We plan to
refactor them to accept configurable filenames while keeping these defaults
for backwards compatibility.


See Also
--------

- :doc:`running` — Command-line phase shift programs
- :doc:`installing_phaseshifts` — Installation instructions
- :doc:`introduction` — Background on LEED and phase shifts
