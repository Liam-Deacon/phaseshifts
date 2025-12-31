.. _running:

*******
Running
*******

The `phsh.py` script (available after installing the package) aims to simplify these
steps with a single command.

The simplest and most reliable cross-platform way to run `phsh.py` is through docker::

  # obtain the image
  docker pull ghcr.io/Liam-Deacon/phaseshifts:latest  # should only need to do this once

  # run phsh.py via the docker image
  docker run ghcr.io/Liam-Deacon/phaseshifts:latest  # will display usage

  # or more generally (adjust as needed)
  docker run ghcr.io/Liam-Deacon//phaseshifts:latest -v /path/to/host/input/data:/data [<phsh-args> ...]


.. tip:: Development docker images can be built locally, e.g.
         :code:`DOCKER_TAG=dev make docker`

.. warning:: There is a `known possible bug <https://github.com/Liam-Deacon/phaseshifts/issues/6>`_
             where the compiled ``libphsh.f`` is not thread-safe (as ascertained by the fortran compiler),
             as such if you anticipate using this library in concurrent environments then it is advised to
             run ``phsh.py`` via :code:`docker run ghcr.io/Liam-Deacon/phaseshifts:latest` as this works around
             this limitation due to the emphereal nature of container instances created using ``docker run``.


Browser-Based Calculator
========================

For users who prefer not to install any software, a browser-based phase shift 
calculator is available. This uses WebAssembly to run the Fortran calculation 
code directly in your web browser.

Features
--------

* **No installation required** - runs entirely in the browser
* **Three calculation methods** - Relativistic, Cavity LEED, and Williams algorithms
* **Interactive visualization** - plot phase shifts vs energy with Chart.js
* **Multiple output formats** - download results in CLEED, ViPErLEED, or CSV format
* **Preset configurations** - quick-start with common elements (Cu, Ni, Fe, Si)

Usage
-----

1. Open the web interface (hosted on GitHub Pages or locally)
2. Select an element from the dropdown (or enter atomic number)
3. Choose the calculation method:
   
   - **Relativistic (phsh_rel)**: Full Dirac equation treatment, essential for heavy elements (Z > 30)
   - **Cavity LEED (phsh_cav)**: Traditional cavity method using Loucks grid
   - **Williams (phsh_wil)**: A.R. Williams' method
   
4. Set the muffin-tin radius and energy range
5. Click "Calculate Phase Shifts"
6. View results as a plot, table, or raw output
7. Download results in your preferred format

Running Locally
---------------

To run the browser interface locally (useful for development or offline use)::

  # From the repository root
  make wasm-serve
  
  # Then open http://localhost:8080/web/ in your browser

Building from Source
--------------------

The WebAssembly version requires Emscripten and f2c to build from source::

  # Install prerequisites (macOS)
  brew install emscripten f2c
  
  # Or install Emscripten from source
  git clone https://github.com/emscripten-core/emsdk.git
  cd emsdk && ./emsdk install latest && ./emsdk activate latest
  source ./emsdk_env.sh
  
  # Build WASM
  make wasm

See ``wasm/README.md`` for detailed build instructions.
