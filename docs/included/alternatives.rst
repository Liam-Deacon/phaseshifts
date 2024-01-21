Alternatives
------------

A number of alternatives are available, notably the following:

1. `AQuaLEED <https://physics.mff.cuni.cz/kfpp/povrchy/files/>`_ (with a useful
   `poster overview of phaseshifts calculations <https://physics.mff.cuni.cz/kfpp/povrchy/files/1179-Poster.pdf>`_).
   This is an officially mentioned piece of software on Michel Van Hove's
   `LEED Calculation Homepage <https://www.icts.hkbu.edu.hk/VanHove_files/leed/leedpack.html>`_.
   Furthermore, although the poster mentions that the software is written in python,
   this software is not (currently) distributed on https://PyPI.org (or via alternative means such as a docker image
   on `DockerHub <https://www.docker.com/products/docker-hub/>`_) and therefore harder to
   intergrate with other python LEED-related projects such as `CLEED <https://github.com/Liam-Deacon/CLEED>`_
   and `cleedpy <https://github.com/empa-scientific-it/cleedpy>`_.
2. A fortran program is described in "`McGreevy, E., & Stewart, A.L. (- Apr 1978). <https://inis.iaea.org/search/search.aspx?orig_q=RN:9399501>`_
   A program for calculating elastic scattering phase shifts for an electron colliding with a one-electron target using perturbation theory.
   Computer Physics Communications, 14(1-2), 99-107.", however this code is not publicly available online (pay-walled by journal).

.. note:: Should you know of alternatives, please either
          `open an issue <https://Liam-Deacon/phaseshifts/issues>`_ or
          (better yet) create a PR with changes to this documentation
          to keep this list up to date.
