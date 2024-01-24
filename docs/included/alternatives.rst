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
2. Elastic Electron-Atom Scattering in Solids and Solid Surfaces
   `(EEASiSSS) <https://www.researchgate.net/profile/John-Rundgren-2/publication/235583683_Optimized_surface-slab_excited-state_muffin-tin_potential_and_surface_core_level_shifts/links/5a266f89a6fdcc8e866bd7e5/Optimized-surface-slab-excited-state-muffin-tin-potential-and-surface-core-level-shifts.pdf>`_
   is authored by John Rundgren and first described in the paper: "J. Rundgren Phys. Rev. B 68 125405 (2003)".
   This program takes a different approach to calculating phase shifts by using optimised muffin-tin potentials
   for surface slabs with preassigned surface core-level shifts.
   Whilst the source code is not publicly available online (to this author's best knowledge), John Rundgren
   has been more than happy to assist when approached in the past.

   .. note:: It would be fantastic to include this software (and document it's use) as part of the phaseshifts python package
             allowing the user to choose the backend they wish to use for calculating phase shifts (e.g. ``EEASiSSS`` or ``phshift2007``).
             As such `John Rundgren <https://www.researchgate.net/profile/John-Rundgren-2>`_ should be contacted to see if
             he would be happy to collaborate on making this possible. This is being tracked by
             `this item <https://github.com/Liam-Deacon/phaseshifts/issues/92>`_.

3. A fortran program is described in "`McGreevy, E., & Stewart, A.L. (- Apr 1978). <https://inis.iaea.org/search/search.aspx?orig_q=RN:9399501>`_
   A program for calculating elastic scattering phase shifts for an electron colliding with a one-electron target using perturbation theory.
   Computer Physics Communications, 14(1-2), 99-107.", however this code is not publicly available online (pay-walled by journal).

.. note:: Should you know of alternatives, please either
          `open an issue <https://Liam-Deacon/phaseshifts/issues>`_ or
          (better yet) create a PR with changes to this documentation
          to keep this list up to date.
