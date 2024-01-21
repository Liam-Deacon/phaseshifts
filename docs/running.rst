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
