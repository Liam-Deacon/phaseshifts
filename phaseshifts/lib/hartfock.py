#!/usr/bin/env python
# encoding: utf-8
"""
Hartree-Fock radial atomic solver (Python translation of legacy Fortran).

This module implements the self-consistent, radial Hartree-Fock (and
Dirac-Fock-style) workflow used by phaseshifts. The solver operates on a
logarithmic radial grid, integrates radial equations using a Numerov-style
scheme, and produces orbitals and charge densities suitable for LEED/XPD
phase-shift generation.

The logarithmic grid is defined as:

    r_i = r_min * (r_max / r_min) ** (i / n_r)

where ``i = 1..n_r`` (1-based indexing retained from the original Fortran).

References
----------
.. [1] Hartree-Fock method,
       https://en.wikipedia.org/wiki/Hartree%E2%80%93Fock_method
.. [2] Dirac equation (relativistic corrections),
       https://en.wikipedia.org/wiki/Dirac_equation
.. [3] Numerov integration scheme,
       https://en.wikipedia.org/wiki/Numerov%27s_method
.. [4] Clebsch-Gordan coefficients and angular momentum algebra,
       https://en.wikipedia.org/wiki/Clebsch%E2%80%93Gordan_coefficients
.. [5] J. P. Perdew and A. Zunger, "Self-interaction correction to density-
       functional approximations for many-electron systems",
       Phys. Rev. B 23, 5048 (1981), https://doi.org/10.1103/PhysRevB.23.5048
.. [6] Local density approximation overview,
       https://en.wikipedia.org/wiki/Local_density_approximation

TODO
----
- Consider replacing tight Python loops with NumPy vectorized kernels to
  accelerate grid operations and radial integrations.
- Explore SciPy-based ODE solvers or specialized Numerov implementations to
  validate stability and improve maintainability.
"""

from __future__ import print_function, division

from math import cos, pi, pow, exp, log, sinh, sqrt
from copy import deepcopy
import io
import sys
import os

try:
    import builtins
except ImportError:  # pragma: no cover - python 2 fallback
    import __builtin__ as builtins  # type: ignore

raw_input = getattr(builtins, "raw_input", input)  # type: ignore
file = getattr(builtins, "file", io.IOBase)  # type: ignore


def get_input(prompt):
    """
    Read a line of user input with Python 2/3 compatibility.

    Parameters
    ----------
    prompt : str
        Prompt string displayed to the user.

    Returns
    -------
    str
        The raw input line.

    Notes
    -----
    This helper preserves the historical interactive workflow used by the
    original Fortran code.

    References
    ----------
    See `Hartree-Fock method <https://en.wikipedia.org/wiki/Hartree%E2%80%93Fock_method>`_.
    """
    if sys.hexversion > 0x03000000:
        return input(prompt)
    else:
        return raw_input(prompt)


class hartfock(object):
    """
    Hartree-Fock/Dirac-Fock radial solver.

    Notes
    -----
    The solver operates on a logarithmic radial grid with distances in Bohr
    radii. Orbitals are stored in ``phe`` with indices ``phe[i][j]`` for grid
    point ``i`` and orbital index ``j``. The charge density is computed as:

    .. math::

        \\rho(r_i) = \\sum_{j=1}^{n_{el}} \\frac{\\mathrm{occ}_j \\, \\phi_{ij}^2}
        {4 \\pi r_i^2}

    where ``occ_j`` is the orbital occupancy and ``phi_{ij}`` corresponds to
    ``phe[i][j]`` in the code.

    The Dirac equation is solved for the orbitals, while the density is stored
    as ``sqrt(F^2 + G^2)`` with the sign of ``G`` (large component). This
    Dirac-Fock approximation is typically small for valence states. Breit
    interaction terms are neglected.

    References
    ----------
    See `Hartree-Fock method <https://en.wikipedia.org/wiki/Hartree%E2%80%93Fock_method>`_.

    TODO
    ----
    - Consider NumPy vectorization for orbital storage and density evaluation.

    Authors
    -------
    Eric Shirley - original implementation
    Liam Deacon - Python translation
    """

    def __init__(
        self,
        iorbs=33,
        iside=600,
        lmax=4,
        ihmax=20,
        nrmax=4000,
        ntmax=10,
        npmax=60,
        input_stream="stdin",
    ):
        """
        Initialize and run the Hartree-Fock workflow.

        Parameters
        ----------
        iorbs : int
            Maximum number of orbitals allocated in arrays.
        iside : int
            Angular quadrature side parameter carried from the original code.
        lmax : int
            Maximum angular quantum number used in expansions.
        ihmax : int
            Maximum number of history points used by iterative solvers.
        nrmax : int
            Number of radial grid points allocated in arrays.
        ntmax : int
            Number of tabulation points for intermediate arrays.
        npmax : int
            Maximum number of pseudopotential entries to process.
        input_stream : str
            Either ``"stdin"`` for interactive input or a file path.

        Notes
        -----
        This initializer mirrors the legacy Fortran workflow: it allocates
        arrays, reads control records, and dispatches to subroutines such as
        ``abinitio`` to perform the self-consistent field steps.

        References
        ----------
        See `Hartree-Fock method <https://en.wikipedia.org/wiki/Hartree%E2%80%93Fock_method>`_.

        TODO
        ----
        - Replace manual list-of-lists storage with NumPy arrays to enable
          vectorized operations and reduce memory overhead.
        """

        # initialise class variables
        no = [None] * iorbs
        nl = [None] * iorbs
        nm = [None] * iorbs
        xnj = [None] * iorbs
        iss = [None] * iorbs
        ev = [None] * iorbs
        ek = [None] * iorbs
        occ = [None] * iorbs
        r = [None] * nrmax
        dr = [None] * nrmax
        r2 = [None] * nrmax
        rho = [None] * nrmax
        njrc = [None] * 4
        vi = [[None] * 7 for _ in range(nrmax)]
        phe = [[None] * iorbs for _ in range(nrmax)]
        orb = [[None] * iorbs for _ in range(nrmax)]
        rpower = [[None] * 16 for _ in range(nrmax)]

        vctab = []
        nr = nel = rd = nst = iuflag = int
        e = etot = etot2 = zorig = xntot = rmin = rmax = dl = float

        if not input_stream != "stdin" and not os.path.isfile(input_stream):
            raise IOError("'%s' is not a valid input file" % input_stream)
        elif input_stream == "stdin":
            pass
        elif os.path.isfile(input_stream):
            input_file = input_stream
            f = open(input_file, "r")

        # read input
        rel = 0
        iline = 0
        ichar = str
        while ichar != "q":
            ichar = "!"  # initalise to comment

            while ichar == "C" or ichar == "!":
                if os.path.isfile(input_stream):
                    line = f.next()
                    ichar = line.strip()[0]
                    iline += 1
                else:
                    ichar = get_input("Enter command character ('q' to quit): ")

            if ichar == "d":
                if os.path.isfile(input_stream):
                    rel = f.next().split("!")[0]
                else:
                    rel = int(get_input("Please enter relativity flag: "))

            # get exchange correlation (Alpha)
            elif ichar == "x":
                if os.path.isfile(input_stream):
                    alfa = float("".join([ch for ch in f.next().split("!")[0] if ch.lower() != "d"]))
                else:
                    alfa = float(get_input("Enter exchange correlation method (0=HARTREE-FOCK, >0=LDA, <0=XALPHA): "))

            elif ichar == "a":
                (
                    etot,
                    nst,
                    rel,
                    alfa,
                    nr,
                    r,
                    dr,
                    r2,
                    dl,
                    e,
                    njrc,
                    vi,
                    zorig,
                    xntot,
                    nel,
                    no,
                    nl,
                    xnj,
                    ev,
                    occ,
                    iss,
                    ek,
                    orb,
                    iuflag,
                    rpower,
                    nm,
                    phe,
                    etot2,
                ) = abinitio(
                    etot,
                    nst,
                    rel,
                    alfa,
                    nr,
                    r,
                    dr,
                    r2,
                    dl,
                    e,
                    njrc,
                    vi,
                    zorig,
                    xntot,
                    nel,
                    no,
                    nl,
                    xnj,
                    ev,
                    occ,
                    iss,
                    ek,
                    orb,
                    iuflag,
                    rpower,
                    nm,
                    phe,
                    etot2,
                    input_stream=f,
                )

            elif ichar == "i":
                (zorig, nr, rmin, rmax, r, dr, r2, dl, njrc, xntot, nel) = initiali(
                    zorig,
                    nr,
                    rmin,
                    rmax,
                    r,
                    dr,
                    r2,
                    dl,
                    njrc,
                    xntot,
                    nel,
                    input_stream=f,
                )

            elif ichar == "q":
                return  # quit

            elif ichar == "w":
                ixflag = 1
                iu = -1
                ir = 0
                hfdisk(
                    iu,
                    ir,
                    etot,
                    nst,
                    rel,
                    nr,
                    rmin,
                    rmax,
                    r,
                    rho,
                    zorig,
                    xntot,
                    ixflag,
                    nel,
                    no,
                    nl,
                    xnj,
                    iss,
                    ev,
                    ek,
                    occ,
                    njrc,
                    vi,
                    phe,
                    orb,
                    input_stream=f,
                )

            elif ichar == "r":
                iu = -1
                ir = 1
                hfdisk(
                    iu,
                    ir,
                    etot,
                    nst,
                    rel,
                    nr,
                    rmin,
                    rmax,
                    r,
                    rho,
                    zorig,
                    xntot,
                    ixflag,
                    nel,
                    no,
                    nl,
                    xnj,
                    iss,
                    ev,
                    ek,
                    occ,
                    njrc,
                    vi,
                    phe,
                    orb,
                )
                setgrid(nr, rmin, rmax, r, dr, r2, dl)

            elif ichar == "u":
                if os.path.isfile(input_stream):
                    iuflag = int(f.next().split("!")[0])
                else:
                    iuflag = int(get_input("Please enter IUFLAG (0=U, 1=SU, 2=R): "))

            elif ichar == "c":
                if os.path.isfile(input_stream):
                    (corpol, rs, rp, sd) = f.next().split("!")[0].split()[:3]
                else:
                    while True:
                        try:
                            (corpol, rs, rp, sd) = get_input("enter ALPHA, RS, RP, RD: ").split()[:3]
                            break
                        except Exception:
                            print("Invalid input - please retry...")

                for k in range(0, nr):
                    fs = pow(1.0 - exp(-pow((r[k] / rs), 2.0)), 2.0)
                    fp = pow(1.0 - exp(-pow((r[k] / rp), 2.0)), 2.0)
                    fd = pow(1.0 - exp(-pow((r[k] / rd), 2.0)), 2.0)
                    vctab[k][0] = -corpol / 2.0 * fs * fs / pow(r[k], 4.0)
                    vctab[k][1] = -corpol / 2.0 * fp * fp / pow(r[k], 4.0)
                    vctab[k][2] = -corpol / 2.0 * fd * fd / pow(r[k], 4.0)

            elif ichar == "f":
                if os.path.isfile(input_stream):
                    (iunit, corpol) = [t(s) for t, s, in zip((int, float), f.next().split("!")[0].split()[:1])]
                    (ilev, inum, eold) = [t(s) for t, s in zip((int, int, float), f.next().split("!")[0].split()[:2])]
                else:
                    while True:
                        try:
                            (iunit, corpol) = [
                                t(s)
                                for t, s in zip(
                                    (int, float),
                                    get_input("Please enter IUNIT, CORPOL: ").split()[:2],
                                )
                            ]
                            break
                        except Exception:
                            print("incorrect input - please retry...")

                    while True:
                        try:
                            (ilev, inum, eold) = [
                                t(s)
                                for t, s in zip(
                                    (int, int, float),
                                    get_input("Please enter ILEV, INUM, EOLD: "),
                                )
                            ]
                            break
                        except Exception:
                            print("incorrect input - please retry...")

                xl = nl[ilev]
                if inum == 1:
                    eav = f.next().split("!")[0]
                else:
                    (e1, e2) = f.next().split("!")[0].split()[:1]
                    eav = (e1 * xl + e2 * (xl + 1.0)) / (xl + xl + 1.0)

                if eav < 0.0:
                    eav = eold + eav

                if iunit == 2:
                    eav = eav / 2.0
                elif iunit == 3:
                    eav = eav / 27.2116
                elif iunit == 4:
                    eav = eav * 0.000123985 / 27.2116

                sd = abs(abs(eav) - abs(ev[ilev]))
                rl = sl = sh = 0.0
                rh = 10.0

                sc = abs(1 + sd * 1000000)  # force initial loop
                while abs(sc - sd) > 0.000001:
                    if sl * sh <= 0.00000001:
                        rc = rl + (rh - rl) / 2.0
                    if sl * sh > 0.00000001:
                        rc = rl + (rh - rl) * (sd - sl) / (sh - sl)
                    sc = 0.0
                    for i in range(1, nr + 1):
                        f = pow(1.0 - exp(-pow((r[i] / rc), 2.0)), 2.0)
                        vcpp = corpol / (2.0 * pow(r[i], 4.0)) * f * f
                        sc += dr[i] * phe[i][ilev] * phe[i][ilev] * vcpp

                    if sc > sd:
                        rl = rc
                        sl = sc
                    elif sc < sd:
                        rh = rc
                        sh = sc
                    print("{} {}".format(rc, sc))

            elif ichar == "p":
                pseudo(
                    etot,
                    nst,
                    rel,
                    alfa,
                    nr,
                    rmin,
                    rmax,
                    r,
                    dr,
                    r2,
                    dl,
                    phe,
                    orb,
                    njrc,
                    vi,
                    zorig,
                    xntot,
                    nel,
                    no,
                    nl,
                    xnj,
                    ev,
                    occ,
                    iss,
                    ek,
                    iuflag,
                    vctab,
                )

            elif ichar == "g":
                iu = f.next().split("!")[0]
                jive = f.next().split("!")[0]  # format(1x,1a11)
                jive2 = f.next().split("!")[0]  # format(1x,1a60)
                jive3 = f.next().split("!")[0]  # format(1x,1a70)

                zizv = abs(r[nr - 1] * vi[nr - 1][1])
                print("{}".format(jive))  # write to iu file
                print("{}".format(jive2))
                print("{}".format(jive3))
                print(3, nr, zizv)
                for i in range(1, nr + 1):
                    print(r[i])
                for k in range(1, nr + 1):
                    print(0, vi[k][1])
                    print(1, vi[k][3])
                    print(2, vi[k][5])
                    print(0.0)
                for j in range(1, nr + 1):
                    rh = 0.0
                    for k in range(1, nel + 1):
                        rh += phe[j][k] * phe[j][k] * occ[k]
                    print("{}".format(rh))

            elif ichar == "v":
                for k in range(1, nr + 1):
                    print(r[k], vi[k][1] * r[k])
                    print(r[k], vi[k][3] * r[k])
                    print(r[k], vi[k][5] * r[k])

            elif ichar == "V":
                fourier(nr, r, dr, r2, vi)

            else:
                print(
                    "'%s' is not a valid command - valid characters are:\n"
                    "a: do abinitio\n"
                    "c: corpol values\n"
                    "d: relativistic switch\n"
                    "f: \n"
                    "g: write to iu file\n"
                    "i: initialise\n"
                    "p: pseudo\n"
                    "u: \n"
                    "v: \n"
                    "V: fourier\n"
                    "q: quit\n" % ichar
                )


def abinitio(
    etot,
    nst,
    rel,
    alfa,
    nr,
    r,
    dr,
    r2,
    dl,
    e,
    njrc,
    vi,
    zorig,
    xntot,
    nel,
    no,
    nl,
    xnj,
    ev,
    occ,
    iss,
    ek,
    orb,
    iuflag,
    rpower,
    nm,
    phe,
    etot2,
    input_stream="stdin",
):
    """
    Run the self-consistent-field (SCF) loop for atomic orbitals.

    Parameters
    ----------
    etot : float
        Total energy accumulator.
    nst : int
        Number of shells to process.
    rel : int
        Relativistic flag (0 for non-relativistic).
    alfa : float
        Exchange-correlation mixing parameter.
    nr : int
        Number of radial grid points.
    r, dr, r2 : list of float
        Radial grid, differential spacing, and squared radius arrays.
    dl : float
        Logarithmic grid spacing.
    e : list of float
        Energy grid container.
    njrc : list of int
        Core radius indices per angular momentum channel.
    vi : list of list of float
        Potential array by grid and channel.
    zorig : float
        Atomic number.
    xntot : float
        Total electron count accumulator.
    nel : int
        Number of electrons.
    no, nl, nm : list of int
        Principal, angular, and magnetic quantum numbers.
    xnj : list of float
        Total angular momentum values (j).
    ev : list of float
        Orbital eigenvalues.
    occ : list of float
        Orbital occupancies.
    iss : list of int
        Spin flags.
    ek : list of float
        Kinetic energy contributions.
    orb : list of list of float
        Orbital buffer for intermediate results.
    iuflag : int
        Output control flag.
    rpower : list of list of float
        Precomputed powers of radius.
    nm : list of int
        Magnetic quantum numbers array.
    phe : list of list of float
        Radial wavefunction storage.
    etot2 : float
        Second energy accumulator.
    input_stream : str
        Input source; ``"stdin"`` for interactive use or file handle.

    Returns
    -------
    tuple
        Updated SCF state (energies, grids, orbitals, and bookkeeping).

    Notes
    -----
    The routine initializes radial powers, reads orbital quantum numbers,
    and iteratively calls ``atsolve`` until the energy error falls below the
    tolerance ``etol``.

    References
    ----------
    See `Hartree-Fock method <https://en.wikipedia.org/wiki/Hartree%E2%80%93Fock_method>`_.

    TODO
    ----
    - Replace nested Python loops with NumPy vectorization for rpower
      construction and orbital initialization.
    - Consider SciPy-based mixing schemes for faster SCF convergence.
    """

    xntot = 0.0
    eerror = evi = float

    # this will be good for going up to and including l=3...
    for i in range(8):
        xi = i
        for k in range(len(r)):
            rpower[k][i] = pow(r[k], xi)

    # read in nfc, nel.  - refer to the documentation for their meanings.
    if input == "stdin":
        (nfc, nel, ratio, etol, xnum) = [
            t(s)
            for t, s in zip(
                (int, int, float, float, float),
                get_input("Please enter NFC NEL RATIO ETOL XNUM: "),
            )
        ]
    elif isinstance(input_stream, file):
        f = input_stream
        (nfc, nel, ratio, etol, xnum) = [
            t(s) for t, s in zip((int, int, float, float, float), f.next().split("!")[0].split()[:5])
        ]

    # for all of the electrons, read in the quantum numbers.
    # get the total Hartree-active charge.  initialise eigenvalues.
    try:
        i = 0
        if input == "stdin":
            for i in range(nfc, nel):
                (no[i], nl[i], nm[i], xnj[i], iss[i], occ[i]) = [
                    t(s)
                    for t, s in zip(
                        (int, int, int, float, float, float),
                        get_input("Please enter [%i] N L M J S OCC: " % i),
                    )
                ]
            ev[i] = 0.0
            xntot += occ[i]
            for j in range(1, nr + 1):
                phe[j][i] = 0.0
                orb[j][i] = 0.0

        else:
            for i in range(nfc, nel):
                (no[i], nl[i], nm[i], xnj[i], iss[i], occ[i]) = [
                    t(s)
                    for t, s in zip(
                        (int, int, int, float, float, float),
                        f.next().split("!")[0].split()[:6],
                    )
                ]
                ev[i] = 0.0
                xntot += occ[i]
                for j in range(1, nr + 1):
                    phe[j][i] = 0.0
                    orb[j][i] = 0.0
    except TypeError:
        raise TypeError("Problem loading N L M J S OCC - entry %i")

    # initialise the parameters for self-consistency loop.
    # ratio is the mixture of old and new field mixing.
    while True:
        atsolve(
            etot,
            nst,
            rel,
            alfa,
            eerror,
            nfc,
            nr,
            r,
            dr,
            r2,
            dl,
            phe,
            njrc,
            vi,
            zorig,
            xntot,
            nel,
            no,
            nl,
            nm,
            xnj,
            ev,
            occ,
            iss,
            ek,
            ratio,
            orb,
            rpower,
            xnum,
            etot2,
            iuflag,
            evi,
        )

        eerror *= (1.0 - ratio) / ratio
        print(" %14.6f%14.6f" % (eerror, etot))  # format (1x,3f14.6)
        if eerror <= etol:
            break

    # write out information about the atom.
    for i in range(1, nel + 1):
        nj = xnj[i] * 2
        print("  %4i%2i%4i%4i%s%4i%10.4f%18.6f\n" % (no[i], nl[i], nm[i], nj, "/2", iss[i], occ[i], ev[i]))
        print("Total energy =  %14.6f  %14.6f" % (etot, etot * 27.2116))

    return (
        etot,
        nst,
        rel,
        alfa,
        nr,
        r,
        dr,
        r2,
        dl,
        e,
        njrc,
        vi,
        zorig,
        xntot,
        nel,
        no,
        nl,
        xnj,
        ev,
        occ,
        iss,
        ek,
        orb,
        iuflag,
        rpower,
        nm,
        phe,
        etot2,
    )


def atsolve(
    etot,
    nst,
    rel,
    alfa,
    eerror,
    nfc,
    nr,
    r,
    dr,
    r2,
    dl,
    phe,
    njrc,
    vi,
    zorig,
    xntot,
    nel,
    no,
    nl,
    nm,
    xnj,
    ev,
    occ,
    iss,
    ek,
    ratio,
    orb,
    rpower,
    xnum,
    etot2,
    iuflag,
    evi,
):
    """
    Solve atomic orbitals for one SCF iteration.

    Parameters
    ----------
    etot : float
        Total energy accumulator (updated in-place).
    nst : int
        Number of shells to process.
    rel : int
        Relativistic flag (0 for non-relativistic).
    alfa : float
        Exchange-correlation mixing parameter.
    eerror : float
        Current eigenvalue error estimate.
    nfc : int
        Number of frozen core orbitals.
    nr : int
        Number of radial grid points.
    r, dr, r2 : list of float
        Radial grid, differential spacing, and squared radius arrays.
    dl : float
        Logarithmic grid spacing.
    phe : list of list of float
        Radial wavefunction storage.
    njrc : list of int
        Core radius indices per angular momentum channel.
    vi : list of list of float
        Potential array by grid and channel.
    zorig : float
        Atomic number.
    xntot : float
        Total electron count accumulator.
    nel : int
        Number of electrons.
    no, nl, nm : list of int
        Principal, angular, and magnetic quantum numbers.
    xnj : list of float
        Total angular momentum values (j).
    ev : list of float
        Orbital eigenvalues.
    occ : list of float
        Orbital occupancies.
    iss : list of int
        Spin flags.
    ek : list of float
        Kinetic energy contributions.
    ratio : float
        Mixing ratio for the potential update.
    orb : list of list of float
        Orbital buffer for intermediate results.
    rpower : list of list of float
        Precomputed powers of radius.
    xnum : float
        Mixing parameter for orbital updates.
    etot2 : float
        Second energy accumulator.
    iuflag : int
        Output control flag.
    evi : float
        Current eigenvalue estimate.

    Returns
    -------
    tuple
        Updated SCF state for the current iteration.

    Notes
    -----
    This routine applies ``setqmm`` to construct effective potentials and
    invokes ``elsolve`` for each orbital, updating eigenvalues and orbitals.

    References
    ----------
    See `Hartree-Fock method <https://en.wikipedia.org/wiki/Hartree%E2%80%93Fock_method>`_.

    TODO
    ----
    - Evaluate NumPy-based storage for ``q0``, ``xm1``, and ``xm2`` to reduce
      Python overhead in repeated loops.
    - Consider SciPy root-finding methods to bracket eigenvalues in ``elsolve``.
    """
    # initialise arrays
    q0 = xm1 = xm2 = v = [None] * len(r)

    # initialise eerror, the biggest change in an eigenvalue, and etot.
    eerror = etot = zeff = 0.0

    # run through all the orbitals.  calculate those not in the core.
    for i in range(nel):
        if i > nfc:
            idoflag = 1
            ctx = SetqmmContext(v, zeff, zorig, rel, nr, r, r2, dl, q0, xm1, xm2, njrc, vi)
            setqmm(i, orb, nl[i], idoflag, ctx)
            zeff = ctx.zeff

        xkappa = -1.0
        if abs(xnj[i]) > nl[i] + 0.25:
            xkappa = -nl[i] - 1
        if abs(xnj[i]) < nl[i] - 0.25:
            xkappa = nl[i]

        print(
            i,
            occ[i],
            no[i],
            nl[i],
            xkappa,
            xnj[i],
            zorig,
            zeff,
            evi,
            phe[1][i],
            v,
            q0,
            xm1,
            xm2,
            nr,
            r,
            dr,
            r2,
            dl,
            rel,
        )

        (
            i,
            occ[i],
            no[i],
            nl[i],
            xkappa,
            xnj[i],
            zorig,
            zeff,
            evi,
            phe[1][i],
            v,
            q0,
            xm1,
            xm2,
            nr,
            r,
            dr,
            r2,
            dl,
            rel,
        ) = elsolve(
            i,
            occ[i],
            no[i],
            nl[i],
            xkappa,
            xnj[i],
            zorig,
            zeff,
            evi,
            phe[1][i],
            v,
            q0,
            xm1,
            xm2,
            nr,
            r,
            dr,
            r2,
            dl,
            rel,
        )

        if abs(ev[i] - evi) > eerror:
            eerror = abs(ev[i] - evi)
        ev[i] = evi

        ekk = 0.0
        ll = 2
        for j in range(nr, 1 - 1, -1):
            dq = phe[j][i] * phe[j][i]
            ekk = ekk + (evi - orb[j][i]) * dr[j] * dq * ll / 3.0
            ll = 6 - ll
        ek[i] = ekk

        # add the kinetic to total, including the frozen core kinetic energy.
        etot += ek[i] * occ[i]
        getpot(
            etot,
            nst,
            rel,
            alfa,
            dl,
            nr,
            dr,
            r,
            r2,
            xntot,
            phe,
            ratio,
            orb,
            occ,
            iss,
            nel,
            nl,
            nm,
            no,
            xnj,
            rpower,
            xnum,
            etot2,
            iuflag,
        )

    return (
        etot,
        nst,
        rel,
        alfa,
        eerror,
        nfc,
        nr,
        r,
        dr,
        r2,
        dl,
        phe,
        njrc,
        vi,
        zorig,
        xntot,
        nel,
        no,
        nl,
        nm,
        xnj,
        ev,
        occ,
        iss,
        ek,
        ratio,
        orb,
        rpower,
        xnum,
        etot2,
        iuflag,
        evi,
    )


def getpot(
    etot,
    nst,
    rel,
    alfa,
    dl,
    nr,
    dr,
    r,
    r2,
    xntot,
    phe,
    ratio,
    orb,
    occ,
    iss,
    nel,
    nl,
    nm,
    no,
    xnj,
    rpower,
    xnum,
    etot2,
    iuflag,
):
    """
    Build the self-consistent potential and update orbital energies.

    Parameters
    ----------
    etot : float
        Total energy accumulator.
    nst : int
        Number of shells to process.
    rel : int
        Relativistic flag (0 for non-relativistic).
    alfa : float
        Exchange-correlation mixing parameter.
    dl : float
        Logarithmic grid spacing.
    nr : int
        Number of radial grid points.
    dr, r, r2 : list of float
        Radial grid spacing, coordinates, and squared radius arrays.
    xntot : float
        Total electron count.
    phe : list of list of float
        Radial wavefunction storage.
    ratio : float
        Mixing ratio for new vs. old potentials.
    orb : list of list of float
        Orbital buffer updated with mixed potentials.
    occ : list of float
        Orbital occupancies.
    iss : list of int
        Spin flags.
    nel : int
        Number of electrons.
    nl, nm, no : list of int
        Angular, magnetic, and principal quantum numbers.
    xnj : list of float
        Total angular momentum values (j).
    rpower : list of list of float
        Precomputed powers of radius.
    xnum : float
        Mixing parameter for orbital updates.
    etot2 : float
        Secondary energy accumulator.
    iuflag : int
        Output control flag.

    Returns
    -------
    tuple
        Updated potential, orbital, and energy state.

    Notes
    -----
    The routine builds direct and exchange terms using Clebsch-Gordan
    coefficients and radial integrals, then mixes the potential using
    ``ratio`` to stabilize SCF convergence.

    References
    ----------
    See `Clebsch-Gordan coefficients <https://en.wikipedia.org/wiki/Clebsch%E2%80%93Gordan_coefficients>`_
    for angular coupling.

    TODO
    ----
    - Replace Clebsch-Gordan loops with ``scipy.special.wigner_3j`` where
      available for clarity and maintainability.
    - Use NumPy vectorization for radial integrals and potential mixing.
    """

    cg = [[[[[None] * 13] * 13] * 13] * 6] * 6
    pin = [[[None] * 17] * 9] * 9
    xq1 = xq2 = xq0 = xqi1 = xqi2 = xqi0 = xqj1 = xqj2 = xqj0 = [None] * nr

    # calculate Clebsch-Gordon coefficients
    clebschgordan(nel, nl, cg)  # update cg matrix only
    getillls(pin)  # update pin matrix only

    ratio1 = 1.0 - ratio
    for i in range(1, nel + 1):
        for k in range(1, nr + 1):
            orb[k][i] = ratio1 * orb[k][i]

    for i in range(1, nel + 1):
        li = nl[i]
        mi = nm[i]

        jstart = i + 1
        if xnj[i] < 0.0 or occ[i] > 1.0 or abs(alfa) > 0.001:
            jstart = i

        for j in range(jstart, nel + 1):  # 2990
            if occ[i] == 0.0 and occ[j] == 0.0:
                continue  # goto 2990

            lj = nl[j]
            mj = nm[j]

            # direct coulomb
            lmx = 2 * li
            if li > lj:
                lmx = 2 * lj

            # l=0 is monopole or spherical term for direct coulomb.  Therefore,
            # when we have occ[i] or occ[j] greater than one, set lmx=0.
            if (occ[i] > 1.0) or (occ[j] > 1.0) or (xnj[i] < 0.0) or (xnj[j] < 0.0):
                lmx = 0

            for la in range(lmx, 0 - 1, -2):  # 2550
                lap = la + 1
                coeff = (
                    ((li + li + 1) * (lj + lj + 1))
                    / pow((la + la + 1), 2.0)
                    * cg[li][li][la][mi][-mi]
                    * cg[lj][lj][la][mj][-mj]
                    * cg[li][li][la][0][0]
                    * cg[lj][lj][la][0][0]
                )

            if mi + mj != int(2 * ((mi + mj) / 2.0)):
                coeff = -coeff
            if i == j:
                coeff /= 2.0
            coeffi = occ[i] * coeff
            coeffj = occ[j] * coeff
            ri = ratio * coeffi
            rj = ratio * coeffj
            rc = coeff * occ[i] * occ[j]

            xouti = 0.0
            xoutj = 0.0
            for k in range(1, nr + 1):
                xqi0[k] = dr[k] * phe[k][i] * phe[k][i] / 2.0
                xqi1[k] = xqi0[k] * rpower[k][la]
                if rpower[k][lap] != 0.0:
                    xqi2[k] = xqi0[k] / rpower[k][lap]
                else:
                    xqi2[k] = 0.0

                xouti += xqi2[k]
                xqj0[k] = dr[k] * phe[k][j] * phe[k][j] / 2.0
                xqj1[k] = xqj0[k] * rpower[k][la]
                if rpower[k][lap] != 0.0:
                    xqj2[k] = xqj0[k] / rpower[k][lap]
                else:
                    xqj2[k] = 0.0

                xoutj += xqj2[k]

            xinti = xqi1[1]
            xintj = xqj1[1]
            xouti = 2.0 * xouti - xqi2[1]
            xoutj = 2.0 * xoutj - xqj2[1]

            for k in range(2, nr + 1):
                xinti += xqi1[k] + xqi1[k - 1]
                xouti -= xqi2[k] - xqi2[k - 1]
                vali = xouti * rpower[k][la]
                if rpower[k][lap] != 0.0:
                    vali += xinti / rpower[k][lap]
                orb[k][j] += ri * vali

                xintj = xintj + xqj1[k] + xqj1[k - 1]
                xoutj = xoutj - xqj2[k] - xqj2[k - 1]
                valj = xoutj * rpower[k][la]
                if rpower[k][lap] != 0.0:
                    valj += xintj / rpower[k][lap]
                orb[k][i] += rj * valj

                etot = etot + rc * (xqi0[k] * valj + xqj0[k] * vali)

            if iss[i] != iss[j] and occ[i] <= 1.0 and occ[j] <= 1.0 and xnj[i] >= 0.0 and xnj[j] >= 0.0:
                continue  # goto 2990
            if abs(alfa) >= 0.001:
                continue  # goto 2990

            # exchange interaction
            lmx = li + lj
            lmin = abs(mi - mj)
            if occ[i] > 1.0 or occ[j] > 1.0 or xnj[i] < 0.0 or xnj[j] < 0.0:
                lmin = 0
            for la in range(lmx, lmin - 1, -2):
                lap = la + 1

                coeff = float((li + li + 1) * (lj + lj + 1)) / float(
                    pow(la + la + 1, 2.0) * pow(cg[li][lj][la][-mi][mj] * cg[li][lj][la][0][0], 2.0)
                )
                if occ[i] > 1.0 or occ[j] > 1.0 or xnj[i] < 0.0 or xnj[j] < 0.0:
                    coeff = pin[li][lj][la] / 4.0
                if i == j:
                    coeff /= 2.0
                coeffi = occ[i] * coeff
                coeffj = occ[j] * coeff
                ri = ratio * coeffi
                rj = ratio * coeffj
                rc = coeff * occ[i] * occ[j]
                xnum2 = xnum * xnum

                xout = 0.0
                for k in range(1, nr + 1):
                    xq0[k] = dr[k] * phe[k][i] * phe[k][j] / 2.0
                    xq1[k] = xq0[k] * rpower[k][la]
                    if rpower[k][lap] != 0.0:
                        xq2[k] = xq0[k] / rpower[k][lap]
                    else:
                        xq2[k] = 0.0

                    xout += xq2[k]

                xint = xq1[1]
                xout = 2.0 * xout - xq2[1]
                for k in range(2, nr + 1):
                    xint += xq1[k] + xq1[k - 1]
                    xout -= xq2[k] - xq2[k - 1]
                    if xq0[k] != 0.0:
                        val = xout * rpower[k][la]
                    if rpower[k][lap] != 0.0:
                        val += xint / rpower[k][lap]
                    etot -= 2.0 * xq0[k] * rc * val
                    xx = phe[k][j] / phe[k][i]
                    if abs(xx) / xnum > 1.0:
                        orb[k][i] -= rj * xnum2 / xx * val
                    else:
                        orb[k][i] -= rj * xx * val

                    xx = phe[k][i] / phe[k][j]
                    if abs(xx) / xnum > 1.0:
                        orb[k][j] -= ri * xnum2 / xx * val
                    else:
                        orb[k][j] -= ri * xx * val

    # here we compute the charge density, if needed, for treating
    # exchange/correlation in a local fashion...
    if abs(alfa) >= 0.001:
        if alfa > 0.0:
            fx = 1.0
            fc = 1.0
        else:
            fx = 1.5 * abs(alfa)
            fc = 0.0

    # note: we don't deal with spin-polarization in local exchange
    # picture, since local exchange is totally wrong for such
    # effects, anyway.  local exchange pretends charge density
    # is paramagnetic.  also, local exchange treats everything
    # as spherical.

    ex = ec = ux1 = ux2 = uc1 = uc2 = 0.0  # initialise

    for i in range(1, nr + 1):
        xn = 0.0
        for j in range(1, nel + 1):
            xn += phe[i][j] * phe[i][j] * occ[j]
        xn1 = xn / 2.0
        xn2 = xn / 2.0
        nst = 2
        (nst, rel, r2[i], xn1, xn2, ex, ec, ux1, ux2, uc1, uc2) = exchcorr(
            nst, rel, r2[i], xn1, xn2, ex, ec, ux1, ux2, uc1, uc2
        )
        exc = fx * ex + fc * ec
        uxc = fx * ux1 + fc * uc1
        etot = etot + dr[i] * xn * exc
        for j in range(1, nel + 1):
            orb[i][j] += uxc * ratio

    for i in range(1, nr + 1):
        if iuflag:
            jj = 1
        ii = jj  # 8960
        icond = True
        while icond:  # 8965
            if ii != nel:
                # goto 8970
                icond = False
                if no[jj] == no[ii + 1] and nl[jj] == nl[ii + 1] and iuflag == 2:
                    icond = True
                if no[jj] == no[ii + 1] and nl[jj] == nl[ii + 1] and iss[jj] == iss[ii + 1] and iuflag == 1:
                    icond = True
                if icond:
                    ii += 1

        orba = 0.0  # 8970
        div = 0.0
        for k in range(jj, ii + 1):
            div += occ[k]
            orba += orb[i][k] * occ[k]

        if div != 0.0:
            orba /= div
            for k in range(jj, ii + 1):
                orb[i][k] = orba

        if ii != nel:
            jj = ii + 1
            continue  # goto 8960

    return (
        etot,
        nst,
        rel,
        alfa,
        dl,
        nr,
        dr,
        r,
        r2,
        xntot,
        phe,
        ratio,
        orb,
        occ,
        iss,
        nel,
        nl,
        nm,
        no,
        xnj,
        rpower,
        xnum,
        etot2,
        iuflag,
    )


def elsolve(
    i,
    occ,
    n,
    l,
    xkappa,
    xj,
    zorig,
    zeff,
    e,
    phi,
    v,
    q0,
    xm1,
    xm2,
    nr,
    r,
    dr,
    r2,
    dl,
    rel,
):
    """
    Solve a radial orbital eigenvalue using bracketing and Numerov steps.

    Parameters
    ----------
    i : int
        Orbital index.
    occ : float
        Orbital occupancy.
    n : int
        Principal quantum number.
    l : int
        Angular momentum quantum number.
    xkappa : float
        Relativistic kappa parameter.
    xj : float
        Total angular momentum value.
    zorig : float
        Atomic number.
    zeff : float
        Effective nuclear charge.
    e : float
        Trial eigenvalue (updated internally).
    phi : list of float
        Radial wavefunction (updated in-place).
    v : list of float
        Effective potential.
    q0, xm1, xm2 : list of float
        Auxiliary arrays for the Numerov scheme.
    nr : int
        Number of radial grid points.
    r, dr, r2 : list of float
        Radial grid, spacing, and squared radius arrays.
    dl : float
        Logarithmic grid spacing.
    rel : float
        Relativistic switch (0.0 for non-relativistic).

    Returns
    -------
    tuple
        Updated orbital state, including ``phi`` and eigenvalue ``e``.

    Notes
    -----
    The routine uses a bisection-like update on the energy bracket and relies
    on ``integ`` to count nodes and detect over/under-shooting.

    References
    ----------
    See `Numerov's method <https://en.wikipedia.org/wiki/Numerov%27s_method>`_.

    TODO
    ----
    - Consider SciPy root-finding (``scipy.optimize.bisect``) to improve
      readability and convergence diagnostics.
    """
    el = -zorig * zorig / float(n * n)
    eh = 0.0
    etol = 0.0000000001

    integ_ctx = IntegContext(zeff, v, xm1, xm2, nr, r, r2, dl, rel)

    while True:
        e = (el + eh) / 2.0  # label 155
        istop = 0
        integ_result = integ(e, l, xkappa, n, istop, phi, integ_ctx)
        nn = integ_result[0]
        ief = integ_result[2]
        if nn < n - l - 1:
            ief = -1
            if ief != 1:  # label 200
                el = e
            if el > -0.001:
                print("Mixing too strong for level : %i" % i)
                return

            if ief != -1:
                eh = e
            if eh - el > etol:
                continue  # goto 155
            if abs(abs(xj) - abs(float(l))) > 0.25:
                augment(e, l, xj, phi, v, nr, r, dl)  # adjust phi array only
            aa = 0.0
            for j in range(1, nr + 1):
                aa += phi[j] * phi[j] * dr[j]

            xnorm = sqrt(aa)
            for j in range(1, nr + 1):
                phi[j] /= xnorm

        break

    return (
        i,
        occ,
        n,
        l,
        xkappa,
        xj,
        zorig,
        zeff,
        e,
        phi,
        v,
        q0,
        xm1,
        xm2,
        nr,
        r,
        dr,
        r2,
        dl,
        rel,
    )


def augment(e, l, xj, phi, v, nr, r, dl):
    """
    Apply relativistic augmentation to the radial wavefunction.

    Parameters
    ----------
    e : float
        Orbital eigenvalue.
    l : int
        Angular momentum quantum number.
    xj : float
        Total angular momentum value.
    phi : list of float
        Radial wavefunction (updated in-place).
    v : list of float
        Effective potential.
    nr : int
        Number of radial grid points.
    r : list of float
        Radial grid.
    dl : float
        Logarithmic grid spacing.

    Returns
    -------
    None
        Updates ``phi`` in-place.

    Notes
    -----
    This routine computes the small-component correction using finite
    differences on the logarithmic grid.

    References
    ----------
    See `Hartree-Fock method <https://en.wikipedia.org/wiki/Hartree%E2%80%93Fock_method>`_.

    TODO
    ----
    - Replace finite-difference stencils with NumPy vectorization to reduce
      per-grid-point overhead.
    """
    phi2 = [None] * len(phi)

    c = 137.038
    cc = c * c
    c2 = cc + cc
    xkappa = -1
    if abs(xj) > l + 0.25:
        xkappa = -l - 1
    if abs(xj) < l - 0.25:
        xkappa = l
    for j in range(4, nr - 3 + 1):
        if phi[j] != 0.0:
            g0 = phi[j]
            ga = phi[j + 1] - phi[j - 1]
            gb = (phi[j + 2] - phi[j - 2]) / 2.0
            gc = (phi[j + 3] - phi[j - 3]) / 3.0
            gg = ((1.5 * ga - 0.6 * gb + 0.1 * gc) / (2.0 * dl) + xkappa * g0) / r[j]
            f0 = c * gg / (e - v[j] + c2)
            phi2[j] = sqrt(g0 * g0 + f0 * f0)
            if g0 < 0.0:
                phi2[j] = -phi2[j]
            else:
                phi2[j] = phi[j]

    for j in range(1, 3 + 1):
        phi2[j] *= phi[4] / phi2[4]

    phi = phi2

    return


class SetqmmContext(object):
    """
    Context container for the ``setqmm`` routine.

    Parameters
    ----------
    v : list of float
        Effective potential array.
    zeff : float
        Effective nuclear charge (updated in-place).
    zorig : float
        Atomic number.
    rel : float
        Relativistic flag value.
    nr : int
        Number of radial grid points.
    r, r2 : list of float
        Radial grid and squared radius arrays.
    dl : float
        Logarithmic grid spacing.
    q0, xm1, xm2 : list of float
        Auxiliary Numerov arrays.
    njrc : list of int
        Core radius indices per angular momentum channel.
    vi : list of list of float, optional
        Potential table by grid and channel.

    Notes
    -----
    This lightweight container reduces the argument count for ``setqmm``.

    References
    ----------
    See `Hartree-Fock method <https://en.wikipedia.org/wiki/Hartree%E2%80%93Fock_method>`_.

    TODO
    ----
    - Consider converting this structure to a small ``dataclass`` once Python
      2.7 support is dropped.
    """

    def __init__(self, v, zeff, zorig, rel, nr, r, r2, dl, q0, xm1, xm2, njrc, vi=None):
        self.v = v
        self.zeff = zeff
        self.zorig = zorig
        self.rel = rel
        self.nr = nr
        self.r = r
        self.r2 = r2
        self.dl = dl
        self.q0 = q0
        self.xm1 = xm1
        self.xm2 = xm2
        self.njrc = njrc
        self.vi = vi


class IntegContext(object):
    """
    Context container for ``integ`` radial integration.

    Parameters
    ----------
    z : float
        Effective nuclear charge.
    v : list of float
        Effective potential array.
    xm1, xm2 : list of float
        Auxiliary Numerov arrays.
    nr : int
        Number of radial grid points.
    r, r2 : list of float
        Radial grid and squared radius arrays.
    dl : float
        Logarithmic grid spacing.
    rel : float
        Relativistic flag value.

    References
    ----------
    See `Numerov's method <https://en.wikipedia.org/wiki/Numerov%27s_method>`_.

    TODO
    ----
    - Consider converting this structure to a ``dataclass`` after Python 2.7
      deprecation to simplify typing and defaults.
    """

    def __init__(self, z, v, xm1, xm2, nr, r, r2, dl, rel):
        self.z = z
        self.v = v
        self.xm1 = xm1
        self.xm2 = xm2
        self.nr = nr
        self.r = r
        self.r2 = r2
        self.dl = dl
        self.rel = rel


def _setqmm_fill_v_from_orb(orb, i, nr, r, zeff, v):
    """
    Fill the potential array from orbital data for a single channel.

    Parameters
    ----------
    orb : list of list of float
        Orbital data by grid and orbital index.
    i : int
        Orbital index.
    nr : int
        Number of radial grid points.
    r : list of float
        Radial grid.
    zeff : float
        Effective nuclear charge.
    v : list of float
        Potential array to update.

    Returns
    -------
    None
        Updates ``v`` in-place.

    References
    ----------
    See `Hartree-Fock method <https://en.wikipedia.org/wiki/Hartree%E2%80%93Fock_method>`_.

    TODO
    ----
    - Replace the loop with NumPy vectorized operations once available.
    """
    for j in range(1, nr + 1):
        v[j] = -zeff / r[j] + orb[j][i]


def _setqmm_fill_v_from_vi(vi, orb, lp2, i, nr, v):
    """
    Fill the potential array using pre-tabulated ``vi`` values.

    Parameters
    ----------
    vi : list of list of float
        Potential table by grid and channel.
    orb : list of list of float
        Orbital data by grid and orbital index.
    lp2 : int
        Angular channel index.
    i : int
        Orbital index.
    nr : int
        Number of radial grid points.
    v : list of float
        Potential array to update.

    Returns
    -------
    None
        Updates ``v`` in-place.

    References
    ----------
    See `Hartree-Fock method <https://en.wikipedia.org/wiki/Hartree%E2%80%93Fock_method>`_.

    TODO
    ----
    - Vectorize with NumPy for faster potential assembly.
    """
    for j in range(1, nr + 1):
        v[j] = vi[j][lp2] + orb[j][i]


def _setqmm_update_xm_from_orb(orb, i, nr, dl, r, r2, a2, za2, zaa, xm1, xm2):
    """
    Update ``xm1`` and ``xm2`` from orbital data via finite differences.

    Parameters
    ----------
    orb : list of list of float
        Orbital data by grid and orbital index.
    i : int
        Orbital index.
    nr : int
        Number of radial grid points.
    dl : float
        Logarithmic grid spacing.
    r, r2 : list of float
        Radial grid and squared radius arrays.
    a2, za2, zaa : float
        Relativistic coefficients for potential derivatives.
    xm1, xm2 : list of float
        Auxiliary Numerov arrays (updated in-place).

    Returns
    -------
    None

    References
    ----------
    See `Hartree-Fock method <https://en.wikipedia.org/wiki/Hartree%E2%80%93Fock_method>`_.

    TODO
    ----
    - Consider NumPy-based finite differences for improved readability.
    """
    for j in range(2, nr - 1):
        dvdl = (orb[j + 1][i] - orb[j - 1][i]) / (2.0 * dl)
        ddvdrr = ((orb[j + 1][i] + orb[j - 1][i] - 2.0 * orb[j][i]) / (dl * dl) - dvdl) / r2[j]
        xm1[j] = -a2 * dvdl / r[j] - za2 / r2[j]
        xm2[j] = -a2 * ddvdrr + zaa / r2[j] / r[j]
    xm1[nr] = xm1[nr - 1]
    xm2[nr] = xm2[nr - 1]
    xm1[1] = xm1[2] + za2 / r2[2] - za2 / r2[1]
    xm2[1] = xm2[2] - zaa / r2[2] / r[2] + zaa / r2[1] / r[1]


def _setqmm_update_xm_from_v(nr, dl, r, r2, a2, v, xm1, xm2):
    """
    Update ``xm1`` and ``xm2`` from the potential array.

    Parameters
    ----------
    nr : int
        Number of radial grid points.
    dl : float
        Logarithmic grid spacing.
    r, r2 : list of float
        Radial grid and squared radius arrays.
    a2 : float
        Relativistic coefficient.
    v : list of float
        Potential array.
    xm1, xm2 : list of float
        Auxiliary Numerov arrays (updated in-place).

    Returns
    -------
    None

    References
    ----------
    See `Hartree-Fock method <https://en.wikipedia.org/wiki/Hartree%E2%80%93Fock_method>`_.

    TODO
    ----
    - Replace explicit loops with NumPy gradients where possible.
    """
    for j in range(2, nr):
        dvdl = (v[j + 1] - v[j - 1]) / (2.0 * dl)
        ddvdrr = ((v[j + 1] + v[j - 1] - 2.0 * v[j]) / (dl * dl) - dvdl) / r2[j]
        xm1[j] = -a2 * dvdl / r[j]
        xm2[j] = -a2 * ddvdrr
    xm1[nr] = xm1[nr - 1]
    xm2[nr] = xm2[nr - 1]
    xm1[1] = xm1[2]
    xm2[1] = xm2[2]


def setqmm(i, orb, l, idoflag, ctx):  # noqa: E741
    """
    Construct effective potentials and Numerov coefficients for an orbital.

    Parameters
    ----------
    i : int
        Orbital index.
    orb : list of list of float
        Orbital data by grid and orbital index.
    l : int
        Angular momentum quantum number.
    idoflag : int
        Update flag controlling whether potentials are rebuilt.
    ctx : SetqmmContext
        Context holding potentials, grid, and auxiliary arrays.

    Returns
    -------
    None
        Updates context fields in-place (``v``, ``xm1``, ``xm2``, ``zeff``).

    Notes
    -----
    The routine assembles the effective potential and Numerov coefficients,
    optionally incorporating tabulated ``vi`` values when provided.

    References
    ----------
    See `Numerov's method <https://en.wikipedia.org/wiki/Numerov%27s_method>`_.

    TODO
    ----
    - Convert the internal loops to NumPy vector operations to reduce
      per-orbital overhead.
    """
    # ns parameter removed; it was unused in the original signature.
    c = 137.038
    alpha = ctx.rel / c
    aa = alpha * alpha
    a2 = aa / 2.0

    lp = l + 1
    lpx = lp
    if lp > 4:
        lpx = 4
    lp2 = l + l + 1
    if lp2 > 7:
        lp2 = 7
    zeff_value = ctx.zorig
    if ctx.njrc[lpx] > 0:
        zeff_value = 0.0
    zaa = zeff_value * aa
    za2 = zeff_value * a2
    v = ctx.v

    if idoflag:
        if not ctx.njrc[lpx]:
            if idoflag == 1:
                _setqmm_fill_v_from_orb(orb, i, ctx.nr, ctx.r, zeff_value, v)
            _setqmm_update_xm_from_orb(
                orb,
                i,
                ctx.nr,
                ctx.dl,
                ctx.r,
                ctx.r2,
                a2,
                za2,
                zaa,
                ctx.xm1,
                ctx.xm2,
            )
    else:
        if idoflag == 1 and ctx.vi is not None:
            _setqmm_fill_v_from_vi(ctx.vi, orb, lp2, i, ctx.nr, v)
        _setqmm_update_xm_from_v(ctx.nr, ctx.dl, ctx.r, ctx.r2, a2, v, ctx.xm1, ctx.xm2)

    # figure out the (Desclaux-Numerov) effective potential.
    xlb = l + pow(0.5, 2.0) / 2.0
    for j in range(1, ctx.nr + 1):
        vj = v[j]
        ctx.q0[j] = vj * (1.0 - a2 * vj) + xlb / ctx.r2[j]

    ctx.zeff = zeff_value


def initiali(
    zorig,
    nr,
    rmin,
    rmax,
    r,
    dr,
    r2,
    dl,
    njrc=[None] * 4,
    xntot=0.0,
    nel=0,
    input_stream="stdin",
):
    """
    Initialize the radial grid and core-radius indices.

    Parameters
    ----------
    zorig : float
        Atomic number.
    nr : int
        Number of radial grid points.
    rmin, rmax : float
        Minimum and maximum radial bounds.
    r, dr, r2 : list of float
        Radial grid, spacing, and squared radius arrays (updated in-place).
    dl : float
        Logarithmic grid spacing (updated in-place).
    njrc : list of int, optional
        Core radius indices per angular momentum channel.
    xntot : float, optional
        Total electron count accumulator.
    nel : int, optional
        Number of electrons.
    input_stream : str
        Input source; ``"stdin"`` or a file handle.

    Returns
    -------
    tuple
        ``(zorig, nr, rmin, rmax, r, dr, r2, dl, njrc, xntot, nel)``.

    Notes
    -----
    The grid spacing is computed by ``setgrid`` using a logarithmic mesh.

    References
    ----------
    See `Hartree-Fock method <https://en.wikipedia.org/wiki/Hartree%E2%80%93Fock_method>`_.

    TODO
    ----
    - Replace the list-based grid with NumPy arrays for faster allocation.
    """

    if input_stream == "stdin":
        (zorig, nr) = [t(s) for t, s in zip((float, int), get_input("Enter Z, NR: ").split())]

    elif isinstance(input_stream, file):
        (zorig, nr) = [t(s) for t, s in zip((float, int), input_stream.next().split("!")[0].split())]

    else:
        raise IOError("input stream is not a file handle or 'stdin'")

    rmin = 0.0001 / zorig
    rmax = 800.0 / sqrt(zorig)

    nr, rmin, rmax, r, dr, r2, dl = setgrid(nr, rmin, rmax, r, dr, r2, dl)
    for idx, _ in enumerate(njrc):
        njrc[idx] = 0

    return (zorig, nr, rmin, rmax, r, dr, r2, dl, njrc, xntot, nel)


def setgrid(nr, rmin, rmax, r, dr, r2, dl):
    """
    Construct a logarithmic radial grid and its derived arrays.

    Parameters
    ----------
    nr : int
        Number of radial grid points.
    rmin, rmax : float
        Minimum and maximum radii.
    r, dr, r2 : list of float
        Radial grid, spacing, and squared radius arrays (updated in-place).
    dl : float
        Logarithmic grid spacing (updated in-place).

    Returns
    -------
    tuple
        ``(nr, rmin, rmax, r, dr, r2, dl)``.

    Notes
    -----
    The grid uses the logarithmic spacing described in the module docstring.

    References
    ----------
    See `Hartree-Fock method <https://en.wikipedia.org/wiki/Hartree%E2%80%93Fock_method>`_.

    TODO
    ----
    - Replace the loop with NumPy vectorized expressions for ``r`` and ``dr``.
    """
    ratio = rmax / rmin
    dl = log(ratio) / float(nr)
    xratio = exp(dl)
    xr1 = sqrt(xratio) - sqrt(1.0 / xratio)
    for i in range(len(r)):
        r[i] = pow(rmin * xratio, i)
        dr[i] = r[i] * xr1
        r2[i] = r[i] * r[i]

    return (nr, rmin, rmax, r, dr, r2, dl)


def integ(e, l, xkappa, n, istop, phi, ctx):  # noqa: E741, C901, MC0001
    """
    Integrate the radial equation and count nodes.

    Parameters
    ----------
    e : float
        Trial eigenvalue.
    l : int
        Angular momentum quantum number.
    xkappa : float
        Relativistic kappa parameter.
    n : int
        Principal quantum number.
    istop : int
        Turning-point index; 0 triggers automatic detection.
    phi : list of float
        Radial wavefunction buffer (updated in-place).
    ctx : IntegContext
        Context containing grid, potential, and auxiliary arrays.

    Returns
    -------
    tuple
        ``(node_count, istop, energy_flag, x0)`` where ``x0`` is the
        log-derivative at the turning point when available.

    Notes
    -----
    The Numerov-style scheme integrates outward and counts nodes to detect
    whether the trial energy is too high or too low.

    References
    ----------
    See `Numerov's method <https://en.wikipedia.org/wiki/Numerov%27s_method>`_.

    TODO
    ----
    - Consider an explicit ODE solver (SciPy) for validation of edge cases.
    """
    # x0 is returned as the log-derivative at the turning point when available.
    dl = ctx.dl
    rel = ctx.rel
    z = ctx.z
    v = ctx.v
    xm1 = ctx.xm1
    xm2 = ctx.xm2
    nr = ctx.nr
    r = ctx.r
    r2 = ctx.r2

    dl2 = dl * dl / 12.0
    dl5 = 10.0 * dl2
    c = 137.038
    alpha = rel / c
    za2 = z * z * alpha * alpha
    a2 = alpha * alpha / 2.0
    xl = l
    xlp = l + 1
    xl2 = 0.5 + xl
    xl4 = xl2 * xl2

    # we set up the leading power.
    # adjust for Desclaux's implementation of Numerov.
    if rel == 0.0:
        ss = xlp
    else:
        rtest = 1.0 - za2
        if rtest < 0.0:
            print("Z>137 is too big.")
            sys.exit(1)

        ss = sqrt(rtest)

    ss2 = ss - 0.5

    # we shall set ief to -1 if energy is too low, +1 if too high.
    energy_flag = 0
    node_count = 0
    x0 = None

    # see Desclaux and documentation to see the origin of the below equations.
    # here, we set up the first two points.
    t = e - v[1]
    xm0 = 1.0 + a2 * t
    tm = xm0 + xm0
    xmx = xm1[1] / xm0
    xk0 = r2[1] * (tm * t - xmx * (xkappa / r[1] + 0.75 * xmx) + xm2[1] / tm) - xl4
    dk0 = 1.0 + dl2 * xk0
    p0 = dk0
    phi[1] = p0 * sqrt(xm0 * r[1]) / dk0

    t = e - v[2]
    xm = 1.0 + a2 * t
    tm = xm + xm
    xmx = xm1[2] / xm
    xk2 = r2[2] * (tm * t - xmx * (xkappa / r[2] + 0.75 * xmx) + xm2[2] / tm) - xl4
    dk2 = 1.0 + dl2 * xk2
    p1 = dk2 * (pow(r[2] / r[1], ss2) - (r[2] - r[1]) * z / xlp) * sqrt(xm0 / xm)
    phi[2] = p1 * sqrt(xm * r[2]) / dk2

    # if istop is set, the we know to stop there.  If it is zero, it shall
    # be set to the classical turning point.
    is0 = istop
    if not istop:
        for j in range(nr - 1, 2 - 1, -1):
            if e > v[j]:
                break
            energy_flag = -1
            return (node_count, istop, energy_flag, x0)
        istop = j

    # initialize number of nodes, and determine the ideal number.
    nnideal = n - l - 1

    # integrate out count nodes, and stop along the way if there are too many
    for i in range(3, istop + 2 + 1):
        t = e - v[i]
        xm = 1.0 + a2 * t
        tm = xm + xm
        xmx = xm1[i] / xm
        p2 = (2.0 - dl5 * xk2) * p1 / dk2 - p0
        xk2 = r2[i] * (tm * t - xmx * (xkappa / r[i] + 0.75 * xmx) + xm2[i] / tm) - xl4
        dk2 = 1.0 + dl2 * xk2
        phi[i] = p2 * sqrt(xm * r[i]) / dk2
        if abs(p2) > 10000000000.0:
            for j in range(1, i + 1):
                phi[j] /= p2

            p0 /= p2
            p1 /= p2
            p2 /= p2

        if p2 * p1 < 0.0:
            node_count += 1
            if node_count > nnideal:
                energy_flag = 1
                return (node_count, istop, energy_flag, x0)

        p0 = p1
        p1 = p2

    if istop > 0:
        psip2 = phi[istop + 2] - phi[istop - 2]
        psip1 = phi[istop + 1] - phi[istop - 1]
        psip = (8.0 * psip1 - psip2) / (12.0 * dl * r[istop])
        x0 = psip / phi[istop]

    if is0:
        return (node_count, istop, energy_flag, x0)

    for i in range(istop + 3, nr - 1 + 1):
        t = e - v[i]
        xm = 1.0 + a2 * t
        tm = xm + xm
        xmx = xm1[i] / xm
        p2 = (2.0 - dl5 * xk2) * p1 / dk2 - p0
        if p2 / p1 > 1.0:
            energy_flag = -1
            return (node_count, istop, energy_flag, x0)

        xk2 = r2[i] * (tm * t - xmx * (xkappa / r[i] + 0.75 * xmx) + xm2[i] / tm) - xl4
        dk2 = 1.0 + dl2 * xk2
        phi[i] = p2 * sqrt(xm * r[i]) / dk2
        if abs(p2) > 10000000000.0:
            for j in range(1, i + 1):
                phi[j] /= p2

                p0 = p2
                p1 /= p2
                p2 /= p2

                if p2 * p1 < 0.0:
                    node_count += 1
                    if node_count > nnideal:
                        energy_flag = 1
                        return (node_count, istop, energy_flag, x0)

            p0 = p1
            p1 = p2

    return (node_count, istop, energy_flag, x0)


def clebschgordan(nel, nl, cg, si, fa):
    """
    Compute Clebsch-Gordan coefficients for angular momentum coupling.

    Parameters
    ----------
    nel : int
        Number of electron orbitals.
    nl : list of int
        Angular quantum numbers for each orbital.
    cg : list
        Preallocated coefficient tensor to update in-place.
    si : list of float
        Sign array (updated in-place).
    fa : list of float
        Factorial array (updated in-place).

    Returns
    -------
    None
        Updates ``cg``, ``si``, and ``fa`` in-place.

    Notes
    -----
    The implementation follows Wigner's formula as described in Rose,
    "Elementary Theory of Angular Momentum".

    References
    ----------
    .. [1] Clebsch-Gordan coefficients,
       https://en.wikipedia.org/wiki/Clebsch%E2%80%93Gordan_coefficients
    .. [2] M. E. Rose, *Elementary Theory of Angular Momentum*, p. 39.

    TODO
    ----
    - Replace this implementation with ``scipy.special.wigner_3j`` and related
      routines once SciPy is a supported dependency.
    """
    lmx = 0
    for i in range(len(nl)):
        if nl[i] > lmx:
            lmx = nl[i]

    si[0] = fa[0] = 1.0

    for i in range(1, len(si)):
        si[i] = -si[i - 1]
        fa[i] = i * fa[i - 1]

    for l1 in range(0, lmx + 1):
        for l2 in range(0, l1 + 1):
            # 52     format (1x,i3,a3,i3)
            for m1 in range(-l1, l1 + 1):
                for m2 in range(-l2, l2 + 1):
                    m3 = m1 + m2
                    lmin = abs(l1 - l2)
                    if lmin < abs(m3):
                        lmin = abs(m3)

                    for l3 in range(lmin, l1 + l2 + 1):
                        prefactor = float(2 * l3 + 1)
                        prefactor *= fa[l3 + l1 - l2] / fa[l1 + l2 + l3 + 1]
                        prefactor *= fa[l3 - l1 + l2] / fa[l1 - m1]
                        prefactor *= fa[l1 + l2 - l3] / fa[l1 + m1]
                        prefactor *= fa[l3 + m3] / fa[l2 - m2]
                        prefactor *= fa[l3 - m3] / fa[l2 + m2]
                        prefactor = sqrt(prefactor)
                        sum1 = 0.0
                        numax = l3 - l1 + l2
                        if l3 + m3 < numax:
                            numax = l3 + m3
                        numin = 0
                        if l1 - l2 - m3 < numin:
                            numin = -l1 - l2 - m3
                        for nu in range(numin, numax + 1):
                            sum1 += (
                                (si[nu + l2 + m2] / fa[nu])
                                * fa[l2 + l3 + m1 - nu]
                                * fa[l1 - m1 + nu]
                                / fa[l3 - l1 + l2 - nu]
                                / fa[l3 + m3 - nu]
                                / fa[nu + l1 - l2 - m3]
                            )

                        cg[l1][l2][l3][m1][m2] = prefactor * sum1
                        cg[l2][l1][l3][m2][m1] = si[l1 + l2 + l3] * prefactor * sum1

    return


def pseudo(
    etot,
    nst,
    rel,
    alfa,
    nr,
    rmin,
    rmax,
    r,
    dr,
    r2,
    dl,
    phe,
    orb,
    njrc,
    vi,
    zorig,
    xntot,
    nel,
    no,
    nl,
    xnj,
    ev,
    occ,
    iss,
    ek,
    iuflag,
    vctab,
    nm,
    input_stream="stdin",
):
    """
    Generate pseudopotentials using the current SCF solution.

    Parameters
    ----------
    etot : float
        Total energy accumulator.
    nst : int
        Number of shells to process.
    rel : int
        Relativistic flag.
    alfa : float
        Exchange-correlation mixing parameter.
    nr : int
        Number of radial grid points.
    rmin, rmax : float
        Minimum and maximum radial bounds.
    r, dr, r2 : list of float
        Radial grid, spacing, and squared radius arrays.
    dl : float
        Logarithmic grid spacing.
    phe, orb : list of list of float
        Radial wavefunction and orbital buffers.
    njrc : list of int
        Core radius indices per angular momentum channel.
    vi : list of list of float
        Potential table by grid and channel.
    zorig : float
        Atomic number.
    xntot : float
        Total electron count accumulator.
    nel : int
        Number of electrons.
    no, nl : list of int
        Principal and angular quantum numbers.
    xnj : list of float
        Total angular momentum values (j).
    ev : list of float
        Orbital eigenvalues.
    occ : list of float
        Orbital occupancies.
    iss : list of int
        Spin flags.
    ek : list of float
        Kinetic energy contributions.
    iuflag : int
        Output control flag.
    vctab : list of list of float
        Core potential corrections.
    nm : list of int
        Magnetic quantum numbers.
    input_stream : str
        Input source; ``"stdin"`` or a file handle.

    Returns
    -------
    None
        Updates orbital and pseudopotential data in-place.

    Notes
    -----
    The routine calls ``pseudize`` for each eligible orbital and can be
    adapted to alternative pseudopotential generation schemes.

    References
    ----------
    See `Hartree-Fock method <https://en.wikipedia.org/wiki/Hartree%E2%80%93Fock_method>`_.

    TODO
    ----
    - Consider modularizing pseudopotential generation for alternative
      approaches (e.g., Kleinman-Bylander forms).
    """

    # initialise
    nm = [0] * nel
    njrcdummy = deepcopy(njrc)
    q0 = xm1 = xm2 = [None] * len(r)
    rpower = [[None] * 7] * len(r)
    zeff = etot2 = float

    # read input
    if input_stream == "stdin":
        (np, corpol, rnorm) = [
            t(s) for t, s in zip((int, float, float), get_input("Please enter NP CORPOL RNORM: ").split())
        ]

    elif isinstance(input_stream, file):
        (np, corpol, rnorm) = [t(s) for t, s in zip((int, float, float), input_stream.readline().split("!")[0].split())]
    else:
        raise IOError("input_stream is not valid!")

    xntot = 0.0

    while True:
        for i in range(np, nel + 1):
            print("l={} ...".format(nl[i]))
            lp2 = nl[i] + nl[i] + 1
            e = ev[i]
            for j in range(1, nr + 1):
                orb[j][i] += vctab[j][nl[i]]
            idoflag = 1
            ctx = SetqmmContext(vi[1][lp2], zeff, zorig, rel, nr, r, r2, dl, q0, xm1, xm2, njrcdummy, vi)
            setqmm(i, orb, nl[i], idoflag, ctx)
            zeff = ctx.zeff
            for j in range(1, nr + 1):
                orb[j][i] = 0.0

            # you can replace the pseudize subroutine with any type of PP
            # generation you want, however, Kleinman-Bylanderization would
            # take more coding...
            pseudize(
                i,
                orb,
                e,
                nl[i],
                xnj[i],
                no[i],
                njrc,
                zeff,
                vi[1][lp2],
                q0,
                xm1,
                xm2,
                nr,
                rmin,
                rmax,
                r,
                dr,
                r2,
                dl,
                rel,
            )
            print("Doing pseudo PP generation...")
            no[i] = nl[i] + 1
            ruse = 0.0
            xkappa = -1.0
            elsolve(
                i,
                occ[i],
                no[i],
                nl[i],
                xkappa,
                xnj[i],
                zorig,
                zeff,
                ev[i],
                phe[1][i],
                vi[1][lp2],
                q0,
                xm1,
                xm2,
                nr,
                r,
                dr,
                r2,
                dl,
                ruse,
            )
            print(nl[i], ev[i])
            xntot += occ[i]
            if lp2 != 1:
                for j in range(1, nr + 1):
                    vi[j][lp2 - 1] = vi[j][lp2]
                break

        print("everything is pseudized")
        for i in range(np, nel + 1):
            inew = 1 + i - np
            no[inew] = no[i]
            nl[inew] = nl[i]
            nm[inew] = nm[i]
            xnj[inew] = xnj[i]
            iss[inew] = 1
            ev[inew] = ev[i]
            occ[inew] = occ[i]
            for j in range(1, nr + 1):
                phe[j][inew] = phe[j][i]

        nel += 1 - np
        for i in range(0, 7 + 1):
            xi = i
            for k in range(1, nr + 1):
                rpower[k][i] = pow(r[k], xi)

        print("everything is scaled down...ready for unscreening")
        xnum = 100.0
        ratio = 1.0
        getpot(
            etot,
            nst,
            rel,
            alfa,
            dl,
            nr,
            dr,
            r,
            r2,
            xntot,
            phe,
            ratio,
            orb,
            occ,
            iss,
            nel,
            nl,
            nm,
            no,
            xnj,
            rpower,
            xnum,
            etot2,
            iuflag,
        )
        print("screening effects in pseudo atom computed...")
        for k in range(1, nel + 1):
            lp2 = nl[k] + nl[k] + 1
            for j in range(1, nr + 1):
                vi[j][lp2] -= orb[j][k]
                if lp2 > 1:
                    vi[j][lp2 - 1] -= orb[j][k]

        print("we got past the unscreening...")
        for j in range(1, nr + 1):
            vl = (vi[j][2] + 2.0 * vi[j][3]) / 3.0
            vso = 2.0 * (vi[j][3] - vi[j][2]) / 3.0
            vi[j][2] = vso
            vi[j][3] = vl
            vl = (2.0 * vi[j][4] + 3.0 * vi[j][5]) / 5.0
            vso = 2.0 * (vi[j][5] - vi[j][4]) / 5.0
            vi[j][4] = vso
            vi[j][5] = vl
            # 2222   format (5f8.4)
            vl = (3.0 * vi[j][6] + 4.0 * vi[j][7]) / 7.0
            vso = 2.0 * (vi[j][7] - vi[j][6]) / 7.0
            vi[j][6] = vso
            vi[j][7] = vl

        rel = 0.0
        print("we got past the spin-orbit jazz")
        izuse = abs(vi[nr - 2][1] * r[nr - 2]) + 0.5
        zuse = izuse

        # check block
        for k in range(1, nr + 1):
            if r[k] > rnorm:
                videal = -zuse / r[k] - corpol / (2.0 * pow(r[k], 4.0))
                vi[k][1] = videal
                vi[k][3] = videal
                vi[k][5] = videal
                vi[k][7] = videal
                vi[k][2] = 0.0
                vi[k][4] = 0.0
                vi[k][6] = 0.0

        print("we got to the end")
        return


def parabreg(f, fp, fpp, rf, vf):
    """
    Fit a parabola through three points to approximate derivatives.

    Parameters
    ----------
    f, fp, fpp : float
        Function value, first derivative, and second derivative estimates.
    rf : list of float
        Radius samples.
    vf : list of float
        Function samples at ``rf``.

    Returns
    -------
    tuple
        ``(f, fp, fpp, rf, vf)`` updated with the fitted values.

    Notes
    -----
    This helper is a local quadratic regression used in pseudization.

    References
    ----------
    See `Hartree-Fock method <https://en.wikipedia.org/wiki/Hartree%E2%80%93Fock_method>`_.

    TODO
    ----
    - Replace with NumPy polynomial fitting or SciPy interpolation routines.
    """
    f = vf[2]
    r21 = rf[2] - rf[1]
    r32 = rf[3] - rf[2]
    v21 = vf[2] - vf[1]
    v32 = vf[3] - vf[2]
    fp = (v21 + v32) / (r21 + r32)
    fpp = (v32 / r32 - v21 / r21) / ((r21 + r32) / 2.0)
    return f, fp, fpp, rf, vf


def hb(x, factor):
    """
    Evaluate the smooth cutoff function used in pseudization.

    Parameters
    ----------
    x : float
        Normalized radius (``r / r_cut``).
    factor : float
        Smoothing factor controlling the cutoff sharpness.

    Returns
    -------
    float
        Smooth cutoff value in ``[0, 1]``.

    Notes
    -----
    The function is based on a hyperbolic-sine shaping factor.

    References
    ----------
    See `Hartree-Fock method <https://en.wikipedia.org/wiki/Hartree%E2%80%93Fock_method>`_.

    TODO
    ----
    - Vectorize via NumPy for bulk evaluations.
    """
    if x > 3.0:
        hb = 0
    if x <= 3.0:
        hb = pow(0.01, pow((sinh(x / factor) / 1.1752), 2.0))
    return hb


def fitx0(
    i,
    orb,
    rcut,
    njrc,
    e,
    l,
    xj,
    n,
    jrt,
    xideal,
    phi,
    zeff,
    v,
    q0,
    xm1,
    xm2,
    nr,
    r,
    dr,
    r2,
    dl,
    rel,
    factor,
):
    """
    Fit the log-derivative target for a pseudized orbital.

    Parameters
    ----------
    i : int
        Orbital index.
    orb : list of list of float
        Orbital data by grid and orbital index.
    rcut : float
        Cutoff radius for pseudization.
    njrc : list of int
        Core radius indices per angular momentum channel.
    e : float
        Trial eigenvalue.
    l : int
        Angular momentum quantum number.
    xj : float
        Total angular momentum value.
    n : int
        Principal quantum number.
    jrt : int
        Index of the matching radius.
    xideal : float
        Target log-derivative value.
    phi : list of float
        Radial wavefunction (updated in-place).
    zeff : float
        Effective nuclear charge.
    v : list of float
        Effective potential array.
    q0, xm1, xm2 : list of float
        Auxiliary Numerov arrays.
    nr : int
        Number of radial grid points.
    r, dr, r2 : list of float
        Radial grid, spacing, and squared radius arrays.
    dl : float
        Logarithmic grid spacing.
    rel : float
        Relativistic flag value.
    factor : float
        Smoothing factor used by ``hb``.

    Returns
    -------
    None
        Updates ``v`` and ``phi`` in-place.

    Notes
    -----
    This routine iteratively adjusts the local potential to match the target
    log-derivative at the matching radius.

    References
    ----------
    See `Hartree-Fock method <https://en.wikipedia.org/wiki/Hartree%E2%80%93Fock_method>`_.

    TODO
    ----
    - Replace the manual update with ``scipy.optimize`` root finding for
      clarity and error control.
    """

    vl = -1000000.0
    vh = 1000000.0
    dummy = nn = xactual = None  # initialise
    while True:
        idoflag = 2  # label 115
        xkappa = -1.0

        ctx = SetqmmContext(v, zeff, dummy, rel, nr, r, r2, dl, q0, xm1, xm2, njrc)
        setqmm(i, orb, l, idoflag, ctx)
        zeff = ctx.zeff

        integ_ctx = IntegContext(zeff, v, xm1, xm2, nr, r, r2, dl, rel)
        nn, jrt, _unused_ief, xactual = integ(e, l, xkappa, n, jrt, phi, integ_ctx)

        if int(nn):
            vl = v[1]
            xla = 1.0
        else:
            if xactual > xideal:
                vh = v[1]
            else:
                vl = v[1]

        xerror = xideal - xactual
        if abs(xerror) < 0.000000001:
            return
        dxdla = 0.0
        for ii in range(1, jrt + 1):
            dxdla += dr[ii] * phi[ii] * phi[ii] * hb(r[ii] / rcut, factor)
        dxdla = 2.0 * dxdla / (phi[jrt] * phi[jrt])
        xla = xerror / dxdla

        vmaybe = v[1] + xla
        if vmaybe > vh or vmaybe < vl:
            xla = (vl + vh) / (2.0 - v[1])
    for ii in range(1, jrt - 1 + 1):
        v[ii] += xla * hb(r[ii] / rcut, factor)


def _pseudize_prompt_cutoff(rcut, factor):
    """
    Get cutoff radius and smoothing factor from input if missing.

    Parameters
    ----------
    rcut : float or None
        Cutoff radius; when ``None`` prompts for input.
    factor : float or None
        Smoothing factor; when ``None`` prompts for input.

    Returns
    -------
    tuple
        ``(rcut, factor)`` populated with floats.
    """
    if rcut is None or factor is None:
        while True:
            try:
                values = get_input("Please enter the cutoff radius, and factor: ").split()[:2]
                rcut = float(values[0])
                factor = float(values[1])
                break
            except (ValueError, IndexError):
                print("Invalid input - please retry...")
    return rcut, factor


def _pseudize_select_cutoff(rcut, phi, r, istop, n, l):
    """
    Select a cutoff radius based on node position when requested.

    Parameters
    ----------
    rcut : float
        Input cutoff value; negative values encode node fraction selection.
    phi : list of float
        Radial wavefunction values.
    r : list of float
        Radial grid.
    istop : int
        Index used for the node search.
    n, l : int
        Principal and angular momentum quantum numbers.

    Returns
    -------
    float
        Updated cutoff radius.
    """
    if rcut < 0.0:
        xnodefrac = -rcut
        j = istop
        while phi[j - 1] / phi[j] <= 1.0:
            j -= 1
        if n > l + 1:
            k = j
        while phi[k - 1] / phi[k] <= 1.0:
            k -= 1
        rcut = r[k] + xnodefrac * (r[j] - r[k])
    return rcut


def _pseudize_indices(rcut, rmin, rmax, nr, r, njrc, lp):
    """
    Compute matching indices for the cutoff and test radii.

    Parameters
    ----------
    rcut : float
        Cutoff radius.
    rmin, rmax : float
        Radial grid bounds.
    nr : int
        Number of grid points.
    r : list of float
        Radial grid.
    njrc : list of int
        Core radius indices per angular momentum channel.
    lp : int
        Angular momentum index in ``njrc`` to update.

    Returns
    -------
    tuple
        ``(rcut, jrc, jrt, rtest)`` updated cutoff, indices, and test radius.
    """
    jrc = 1.0 + float(nr - 1) * log(rcut / rmin) / log(rmax / rmin)
    rcut = r[jrc]
    rtest = 2.0 * rcut
    jrt = 1.0 + float(nr - 1) * log(rtest / rmin) / log(rmax / rmin)
    njrc[lp] = jrt
    return rcut, jrc, jrt, r[jrt]


def _pseudize_normalize_phi(phi, jrt):
    """
    Normalize ``phi`` so that ``phi[jrt]`` equals one.

    Parameters
    ----------
    phi : list of float
        Radial wavefunction (updated in-place).
    jrt : int
        Normalization index.
    """
    for ii in range(len(phi)):
        phi[ii] /= phi[jrt]


def _pseudize_integral_norm(phi, dr, jrt):
    """
    Compute the trapezoidal norm integral up to ``jrt``.

    Parameters
    ----------
    phi : list of float
        Radial wavefunction values.
    dr : list of float
        Radial spacing values.
    jrt : int
        Matching index.

    Returns
    -------
    float
        Approximate norm integral up to ``jrt``.
    """
    total = 0.0
    for ii in range(1, jrt - 1 + 1):
        total += dr[ii] * phi[ii] * phi[ii]
    total += dr[jrt] * phi[jrt] * phi[jrt] / 2.0
    return total


def _pseudize_reference_values(ev, l, xkappa, n, jrt, phi, integ_ctx, dr):
    """
    Compute reference log-derivative and norm for the all-electron orbital.

    Parameters
    ----------
    ev : float
        Eigenvalue.
    l : int
        Angular momentum quantum number.
    xkappa : float
        Relativistic kappa parameter.
    n : int
        Principal quantum number.
    jrt : int
        Matching index.
    phi : list of float
        Radial wavefunction (updated in-place).
    integ_ctx : IntegContext
        Integration context for ``integ``.
    dr : list of float
        Radial spacing values.

    Returns
    -------
    tuple
        ``(x00, xn00, c00)`` reference log-derivative, norm, and slope.
    """
    x00 = integ(ev, l, xkappa, n, jrt, phi, integ_ctx)[3]
    _pseudize_normalize_phi(phi, jrt)
    xn00 = _pseudize_integral_norm(phi, dr, jrt)
    de = 0.0001
    xp = integ(ev + de / 2.0, l, xkappa, n, jrt, phi, integ_ctx)[3]
    xm = integ(ev - de / 2.0, l, xkappa, n, jrt, phi, integ_ctx)[3]
    c00 = (xm - xp) / (2.0 * de)
    return x00, xn00, c00


def _pseudize_apply_polynomial_potential(v, r, r2, dl, rcut, jrc):
    """
    Replace the core potential with a smooth polynomial inside ``rcut``.

    Parameters
    ----------
    v : list of float
        Potential array (updated in-place).
    r, r2 : list of float
        Radial grid and squared radii.
    dl : float
        Logarithmic grid spacing.
    rcut : float
        Cutoff radius.
    jrc : int
        Cutoff index on the radial grid.
    """
    v0 = v[jrc]
    dvdl = (8.0 * (v[jrc + 1] - v[jrc - 1]) - (v[jrc + 2] - v[jrc - 2])) / (12.0 * dl)
    ddvdll = (16.0 * (v[jrc + 1] + v[jrc - 1]) - 30.0 * v[jrc] - v[jrc + 2] - v[jrc - 2]) / (12.0 * dl * dl)
    dldr = 1.0 / r[jrc]
    ddldrr = -1.0 / r2[jrc]
    v1 = dvdl * dldr
    v2 = dvdl * ddldrr + ddvdll * dldr * dldr
    b4 = (v2 * rcut - v1) / (8.0 * pow(rcut, 3.0))
    b2 = (v1 - 4.0 * b4 * pow(rcut, 3.0)) / (2.0 * rcut)
    b0 = v0 - b4 * pow(rcut, 4.0) - b2 * pow(rcut, 2.0)
    for ii in range(1, jrc + 1):
        rr = r[ii]
        v[ii] = b0 + b2 * pow(rr, 2.0) + b4 * pow(rr, 4.0)


def _pseudize_compute_xi(phi0, dr, r, rcut, factor, jrt):
    """
    Compute normalization integrals used for the delta-l update.

    Parameters
    ----------
    phi0 : list of float
        Reference radial wavefunction.
    dr : list of float
        Radial spacing values.
    r : list of float
        Radial grid.
    rcut : float
        Cutoff radius.
    factor : float
        Smoothing factor.
    jrt : int
        Matching index.

    Returns
    -------
    tuple
        ``(xi0, xi1, xi2)`` integral coefficients.
    """
    xi0 = xi1 = xi2 = 0.0
    for ii in range(1, jrt + 1):
        f = hb(r[ii] / rcut, factor)
        ph2 = dr[ii] * phi0[ii] * phi0[ii]
        xi0 += ph2
        xi1 += ph2 * f
        xi2 += ph2 * f * f
    ph2 = phi0[jrt] * phi0[jrt]
    return xi0 / ph2, xi1 / ph2, xi2 / ph2


def _pseudize_compute_deltal(xi0, xi1, xi2, c00):
    """
    Compute the delta-l correction for the norm-conserving update.

    Parameters
    ----------
    xi0, xi1, xi2 : float
        Integral coefficients from ``_pseudize_compute_xi``.
    c00 : float
        Target slope value from the all-electron solution.

    Returns
    -------
    float
        Delta-l update value.
    """
    quant = xi1 * xi1 + xi2 * (c00 - xi0)
    if quant > 0.0:
        return (sqrt(xi1 * xi1 + xi2 * (c00 - xi0)) - xi1) / xi2
    return (c00 - xi0) / (2.0 * xi1)


def _pseudize_update_phi(phi, phi0, yl, r, rcut, factor, deltal, jrt):
    """
    Update the pseudo-orbital using the current delta-l correction.

    Parameters
    ----------
    phi : list of float
        Pseudo-orbital values (updated in-place).
    phi0 : list of float
        Reference orbital values.
    yl : list of float
        Work buffer for the cutoff function.
    r : list of float
        Radial grid.
    rcut : float
        Cutoff radius.
    factor : float
        Smoothing factor.
    deltal : float
        Delta-l correction.
    jrt : int
        Matching index.
    """
    for ii in range(1, jrt):
        yl[ii] = phi0[ii] * hb(r[ii] / rcut, factor)
        phi[ii] = phi0[ii] + deltal * yl[ii]
        if phi[ii] < 0.0:
            print("Big trouble# # #  cross axis# # # ")
            sys.exit(1)


def _pseudize_update_potential(
    v,
    vraw,
    phi,
    phi0,
    yl,
    rf,
    vf,
    r,
    rcut,
    factor,
    jrt,
    fitx0_args,
    integ_args,
):
    """
    Update the local potential for the pseudo-orbital iteration.

    Parameters
    ----------
    v, vraw : list of float
        Working and reference potentials (``v`` updated in-place).
    phi, phi0 : list of float
        Current and reference orbitals.
    yl : list of float
        Cutoff buffer values.
    rf, vf : list of float
        Local interpolation buffers.
    r : list of float
        Radial grid.
    rcut : float
        Cutoff radius.
    factor : float
        Smoothing factor.
    jrt : int
        Matching index.
    fitx0_args : tuple
        Arguments for ``fitx0``.
    integ_args : tuple
        Arguments for ``integ``.

    Returns
    -------
    float
        Updated log-derivative value ``x0``.
    """
    (f, fp, fpp, psi, psip, psipp) = (None, None, None, None, None, None)
    for ii in range(1, jrt - 1 + 1):
        if phi[ii] == 0.0 or yl[ii] == 0.0:
            break
        jj = 2 if ii == 1 else ii
        for j in range(jj - 1, jj + 1 + 1):
            rf[2 + j - jj] = r[j]
            vf[2 + j - jj] = hb(r[j] / rcut, factor)
        (f, fp, fpp, rf, vf) = parabreg(f, fp, fpp, rf, vf)
        for j in range(jj - 1, jj + 1):
            vf[2 + j - jj] = phi0[j]
        (psi, psip, psipp, rf, vf) = parabreg(psi, psip, psipp, rf, vf)
        v[ii] = vraw[ii] + (1.0 - phi0[ii] / phi[ii]) * (2.0 * psip / psi * fp / f + fpp / f) / 2.0
    fitx0(*fitx0_args)
    return integ(*integ_args)[3]


def _pseudize_iteration(
    phi,
    phi0,
    v,
    vraw,
    yl,
    rf,
    vf,
    r,
    dr,
    rcut,
    factor,
    jrt,
    ev,
    l,
    xkappa,
    n,
    integ_ruse_ctx,
    fitx0_args,
    xi0,
    xi1,
    xi2,
    c00,
):
    """
    Iterate the pseudo-orbital until the log-derivative matches the target.

    Parameters
    ----------
    phi, phi0 : list of float
        Working and reference orbitals.
    v, vraw : list of float
        Working and reference potentials.
    yl, rf, vf : list of float
        Work buffers used in the update.
    r, dr : list of float
        Radial grid and spacing.
    rcut : float
        Cutoff radius.
    factor : float
        Smoothing factor.
    jrt : int
        Matching index.
    ev, l, xkappa, n : float or int
        Orbital parameters for ``integ``.
    integ_ruse_ctx : IntegContext
        Integration context with ``ruse`` enabled.
    fitx0_args : tuple
        Arguments for ``fitx0``.
    xi0, xi1, xi2 : float
        Integral coefficients.
    c00 : float
        Target slope from the all-electron solution.

    Returns
    -------
    tuple
        ``(c0, x0)`` final slope and log-derivative values.
    """
    deltal = _pseudize_compute_deltal(xi0, xi1, xi2, c00)
    print("DELTAL = %11.8f" % deltal)  # format (1x,1a9,1f11.8)
    while True:
        _pseudize_update_phi(phi, phi0, yl, r, rcut, factor, deltal, jrt)
        x0 = _pseudize_update_potential(
            v,
            vraw,
            phi,
            phi0,
            yl,
            rf,
            vf,
            r,
            rcut,
            factor,
            jrt,
            fitx0_args,
            (ev, l, xkappa, n, jrt, phi, integ_ruse_ctx),
        )
        _pseudize_normalize_phi(phi, jrt)
        xn0 = _pseudize_integral_norm(phi, dr, jrt)
        de = 0.0001
        xp = integ(ev + de / 2.0, l, xkappa, n, jrt, phi, integ_ruse_ctx)[3]
        xm = integ(ev - de / 2.0, l, xkappa, n, jrt, phi, integ_ruse_ctx)[3]
        c0 = (xm - xp) / (2.0 * de)
        print(c0, x0)
        print(xn0)
        if abs(c0 - c00) > 0.000000001:
            dqddel = 2.0 * (xi1 + deltal * xi2)
            deltal = deltal + (c00 - c0) / dqddel
        else:
            return c0, x0


def pseudize(
    i,
    orb,
    ev,
    l,
    xj,
    n,
    njrc,
    zeff,
    v,
    q0,
    xm1,
    xm2,
    nr,
    rmin,
    rmax,
    r,
    dr,
    r2,
    dl,
    rel,
    phi,
    rcut=None,
    factor=None,
):
    """
    Construct a norm-conserving pseudized orbital.

    Parameters
    ----------
    i : int
        Orbital index.
    orb : list of list of float
        Orbital data by grid and orbital index.
    ev : float
        Orbital eigenvalue.
    l : int
        Angular momentum quantum number.
    xj : float
        Total angular momentum value.
    n : int
        Principal quantum number.
    njrc : list of int
        Core radius indices per angular momentum channel.
    zeff : float
        Effective nuclear charge.
    v : list of float
        Effective potential array.
    q0, xm1, xm2 : list of float
        Auxiliary Numerov arrays.
    nr : int
        Number of radial grid points.
    rmin, rmax : float
        Minimum and maximum radial bounds.
    r, dr, r2 : list of float
        Radial grid, spacing, and squared radius arrays.
    dl : float
        Logarithmic grid spacing.
    rel : float
        Relativistic flag value.
    phi : list of float
        Radial wavefunction buffer (updated in-place).
    rcut : float, optional
        Cutoff radius; if ``None`` it is read from input.
    factor : float, optional
        Smoothing factor for the cutoff function.

    Returns
    -------
    None
        Updates ``phi`` and ``v`` in-place.

    Notes
    -----
    The method constructs a smooth, norm-conserving pseudo-orbital and iterates
    until the log-derivative matches the all-electron solution.

    References
    ----------
    See `Hartree-Fock method <https://en.wikipedia.org/wiki/Hartree%E2%80%93Fock_method>`_.

    TODO
    ----
    - Factor this routine into smaller NumPy-friendly kernels to simplify
      testing and maintenance.
    """
    x0 = x00 = None
    rf = [None] * len(njrc)
    vf = [None] * len(njrc)

    lp = l + 1
    xkappa = -1.0
    istop = nr
    while ev > q0[istop]:
        istop -= 1

    integ_ctx = IntegContext(zeff, v, xm1, xm2, nr, r, r2, dl, rel)
    istop = integ(ev, l, xkappa, n, istop, phi, integ_ctx)[1]

    rcut, factor = _pseudize_prompt_cutoff(rcut, factor)
    rcut = _pseudize_select_cutoff(rcut, phi, r, istop, n, l)
    rcut, jrc, jrt, rtest = _pseudize_indices(rcut, rmin, rmax, nr, r, njrc, lp)
    print("RCUTOFF = %8.4f  JRC = %5i" % (rcut, jrc))  # 94 format(1x,2d15.8)
    print("RTEST   = %8.4f  JRT = %5i" % (rtest, jrt))  # 1x,1a10,1f8.4,1a8,1i5
    x00, xn00, c00 = _pseudize_reference_values(ev, l, xkappa, n, jrt, phi, integ_ctx, dr)
    print(c00, x00)  # format 94
    print(xn00)  # format 94
    ruse = 0.0
    integ_ruse_ctx = IntegContext(zeff, v, xm1, xm2, nr, r, r2, dl, ruse)
    _pseudize_apply_polynomial_potential(v, r, r2, dl, rcut, jrc)
    fitx0_args = (
        i,
        orb,
        rcut,
        njrc,
        ev,
        l,
        xj,
        lp,
        jrt,
        x00,
        phi,
        zeff,
        v,
        q0,
        xm1,
        xm2,
        nr,
        r,
        dr,
        r2,
        dl,
        ruse,
        factor,
    )
    fitx0(*fitx0_args)

    phi0 = deepcopy(phi)
    vraw = deepcopy(v)

    xi0, xi1, xi2 = _pseudize_compute_xi(phi0, dr, r, rcut, factor, jrt)
    yl = [None] * len(phi0)
    c0, x0 = _pseudize_iteration(
        phi,
        phi0,
        v,
        vraw,
        yl,
        rf,
        vf,
        r,
        dr,
        rcut,
        factor,
        jrt,
        ev,
        l,
        xkappa,
        n,
        integ_ruse_ctx,
        fitx0_args,
        xi0,
        xi1,
        xi2,
        c00,
    )

    print(c0, x0)
    print("NCPP achieved # # # ")
    return


def fourier(nr, r, dr, r2, vi, output="stdout"):
    """
    Compute a radial Fourier-like transform of the potential.

    Parameters
    ----------
    nr : int
        Number of radial grid points.
    r, dr, r2 : list of float
        Radial grid, spacing, and squared radius arrays.
    vi : list of list of float
        Potential table by grid and channel.
    output : str, optional
        Output destination; ``"stdout"`` prints to console.

    Returns
    -------
    None
        Writes the transformed data to ``output`` or stdout.

    Notes
    -----
    This is a direct quadrature over the logarithmic grid rather than a fast
    Fourier transform.

    References
    ----------
    See `Hartree-Fock method <https://en.wikipedia.org/wiki/Hartree%E2%80%93Fock_method>`_.

    TODO
    ----
    - Consider NumPy vectorization for the quadrature loop.
    - Evaluate FFT-based alternatives when a uniform grid representation is
      appropriate.
    """

    a = [None] * len(r)
    v1 = [None] * len(r)

    for l in range(0, 3 + 1):
        lp2 = l + l + 1
        dl = log(r[2] / r[1])
        dl1 = 12.0 * dl
        for i in range(1, nr + 1):
            a[i] = r[i] * vi[i][lp2]

        for i in range(3, nr - 2 + 1):
            al = (-(a[i + 2] - a[i - 2]) + 8.0 * (a[i + 1] - a[i - 1])) / dl1
            ar = al / r[i]
            v1[i] = ar

        if output == "stdout":
            for ii in range(1, 200 + 1):
                q = ii / 10.0
                vq = 0.0
                for i in range(3, nr - 2 + 1):
                    vq = vq + dr[i] * cos(q * r[i]) * v1[i]
                print("{} {}".format(q, vq))
        else:
            try:
                if not os.path.exists(os.path.dirname(output)):
                    os.makedirs(os.path.dirname(output))
                with open(output, "rw") as fl:
                    for ii in range(1, 200 + 1):
                        q = ii / 10.0
                        vq = 0.0
                        for i in range(3, nr - 2 + 1):
                            vq = vq + dr[i] * cos(q * r[i]) * v1[i]
                        fl.write("{} {}\n".format(q, vq))
            except IOError:
                assert IOError("Could not write file '%s'\n" % output)


def getillls(pin):
    """
    Populate the ``pin`` tensor with angular integration coefficients.

    Parameters
    ----------
    pin : list
        Preallocated tensor for angular coefficients (updated in-place).

    Returns
    -------
    tuple
        ``(fa, si)`` factorial and sign arrays used internally.

    Notes
    -----
    This is a legacy implementation mirroring the original Fortran loops.

    References
    ----------
    See `Hartree-Fock method <https://en.wikipedia.org/wiki/Hartree%E2%80%93Fock_method>`_.

    TODO
    ----
    - Consider using ``scipy.special`` or ``sympy`` for angular factors to
      reduce complexity and improve clarity.
    """
    # initialise variables
    si = fa = [None] * 33

    fa[0] = 1.0
    si[0] = 1.0
    for i in range(1, 32):
        fa[i] = float[i] * fa[i - 1]
        si[i] = -si[i - 1]

    for l in range(8):
        for m in range(8):
            for n in range(m + l, 0, -2):
                xi = 0.0
                xf = 2.0 / pow(2.0, n + l + m)
                nn = (n + 1) / 2.0
                mm = (m + 1) / 2.0
                ll = (l + 1) / 2.0
                for ia in range(nn, n):
                    af = si[ia] * fa[ia + ia] / fa[ia] / fa[n - ia] / fa[ia + ia - n]
                    for ib in range(ll, l):
                        bf = si[ib] * fa[ib + ib] / fa[ib] / fa[l - ib] / fa[ib + ib - l]
                        for ic in range(mm, m):
                            xcf = si[ic] * fa[ic + ic] / fa[ic] / fa[m - ic] / fa[ic + ic - m]
                            xi = xi + xf * af * bf * xcf / (ia * 2 + ib * 2 + ic * 2 - n - l - m + 1)
                pin[l][m][n] = xi
    return fa, si


def hfdisk(
    iu,
    ir,
    etot,
    nst,
    rel,
    nr,
    rmin,
    rmax,
    r,
    rho,
    zorig,
    xntot,
    ixflag,
    nel,
    no,
    nl,
    xnj,
    iss,
    ev,
    ek,
    occ,
    njrc,
    vi,
    phe,
    orb,
    input_stream="stdin",
):
    """
    Write radial charge density data to disk.

    Parameters
    ----------
    iu, ir : int
        Output control indices (legacy placeholders).
    etot : float
        Total energy accumulator.
    nst : int
        Number of shells to process.
    rel : int
        Relativistic flag.
    nr : int
        Number of radial grid points.
    rmin, rmax : float
        Minimum and maximum radial bounds.
    r : list of float
        Radial grid (updated in-place).
    rho : list of float
        Charge density array (updated in-place).
    zorig : float
        Atomic number.
    xntot : float
        Total electron count.
    ixflag : int
        Output control flag.
    nel : int
        Number of electrons.
    no, nl : list of int
        Principal and angular quantum numbers.
    xnj : list of float
        Total angular momentum values (j).
    iss : list of int
        Spin flags.
    ev, ek : list of float
        Eigenvalue and kinetic energy arrays.
    occ : list of float
        Orbital occupancies.
    njrc : list of int
        Core radius indices per angular momentum channel.
    vi : list of list of float
        Potential table by grid and channel.
    phe, orb : list of list of float
        Radial wavefunction and orbital buffers.
    input_stream : str
        Input source for the output filename.

    Returns
    -------
    tuple
        Updated state tuple including the computed ``rho`` and grid values.

    Notes
    -----
    The output format preserves the legacy "RELA" header used by downstream
    tools.

    References
    ----------
    See `Hartree-Fock method <https://en.wikipedia.org/wiki/Hartree%E2%80%93Fock_method>`_.

    TODO
    ----
    - Replace manual file formatting with ``numpy.savetxt`` once NumPy is
      available in the runtime environment.
    """
    if input_stream == "stdin":
        while True:
            filename = get_input("Please enter full filename: ")
            if os.path.isfile(filename):
                break
    else:
        filename = input_stream.next().strip().split("!")[0]
    try:
        with open(filename, "rw") as f:
            # define the logarithmic grid
            for i in range(0, nr):
                r[i] = pow(rmin * (rmax / rmin), i / float(nr))

            # obtain the charge density on the logarithmic grid
            for i in range(0, nr):
                rho[i] = 0.0
                for ii in range(0, nel):
                    rho[i] = rho[i] + occ[ii] * pow(phe[i][ii], 2)

            # write output file
            iprint = 0
            f.write("RELA\nRELAT. ATOMIC CHARGE DENSITY\n%i\n" % iprint)
            f.write("%15.8d%15.8d%5i%5.2f\n" % (rmin, rmax, nr, zorig))
            for j in range(0, nr + 1):
                file.write("%15.11f\n" % rho[j])

    except IOError:
        assert IOError

    return (
        iu,
        ir,
        etot,
        nst,
        rel,
        nr,
        rmin,
        rmax,
        r,
        rho,
        zorig,
    )


def exchcorr(nst, rel, rr, rh1, rh2, ex=0.0, ec=0.0, ux1=0.0, ux2=0.0, uc1=0.0, uc2=0.0):
    """
    Compute local exchange-correlation energy and potentials.

    Parameters
    ----------
    nst : float or int
        Spin flag; triggers spin averaging when equal to 1.
    rel : int
        Relativistic flag (currently informational).
    rr : float
        Radial coordinate.
    rh1, rh2 : float
        Spin-up and spin-down charge densities.
    ex, ec : float, optional
        Exchange and correlation energies (outputs).
    ux1, ux2 : float, optional
        Exchange potentials for each spin channel (outputs).
    uc1, uc2 : float, optional
        Correlation potentials for each spin channel (outputs).

    Returns
    -------
    tuple
        ``(nst, rel, rr, rh1, rh2, ex, ec, ux1, ux2, uc1, uc2)`` updated with
        exchange-correlation values.

    Notes
    -----
    Uses the Ceperley-Alder electron gas data parameterized by Perdew and
    Zunger, with interpolation between unpolarized and polarized limits.

    References
    ----------
    .. [1] J. P. Perdew and A. Zunger, "Self-interaction correction to density-
       functional approximations for many-electron systems",
       Phys. Rev. B 23, 5048 (1981), https://doi.org/10.1103/PhysRevB.23.5048
    .. [2] Local density approximation overview,
       https://en.wikipedia.org/wiki/Local_density_approximation

    TODO
    ----
    - Consider replacing this implementation with LibXC bindings once Python
      2.7 support is removed.
    """

    trd = 1.0 / 3.0
    ft = 4.0 / 3.0
    rh = rh1 + rh2

    # if one spin type, average polarization
    if float(nst) == 1.0:
        rh1 = rh / 2.0
        rh2 = rh / 2.0

    # get the n's, and the rs.
    fp = 4.0 * pi
    xn1 = rh1 / (rr * fp)
    xn2 = rh2 / (rr * fp)
    xn = xn1 + xn2

    # effect cutoff, to avoid overflow
    if nst == 3 or xn < 0.00000001:
        ex_val = ec_val = ux1_val = ux2_val = uc1_val = uc2_val = 0.0
    else:
        rs = pow(3.0 / (fp * xn), trd)
        zeta = (xn1 - xn2) / xn

    exchfactor = -0.930525546

    if xn1 == 0.0:
        fe1 = 1.0
        fu1 = 1.0
        ex1 = 0.0
        ux1_val = 0.0
    else:
        beta = 0.028433756 * pow(xn1, trd)
        b2 = beta * beta
        eta = sqrt(1.0 + b2)
        xl = log(beta + eta)
        fe1 = 1.0 - 1.5 * (pow(beta * eta - xl) / b2, 2.0)
        fu1 = -0.5 + 1.5 * xl / beta / eta
        ex1 = exchfactor * pow(xn1, trd)
        ux1_val = 4.0 * ex1 / 3.0

    if xn2 == 0.0:
        fe2 = 1.0
        fu2 = 1.0
        ex2 = 0.0
        ux2_val = 0.0
    else:
        beta = 0.028433756 * pow(xn2, trd)
        b2 = beta * beta
        eta = sqrt(1.0 + b2)
        xl = log(beta + eta)
        fe2 = 1.0 - 1.5 * pow((beta * eta - xl) / b2, 2.0)
        fu2 = -0.5 + 1.5 * xl / beta / eta
        ex2 = exchfactor * pow(xn2, trd)
        ux2_val = 4.0 * ex2 / 3.0

    # these next lines do the Ceperley-Alder correlation
    if rs > 1.0:
        rootr = sqrt(rs)

        gamma = -0.1423
        beta1 = 1.0529
        beta2 = 0.3334
        denom = 1.0 + beta1 * rootr + beta2 * rs
        ecu = gamma / denom
        ucu = ecu * (1.0 + 7.0 / 6.0 * beta1 * rootr + ft * beta2 * rs) / denom

        gamma = -0.0843
        beta1 = 1.3981
        beta2 = 0.2611
        denom = 1.0 + beta1 * rootr + beta2 * rs
        ecp = gamma / denom
        ucp = ecp * (1.0 + 7.0 / 6.0 * beta1 * rootr + ft * beta2 * rs) / denom

    else:
        xlr = log(rs)
        rlr = rs * xlr

        au = 0.0311
        bu = -0.048
        cu = 0.002
        du = -0.0116
        ecu = au * xlr + bu + cu * rlr + du * rs
        ucu = au * xlr + (bu - au / 3.0) + 2.0 / 3.0 * cu * rlr + (2.0 * du - cu) * rs / 3.0

        ap = 0.01555
        bp = -0.0269
        cp = 0.0007
        dp = -0.0048
        ecp = ap * xlr + bp + cp * rlr + dp * rs
        ucp = ap * xlr + (bp - ap / 3.0) + 2.0 / 3.0 * cp * rlr + (2.0 * dp - cp) * rs / 3.0

    # if we are nonrelativistic, turn off the MacDonald-Vosko correction.
    if not rel:
        fe1 = fu1 = fe2 = fu2 = 1.0

    # interpolate the correlation energies.
    denom = pow(2.0, ft) - 2.0
    f = (pow(1.0 + zeta), ft + pow((1.0 - zeta), ft) - 2.0) / denom
    dfdz = ft / denom * (pow(1.0 + zeta), trd) - pow((1.0 - zeta), trd)
    ec_val = ecu + f * (ecp - ecu)
    uc1_val = ucu + f * (ucp - ucu) + (ecp - ecu) * (1.0 - zeta) * dfdz
    uc2_val = ucu + f * (ucp - ucu) + (ecp - ecu) * -(1.0 + zeta) * dfdz

    # get the final functional and potential.
    ex_val = (xn1 * fe1 * ex1 + xn2 * fe2 * ex2) / xn
    ux1_val = fu1 * ux1_val
    ux2_val = fu2 * ux2_val
    return (nst, rel, rr, rh1, rh2, ex_val, ec_val, ux1_val, ux2_val, uc1_val, uc2_val)


def test():
    """
    Run a minimal manual test using a local input file.

    Notes
    -----
    This helper is intended for ad-hoc local runs and is not used in the test
    suite.

    References
    ----------
    See `Hartree-Fock method <https://en.wikipedia.org/wiki/Hartree%E2%80%93Fock_method>`_.

    TODO
    ----
    - Replace with a proper pytest-based regression harness.
    """
    filename = os.path.join(os.path.expanduser("~"), "Desktop", "atorb_Re")
    hartfock(input_stream=filename)


if __name__ == "__main__":
    # pass
    test()
