#!/usr/bin/env python
# encoding: utf-8

from __future__ import print_function, division

from math import cos, pi, pow, exp, log, sinh, sqrt
from copy import deepcopy
import sys
import os


def get_input(prompt):
    if sys.hexversion > 0x03000000:
        return input(prompt)
    else:
        return raw_input(prompt)


class hartfock(object):
    """
    hartfock class
    -------------------

    there are nr grid points, and distances are in bohr radii...

    r(i)=rmin*(rmax/rmin)**(float(i)/float(nr)) , i=1,2,3,...nr-1,nr



    The orbitals are store in phe(), first index goes 1...nr, the
    second index is the orbital index (i...nel)

    Look at the atomic files after printing this out to see everything...

    Suffice it to say, that the charge density at radius r(i)
    in units of electrons per cubic Bohr radius is given by

    sum of j=1...nel,
    occ[j]*phe(i,j)*phe(i,j)/(4.*3.14159265....*r(i)*r(i))...

    Think of the phe functions as plotting the radial wave-functions
    as a function of radius...on our logarithmic mesh...

    Final note
    ----------

    the Dirac equation is solved for the orbitals, whereas their density
    is treated by setting phe to the square root of Dirac's F*F+G*G
    times the sign of G...

    so we are doing Dirac-Fock, except that we are not treating exchange
    exactly, in terms of working with major and minor components of the
    orbitals, and the phe's give the CORRECT CHARGE DENSITY...

    the above approximation ought to be very small for valence states,
    so you need not worry about it...

    the Breit interaction has been neglected altogether...it should not
    have a huge effect on the charge density you are concerned with...

    authors
    -------
    Eric Shirley - original implementation
    Liam Deacon - python translation

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
        Description
        -----------
        Performs calculation of Hartfock subroutine

        Parameters
        ----------
        iorbs : int
            number of orbitals
        iside : int
            number of sides
        lmax : int
            maximum angular quantum number
        ihmax : int
            height?
        nrmax : int
            number of radii in grid
        ntmax : int
            ?
        npmax : int
            ?
        input_stream : str
            may be either 'stdin' for user input or else a path to input file
        """

        # initialise class variables
        io2 = iorbs * (iorbs + 1) / 2.0
        ijive = io2 * (io2 + 1) / 2.0
        no = nl = nm = xnj = iss = ev = ek = occ = [None] * iorbs
        r = dr = r2 = v = rho = [None] * nrmax
        njrc = [None] * 4
        vi = [[None] * 7] * nrmax
        phe = orb = [[None] * iorbs] * nrmax
        rpower = [[None] * 16] * nrmax

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
                    alfa = float(
                        "".join(
                            [
                                ch
                                for ch in f.next().split("!")[0]
                                if ch is not "d" and ch is not "D"
                            ]
                        )
                    )
                else:
                    alfa = float(
                        get_input(
                            "Enter exchange correlation method"
                            " (0=HARTREE-FOCK, >0=LDA, <0=XALPHA): "
                        )
                    )

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
                            (corpol, rs, rp, sd) = get_input(
                                "enter ALPHA, RS, RP, RD: "
                            ).split()[:3]
                            break
                        except:
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
                    (iunit, corpol) = [
                        t(s)
                        for t, s, in zip(
                            (int, float), f.next().split("!")[0].split()[:1]
                        )
                    ]
                    (ilev, inum, eold) = [
                        t(s)
                        for t, s in zip(
                            (int, int, float), f.next().split("!")[0].split()[:2]
                        )
                    ]
                else:
                    while True:
                        try:
                            (iunit, corpol) = [
                                t(s)
                                for t, s in zip(
                                    (int, float),
                                    get_input("Please enter IUNIT, CORPOL: ").split()[
                                        :2
                                    ],
                                )
                            ]
                            break
                        except:
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
                        except:
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
    """abinitio subroutine"""

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
            t(s)
            for t, s in zip(
                (int, int, float, float, float), f.next().split("!")[0].split()[:5]
            )
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
        print(
            "  %4i%2i%4i%10.4f%18.6f\n"
            % (no[i], nl[i], nm[i], nj, "/2", iss[i], occ[i], ev[i])
        )
        print("Total energy =  %14.6f  14.6f" % (etot, etot * 27.2116))

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
    """atsolve subroutine"""
    # initialise arrays
    q0 = xm1 = xm2 = v = [None] * len(r)

    # initialise eerror, the biggest change in an eigenvalue, and etot.
    eerror = etot = zeff = 0.0

    # run through all the orbitals.  calculate those not in the core.
    for i in range(nel):
        if i > nfc:
            idoflag = 1
            (
                i,
                orb,
                nl[i],
                iss[i],
                idoflag,
                v,
                zeff,
                zorig,
                rel,
                nr,
                r,
                r2,
                dl,
                q0,
                xm1,
                xm2,
                njrc,
                vi,
            ) = setqmm(
                i,
                orb,
                nl[i],
                iss[i],
                idoflag,
                v,
                zeff,
                zorig,
                rel,
                nr,
                r,
                r2,
                dl,
                q0,
                xm1,
                xm2,
                njrc,
                vi,
            )

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
    """getpot subroutine"""

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

            if (
                iss[i] != iss[j]
                and occ[i] <= 1.0
                and occ[j] <= 1.0
                and xnj[i] >= 0.0
                and xnj[j] >= 0.0
            ):
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
                    pow(la + la + 1, 2.0)
                    * pow(cg[li][lj][la][-mi][mj] * cg[li][lj][la][0][0], 2.0)
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
                if (
                    no[jj] == no[ii + 1]
                    and nl[jj] == nl[ii + 1]
                    and iss[jj] == iss[ii + 1]
                    and iuflag == 1
                ):
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
    """elsolve subroutine"""
    el = -zorig * zorig / float(n * n)
    eh = 0.0
    etol = 0.0000000001

    ief = x0 = float
    nn = int

    while True:
        e = (el + eh) / 2.0  # label 155
        istop = 0
        (
            e,
            l,
            xkappa,
            n,
            nn,
            istop,
            ief,
            x0,
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
        ) = integ(
            e,
            l,
            xkappa,
            n,
            nn,
            istop,
            ief,
            x0,
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
        )
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
    """augment subroutine"""
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


def setqmm(
    i, orb, l, ns, idoflag, v, zeff, zorig, rel, nr, r, r2, dl, q0, xm1, xm2, njrc, vi
):
    """setqmm subroutine"""
    c = 137.038
    alpha = rel / c
    aa = alpha * alpha
    a2 = aa / 2.0

    lp = l + 1
    lpx = lp
    if lp > 4:
        lpx = 4
    lp2 = l + l + 1
    if lp2 > 7:
        lp2 = 7
    zeff = zorig
    if njrc[lpx] > 0:
        zeff = 0.0
    zaa = zeff * aa
    za2 = zeff * a2

    if idoflag:
        if not njrc[lpx]:
            if idoflag == 1:
                for j in range(1, nr + 1):
                    v[j] = -zeff / r[j] + orb[j][i]

            for j in range(2, nr - 1):
                dvdl = (orb[j + 1][i] - orb[j - 1][i]) / (2.0 * dl)
                ddvdrr = (
                    (orb[j + 1][i] + orb[j - 1][i] - 2.0 * orb[j][i]) / (dl * dl) - dvdl
                ) / r2[j]
                xm1[j] = -a2 * dvdl / r[j] - za2 / r2[j]
                xm2[j] = -a2 * ddvdrr + zaa / r2[j] / r[j]

            xm1[nr] = xm1[nr - 1]
            xm2[nr] = xm2[nr - 1]
            xm1[1] = xm1[2] + za2 / r2[2] - za2 / r2[1]
            xm2[1] = xm2[2] - zaa / r2[2] / r[2] + zaa / r2[1] / r[1]
    else:
        if idoflag == 1:
            for j in range(1, nr + 1):
                v[j] = vi[j][lp2] + orb[j][i]

        for j in range(2, nr - 1 + 1):
            dvdl = (v[j + 1] - v[j - 1]) / (2.0 * dl)
            ddvdrr = ((v[j + 1] + v[j - 1] - 2.0 * v[j]) / (dl * dl) - dvdl) / r2[j]
            xm1[j] = -a2 * dvdl / r[j]
            xm2[j] = -a2 * ddvdrr

        xm1[nr] = xm1[nr - 1]
        xm2[nr] = xm2[nr - 1]
        xm1[1] = xm1[2]
        xm2[1] = xm2[2]

    # figure out the (Desclaux-Numerov) effective potential.
    xlb = l + pow(0.5, 2.0) / 2.0
    for j in range(1, nr + 1):
        vj = v[j]
        q0[j] = vj * (1.0 - a2 * vj) + xlb / r2[j]

    return (
        i,
        orb,
        l,
        ns,
        idoflag,
        v,
        zeff,
        zorig,
        rel,
        nr,
        r,
        r2,
        dl,
        q0,
        xm1,
        xm2,
        njrc,
        vi,
    )


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
    Description
    -----------
    Initialise the radial charge grid

    Parameters
    ----------
    zorig : float

    nr : int
        Number of radial grid points

    rmin : float
        Minimum radius

    rmax : float
        Maximum radius

    r : list
        Dummy list of radii

    dr : list
        Dummy list to be populated

    r2 : list
        Dummy list to be populated

    dl : float
        Dummy float

    njrc : list
        Dummy list to be populated

    xntot :

    nel : int

    Returns
    -------
    tuple : (zorig, nr, rmin, rmax, r, dr, r2, dl, njrc, xntot, nel)

    """

    if input_stream == "stdin":
        (zorig, nr) = [
            t(s) for t, s in zip((float, int), get_input("Enter Z, NR: ").split())
        ]

    elif isinstance(input_stream, file):
        (zorig, nr) = [
            t(s)
            for t, s in zip((float, int), input_stream.next().split("!")[0].split())
        ]

    else:
        raise IOError("input stream is not a file handle or 'stdin'")

    rmin = 0.0001 / zorig
    rmax = 800.0 / sqrt(zorig)

    nr, rmin, rmax, r, dr, r2, dl = setgrid(nr, rmin, rmax, r, dr, r2, dl)
    njrc[j] = [0 for j in range(len(njrc))]

    return (zorig, nr, rmin, rmax, r, dr, r2, dl, njrc, xntot, nel)


def setgrid(nr, rmin, rmax, r, dr, r2, dl):
    """Set the radial grid values"""
    ratio = rmax / rmin
    dl = log(ratio) / float(nr)
    xratio = exp(dl)
    xr1 = sqrt(xratio) - sqrt(1.0 / xratio)
    for i in range(len(r)):
        r[i] = pow(rmin * xratio, i)
        dr[i] = r[i] * xr1
        r2[i] = r[i] * r[i]

    return (nr, rmin, rmax, r, dr, r2, dl)


def integ(
    e, l, xkappa, n, nn, istop, ief, x0, phi, z, v, q0, xm1, xm2, nr, r, dr, r2, dl, rel
):
    """integrate out count nodes"""
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
    ief = 0

    # see Desclaux and documentation to see the origin of the below equations.
    # here, we set up the first two points.
    t = e - v(1)
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
    p1 = dk2 * pow((r[2] / r[1], ss2) - (r[2] - r[1]) * z / xlp) * sqrt(xm0 / xm)
    phi[2] = p1 * sqrt(xm * r[2]) / dk2

    # if istop is set, the we know to stop there.  If it is zero, it shall
    # be set to the classical turning point.
    is0 = istop
    if not istop:
        for j in range(nr - 1, 2 - 1, -1):
            if e > v[j]:
                break
            ief = -1
            return ief
        istop = j

    # initialize number of nodes, and determine the ideal number.
    nn = 0
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
            nn += 1
            if nn > nnideal:
                ief = 1
                return ief

        p0 = p1
        p1 = p2

    if istop > 0:
        psip2 = phi[istop + 2] - phi[istop - 2]
        psip1 = phi[istop + 1] - phi[istop - 1]
        psip = (8.0 * psip1 - psip2) / (12.0 * dl * r[istop])
        x0 = psip / phi[istop]

    if not is0:
        return

    for i in range(istop + 3, nr - 1 + 1):
        t = e - v[i]
        xm = 1.0 + a2 * t
        tm = xm + xm
        xmx = xm1[i] / xm
        p2 = (2.0 - dl5 * xk2) * p1 / dk2 - p0
        if p2 / p1 > 1.0:
            ief = -1
            return ief

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
                    nn += 1
                    if nn > nnideal:
                        ief = 1
                        return

            p0 = p1
            p1 = p2

    return ief


def clebschgordan(nel, nl, cg, si, fa):
    """
    routine to generate Clebsch-Gordan coefficients, in the form of
    cg(l1,l2,L,m1,m2) = <l1,m1;l2,m2|L,m1+m2>, according to Rose's
    'Elementary Theory of Angular Momentum', p. 39, Wigner's formula.
    those coefficients listed are only those for which l1.ge.l2.
    coefficients known to be zero because of either the L or M
    selection rules are not computed, and should not be sought.
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
    """pseudo subroutine"""

    # initialise
    nm = [0] * nel
    njrcdummy = deepcopy(njrc)
    q0 = xm1 = xm2 = [None] * len(r)
    rpower = [[None] * 7] * len(r)
    zeff = etot2 = float

    # read input
    if input_stream == "stdin":
        (np, corpol, rnorm) = [
            t(s)
            for t, s in zip(
                (int, float, float), get_input("Please enter NP CORPOL RNORM: ").split()
            )
        ]

    elif isinstance(input_stream, file):
        (np, corpol, rnorm) = [
            t(s)
            for t, s in zip(
                (int, float, float), input_stream.readline().split("!")[0].split()
            )
        ]
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
            ns = 1
            setqmm(
                i,
                orb,
                nl[i],
                ns,
                idoflag,
                vi[1][lp2],
                zeff,
                zorig,
                rel,
                nr,
                r,
                r2,
                dl,
                q0,
                xm1,
                xm2,
                njrcdummy,
                vi,
            )
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
    f = vf[2]
    r21 = rf[2] - rf[1]
    r32 = rf[3] - rf[2]
    v21 = vf[2] - vf[1]
    v32 = vf[3] - vf[2]
    fp = (v21 + v32) / (r21 + r32)
    fpp = (v32 / r32 - v21 / r21) / ((r21 + r32) / 2.0)
    return f, fp, fpp, rf, vf


def hb(x, factor):
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
    """fitx0 subroutine"""

    vl = -1000000.0
    vh = 1000000.0
    dummy = nn = ief = xactual = None  # initialise
    while True:
        idoflag = 2  # label 115
        ns = 1
        xkappa = -1.0

        (
            i,
            orb,
            l,
            ns,
            idoflag,
            v,
            zeff,
            dummy,
            rel,
            nr,
            r,
            r2,
            dl,
            q0,
            xm1,
            xm2,
            njrc,
        ) = setqmm(
            i,
            orb,
            l,
            ns,
            idoflag,
            v,
            zeff,
            dummy,
            rel,
            nr,
            r,
            r2,
            dl,
            q0,
            xm1,
            xm2,
            njrc,
        )

        (
            e,
            l,
            xkappa,
            n,
            nn,
            jrt,
            ief,
            xactual,
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
        ) = integ(
            e,
            l,
            xkappa,
            n,
            nn,
            jrt,
            ief,
            xactual,
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
        )

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
    """pseudize subroutine"""
    nn = ief = x0 = x00 = xm = xp = fp = fpp = psi = psip = psipp = None
    xdummy = [None] * len(phi)
    rf = vf = [None] * len(njrc) - 1

    lp = l + 1
    xkappa = -1.0
    istop = nr
    while ev > q0[istop]:
        istop -= 1

    (
        ev,
        l,
        xkappa,
        n,
        nn,
        istop,
        ief,
        xdummy,
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
    ) = integ(
        ev,
        l,
        xkappa,
        n,
        nn,
        istop,
        ief,
        xdummy,
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
    )

    if rcut is None or factor is None:
        (rcut, factor) = get_input(
            "Please enter the cutoff radius, and factor: "
        ).split()

    if rcut < 0.0:
        xnodefrac = -rcut
        j = istop
        while phi[j - 1] / phi[j] <= 1.0:
            j -= 1
        if n > l + 1:
            k = j
        while phi[k - 1] / phi[k] <= 1.0:
            k -= 1
    else:
        k = 1

    rcut = r[k] + xnodefrac * (r[j] - r[k])

    jrc = 1.0 + float(nr - 1) * log(rcut / rmin) / log(rmax / rmin)
    rcut = r[jrc]
    rtest = 2.0 * rcut
    jrt = 1.0 + float(nr - 1) * log(rtest / rmin) / log(rmax / rmin)
    njrc[lp] = jrt
    rtest = r[jrt]
    # switch = phi[jrt] / abs(phi[jrt])
    print("RCUTOFF = %8.4f  JRC = %5i" % (rcut, jrc))  # 94 format(1x,2d15.8)
    print("RTEST   = %8.4f  JRT = %5i" % (rtest, jrt))  # 1x,1a10,1f8.4,1a8,1i5
    integ(
        ev,
        l,
        xkappa,
        n,
        nn,
        jrt,
        ief,
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
        rel,
    )
    for ii in range(len(phi)):
        phi[ii] /= phi[jrt]

    xn00 = 0.0

    for ii in range(1, jrt - 1 + 1):
        xn00 += dr[ii] * phi[ii] * phi[ii]

    xn00 += dr[jrt] * phi[jrt] * phi[jrt] / 2.0
    de = 0.0001
    ee = ev + de / 2.0
    integ(
        ee,
        l,
        xkappa,
        n,
        nn,
        jrt,
        ief,
        xp,
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
    )
    ee = ev - de / 2.0
    integ(
        ee,
        l,
        xkappa,
        n,
        nn,
        jrt,
        ief,
        xm,
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
    )
    c00 = (xm - xp) / (2.0 * de)
    print(c00, x00)  # format 94
    print(xn00)  # format 94
    ruse = 0.0
    v0 = v[jrc]
    dvdl = (8.0 * (v[jrc + 1] - v[jrc - 1]) - (v[jrc + 2] - v[jrc - 2])) / (12.0 * dl)
    ddvdll = (
        16.0 * (v[jrc + 1] + v[jrc - 1]) - 30.0 * v[jrc] - v[jrc + 2] - v[jrc - 2]
    ) / (12.0 * dl * dl)
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

        fitx0(
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

    phi0 = deepcopy(phi)
    vraw = deepcopy(v)

    xi0 = 0.0
    xi1 = 0.0
    xi2 = 0.0
    for ii in range(1, jrt + 1):
        f = hb(r[ii] / rcut, factor)
        ph2 = dr[ii] * phi0[ii] * phi0[ii]
        xi0 = xi0 + ph2
        if ii <= jrt:
            xi1 += ph2 * f
            xi2 += ph2 * f * f

    ph2 = phi0[jrt] * phi0[jrt]
    xi0 /= ph2
    xi1 /= ph2
    xi2 /= ph2
    quant = xi1 * xi1 + xi2 * (c00 - xi0)
    if quant > 0.0:
        deltal = (sqrt(xi1 * xi1 + xi2 * (c00 - xi0)) - xi1) / xi2
    else:
        deltal = (c00 - xi0) / (2.0 * xi1)

    print("DELTAL = %11.8f" % deltal)  # format (1x,1a9,1f11.8)

    yl = [None] * len(phi0)

    while True:  # 225
        for ii in range(1, jrt):
            yl[ii] = phi0[ii] * hb(r[ii] / rcut, factor)
            phi[ii] = phi0[ii] + deltal * yl[ii]
            if phi[ii] < 0.0:
                print("Big trouble# # #  cross axis# # # ")
                sys.exit(1)

        for ii in range(1, jrt - 1 + 1):
            if phi[ii] == 0.0 or yl[ii] == 0.0:
                break
            jj = ii
            if ii == 1:
                jj = 2
            for j in range(jj - 1, jj + 1 + 1):
                rf[2 + j - jj] = r[j]
                vf[2 + j - jj] = hb(r[j] / rcut, factor)
            (f, fp, fpp, rf, vf) = parabreg(f, fp, fpp, rf, vf)
            for j in range(jj - 1, jj + 1):
                vf[2 + j - jj] = phi0[j]
            (psi, psip, psipp, rf, vf) = parabreg(psi, psip, psipp, rf, vf)
            v[ii] = (
                vraw[ii]
                + (1.0 - phi0[ii] / phi[ii])
                * (2.0 * psip / psi * fp / f + fpp / f)
                / 2.0
            )

            fitx0(
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
            integ(
                ev,
                l,
                xkappa,
                n,
                nn,
                jrt,
                ief,
                x0,
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
            )

        for ii in range(1, jrt + 1):
            phi[ii] = phi[ii] / phi[jrt]

        xn0 = 0.0

        for ii in range(1, jrt - 1 + 1):
            xn0 = xn0 + dr[ii] * phi[ii] * phi[ii]

        xn0 = xn0 + dr[jrt] * phi[jrt] * phi[jrt] / 2.0
        de = 0.0001
        ee = ev + de / 2.0
        integ(
            ee,
            l,
            xkappa,
            n,
            nn,
            jrt,
            ief,
            xp,
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
        )
        ee = ev - de / 2.0
        integ(
            ee,
            l,
            xkappa,
            n,
            nn,
            jrt,
            ief,
            xm,
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
        )
        c0 = (xm - xp) / (2.0 * de)
        print(c0, x0)
        print(xn0)
        if abs(c0 - c00) > 0.000000001:
            dqddel = 2.0 * (xi1 + deltal * xi2)
            deltal = deltal + (c00 - c0) / dqddel
        else:
            break

    print(c0, x0)
    print("NCPP achieved # # # ")
    return


def fourier(nr, r, dr, r2, vi, output="stdout"):
    """fourier subroutine"""

    a = [None] * len(r)
    v1 = [None] * len(r)

    for l in range(0, 3 + 1):
        lp2 = l + l + 1
        dl = log(r[2] / r[1])
        dl1 = 12.0 * dl
        # dl2 = 12. * dl * dl

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
    """getills subroutine"""
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
                        bf = (
                            si[ib] * fa[ib + ib] / fa[ib] / fa[l - ib] / fa[ib + ib - l]
                        )
                        for ic in range(mm, m):
                            xcf = (
                                si[ic]
                                * fa[ic + ic]
                                / fa[ic]
                                / fa[m - ic]
                                / fa[ic + ic - m]
                            )
                            xi = xi + xf * af * bf * xcf / (
                                ia * 2 + ib * 2 + ic * 2 - n - l - m + 1
                            )
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
    """ """
    # rden = 3.0
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
            # nden = nr * (log(rden / rmin) / log(rmax / rmin))
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


def exchcorr(
    nst, rel, rr, rh1, rh2, ex=0.0, ec=0.0, ux1=0.0, ux2=0.0, uc1=0.0, uc2=0.0
):
    """
    Description
    -----------
    Exchange correlation routine, via Ceperley-Alder, as parametrized by
    Perdew and Zunger, Phys. Rev. B 23, 5048.  we use their interpolation
    between the unpolarized and polarized gas for the correlation part.

    Parameters
    ----------
    nst : float, int or bool
        Will cause spin averaging if polarisation is equal to 1

    rel : int
        Flag to signify whether to account for relativistic effects.

    rr :

    rh1 :

    rh2 :

    ex :

    ec :

    ux1 :

    ux2 :

    uc1 :

    uc2 :


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
        ex = ec = ux1 = ux2 = uc1 = uc2 = 0.0
    else:
        rs = pow(3.0 / (fp * xn), trd)
        zeta = (xn1 - xn2) / xn

    exchfactor = -0.930525546

    if xn1 == 0.0:
        fe1 = 1.0
        fu1 = 1.0
        ex1 = 0.0
        ux1 = 0.0
    else:
        beta = 0.028433756 * pow(xn1, trd)
        b2 = beta * beta
        eta = sqrt(1.0 + b2)
        xl = log(beta + eta)
        fe1 = 1.0 - 1.5 * (pow(beta * eta - xl) / b2, 2.0)
        fu1 = -0.5 + 1.5 * xl / beta / eta
        ex1 = exchfactor * pow(xn1, trd)
        ux1 = 4.0 * ex1 / 3.0

    if xn2 == 0.0:
        fe2 = 1.0
        fu2 = 1.0
        ex2 = 0.0
        ux2 = 0.0
    else:
        beta = 0.028433756 * pow(xn2, trd)
        b2 = beta * beta
        eta = sqrt(1.0 + b2)
        xl = log(beta + eta)
        fe2 = 1.0 - 1.5 * pow((beta * eta - xl) / b2, 2.0)
        fu2 = -0.5 + 1.5 * xl / beta / eta
        ex2 = exchfactor * pow(xn2, trd)
        ux2 = 4.0 * ex2 / 3.0

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
        ucu = (
            au * xlr
            + (bu - au / 3.0)
            + 2.0 / 3.0 * cu * rlr
            + (2.0 * du - cu) * rs / 3.0
        )

        ap = 0.01555
        bp = -0.0269
        cp = 0.0007
        dp = -0.0048
        ecp = ap * xlr + bp + cp * rlr + dp * rs
        ucp = (
            ap * xlr
            + (bp - ap / 3.0)
            + 2.0 / 3.0 * cp * rlr
            + (2.0 * dp - cp) * rs / 3.0
        )

    # if we are nonrelativistic, turn off the MacDonald-Vosko correction.
    if not rel:
        fe1 = fu1 = fe2 = fu2 = 1.0

    # interpolate the correlation energies.
    denom = pow(2.0, ft) - 2.0
    f = (pow(1.0 + zeta), ft + pow((1.0 - zeta), ft) - 2.0) / denom
    dfdz = ft / denom * (pow(1.0 + zeta), trd) - pow((1.0 - zeta), trd)
    ec = ecu + f * (ecp - ecu)
    uc1 = ucu + f * (ucp - ucu) + (ecp - ecu) * (1.0 - zeta) * dfdz
    uc2 = ucu + f * (ucp - ucu) + (ecp - ecu) * -(1.0 + zeta) * dfdz

    # get the final functional and potential.
    ex = (xn1 * fe1 * ex1 + xn2 * fe2 * ex2) / xn
    ux1 = fu1 * ux1
    ux2 = fu2 * ux2
    uc1 = uc1
    uc2 = uc2

    return (nst, rel, rr, rh1, rh2, ex, ec, ux1, ux2, uc1, uc2)


def test():
    """for testing"""
    filename = os.path.join(os.path.expanduser("~"), "Desktop", "atorb_Re")
    hartfock(input_stream=filename)


if __name__ == "__main__":
    # pass
    test()
