from __future__ import print_function, division

from math import copysign, log, exp, sqrt, pi, atan, cos, sin, log10, cosh, sinh
import sys


def slices(s, *args):
    """return list of slices"""
    position = 0
    for length in args:
        yield s[position : position + length]
        position += length


def phsh_cav(mufftin_file, phasout_file, dataph_file):
    """
    Description
    -----------
    Potential-to-phase-shift calculation (cavleed package)
    uses loucks grid (e.g. as supplied by the muffin-tin potential
    program). Energies input in hartrees.
    """
    v = rx = [None] * 250
    phs = [None] * 20
    delstore = [[None] * 15] * 401
    estore = [None] * 401

    # first input channels
    mtz = open(mufftin_file, "r")  # unit=4

    # now output channels
    zph = open("zph.o", "rw")  # unit=6
    phasout = open(phasout_file, "rw")  # unit=9
    dataph = open(dataph_file, "rw")  # unit=8

    # standard values for phase shifts calculation
    dataph.write("\n")  # 110
    emin = 1.0
    emax = 12.0
    estep = 0.25
    ianz = (emax - emin) / estep + 1.01
    nl = 12
    nr = mtz.readline()  # 103
    for kkk in range(1, nr + 1):  # 2
        name = mtz.readline()  # 100
        (z, rmt, mtz) = mtz.readline()  # 101
        ntab = mtz.readline()  # 103
        mtz /= 2.0
        for ix in range(1, ntab + 1):  # 19
            rx[ix], v[ix] = mtz.readline().split()  # 219
        zph.write("{} {} {} {}".format(name, z, rmt, mtz))  # 200
        phasout.write("non-relativistic phase shifts for %s\n" % "".join(name))
        phasout.write("  %9.4f  %9.4f  %3i  %3i\n" % (emin, estep, ianz, nl))
        e = emin
        ncount = 0
        while True:
            e *= 2.0
            ncount += 1
            ps(v, rx, ntab, rmt, e, phs, nl)
            e *= 0.5
            for l in range(1, nl + 1):
                zph.write("{} {}".format(e, phs[l]))  # 201
                phasout.write("%9.4f        %8.4f" % (e * 27.21, phs[l]))
            # store phase shifts
            for kk in range(1, nl + 1):
                delstore[ncount][kk] = phs[kk]
            estore[ncount] = e
            e += estep
            if e > emax:
                break

        # write phase shifts as function of energy for plotting
        for kk in range(1, nl + 1):  # 145
            dataph.write("%i2\n" % kk - 1)
            for ii in range(1, ncount + 1):  # 146
                dataph.write("{} {}\n".format(estore[ii], delstore[ii][kk]))
            dataph.write("\n")

        # explicitly close file handles
        zph.close()
        phasout.close()
        dataph.close()


def ps(v, rx, ngrid, rad, e, phs, nl, zph=sys.stdout):
    """
    Desciption
    ----------
    Calculate nl phase shifts(l=0, nl - 1) for an
    atomic potential tabulated on the loucks radial grid.

    Parameters
    ----------
    v : list(float)
        atomic potential (Rydbergs)
    rx : list(float)
        loucks' exponential grid  rx[i] = exp(-8.8+0.05(i-1))
    ngrid : int
        number of tabulated points
    rad : float
        limit of integration of Schrodinger equation (a.u.)
    e : float
        energy (rydbergs)
    phs : list(float)
        phase shifts for l = 0 to nl-1

    Notes
    -----
    reference - Loucks TL., (1967), a.p.w. method, Benjamin, NV.
    """
    # check initialisations
    wf = [None] * 250
    bj = bn = [None] * 25
    xr = [None] * 10
    rf = [None] * 5

    dx = 0.05
    dac = 2.083333e-04
    db = 2.083333e-03
    index = [20.0 * (log(x) + 8.8) + 2.0 for x in range(1)]

    # tabulation of spherical Bessel functions in bj and bn
    es = sqrt(e)
    x = es * rad
    ll = nl + 1
    calcbf(bj, bn, ll, x)

    # integration of the radial Schrodinger equation by the Numerov
    # method (see Loucks p56 ff).  wf contains wave function x radius,
    # on the Loucks grid
    x1 = exp(-8.8)
    for l1 in range(nl):  # 6
        fl = l1 - 1
        fl2 = fl + 0.5
        ff = fl2 * fl2
        y1 = pow(x1, fl2)
        y2 = exp(dx * fl2) * y1
        zph.write("0l?%5.1fy1,y2?%14.5f%14.5f\n" % (fl, y1, y2))  # 60
        gam1 = ff + rx[0] * rx[0] * (v[0] - e)
        gam2 = ff + rx[1] * rx[1] * (v[1] - e)
        wf[0] = y1 * sqrt(rx[0])
        wf[1] = y2 * sqrt(rx[1])
        for ix in range(3, ngrid + 1):  # 2
            gam = ff + rx[ix] * rx[ix] * (v[ix] - e)
            a = 1.0 - dac * gam
            b = -2.0 - db * gam2
            c = 1.0 - dac * gam1
            yn = -(b * y2 + c * y1) / a
            wf[ix] = yn * sqrt(rx[ix])
            y1 = y2
            y2 = yn
            gam1 = gam2
            gam2 = gam

        # Lagrangian interpolation for wavefunction and derivative at radius x.
        # wfn holds wavefunction x radius, and dwfn derivative x radius
        x = rad
        jr = index[rad]

        for j in range(5):  # 3
            xr[j] = rx[jr - 5 + j]
            xr[j + 5] = xr[j]
            rf[j] = wf[jr - 5 + j]

        wfn = 0.0
        dwfn = 0.0
        a = (x - xr[0]) * (x - xr[1]) * (x - xr[2]) * (x - xr[3]) * (x - xr[4])

        for i in range(5):  # 5
            term = (
                a / (x - xr[i]) / (xr[i] - xr[i + 1]) / (xr[i] - xr[i + 2]) / (xr[i] - xr[i + 3]) / (xr[i] - xr[i + 4])
            )
            sum = 0.0

            for j in range(5):
                if i != j:  # goto 4
                    sum += term / (x - xr[j])

            wfn += term * rf[i]
            dwfn += sum * rf[i]

        # logarithmic derivative
        loga = dwfn / wfn - 1.0 / float(rad)

        # phase shifts
        x = es * rad
        a = es * (fl * bj[l1] / x - bj[l1 + 1]) - loga * bj[l1]
        b = es * (fl * bn[l1] / x - bn[l1 + 1]) - loga * bn[l1]
        phs[l1] = pi / 2.0

        if abs(b) > 1.0e-8:
            phs[l1] = atan(a / float(b))
        zph.write("0phase shift?%10.4f" % phs[l1])

    return


def calcbf(bj, bn, nl, x, zph=sys.stdout):
    if abs(x) < 1.0e-6:
        zph.write("argument {} too small for routine calcbf".format(x))
        return

    bj[0] = sin(x) / x
    bn[0] = -cos(x) / x
    if nl == 1:
        return

    bj[1] = (bj[0] - cos(x)) / x
    bn[1] = (bn[0] - sin(x)) / x
    if nl == 2:
        return

    if nl * (nl + 1) > x * x:
        # goto 2
        pass

    # forward recurrence for bj's
    fl = 3.0
    for l in range(2, nl):  # 1
        bj[l] = fl * bj[l - 1] / x - bj[l - 2]
        fl += 2.0
    # goto 5

    # backward recurrence for bj's
    bj0 = bj[0]  # label 2
    bj1 = bj[1]
    nn = max(10, 2 * nl)
    a = 0.0
    b = 1.0
    fl = float(2 * nn + 1)

    for i in range(nn):
        l = nn - i + 1
        c = fl * b / x - a
        if l <= nl:
            bj[l] = c
        a = b
        b = c
        fl -= 2.0

    # normalisation
    b = bj0 / bj[0]

    if abs(bj0) < 0.01:
        b = bj1 / bj[1]
    for l in range(nl):
        bj[l] = b * bj[l]

    # forward recurrence for bn's
    fl = 3.0
    for l in range(2, nl):
        bn[l] = fl * bn[l - 1] / x - bn[l - 2]
        fl += 2.0

    return


def phsh_wil(
    mufftin_file="mtz.i",
    phasout_file="phasout.o",
    dataph_file="dataph.o",
    zph_file="zph.o",
):
    """
    A.R. Williams phase shift program (given a muffin-tin potential)
    """
    # initialise arrays
    e = [None] * 402
    s = c = [[None] * 16] * 402
    delta = delold = [None] * 16
    dell = [None] * 10
    delstore = [[[None] * 16] * 402] * 9

    # common
    e1 = e2 = ne = neuo = nr = nl = float
    r = y = f = ilst = nrr = float

    # input channels
    mtz = open(mufftin_file, "r")  # unit=5

    # now output channels
    zph = open(zph_file, "w")  # unit=6
    phasout = open(phasout_file, "w")  # unit=7
    dataph = open(dataph_file, "w")  # unit=8

    # read ip
    #   ip = 0: only radial wavedef
    #   ip = 1: phase shifts in addition
    #   ip = 2: s and c
    #   ip = 3: produce logarithm of phase shifts
    #   nrr =  number of inequivalent atoms for which we want phase shifts
    ip = 1
    nl2 = mtz.readline()
    zph.write("{}\n".format(nl2))
    dataph.write("titletext: {}".format("?"))  # (11htitletext: ,8hdelta(e))

    # input
    for kkk in range(nrr):
        s16()
        tx = 2.0 * r[nr] / (nr - 1.0)
        de = (e2 - e1) / (max(ne - 1.0, 1.0))
        for i in range(1, ne + 1):
            e[i] = e1 + (i - 1.0) * de

            # radial integration
            s10(e[i])
            t3 = r[nr] * e[i]
            t4 = r[nr] * t3

            for lp1 in range(1, nl + 1):
                l = lp1 - 1
                tlp1 = 2 * l + 1
                t5 = pow(r[nr], lp1)
                ut = f[tlp1][ilst] / tx + float(l) * y[tlp1][ilst] / r[nr]
                t1 = f44(l, t4) * y[2 * lp1][ilst] + t3 * f44(lp1, t4) * y[tlp1][ilst] * t5
                t2 = (f45(l, t4) * ut - t3 * f45(l - 1, t4) * y[tlp1][ilst]) * r[nr] / t5
                s[i][lp1] = t1
                c[i][lp1] = t2

            # produce phase shifts
            for lp in range(1, nl + 1):  # 8
                delold[lp] = 0.0

            for i in range(1, ne + 1):  # 10
                for lp in range(1, nl + 1):  # 11
                    delta[lp] = atan(-pow(abs(e[i]), lp - 0.5) * s[i][lp] / c[i][lp])

                # remove discontinuities by multiples of pi
                for lp in range(1, nl + 1):
                    ls = 0
                    deldif = delta[lp] - delold[lp]  # label 111

                ls += 1
                delta[lp] = delta[lp] - (copysign(pi, deldif) * abs(pi))

                zph.write(
                    "too large change in phase shift [l = {l}]"
                    " since last energy\ndiscontinuity by mulitple"
                    " of pi possible\n".format(l=lp)
                )
                delold[lp] = delta[lp]
            if neuo == 2:
                e[i] = 0.5 * e[i]

        # print phase shifts
        for lp in range(1, nl + 1):
            zph.write("%14.7e%14.7e\n" % (e[i], delta[lp]))  # format(1p8e14.7,  / , 14x, 1p7e14.7,  / )
            # write phase shifts in format used by leed program
            # 71    format(1f7.4)
            # 72    format(10f7.4)

            # store phase shifts
            for kk in range(1, nl + 1):  # 144
                delstore[kkk][i][kk] = delta[kk]
                if ip < 3:
                    continue
                for j in range(1, 9 + 1):  # 14
                    dell[j] = -4.0
                    if delta[j] < 1.0e-4:
                        continue
                    dell[j] = log10(delta[j])

            # write phase shifts as function of energy for plotting
            for kk in range(1, nl + 1):  # 145
                dataph.write('"l = %2i\n' % kk - 1)  # format(3h"l = ,i2)
                for ii in range(1, ne + 1):  # 146
                    dataph.write("{e} {delta}\n".format(e=e[ii], delta=delstore[kkk][ii][kk]))
                dataph.write("\n")  # 15 continue
        # 2 continue

        print("be careful about the order of the elements")

        for ii in range(1, ne + 1):  # 147
            phasout.write("%7.4f" % e[ii])  # format 1f7.4

            for i in range(1, nrr + 1):
                for lp in range(1, nl + 1):
                    phasout.write("%7.4f\n" % delstore[i][ii][lp])

        if ip < 2:
            continue

    for lp1 in range(1, nl + 1):  # 9
        s41(e, s[1][lp1], ne)
        s41(e, c[1][lp1], ne)

    return


def s16(r, rt, v, z, mtz=sys.stdin, zph=sys.stdout):
    """
    Description
    -----------
    s16 inputs data

    cs: core shift (position of zero of energy)
    z: atomic number
    e1,e2: first and last energies desired (in rydbergs or
             hartrees, cf. neui)
    ne: number of energies desired
    nl: number of phase shifts desired ( = lmax+1)
    nr: number of radial grid points used in calculation
    ix = 0: signal to stop
    ix = 1: signal to expect a (new) potential to be read in
    rt: muffin-tin radius (in bohr radii)
    neui,neuo
    if = 1: rydberg unit used for input (neui) and output (neuo) of
            energies and potential
    if = 2: hartree unit (double rydberg) used instead of rydberg
            unit for input (neui) and output (neuo)
    potyp = 1: radial input as v[r)
    potyp = 2: radial input as r*v[r)
    """
    # dimension r(201), v[201, 15)
    rs = zs = ztt = [None] * 202
    # set default values of variables in namelist  / nl16 /
    ix = 1
    e1 = 4.0
    e2 = 24.0
    nl = 9
    nr = 101
    neui = 1
    potyp = 2
    if mtz == sys.stdin:  # read from standard input
        nl16 = mtz.readline()
        mtz.flush()
    else:
        nl16 = mtz.readline()
    cs = 0.0

    if ix < 1:
        return

    if neui != 1:
        cs = 2.0 * cs
        e1 = 2.0 * e1
        e2 = 2.0 * e2

    zph.write("{}\n".format(nl16))
    drdn2 = pow(nr - 1.0, 2.0) / rt

    # read format used for input of r vs. v[r) or r vs. r*v[r)
    # (v is assumed positive)
    for i in range(1, 200 + 1):  # 16
        rs[i], zs[i] = [t(s) for t, s in zip((float, float), mtz.readline().split()[:1])]
        # the next lines assume that the input potential and cs are negative
        zs[i] = -zs[i]
        if rs[i] < 0:
            break

    nrs = i - 1
    if neui != 1:
        for i in range(1, nrs + 1):
            zs[i] = 2.0 * zs[i]

    if potyp != 2:
        pass  # should be: continue ?

    for i in range(1, nrs):
        zs[i] = (zs[i] - cs) * rs[i]
    # go to 21

    for i in range(2, nrs + 1):
        zs[i] = (zs[i] / rs[i] - cs) * rs[i]
        iv = 1  # label 21
        r[1] = 0.0
        ztt[1] = z + z

    for i in range(2, nr + 1):
        r[i] = pow(i - 1.0, 2.0) / drdn2

        while r[i] > rs[iv + 2] and iv + 3 < nrs:
            iv = iv + 1

        ztt[i] = f12(rs[iv], zs[iv], r[i], 4)

        for lp1 in range(1, nl + 1):
            v[i][lp1] = -ztt[i] / r[i]

    return


def f12(x, y, z, n):
    """
    f12 performs iterative interpolation in a table of n values of
    x and y to find the value of y at z
    """
    w = [None] * 11
    w[1] = y[1]

    for i in range(2, n + 1):
        w[i] = y[i]
        u = z - x[i]
        ip1 = i + 1

        for j in range(2, i + 1):
            k = ip1 - j
            w[k] = w[k + 1] + u * (w[k] - w[k + 1]) / (x[k] - x[i])

    return w[1]


def s5(e, nl, nr, v, r, f, y):
    """
    s5 -- Hamming's method for the integration of systems of first
    order differential equations
    """
    eest = [None] * 31
    vme = [None] * 15
    # common / cmrv / r(201), v[201, 15), nr, nl, z
    # common / cm5 / y(30, 4), f(30, 4), ip1
    nj = 2 * nl

    for j in range(1, nj + 1):
        eest[j] = 0.0

    for i in range(5, nr + 1):
        for lp1 in range(1, nl + 1):
            vme[lp1] = (v[i, lp1] - e) * r[i]
        t1 = 2.0 / (i - 1.0)
        ip1 = (i - 1 % 4) + 1
        im2 = (ip1 % 4) + 1
        im1 = (im2 % 4) + 1
        ip0 = (im1 % 4) + 1

        for j in range(1, nj + 1):
            f[j][im2] = y[j][ip1] + (2.0 * (f[j][ip0] + f[j][im2]) - f[j][im1]) / 0.75
            y[j][ip1] = f[j][im2] - 0.925619835 * eest[j]

        for j in range(1, nj + 1, 2):
            jp1 = j + 1
            lp1 = jp1 / 2.0
            flp1 = lp1
            f[j][ip1] = (flp1 * y[j][ip1] + r[i] * y[jp1][ip1]) * t1
            f[jp1][ip1] = (vme[lp1] * y[j][ip1] - flp1 * y[jp1][ip1]) * t1

        for j in range(1, nj + 1):
            y[j][ip1] = y[j][ip0] + (y[j][ip0] - y[j][im2] + 3.0 * (f[j][ip1] + 2.0 * f[j][ip0] - f[j][im1])) / 8.0
            eest[j] = f(j, im2) - y(j, ip1)
            y[j][ip1] += 0.743801653e-1 * eest[j]

        for j in range(1, nj + 1, 2):
            jp1 = j + 1
            lp1 = jp1 / 2.0
            flp1 = lp1
            f[j][ip1] = (flp1 * y[j][ip1] + r[i] * y[jp1][ip1]) * t1
            f[jp1][ip1] = (vme[lp1] * y[j][ip1] - flp1 * y[jp1][ip1]) * t1
    return


def s10(e, cmrv, r, v, nr, nl, z, cm5, y, f, ilst, zph=sys.stdout):
    """
    power series expansion of the solution about the origin
    and radial integration in s5
    """
    tlp1 = int
    a = b = tr = [None] * 10

    ni = 2 * nl
    tz = 2.0 * z
    a[1] = 1.0

    for i in range(1, ni + 1, 2):  # 5
        lp1 = (i + 1) / 2.0
        tlp1 = 2 * lp1
        ep = e - v[4][lp1] - tz / r[4]
        y[i][1] = 0.0
        y[i + 1][1] = 0.0
        a[1] = a[1] / (2.0 * lp1 - 1)
        b[1] = -z * a[1] / float(lp1)

        for j in range(2, 4 + 1):
            tr[j] = pow(r[j], lp1)
            y[i][j] = a[1] * tr[j]
            y[i + 1][j] = b[1] * tr[j]  # 1

        for k in range(1, 9 + 1):  # 3
            a[k + 1] = b[k] / float(k)
            b[k + 1] = -(ep * a[k] + tz * a[k + 1]) / float(tlp1 + k)

            for j in range(2, 4 + 1):  # 2
                tr[j] = tr[j] * r[j]
                y[i][j] = y[i][j] + tr[j] * a[k + 1]
                y[i + 1][j] = y[i + 1][j] + tr[j] * b[k + 1]

            if abs(tr[4] * a[k + 1] / y[i][4]) < 1.0e-4:
                # go to 5
                break

        for k in range(1, 10):
            zph.write(" %10.2e%10i %11e%10.2e" % (e, lp1, r[4], a[k]))
            # 4  format(1pe10.2, i10, 11e10.2)

    for j in range(2, 4 + 1):  # 6
        t1 = 2.0 / (j - 1.0)

        for i in range(1, ni + 1, 2):
            ip1 = i + 1
            lp1 = ip1 / 2.0
            flp1 = lp1
            f[i][j] = (flp1 * y[i][j] + r[j] * y[ip1][j]) * t1
            f[ip1][j] = t1 * ((v[j][lp1] - e) * r[j] * y[i][j] - flp1 * y[ip1][j])
    s5(e)
    return


def f44(l, x):
    """Evaluates the special version of the spherical Bessel function"""
    s = [None] * 20
    js = l + l + 1

    if abs(x / float(js)) > 10.0:
        pass
    else:
        fi = 1.0
        if l >= 1:
            for k in range(2, js, 2):
                fi *= float(k)

        t1 = 1.0 / fi
        dt = 1.0
        t = 1.0
        i = 0
        for k in range(100):
            i += 2
            dt *= -x / float(i * (i + js))
            t += dt
            if abs(dt) < 1.0e-8:
                break

        t1 *= t
        return t1

    t = sqrt(abs(x))
    if x < 0.0:
        s[1] = sinh(t) / t
        if l < 1:
            return s[1]

        s[1] = cosh(t)

    else:
        s[1] = sin(t) / t
        if l <= 0:
            return s[2]

        s[0] = cos(t)

    iss = l + 2

    for i in range(3, iss + 1):
        s[i] = (s[i - 1] * float(2 * i - 5) - s[i - 2]) / x

    return s[iss]


def f45(l, x):
    """Evaluates special version of the spherical Neumann function"""

    iss = l  # check this

    s = [None] * 20
    if l < 0:
        return -f44(l + 1, x)

    js = l + l + 1

    if abs(x / float(js)) <= 10.0:
        fi = 1.0
        if l < 1:
            t1 = fi / float(js)

        else:
            for k in range(3, js + 1, 2):
                fi *= float(k)
                t1 = fi / float(js)

        dt = 1.0
        t = 1.0
        i = 0

        for k in range(1, 100 + 1):
            i += 2
            dt = -dt * x / float(i * (i - js))
            t += dt
            if abs(dt) < 1.0e-8:
                break
            t1 *= t

        return t1

    t = sqrt(abs(x))
    if x < 0.0:
        s[2] = cosh(t)
        if l < 1:
            return s[iss]
        s[1] = -sinh(t) / t

    else:
        s[2] = cos(t)
        if l <= 0:
            return s[2]

        s[1] = -sin(t) / t

    iss = l + 2

    for i in range(3, iss):
        s[i] = s[i - 1] * (2.0 * i - 5.0) - x * s[i - 2]

    return s[iss]


def s41(x, y, n, c, o, d, b, zph=sys.stdout):
    """produces xy data"""
    x = y = [None] * 100
    p = [None] * 97
    y1 = 0
    y2 = 0

    for i in range(n):
        y1 = min(y1, y[i])
        y2 = max(y2, y[i])

    for i in range(len(p)):
        p[i] = b

    t = (len(p) - 1) / (y2 - y1)
    j0 = -y1 * t + 1.5
    p[j0] = 0

    if n >= 30:
        zph.write("{}\n".format(p))  # format(1h1, 34x, 97a1, / / )

    if n < 30:
        zph.write("{}\n".format(p))  # format( / / /  / , 35x, 97a1, / / )

    p[j0] = d

    for i in range(n):
        j = t * (y[i] - y1) + 1.5
        p[j] = c
        zph.write("%16.6e%16.6e  %s\n" % (x[i], y[i], p))  # (1h , 1p2e16.6, 2x, 97a1)
        p[j] = b
        p[j0] = d

    return


def phsh_rel(
    mufftin_file="mtz.i",
    phasout_file="phasout.o",
    dataph_file="dataph.o",
    inpdat_file="inpdat.o",
):
    """
    Description
    -----------
    Calculates relativistic phase shifts

    Parameters
    ----------
    mufftin_file : str
        path to muffin-tin potential output file (to be read into function)
    phasout_file : str
        path to output file containing phase shifts
    dataph_file : str
        path to phase data output file
    inpdat_file : str
        path to redundant input data file (not required)
    """
    # initialise arrays and variables
    jf = jfs = jfd = jfu = float
    name = [float] * 4
    jf = [[None] * 19] * 251
    energ = [None] * 251

    adata = [None] * 7

    zero = 0.0
    one = 1.0
    half = 0.5
    des = 0.025e0
    sub = "sub"
    # record = "nos"

    opt = opt1 = opts = None  # dummy variables
    rmaxi = lsm1 = 1  # dummy

    """
    1 format (3d12.4,4x,i3,4x,d12.4)
    2 format (5e14.6)
 3    format    ( / / / ,10x,a6, / / /  / ,10x,18hmuffin-tin radius = ,
     1f10.6,4x,14hconstant pot. = ,f10.6, /  / ,10x,'atomic data set for z = ',
     2i3,' and l = ',i3,4x,'     ',a2, /  / ,14x,5he(ev),12x,7hl - 0.5,13x,
     37hl + 0.5,12x,10hs-averaged, / )
  4   format (10x,f10.6,3f20.8)
  9   format (1x,a1,i3,3d15.7,22x,a3)
   12 format (5d14.6)
   15 format (a28,t29,a2,t35,a30,t65,f10.7,t76,a3)
   17 format (1h1, / ,10x,'relativistic phase shifts for ',a30, /  / ,10x,
     'exca  = ',f10.6,4x,'excb  = ',f10.6,4x,'exco  = ',f10.6 /  / ,10x,
     'lattice constant  = ',f10.6,' ,',f10.6,' ,',f10.6)
    """

    # sort input and output streams
    inpdat = open(inpdat_file, "w")  # unit=4
    mtz = open(mufftin_file, "r")  # unit=5
    phasout = open(phasout_file, "w")  # unit=7
    dataph = open(dataph_file, "w")  # unit=8

    inpdat.write("input data\n")  # 13 format (1h1, /  / ,t61,'input data',)
    name = mtz.readline().split()  # 10
    (es, de, ue, lsm, _vc) = [t(s) for t, s in zip((float, float, float, int, float), mtz.readline().split())]

    nl = 8  # nl is the number of plotted phase shifts
    inpdat.write("%12.4d%12.4d%12.4d    %s%s%3i\n" % (es, de, ue, opt, opt1, lsm))

    # 75 format (a28, a2, 4x, a30, f10.7, 1x, a3)
    (nz, adata[1], jri, alc, blc, clc, exca, excb, exco) = slices(mtz.readline(), 4, 10, 4, 21, 10, 10, 10, 10, 10, 10)
    # 16 format (i4,f10.6,i4,t21,6f10.6)
    inpdat.write(
        "%4i%10.6f%4i  %10.6f%10.6f%10.6f%10.6f%10.6f%10.6f" % (nz, adata[1], jri, alc, blc, clc, exca, excb, exco)
    )  # 4,76
    # 76 format (i4,f10.6,i4,2x,6f10.6)

    if jri <= 0 or jri > 340:
        inpdat.write("# end of input data\n")
        return

    zp = mtz.readline().split()
    inpdat.write("{}\n".format("".join([zp[j] for j in range(1, jri + 1)])))
    rhoz = -0.90306514e01
    delrho = 0.3125e-01
    rm = exp(rhoz)
    xrx = exp(delrho)

    for j in range(1, jri + 1):
        if zp[j] < 0.0:
            zp[j] = -zp[j]
        if j == jri:
            break
        rm *= xrx

    if de <= zero:
        es = -half
        de = des
        ue = one
    n = (ue - es) / de + half
    n += 1

    phasout.write("relativistic phase shifts for {}\n".format("".join(name)))
    phasout.write("%10.4f%9.4f%5i%5i\n" % (es, de, n, lsm))

    es /= 13.6
    de /= 13.6
    ue /= 13.6

    if n > 250:
        n = 250

    l = 0
    e = es

    kap = -1
    l = 1

    for j in range(1, n + 1):  # label 30
        ttr = dlgkap[e][kap] / 12.5663706
        sbfit(ttr, e, l - 1, rmaxi, jfs)
        jf[j][l] = jfs
        e *= 13.6
        energ[j] = e
        e /= 13.6
        e += de

    kap = -(l + 1)
    lind = l + 1
    e = es
    lvor = 0

    for j in range(1, n + 1):  # label 50
        dlu = dlgkap[e][kap] / 12.5663706
        dld = dlgkap[e][l] / 12.5663706
        sbfit(dld, e, l, rmaxi, jfd)
        sbfit(dlu, e, l, rmaxi, jfu)
        lk = 0

        zfdiff = -(jfd - jfu) * lvor
        if zfdiff > 2.5:
            lk = l
        if zfdiff < -2.5:
            lk = l + 1

        jfs = l * jfd - kap * jfu + lvor * lk * pi
        jfs /= 2 * l + 1

        if jfs > 0.5 * pi:
            jfs -= pi
        if jfs < -0.5 * pi:
            jfs += pi

        jf[j][lind] = jfs
        if lk == 0:
            lvor = copysign(1.0, jfs) * 1.0
            e *= 13.6
            e /= 13.6
            e += de

    l += 1
    for i in range(1, n + 1):  # 90 or 80
        lsm1 += 1
        phasout.write("%9.4f%s\n" % (energ[i], "".join(["%8.5f" % jf[i][l] for l in range(1, lsm1 + 1)])))  # noqa: E741

    for kk in range(1, nl + 1):  # label 145
        dataph.write('"l = %2i\n' % kk - 1)

        for ii in range(1, n + 1):  # label 146
            dataph.write("{} {}\n".format(energ[ii], jf[ii][kk]))

        dataph.write("\n")

        es *= 13.6
        de *= 13.6
        ue *= 13.6
        continue  # 999

    inpdat.write("# end of input data\n")

    return jfs


def dlgkap(e, kappa, z, t, radfun, u, w, zzzz, pot, vcz, ipt, jri):
    """
    Description
    -----------
    dlgkap calculates the logrithmic derivative of the large
    component using the procedure described by loucks in app# endix 7.
    the small multiplicative factor is included.
    potential  in the form of 2zp is to be passed in the common  / zzzz /
    the radial functions are made available in the common  / radfun /
    waber mesh (xnot = -9.03065 and dx = 1 / 32) is used
    jri is the mesh point of the apw sphere radius
    e is the energy to be used (in Rydbergs)
    """

    # initialise parameters
    up = wp = [None] * len(pot)
    sxm = sxk = [None] * 4

    ustart = 1.0e-25
    zilch = 1.0e-30
    test = 1.0e6
    xs = -9.03065133e00
    dx = 3.125e-2
    c = 2.740746e2
    cin = 1.3312581146e-5
    hf = 0.5e0
    th = 0.3333333333e0
    t2 = 2.0e0
    t7 = 7.0e0
    t11 = 11.0e0
    t12 = 12.0e0
    t14 = 14.0e0
    t26 = 26.0e0
    t32 = 32.0e0
    # zero = 0.1e+0

    # 83 format (10h hard test,d20.8, i5, 4d20.8 )
    # **********set up for relativistic or no relativistic effect************
    if ipt <= 0:
        cin = 0.00
    # ******************* set up starting values  ***************************
    dx2 = hf * dx
    x20 = 0.3 * dx
    xmft = 4.4444444444444e-2 * dx
    ts = exp(xs)
    tdx = exp(dx)
    hoc = (vcz * ts + pot[1]) / c
    xk = kappa
    u[1] = ustart  # label 49

    if abs(hoc / xk) > 0.05:
        p = (xk + sqrt(xk * xk - hoc * hoc)) / hoc
    else:
        p = (xk + abs(xk)) / hoc - hf * hoc / abs(xk)
        tc = (e + vcz) * ts + pot[1]  # label 7
        vc = cin * tc
        w[1] = c * p * ustart

    # ***************** start Runge-Kutte procedure  ************************
    x = xs  # label 11
    n = 1
    while True:
        ik = 0  # label 25
        np1 = n + 1
        xc = x
        bgc = pot[n]
        wc = w[n]
        uc = u[n]
        ik += 1  # label 20
        t = exp(xc)
        tc = (e + vcz) * t + bgc
        vc = cin * tc
        sxk[ik] = dx2 * (wc * (vc + t) - xk * uc)  # label 12
        sxm[ik] = dx2 * (xk * wc - tc * uc)
        if ik == 1:
            xc += dx2
            uc += sxk[1]
            wc += sxm[1]
            bgc = hf * (bgc + pot[np1])
        elif ik == 2:
            uc += sxk[2] - sxk[1]
            wc += sxm[2] - sxm[1]
        elif ik == 3:
            xc += dx2
            uc += t2 * sxk[3] - sxk[2]
            wc += t2 * sxm[3] - sxm[2]
            bgc = pot[np1]
        elif ik == 4:
            w[np1] = w[n] + (sxm[1] + sxm[4] + t2 * (sxm[2] + sxm[3])) * th

        u[np1] = u[n] + (sxk[1] + sxk[4] + t2 * (sxk[2] + sxk[3])) * th
        up[np1] = (vc + t) * w[np1] - xk * u[np1]
        wp[np1] = xk * w[np1] - tc * u[np1]
        x += dx  # label 24
        n = np1
        if n >= 6:
            break

    # end of starting integration.  begin milne procedure.
    t = exp(x)
    while n < jri:
        t *= tdx  # label 26
        np1 = n + 1
        nm1 = n - 1
        nm2 = n - 2
        nm3 = n - 3
        nm4 = n - 4
        nm5 = n - 5
        tc = (e + vcz) * t + pot[np1]
        vc = cin * tc
        unp = u[nm5] + x20 * (t11 * (up[n] + up[nm4]) + t26 * up[nm2] - c - t14 * (up[nm1] + up[nm3]))  # label 27
        wnp = w[nm5] + x20 * (t11 * (wp[n] + wp[nm4]) + t26 * wp[nm2] - c - t14 * (wp[nm1] + wp[nm3]))
        nit = 0

        while n < jri:
            up[np1] = (vc + t) * wnp - xk * unp  # label 33
            wp[np1] = xk * wnp - tc * unp
            unp2 = u[nm3] + (t7 * (up[np1] + up[nm3]) + t32 * (up[nm2] + up[n]) + t12 * up[nm1]) - c * xmft
            wnp2 = w[nm3] + (t7 * (wp[np1] + wp[nm3]) + t32 * (wp[nm2] + wp[n]) + t12 * wp[nm1]) - c * xmft

            # compare predictor with corrector
            if abs(test * (unp2 - unp)) > abs(unp2):
                if nit < 5:
                    nit += 1  # label 81
                    wnp = wnp2
                    unp = unp2

            if abs(test * (wnp2 - wnp)) <= abs(wnp2):
                w[np1] = wnp2  # label 32
                u[np1] = unp2
                n = np1
                break

    # end of milne procedure
    if abs(u[jri]) > zilch:
        u[jri] = copysign(zilch, u[jri]) * abs(zilch)

    p = (t + vc) / t  # 37
    wnp = p * w[jri] / u[jri]
    unp = wnp - (kappa + 1) / t
    return 12.5663706 * t * t * unp


def sbfit(t, e, l, r, jfs):
    se = int(e)
    sr = int(r)
    st = int(t)
    kappa = sqrt(se)
    x = kappa * sr
    bj1 = sin(x) / x
    bn1 = -cos(x) / x
    bj2 = bj1 / x + bn1
    bn2 = bn1 / x - bj1

    if l > 0:
        ls = 1

        while True:
            ls += 1
            bjt = (2 * ls - 1) * bj2 / x - bj1
            bnt = (2 * ls - 1) * bn2 / x - bn1
            bj1 = bj2
            bj2 = bjt
            bn1 = bn2
            bn2 = bnt
            if l + 1 - ls > 0:
                continue
            break

    dl = st / (sr * sr)
    dl = dl - l / sr
    an = dl * bj1 + kappa * bj2
    ad = dl * bn1 + kappa * bn2
    jfs = 3.141592654 / 2.0

    if abs(ad) - 1.0e-8 > 0:
        jfs = atan(an / ad)

    return jfs
