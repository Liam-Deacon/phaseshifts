from math import pi, min, max, log, sin, cos, tan, asin, acos, atan

def cavpot(mtz_string, slab_flag, atomic_file, cluster_file, mufftin_file, output_file, info_file):
       
    # initialize variables
    nieq = ntot = int
    rx = rs = pot = [float] * 550
    rc = [[float] *3] * 3
    title = str
    sum = float
    rk = trk = sig = rho = vh = vs = vmad = ad = [[float]]
    zm = tzm = z = zc = rmt = [float]
    jrmt = jrmt2 = nrr = [int]
    ncon = nx = na = [int]
    ia = [[int]]
    name = [str]
    wfn = str
    wfn0 = "rela"
    wfn1 = "herm"
    wfn2 = "clem"
    wfn3 = "pote"
    #common /wk/ wk1(250),wk2(250)
    #common /wf/ wf2(250,14),wc(14),lc(14)
    ngrid = 250
    mc = 30
    pos1 = pos2 = int

    index[x] = 20.0*(log(x) + 8.8) + 2.0

    # first input channels
    atomic_f = open(atomic_file, 'r') # unit 4
    cluster_f = open(cluster_file, 'r') # unit 7

    # now output channels
    info_f = open(info_file, 'w') # unit 11
    mufftin_f = open(mufftin_file, 'w') # unit 9

    # initialisation of Loucks' exponential mesh
    x = -8.8

    for ix in range(ngrid+1):
        rx[ix] = exp(x)
        x += 0.05

    title = cluster_f.read_line() 
    info_f.write("%s\n" % title)

    ''' Input of crystallographic data:
        spa = lattice constant in a.u.
        rc(i,j) = i'th coordinate of the j'th axis of unit cell,
                  in units of spa
        rk(i,j) = i'th coordinate of the j'th atom in unit cell,
                  in units of spa
        nieq = number of inequivalent atoms in unit cell
    
        for an atom of type ir:
           nrr(ir) = number in unit cell
           z(ir) = atomic number
           zc(ir) = valence charge
           rmt(ir) = muffin-tin radius in bohr radii
    '''
    spa = [t(s) for s in zip(
                        (float, float, float), cluster_f.readline().split())]
    for i in range(1, 3+1):
        for j in range(1, 3+1):
            rc[i][j] = cluster_f.readline()  # read(7,101)((rc(i,j),i=1,3),j=1,3)

    for i in range(1, 3+1):
        for j in range(1, 3+1): 
            rc[i][j] *= spa

    # read number of inequivalent atoms in unit cell
    nieq = int(cluster_f.readline().split('!')[0])  # read(7,102) 

    # allocate arrays for nieq atoms
    z = zc = [0.] * nieq
    vh = sig = rho = vs = [[0.] * 250] * nieq
    vmad = [[0.] * 550] * nieq
    ad = [[0.]* 30 ] * nieq
    rmt = [0.] * nieq
    jrmt = [int] * nieq
    jrmt2 = [int] * nieq
    nrr = [int] * nieq
    ncon = [int] * nieq
    nx = [int] * nieq
    na = [[int] * 30] * nieq
    ia = [[int] * 30] * nieq
    name = [[int] * 8] * nieq

    # initialise values
    for ir in range(1, nieq+1):
        for i in range(1, ngrid+1):
            vh[i][ir] = 0.0
            vs[i][ir] = 0.0
            vmad[i][ir] = 0.0
            sig[i][ir] = 0.0
            rho[i][ir] = 0.0

    vhar = 0.0
    vex = 0.0

    jj = 0
    zz = 0.0
    inrr = 0

    # loop through each ineq atom
    for ir in range(1, nieq+1):
        # read name and split character string
        name[ir] = cluster_f.readline().split('!')[0]

        # read parameters num_atoms, z, valence, rydberg_mufftin_radius
        nrr[ir], z[ir], zc[ir], rmt[ir], = [t(s) for s in zip(
            (int, float, float, float), 
                cluster_f.readline().split('!')[0].split())]  # format 103

        # reallocate enough space for extra atoms
        if rk is None:
            rk = [[float] * 3] * nrr[ir]
        else:
            pass
            #allocate(trk(3,inrr+nrr(ir)))
            #trk(1:3,1:inrr) = rk
            #deallocate(rk)
            #call move_alloc(trk,rk)

        if zm is None:
           zm = [float] * nrr[ir]
        else:
           pass
           #allocate(tzm(inrr+nrr(ir)))
           #tzm(1:size(zm)) = zm
           #deallocate(zm)
           #call move_alloc(tzm,zm)
         
        inrr += nrr[ir]

        zz += abs(zc[ir])
        jrmt[ir] = index[rmt[ir]]
        n = nrr[ir]

        for j in range(1, n+1):
            jj += 1
            zm[jj] = zc[ir]
            for i in range(1, 3+1):
                rk[jj] = [t(s) for s in zip(
                          (float, float, float), 
                          cluster_f.read_line().split('!')[0].split())] 

            for i in range(1, 3+1):
                rk[i][jj] *= spa


        # n = total number of atoms in unit cell
        # av = total volume of unit cell
        # oma = atomic volume
        # rws = wigner-seitz radius
        n = jj
        rcc1 = rc[2][2]*rc[3][3] - rc[3][2]*rc[2][3]
        rcc2 = rc[3][2]*rc[1][3] - rc[1][2]*rc[3][3]
        rcc3 = rc[1][2]*rc[2][3] - rc[2][2]*rc[1][3]
        av = abs(rc[1][1]*rcc1 + rc[2][1]*rcc2 + rc[3][1]*rcc3)
        oma = av / float(n)
        rws = (0.75*oma / pi)**(1.0/3.0)
        jrws = index(rws)
        write(11,201)((rc(i,j),i=1,3),j=1,3)
        write(11,202)av,oma,rws
        jj = 0
        
        for ir in range(1, nieq+1):
            write(11,203)ir,(name(i,ir),i=1,4),nrr(ir)
            inr = nrr[ir]

            for iir in range(1, inr+1):
                jj += 1
                write(11,204)(rk(i,jj),i=1,3)

            write(11,205)z(ir),zc(ir),rmt(ir)

       info_f.write(11,216)(rx[ix],ix=1,ngrid)

       # for each atomic type, read in atomic wavefunctions for neutral
       # atom, in either herman-skillman or clementi form, producing?
       # rho = 4*pi*charge density * radius**2
       mix = 0

       for ir in range(1, nieq+1):
           read(4,100) wfn
           # option 0)  relativistic charge density input
           if wfn == wfn0: 
               rela(rho(1,ir),rx,nx(ir),ngrid)
           # option 1)  herman-skillman input
           if wfn == wfn1: 
               hsin(rho(1,ir),rx,nx(ir),ngrid)
           # option 2)  clementi input
           if wfn == wfn2: 
               clemin(rho(1,ir),rx,nx(ir),ngrid)
           # option 3)  potential input
           if wfn == wfn3: 
               pass  # goto 14

            # rho is normalised using total electronic charge on the atom
            # calculated by the trapezoidal rule
            nix = nx[ir]
            mix = max0[nix][mix]
            sum = 0.0
            w1 = 0.025*rho[1][ir]*rx[1]

            for ix in range(2, nix+1):
                w2 = 0.025*rho[ix][ir]*rx[ix]
                sum += w1 + w2
                w1 = w2

            ze = sum

            # solve poisson's equation:
            # sig = coulomb potential
            # rho = 4*pi*charge density*radius squared
            poison(rho(1,ir),z(ir),nix,sig(1,ir))

            x=-8.8

            for ix in range(1, nix+1):
                ce = exp(-0.5*x)
                sig[ix][ir] = ce*(-2.0*z[ir]*ce + sig[ix][ir])
                rho[ix][ir] /= (rx[ix]**2)
                x += 0.05

            write(11,206)(name(i,ir),i=1,4),ze,rx(nix),nix
            write(11,207)(sig(ix,ir),ix=1,nix)

        # details of neighbouring shells for each atomic type ir?
        # ncon(ir) = number of shells included
        # ia(j,ir) = atomic type in j'th shell
        # na(j,ir) = number of atoms in j'th shell
        # ad(j,ir) = distance to j'th shell
        rmax = rx[mix]

        nbr(ia,na,ad,ncon,nrr,nieq,rc,rk,n,rmax,mc)
        write(11,208)

        for ir in range(1, nieq+1):
            write(11,209)ir
            nc = ncon(ir)
            ic = (nc-1)/12+1
            kc = 0

            for i in range(1, ic):
                jc = kc + 1
                kc = min(nc, kc + 12)
                write(11,210)(ad(j,ir),j=jc,kc)
                write(11,211)(na(j,ir),j=jc,kc)
                write(11,212)(ia(j,ir),j=jc,kc)

        read(7,102) nform

        # calculation of the muffin-tin potential for each neutral
        # atom, following the mattheiss prescription
        # read in alpha for the slater exchange term
        read(7,101)alpha
        info_f.write("%f\n", alpha)
        pd = 6.0 / (pi*pi)

        for ir in range(1, nieq+1):
            jrx = max(jrws, jrmt[ir])
            # summing the potentials from neutral atoms
            # vh = hartree potential
            sumax(vh(1,ir),sig,rx,nx,ncon(ir),ia(1,ir),na(1,ir),
                  ad(1,ir),jrx,ngrid,nieq)

            # summing the charge density about each atomic type
            # vs = total charge density,: slater exchange term
            sumax(vs(1,ir),rho,rx,nx,ncon(ir),ia(1,ir),na(1,ir),
                  ad(1,ir),jrx,ngrid,nieq)

           for ix in range(1, jrx+1):
               vs[ix][ir] = -1.5*alpha*(pd*vs[ix][ir])**(1.0/3.0)

        # calculate the muffin-tin zero
        vint = 0.
        nh = int(cluster_f.read()) # (7,102) nh

        if nh == 0 and nieq == 1:
            mtzm(vh(1,1),vs(1,1),rx,ngrid,
                 rmt(1),rws,jrmt(1),jrws,vhar,vex)

        if nh != 0:
            mtz(sig,rho,rx,ngrid,rmt,nrr,nx,nieq,rc,rk,n,
                vhar,vex,alpha,av,nh)

        # slab or bulk calculation?
        # slab_flag = 1 for slab or 0 for bulk (default = 0)
        if slab_flag == 1:
            # input the mtz value from the bulk substrate calculation
            if (len(mtz_string) >= 1:
                esht = float(mtz_string)
            else:
                sys.stderr.write(
                        "error: mtz input is invalid for slab calculation\n")
                close(4)
                cluster_f.close()
                close(9)
                close(10)
                info_f.close()
           
                return

            esh = esht - (vhar + vex)

        else:
            # check slab_flag is valid
            if (slab_flag != 0):
                sys.stder.write("warning: slag_flag defaulted to 0\n")

        # if you are interested in adatoms on this substrate
        # rerun a slab calculation with the adatoms
        # and use this mtz value as input when asked
        print (vhar + vex)
        cavpot = vhar + vex  # return value of function

        goto 16

        # option 3)  read in potential of neutral atom, vh, on radial
        #            grid, rx, for correction by madelung summation
14      read(4,104) ngrid,(rx(ix),ix=1,ngrid)

        for ir in range(1, nieq+1):
            read(4,104)jrx,(vh(ix,ir),ix=1,jrx)
            jrmt[ir] = jrx

        # the madelung correction for ionic materials. subroutine mad
        # computes the spherically and spatially averaged fields for
        # the lattice of point charges about each atomic type
16      if zz !=0: 
            mad(vmad,rx,ngrid,rmt,nrr,jrmt,nieq,rc,rk,zm,n,av)

        # the total muffin-tin potential is accumulated into sig,
        # referred to the muffin-tin zero
        vint = vhar + vex
        if nform == 0: 
            mufftin_f.write("%4i\n" % nieq)

        for ir in range(1, nieq+1):
            write(11,213) (name(i,ir),i=1,4), vint, rmt(ir)
            jrx = jrmt[ir]

            for ix in range(1, jrx+1):
                vh[ix][ir] = vh[ix][ir] - vhar
                vs[ix][ir] = vs[ix][ir] - vex
                sig[ix][ir] = vh[ix][ir] + vs[ix][ir] + vmad[ix][ir]
                info_f.write("%12.5f%20.4e%20.4e%20.4e%20.4e\n" % (
                             rx[ix], vh[ix][ir], vs[ix][ir], 
                             vmad[ix][ir], sig[ix][ir]))  # format 214

        # write output in a format to be read by williams phase shift
        # program (nform=1), by cavleed phase shift program (nform=0), or
        # by the relativistic phase shift program (nform=2)

        # also prepare to shift the potential by an amount of the order
        # of the bulk muffintin zero.
        # this is needed only if the cluster.i file corresponds to a
        # surface adsorbate
        if nform == 1:
            mufftin_f.write(" &nl2 nrr=%2i &end\n" % nieq)
        elif nform == 2:
            # define german grid rx and save old grid in rs
            rm = 60.0
            dx = 0.03125
            nmx = 421
            rs[1] = rx[1]
            rx[1] = rm*exp(dx*(1 - nmx))
            j = 1
            rm = exp(dx)

            while j < nmx:
                k = j+1
                rs(k) = rx(k)
                rx(k) = rm*rx(j)
                j = k

        for ir in range(1, nieq+1):
            jrx = jrmt(ir)
            if nform == 0:
                mufftin_f.write("%s\n" % name[ir])
                mufftin_f.write(9,218) z[ir], rmt[ir], vint
            elif nform == 1:
                mufftin_f.write(9,221) z[ir], rmt[ir]
            else:
                # es=emin for phase shift calculation (ev)
                # de=delta e for phase shift calculation (ev)
                # ue=emax for phase shift calculation (ev)
                # lsm=maximum number of phase shifts desired
                es = 20.
                de = 5.
                ue = 300.
                lsm = 12
                mufftin_f.write("%s" % name[ir])
                mufftin_f.write("%12.4f%12.4f%12.4f    %3i    %12.4f" %(
                                                        es, de, ue, lsm, vint))

                # interpolation to grid rx
                for k in range(1, jrx+1):
                    sig[k][ir] = (sig[k][ir] - esh)*rs[k]

                    nmxx = nmx
                    chgrid(sig[1][ir], rs, jrx, pot, rx, nmxx)
                    iz = z[ir]
                    mufftin_f.write("%4i%10.6f%4i" % (iz, rmt[ir], nmxx))
                    jrx = nmxx

            if nform == 0:
                mufftin_f.write(9,102)jrx
            elif nform == 1:
                for ix in range(1, jrx+1):
                    mufftin_f.write(9,219) rx[ix], rx[ix]*(sig[ix][ir] - esh)
                rneg = -1.
                mufftin_f.write(9,219) rneg
            if (nform==0):
                for ix in range(1, jrx+1):
                    mufftin_f.write(9,219)rx(ix),(sig(ix,ir)-esh)
            else:
                mufftin_f.write("%14.7e%14.7e%14.7e%14.7e%14.7e", pot)


        # close file handles
        close(4)
        close(7)
        close(9)
        close(10)
        close(11)

        return

101    format(3f8.4)
102    format(i4)
103    format(i4,3f8.4)
104    format(i4/(5e14.5))
200    format("muffin-tin potential program",5x,20a4)
201    format(///" axes of unit cell"/(6x,3f8.4))
202    format(" unit cell volume",f15.4/" atomic volume",f18.4/
     + " wigner-seitz radius",f12.4)
203    format(///" type",i2," atom",2x,4a4/i4," atoms in unit cell")
204    format(6x,3f8.4)
205    format(" atomic number",f15.1/"valence",f21.1/
     + "muffin-tin radius",f14.4)
206    format(///1x,4a4," electronic charge",f12.5/
     + "coulomb potential for isolated atom, out to radius",
     + f12.5,10x,"nx",i4/)
207    format(5(10e12.4/))
208    format("1")
209    format(//"nearest neighbour shells for type",i2," atom")
210    format(" distance",1x,15(f8.4))
211    format(" number",3x,15(i5,3x))
212    format(" type",5x,15(i5,3x))
213    format("1",4a4,5x,
     + "potentials in rydbergs correct to muffin-tin zero",f8.4,
     + "muffin-tin radius",f8.4//5x,"radius",5x,"hartree potential",9x,
     + "exchange",4x,"madelung correction",5x,"total potential")
214    format(f12.5,4e20.6)
215    format(///"statistical exchange parameter, alpha",f10.4)
216    format(///"loucks' radial mesh"//5(10f11.5/))
217    format(4a4)
218    format(3f8.4)
219    format(2e14.5)
220    format(" &nl2 nrr=",i2," &end")
221    format(" &nl16 z=",f7.4,",rt=",f7.4," &end")

      end function

def chgrid(fx,x,nx,fy,y,ny):
    ''' 
    Performs piecewise quadratic interpolation from grid x to grid y,  by
    Aitken's divided difference scheme.   nx,ny are array dimensions?
    note that ny is reset in chgrid
    '''
    iy = 1

    for ix in range(3, nx+1):  # label 2 
        while True:
            if iy > ny:
                break
            
            yy = y[iy]
            
            if yy > x[ix]:
                continue
            
            a1 = x[ix-2] - yy
            a2 = x[ix-1] - yy
            a3 = x[ix] - yy
            a12 = (fx[ix-2]*a2 - fx[ix-1]*a1) / (x[ix-1 - x[ix-2])
            a13 = (fx[ix-2]*a3 - fx[ix]*a1) / (x[ix] - x[ix-2])
            fy[iy] = (a12*a3 - a13*a2) / (x[ix] - x[ix-1])
            
            if iy > ny:
                break
           
            iy += 1

    ny = iy - 1

def clemin(rho,rx,nx,ngrid):
    '''
    routine for input of wavefunctions in the clementi parametrised
    form, and calculation of charge density on radial mesh rx
        rho = 4*pi*sum over states of (modulus(wave fn)**2) *
             radius**2
        nc = number of atomic states
    for each atomic state i?
        lc(i) = angular momentum
        frac = fractional occupation
        wc(i) = number of electrons
        wfc(ix,i) = wavefunction x radius at grid point ix
    '''
    # dimension rho(ngrid),rx(ngrid)
    # real name(4),sum
    # common /wk/ ex(20),fac(20),fnt(20),nt(20)
    # common /wf/ wfc(250,14),wc(14),lc(14)
    
    name = atomic_f.readline().split('!')[0]
    iprint = int(atomic_f.read(4).split('!')[0])
    nc = int(atomic_f.read(4).split('!')[0])
    
    # initialize wfc
    for ic in range(1, nc+1):
        for ig in range(1, ngrid+1):
            wfc[ig][ic] = 0.0
    
    # input of clementi parameters
    ic = 1
    
    while True:
        ns = int(atomic_f.readline())  # outer loop
        if ns <= 0:
            break # outer loop
        
        for i in range(1, ns+1):
            nt[i], ex[i] = [t(s) for s in zip((int, float), atomic_f.readline().split())] 
        
        for j in range(1, ns+1):
            a = 1.0
            b = 2.0
            k = nt[j]
            c = float(k)
            kd = k + k
        
            for i in range(2, kd+1):
                a *= b
                b += 1.0
        
            fnt[j] = exp(-0.5*log(a) + (c + 0.5)*log(2.0*ex[j]))
        
        while True:
            lc[ic] = atomic_f.read(4,101)  # inner loop
            if lc[ic] < 0:
                break # inner loop
                
            fac = [t(s) for s in zip(
                    [float] * ns, atomic_f.readline().split())]
            frac = float(atomic_f.readline())
            wc[ic] = 2.0*float(2*lc[ic] + 1)*frac
            
            for ix in range(1, ngrid+1):
                sum = 0.0
            
                for k in range(1, ns+1):
                    exx = ex[k] * rx[ix]
                    if exx > 80.0:
                        break  # inner loop
                    sum += fac[k]*fnt[k]*(rx[ix]**(nt[k]))*exp(-exx)
            
                wfc[ix][ic] = sum
            
            ic += 1
        
        # outer loop
        if lc[ic] < 0:
            continue
        if exx > 80.:
            break
    
    # calculation of charge density
    for ix in range(1, ngrid+1):
        sum = 0.0
        for ic in range(1, nc+1):
            sum += wc[ic]*wfc[ix][ic]*wfc[ix][ic]
        rho[ix] = sum
    
        if sum < 1.0e-9:
            break  # for loop
    
    nx = ix
    
    if iprint is False:
        return
    
    info_f.write("%s" % name) # fmt 200
    
    for ic in range(1, nc+1):
        info_f.write("l%3i\n" % lc(ic))
        for ix in range(1, ngrid):
            info_f.write("%11.5f", % wfc[ix][ic])
        info_f.write("\n")
        
    info_f.write("charge density out to radius%12.5f%snx%4i" 
                 % (rx(nx), ' ' * 10, nx))
    for ix in range(1, nx+1):
        info_f.write("%12.4e" % rho(ix))
    info_f.write("\n")
    
#---------------------------------------------------------------------
#  routine for input of atomic wavefunctions from herman-skillman
#  tables, and calculation of charge density on the radial mesh rx
#    rho = 4*pi*sum over states of (modulus(wave fn)**2) *
#          radius**2
#    nm ? h-s grid interval doubles every nm mesh points
#    nc = number of atomic states
# for each atomic state i?
#    lc(i) = angular momentum
#    frac = fractional occupation
#    wc(i) = number of electrons
#    wfc(ix,i) = wavefunction x radius at grid point ix
def hsin(rho,rx,nx,ngrid):

    # dimension rho(ngrid),rx(ngrid)
    # real name(4),sum
    # common /wk/ rr(250),rs(250)
    # common /wf/ wfc(250,14),wc(14),lc(14)
    
    name, z = [t(s) for s in zip(str, float), atomic_f.readline().split())]
    read(4,101)iprint
    read(4,101)nm
    read(4,101)nc
    
    for ig in range(1, 250+1):
        rs(ig) = 0.0
    
        for ic in range(1, nc+1):
            wfc[ig][ic] = 0.0
    
    # initialisation of Herman-Skillman mesh
    dr = 0.005*0.88534138 / exp(log(z) / 3.0)
    rr[1] = 0.0
    
    for i in range(2, 250+1):
        if (i % nm == 2:
            dr += dr
        rr[i] = rr[i-1] + dr
    
    ns = 0
    
    for ic in range(1, nc+1):
        read(4,101)lc(ic),n,frac
        ns = max0(ns,n)
        wc(ic) = 2.0*float(2*lc(ic)+1)*frac
        read(4,102)(wfc(ix,ic),ix=1,n)
    
    # calculation of charge density
    for ix in range(1, ns+1):
        sum = 0.0
    
        for ic in range(1, nc+1):)
            sum += wc[ic]*wfc[ix][ic]*wfc[ix][ic]
    
        rs[ix] = sum
    
    # interpolation to grid rx
    nx = ngrid
    chgrid(rs,rr,ns,rho,rx,nx)
    if iprint is False:
        return
    info_f.write(200) name,(rr(ix),ix=1,ns)
    
    for ic in range(1, nc+1):
        write(11,201)lc(ic),(wfc(ix,ic),ix=1,ns)
    
    ix = 1
    while (ix < nx and rho[ix] > 1.0e-9):
        ix += 1
    
    nx = ix
    write(11,202)rx(nx),nx,(rho(ix),ix=1,nx)
    
    return
    
    100    format(4a4/f9.4)
    101    format(2i4,f9.4)
    102    format(5f9.4)
    200    format("1",4a4," atomic wavefunctions (herman-skillman)",
    + " x radius"//" herman-skillman mesh"//5(10f12.5/))
    201    format("l",i3//5(10f11.5/))
    202    format("charge density out to radius",f12.5,10x,
    + "nx",i4//5(10e12.4/))
    
    return

def mad(vmad,rx,ngrid,rmt,nrr,nx,nr,rc,rk,zm,n,av):
    '''
    Calculates the spherically and spatially averaged
    fields from a lattice of point charges, and tabulates them on
    a radial mesh rx, about each atomic type in the unit cell
    ** nb this routine works in hartrees, but converts to rydbergs **
    rc(i,j) = the i'th coordinate of the j'th axis of the unit cell
    rk(i,j) = the i'th coordinate of the j'th atom in the unit cell
    vmad(j,ir) = the j'th tabulated value of the spherically averaged
    potential about a type-ir atom
    zm(k)=charge on the k'th atom
    rmt(ir) = muffin-tin radius of a type-ir atom
    nr = number of inequivalent atoms in the cell
    av = volume of unit cell
    g(i,j) = i'th coordinate of the j'th reciprocal lattice vector
    vmm(ir) = the integral of the potential about a type-ir atom
    out to the muffin-tin radius
    '''
       dimension vmad(ngrid,nr),rx(ngrid),rc(3,3),rk(3,n),zm(n)
       dimension rmt(nr),nrr(nr),nx(nr)
       common /wk/ g(3,3),vmm(5),fr(5),ra(3),ga(3)
       data pi,test/3.1415926536,1.0e-4/

       rad(a1,a2,a3)=sqrt(a1*a1+a2*a2+a3*a3)


       do ir=1,nr
         fr(ir)=0.0
         do j=1,ngrid
           vmad(j,ir)=0.0
         end do
       end do

      # the reciprocal lattice is defined by three vectors, g
       atv=2.0*pi/av
       g(1,1)=(rc(2,1)*rc(3,2)-rc(3,1)*rc(2,2))*atv
       g(2,1)=(rc(3,1)*rc(1,2)-rc(1,1)*rc(3,2))*atv
       g(3,1)=(rc(1,1)*rc(2,2)-rc(2,1)*rc(1,2))*atv

       g(1,2)=(rc(2,2)*rc(3,3)-rc(3,2)*rc(2,3))*atv
       g(2,2)=(rc(3,2)*rc(1,3)-rc(1,2)*rc(3,3))*atv
       g(3,2)=(rc(1,2)*rc(2,3)-rc(2,2)*rc(1,3))*atv

       g(1,3)=(rc(2,3)*rc(3,1)-rc(3,3)*rc(2,1))*atv
       g(2,3)=(rc(3,3)*rc(1,1)-rc(1,3)*rc(3,1))*atv
       g(3,3)=(rc(1,3)*rc(2,1)-rc(2,3)*rc(1,1))*atv

      # maximum value of rk, and minimum values of rc,g - prior to
      # choosing the separation constant al and limits for summations
       rkmax=0.0

       do j=1,n
         rkmax=max(rkmax,rad(rk(1,j),rk(2,j),rk(3,j)))
       end do

       rcmin=1.0e6
       gmin=1.0e6

       do j=1,3
         rcmin=min(rcmin,rad(rc(1,j),rc(2,j),rc(3,j)))
         gmin=min(gmin,rad(g(1,j),g(2,j),g(3,j)))
       end do

      # al is chosen to give equal numbers of terms in real and
      # reciprocal space summations
       fac1=test*log(test)**4
       fac2=(4.0*pi*rcmin**4)/(av*gmin**4)
       al=exp(log(fac1/fac2)/6.0)
       itr=1+ifix((al*rkmax-log(test))/(al*rcmin))
       limr=itr+itr+1
       fac1=4.0*pi*al*al/(av*gmin**4)
       itg=1+ifix(exp(log(fac1/test)/4.0))
       limg=itg+itg+1
       write(11,200)((g(i,j),i=1,3),j=1,3)
       write(11,201)rcmin,gmin,rkmax,test,al

      # real space summation
       write(11,202)itr
      # the prefactors fr from the real space summation are calculated
       as=-float(itr)-1.0
       ax=as

       do 5 jx=1,limr
         ax=ax+1.0
         ay=as

         do 5 jy=1,limr
           ay=ay+1.0
           az=as

           do 5 jz=1,limr
             az=az+1.0

             do i=1,3
               ra(i)=ax*rc(i,1)+ay*rc(i,2)+az*rc(i,3)
             end do

             do 5 j=1,n
               k=1

               do 5 kr=1,nr
                 r=rad(ra(1)+rk(1,j)-rk(1,k),ra(2)+rk(2,j)-rk(2,k),
     +             ra(3)+rk(3,j)-rk(3,k))
                 if (r<1.0e-4) goto 5
                 fr(kr)=fr(kr)+zm(j)*exp(-al*r)/r
5                k=k+nrr(kr)

       k=1

       do kr=1,nr
         x=rmt(kr)
         a=exp(-al*x)
         ai1=((1.0-a)/al-x*a)/al
         ai2=(x*0.5*(1.0/a+a)-0.5*(1.0/a-a)/al)/al/al
         vmm(kr)=4.0*pi*(zm(k)*ai1+ai2*fr(kr))
         nix=nx(kr)

         do j=1,nix
           x=rx(j)
           a=exp(al*x)
           vmad(j,kr)=fr(kr)*0.5*(a-1.0/a)/(al*x)+zm(k)/(a*x)
         end do

         k=k+nrr(kr)

       end do

       write(11,203)(vmm(kr),kr=1,nr)

      # next comes the summation in reciprocal space
       write(11,204)itg
       as=-float(itg)-1.0
       ax=as

       do 13 jx=1,limg
         ax=ax+1.0
         ay=as

         do 13 jy=1,limg
           ay=ay+1.0
           az=as

           do 13 jz=1,limg
             az=az+1.0

             do i=1,3
               ga(i)=ax*g(i,1)+ay*g(i,2)+az*g(i,3)
             end do

             gm=rad(ga(1),ga(2),ga(3))
             gs=gm*gm
             fac1=0.0
             if (gs<1.0e-4) goto 13
             fac1=4.0*pi*al*al/(av*gs*(gs+al*al))
             k=1

             do kr=1,nr
               fac2=0.0

               do j=1,n
                 gr=0.0

                 do i=1,3
                   gr=gr+ga(i)*(rk(i,k)-rk(i,j))
                 end do

                 fac2=fac2+cos(gr)*zm(j)
               end do

               x=rmt(kr)
               ai3=(sin(gm*x)/gm-x*cos(gm*x))/gs
               vmm(kr)=vmm(kr)+4.0*pi*ai3*fac1*fac2
               nix=nx(kr)

               do i=1,nix
                 x=rx(i)
                 vmad(i,kr)=vmad(i,kr)+fac1*fac2*sin(gm*x)/(gm*x)
               end do

               k=k+nrr(kr)
             end do

13     continue

       write(11,203)(vmm(kr),kr=1,nr)

      # refer to muffin-tin zero
       vm=0.0
       amt=0.0

       do ir=1,nr
         vm=vm+float(nrr(ir))*rmt(ir)**3
         amt=amt+float(nrr(ir))*vmm(ir)
       end do

       amt=amt/(av-4.0*pi*vm/3.0)

      # express the final potential in rydbergs
       amt=-2.0*amt
       write(11,205)amt

       do kr=1,nr
         nix=nx(kr)
         do j=1,nix
           vmad(j,kr)=2.0*vmad(j,kr)-amt
         end do
       end do

       return

200    format(///"madelung correction"//"reciprocal lattice"/(6x,3f8.4))
201    format("rcmin",f10.4,10x,"gmin",f10.4,10x,"rkmax",f10.4,
     + 10x,"test",e12.4/" separation constant",e12.4)
202    format("real space summation",11x,"itr",i3)
203    format(" vmm (hartrees) ",5e12.4)
204    format("reciprocal space summation",5x,"itg",i3)
205    format("madelung muffin-tin zero",5e12.4)

      return

#---------------------------------------------------------------------
#  subroutine for calculation of the muffin-tin zero level?
#  the average value of the potential between the muffin-tin
#  spheres in the unit cell
def mtz(sig,rho,rx,ngrid,rmt,nrr,nx,nr,rc,rk,n,vhar,vex,alpha,av,nh):

       dimension sig(ngrid,nr),rho(ngrid,nr),rx(ngrid),rmt(nr),
     + nrr(nr),nx(nr),vg(20),rc(3,3),rk(3,n)
       common /wk/ x(3),rb(3)
       data pi,ng/3.14159265358,20/

      # grid reference for radius on loucks' mesh
       index(y)=20.0*(log(y)+8.8)+1.0
       rad(a1,a2,a3)=sqrt(a1*a1+a2*a2+a3*a3)

       pd=6.0/pi/pi

       do ig=1,ng
         vg(ig)=0.0
       end do

       ig=0
       vhar=0.0
       vex=0.0
       npoint=0
       nint=0
       dh=1.0/float(nh)

1      ah=dh/2.0

       ax=-ah

       do 7 ix=1,nh
         ax=ax+dh
         ay=-ah

         do 7 iy=1,nh
           ay=ay+dh
           az=-ah

           do 7 iz=1,nh
             az=az+dh

             do i=1,3
               x(i)=ax*rc(i,1)+ay*rc(i,2)+az*rc(i,3)
             end do

             npoint=npoint+1
             # gives sample point x inside the unit cell - test whether
             # interstitial
             bx=-1.0

             do 4 jx=1,2
               bx=bx+1.0
               by=-1.0

               do 4 jy=1,2
                 by=by+1.0
                 bz=-1.0

                 do 4 jz=1,2
                   bz=bz+1.0

                   do i=1,3
                     rb(i)=x(i)-bx*rc(i,1)-by*rc(i,2)-bz*rc(i,3)
                   end do

                   i=0

                   do 4 ir=1,nr
                     inr=nrr(ir)

                     do 4 iir=1,inr
                       i=i+1
                       xr=rad(rb(1)-rk(1,i),rb(2)-rk(2,i),rb(3)-rk(3,i))
                       if (xr<rmt(ir)) goto 7

4                      continue

             # we have an interstitial point
             nint=nint+1
             # sum coulomb and exchange energies from atoms within 2 unit
             # cells around this point
             sumc=0.0
             sume=0.0
             bx=-3.0

             do 6 jx=1,5
               bx=bx+1.0
               by=-3.0

               do 6 jy=1,5
                 by=by+1.0
                 bz=-3.0

                 do 6 jz=1,5
                   bz=bz+1.0

                   do i=1,3
                     rb(i)=bx*rc(i,1)+by*rc(i,2)+bz*rc(i,3)-x(i)
                   end do

                   j=0

                   do 6 jr=1,nr
                     jnr=nrr(jr)

                     do 6 jjr=1,jnr
                       j=j+1
                       xr=rad(rb(1)+rk(1,j),rb(2)+rk(2,j),rb(3)+rk(3,j))
                       j2=index(xr)
                       if (j2>=nx(jr)) goto 6
                       j1=j2-1
                       j3=j2+1
                       x1=rx(j1)
                       x2=rx(j2)
                       x3=rx(j3)
                       termc=(xr-x2)*(xr-x3)/(x1-x2)/(x1-x3)*sig(j1,jr)
     +                      +(xr-x1)*(xr-x3)/(x2-x1)/(x2-x3)*sig(j2,jr)
     +                      +(xr-x2)*(xr-x1)/(x3-x2)/(x3-x1)*sig(j3,jr)
                       terme=(xr-x2)*(xr-x3)/(x1-x2)/(x1-x3)*rho(j1,jr)
     +                      +(xr-x1)*(xr-x3)/(x2-x1)/(x2-x3)*rho(j2,jr)
     +                      +(xr-x2)*(xr-x1)/(x3-x2)/(x3-x1)*rho(j3,jr)
                        sumc=sumc+termc
                        sume=sume+terme

6                       continue


           if (sume<=1.e-8):
             sume=.0
           else:
             sume=-1.5*alpha*(pd*sume)**(1./3.)
           endif

           vhar=vhar+sumc
           vex=vex+sume
           jg=mod(ig,20)+1
           vg(jg)=vg(jg)+sumc+sume
           ig=ig+1

7          continue

       dh=dh/2.0
       nh=nh+nh
       if (nint==0) goto 1

       ant=float(nint)
       vhar=vhar/ant
       vex=vex/ant
       vint=vhar+vex

      # estimate standard deviation
       if (nint < ng) ng=nint
       nag=nint/ng
       ag=float(nag)

       do ig=1,ng
         vg(ig)=vg(ig)/ag
       end do

       var=0.0

       do ig=1,ng
         var=var+(vint-vg(ig))**2
       end do

       var=sqrt(var/float(ng*(ng-1)))
      # the current monte-carlo volume for the interstitial region
      # is volc
       volc=ant/float(npoint)*av
      # volt is the true volume of the region between muffin-tin
      # spheres in the unit cell
       vm=0.0

       do ir=1,nr
         vm=vm+float(nrr(ir))*rmt(ir)**3
       end do

       volt=av-4.0*pi*vm/3.0

       write(11,200)nint,npoint,ng,nag,volt,volc
       write(11,201)vhar,vex,vint,var

       return

200    format(///"muffin-tin zero calculation, sampling with",i6,
     + " points from grid of",i6/" variance estimated from",
     + i4," groups of",i5//" true volume of interstitial region",
     + f11.4,5x,"monte-carlo volume",11x,f9.4)
201    format(" average hartree potential",6x,f14.5,5x,
     + "average exchange potential",f12.5/
     + "muffin-tin zero",f12.5,10x,"standard deviation",f12.5)

      return

#---------------------------------------------------------------------
#  subroutine for calculation of the muffin-tin zero level for
#  monoatomic crystals, using a spherical average of the potential
#  between muffin-tin radius and wigner-seitz radius, as in eq 3.31
#  of loucks, transformed to the exponential grid rx?
#                       rx(i)=exp(-8.8+0.05(i-1))
#  integration by trapezium rule.  jrmt,jrws are grid points outside
#  muffin-tin radius and wigner-seitz radius respectively
def mtzm(vh,vs,rx,ngrid,rmt,rws,jrmt,jrws,vhar,vex):

    # dimension vh(ngrid),vs(ngrid),rx(ngrid)
    # double precision sumh,sume
    
    dx = 0.05
    ddx = 0.5*dx
    dxx = exp(3.*dx)
    x = log(rx[jrmt] / rmt)
    rdx = x/dx
    xx = rx[jrmt-1]**3
    xxmt = xx*dxx
    sumh = 0.5*x*(rdx*xx*vh[jrmt-1] + (2.-rdx)*xxmt*vh[jrmt])
    sume = 0.5*x*(rdx*xx*vs[jrmt-1] + (2.-rdx)*xxmt*vs[jrmt])
    xx = xxmt
    jrw = jrws - 1
    
    if (jrmt != jrw):
      vh1 = ddx*xx*vh[jrmt]
      vx1 = ddx*xx*vs[jrmt]
      jrm = jrmt + 1
    
      for j in range(jrm, jrw+1):
          xx *= dxx
          vh2 = ddx*xx*vh[j]
          vx2 = ddx*xx*vs[j]
          sumh += vh1 + vh2
          sume += vx1 + vx2
          vh1=vh2
          vx1=vx2
    
    x = log(rws/rx[jrw])
    
    rdx = x/dx
    xxws = xx*dxx
    sumh += 0.5*x*((2.-rdx)*xx*vh(jrw)+rdx*xxws*vh(jrws))
    sume += 0.5*x*((2.-rdx)*xx*vs(jrw)+rdx*xxws*vs(jrws))
    c = 3. / (rws*rws*rws - rmt*rmt*rmt)
    vhar = c*sumh
    vex = c*sume
    vint = vhar + vex
    
    info_f.write("\n,uffin-tin zero by spherical average\n"
                 " average Hartree potential      %14.5f     "
                 "average exchange potential%12.5f\n"
                 "muffin-tin zero%12.5f\n" % (vhar, vex, vint) )

#---------------------------------------------------------------------
#  routine to supply nearest neighbour data for atoms in
#  a crystal structure, given?
#  rc(i,j)? the i'th coordinate of the j'th axis of the unit cell
#  rk(i,j)? the i'th coordinate of the j'th atom in the unit cell
#  nrr(ir)? the number of type-ir atoms in the unit cell
#  the information returned, for a type-ir atom, is
#  ncon(ir)? the number of nearest neighbour shells of a type-ir
#  atom included, out to a distance of rmax, but <= mc
#  ia(j,ir)? the type of atoms in the j'th neighbouring shell
#  na(j,ir)? the number of atoms in the j'th shell
#  ad(j,ir)? the radius of the j'th shell
def nbr(ia,na,ad,ncon,nrr,nr,rc,rk,n,rmax,mc):

    # dimension ia(mc,nr),na(mc,nr),ad(mc,nr),ncon(nr),nrr(nr), rc(3,3),rk(3,n)
    # common /wk/ rj(3)
    
    # initialisation
    rad[a1][a2][a3] = sqrt(a1*a1 + a2*a2 + a3*a3)
    
    rcmin = 1.0e6
    
    for i in range(1, 3+1):
        rcmin = min(rcmin, rad(rc[1][i], rc[2][i], rc[3][i]))
    
    for ir in range(1, nr+1):
        for ic in range(1, mc+1):
            ia[ic][ir] = 0
            na[ic][ir] = 0
            ad[ic][ir] = 1.0e6
    
    # search over adjacent unit cells to include mc nearest neighbours
    itc = ifix(rmax/rcmin) + 1
    limc = itc + itc + 1
    as = -float(itc + 1)
    ax = as
    
    for jx in range(1, limc+1):  # label 10
        ax += 1.0
        ay = as
    
        for jy in range(1, limc+1):
            ay += 1.0
            az = as
    
        do 10 jz=1,limc
          az=az+1.0
    
          do j=1,3
            rj(j)=ax*rc(j,1)+ay*rc(j,2)+az*rc(j,3)
          end do
    
          # rj is current unit cell origin.
          # for each atom in this unit cell find displacement r
          # from kr-type atom in basic unit cell
          j=0
    
          do 10 jr=1,nr
            jnr=nrr(jr)
    
            do 10 jjr=1,jnr
              j=j+1
              k=1
    
              do kr=1,nr
                r=rad(rj(1)+rk(1,j)-rk(1,k),rj(2)+rk(2,j)-rk(2,k),
    +               rj(3)+rk(3,j)-rk(3,k))
                if (r > rmax) break  # loop
                # compare r with nearest neighbour distances
                # already found
                ic=0
    4                  ic=ic+1
                if (ic>mc) break  # loop
                dr=r-ad(ic,kr)
                if (abs(dr) < 1.0e-4) dr=0.0
                if (dr) 6,5,4
    5                  if (ia(ic,kr) !=jr) goto 4
                na(ic,kr)=na(ic,kr)+1
                break  # loop
    6                  if (ic == mc) goto 8
                iic=ic+1
    
                do jjc=iic,mc
                  jc=mc+iic-jjc
                  ia(jc,kr)=ia(jc-1,kr)
                  na(jc,kr)=na(jc-1,kr)
                  ad(jc,kr)=ad(jc-1,kr)
                end do
    
    8                  ia(ic,kr)=jr
    
                na(ic,kr)=1
                ad(ic,kr)=r
                k=k+nrr(kr)
              end do
    10     continue
    
    for ir in range(1, nr+1):
        ncon[ir] = 0
    
        for ic in range(1, mc+1)
            if na[ic][ir] == 0:
                continue # goto 12
            
            ncon[ir] += 1
    
    return

def poison(psq,z,j,w):
    '''Taken from Loucks' book, appendix 1 '''

    # dimension psq(j),w(j)
    # double precision e(250),f(250),acc,a,b,c,d,c2
    a = 1.0 - 0.0025/48.0
    
    # eq. a1.11
    b = -2.0 - 0.025/48.0
    
    # eq. a1.12
    c = 0.0025 / 6.0
    d = exp(0.025)
    c2 = -b/a
    e[1] = 0.0
    
    # eq. a1.29
    f[1] = d
    
    # eq.a1.30
    x = -8.75
    j1 = j - 1
    
    for i in range(2, j1+1):
        acc = c*exp(0.5*x)*(d*psq[i+1] + 10.0*psq[i] + psq[i-1]/d)
        # eqs. a1.13, a1.6
        f[i] = c2 - 1.0/f[i-1]
        # eq. a1.20
        e[i] = (acc/a + e[i-1])/f(i)
        # eq. a1.21
        x +=0.05
    
    w[j] = 2.0*z*exp(-0.5*x)
    acc = w[j]
    
    # eq.a1.15
    for i in range(1, j1+1):
        jc = j - i
        acc = e[jc] + acc/f[jc]
        w[jc] = acc

def rela(rho,rx,nx,ngrid):
    '''
    Reads input of charge density from relativistic orbitals
    (Eric Shirley program), and calculation of charge density on
    the radial mesh rx
      rho = 4*pi*sum over states of (modulus(wave fn)**2) *
            radius**2
      rmin= minimum radial coordinate defining the logarithmic mesh used
            in relativistic calculation
      rmax= maximum radial coordinate defining the logarithmic mesh used
            in relativistic calculation
      nr  = number of points in the mesh
      
    the mesh is defined as r(i)=rmin*(rmax/rmin)**(float(i)/float(nr))
    for each atomic state i?
    '''    
    # dimension rho(ngrid),rx(ngrid)
    # real name(4),sum
    # common /wk/ rr(2000),rs(2000)
    
    name, iprint = [t(s) for s in zip((str, int), read(4,100).split('!')[0].split())] 
    rmin, rmax, nr, z = [t(s) for s in zip((float, float, int, float), read(4,54))] 
    # 54     format (d15.8,d15.8,i5,f5.2)
    
    # initialization of logarithmic grid
    for i in range(1, nr+1):
        rr[i] = rmin*(rmax/rmin)**(float(i)/float(nr))
    
    ns = nr
    # read in charge density
    read(4,56) (rs(j),j=1,nr)  # 56     format (f15.10)
    
    # interpolation to grid rx
    nx = ngrid
    chgrid(rs,rr,ns,rho,rx,nx)
    if iprint is False:
        return
    print("1%s relat. wavefunctions (Eric Shirley) r radius logarithmic mesh",200) name,(rr(ix),ix=1,ns) 
    
    for ix in range(1, nx+1):
        if rho[ix] < 1.0e-9: break  # loop
    
    nx = ix
    
    print("charge density out to radius%12,5f%snx%4i"202)rx(nx),nx,(rho(ix),ix=1,nx)

def sumax(acc,chi,rx,nx,ncon,ia,na,ad,imax,ngrid,nr):
    ''' 
    Performs the summation of contributions from
    neighbouring atoms (eq. 3.22,3.26,3.28).  integration by
    trapezium rule on radial grid  rx(i)=exp(-8.8+0.05(i-1))
    '''
    # dimension acc(ngrid),chi(ngrid,nr),rx(ngrid),nx(nr)
    # dimension ia(ncon),na(ncon),ad(ncon)
    # double precision sum

    index(x) = 20.*(log(x) + 8.8) + 2.
    dx = 0.05
    ddx = 0.5*dx
    dxx = exp(2.*dx)
    ic = ia[1]
    
    for i in range(1, imax+1):
        acc[i] = chi[i][ic]
    
    while True:
        for ja in range(2, ncon):  # label 4
            ic = ia[ja]
            nix = nx[ic]
    
            for i in range(1, imax+1):
                sum = 0.0
                x1 = abs(rx[i] - ad[ja])
                ix1 = index(x1)
    
                if ix1 > nix:
                    continue # goto 4
    
                dx1 = log(rx(ix1)/x1)
                rdx1 = dx1/dx
                x2 = min((rx[i] + ad[ja]), rx[nix])
                ix2 = min(index[x2],nix)
                dx2 = log(rx[ix2]/x2)
                rdx2 = dx2/dx
                xx = rx[ix2-1]**2
                xx1 = xx*dxx
            
                if ix1 == ix2:
                    break # goto 3
            
                sum += 0.5*dx2*((2.-rdx2)*xx*chi(ix2-1,ic)+rdx2*xx1*chi(ix2,ic))
                xx = rx[ix1-1]**2
                xx1 = xx*dxx
                sum += 0.5*dx1*(rdx1*xx*chi(ix1-1,ic)+(2.-rdx1)*xx1*chi(ix1,ic))
                ix1 += 1
                if ix1 == ix2:
                    # goto 4
                xx = xx1
                t1 = ddx*xx*chi[ix1][ic
                ix2 -= 1
            
                for ix in range(ix1, ix2+1):
                    xx *= dxx
                    t2 = ddx*xx*chi[ix][ic]
                    sum += t1 + t2
                    t1 = t2
            
                break
            continue
    
    sum = 0.5*(dx2-dx1)*((rdx1+rdx2)*xx*chi(ix1-1,ic)+(2.-rdx1-rdx2)*xx1*chi(ix1,ic)) # label 3
    
    acc(i) += 0.5*sum*float(na[ja]) / (ad[ja]*rx[i]) # label 4?
    
    return

#---------------------------------------------------------------------
#  subroutine phsh2cav
#  potential-to-phase-shift calculation(cavleed package)
#  uses loucks grid (e.g. as supplied by the muffin-tin potential
#  program).  energies input in hartrees.
      subroutine phsh_cav(mufftin_file, phasout_file,
     +                    dataph_file, zph_file)

       character(len=*), intent(inout) :: mufftin_file, zph_file
       character(len=*), intent(inout) :: phasout_file, dataph_file
       dimension v(250),rx(250),phs(20)
       double precision name(2),mtz,delstore(401,15),estore(401)
       integer unit_muff, unit_phas, unit_data, unit_zph

       # check for null input strings
       if (len_trim(mufftin_file)<1) mufftin_file = "mufftin.d"
       if (len_trim(phasout_file)<1) phasout_file = "phasout"
       if (len_trim(dataph_file)<1) dataph_file = "dataph"
       if (len_trim(zph_file)<1) zph_file = "zph.o"

       # get free file units
       unit_muff = 5
       unit_data = 6
       unit_phas = 7
       unit_zph = 8

       # first input channels
       open(unit=unit_muff,file=mufftin_file,status='old')

       # now output channels
       open(unit=unit_zph,file=zph_file,status='unknown')
       open(unit=unit_phas,file=phasout_file,status='unknown')
       open(unit=unit_data,file=dataph_file,status='unknown')

      # standard values for phase shifts calculation
       write(unit_data,110)
       emin=1.
       emax=12.
       estep=.25
       ianz=(emax-emin)/estep +1.01
       nl=12
       read(unit_muff,103) nr

       do kkk=1,nr
         read(unit_muff,100) name
         read(unit_muff,101) z,rmt, mtz
         read(unit_muff,103) ntab
         mtz=mtz/2.

         do ix=1,ntab
           read(unit_muff,219) rx(ix), v(ix)
         end do

         write(unit_zph,200)name,z,rmt,mtz
         write(unit_data,181)(name(i),i=1,2)
 181     format('non-relativistic phase shifts for ',2a4)
         write(unit_data,1030) emin,estep,ianz,nl
 1030    format(2f9.4,2(2x,i3))
         e=emin
         ncount=0
1        e=2.*e
         ncount=ncount+1
         call ps(v,rx,ntab,rmt,e,phs,nl,unit_zph)
         e=0.5*e
         write(unit_zph,201)e,(phs(l),l=1,nl)
         write(unit_data,1040)e*27.21,(phs(l),l=1,nl)
 1040    format(f9.4,8f8.4)

         # store phase shifts
         do kk=1,nl
           delstore(ncount,kk)=phs(kk)
         end do

         estore(ncount)=e
         e=e+estep
         if (e <= emax) goto 1

         # write phase shifts as function of energy for plotting
         do kk=1,nl
           write(unit_data,107) kk-1

           do ii=1,ncount
             write(unit_data,*) estore(ii),delstore(ii,kk)
           end do

           write(unit_data,*)
         end do

       end do

       # close file handles
       close(unit_muff)
       close(unit_data)
       close(unit_zph)
       close(unit_phas)

       return

100    format(2a8)
101    format(3f8.4)
103    format(i4)
107    format('"l=',i2)
110    format ("titletext: ","delta(e)")
200    format("phase shifts for ",2a8,2x,"atomic number,",f6.1/
     + "muffin-tin radius",f8.4,6x,
     + " muffin-tin zero level",f8.4," hartrees")
201    format("energy",f8.4," hartrees"/(10f12.5))
219    format(2e14.5)

       return

      return

#***********************************************************************
#  subroutine to calculate nl phase shifts (l=0,nl-1) for an
#  atomic potential tabulated on the loucks radial grid.
#     v?  atomic potential (rydbergs)
#     rx?  loucks' exponential grid  rx(i)=exp(-8.8+0.05(i-1))
#     ngrid?  number of tabulated points
#     rad?  limit of integration of schrodinger equation (a.u.)
#     e?  energy (rydbergs)
#     phs?  phase shifts for l=0 to nl-1
#  reference? loucks t l, (1967), a.p.w. method, benjamin, ny.
#***********************************************************************
      subroutine ps(v,rx,ngrid,rad,e,phs,nl,file_unit)

       integer, intent(in)    :: file_unit
       dimension v(ngrid),rx(ngrid),phs(nl)
       dimension wf(250),bj(25),bn(25),xr(10),fr(5)
       data pi,dx,dac,db/3.141592653589,0.05,2.083333e-04,2.083333e-03/

       index(x)=20.*(log(x)+8.8)+2.
       # tabulation of spherical bessel functions in bj and bn
       es=sqrt(e)
       x=es*rad
       z=x
       ll=nl+1
       call calcbf(bj,bn,ll,x,file_unit)

       # integration of the radial schrodinger equation by the numerov
       # method (see loucks p56 ff). wf contains wave function x radius,
       # on the loucks grid
       x1=exp(-8.8)

       do l1=1,nl
         fl=l1-1
         fl2=fl+0.5
         ff=fl2*fl2
         y1=x1**fl2
         y2=exp(dx*fl2)*y1

         write(file_unit,60)fl,y1,y2
60       format('0l',f5.1,5x,'y1,y2',2e14.5)

         gam1=ff+rx(1)*rx(1)*(v(1)-e)
         gam2=ff+rx(2)*rx(2)*(v(2)-e)
         wf(1)=y1*sqrt(rx(1))
         wf(2)=y2*sqrt(rx(2))

         do ix=3,ngrid
           gam=ff+rx(ix)*rx(ix)*(v(ix)-e)
           a=1.-dac*gam
           b=-2.-db*gam2
           c=1.-dac*gam1
           yn=-(b*y2+c*y1)/a
           wf(ix)=yn*sqrt(rx(ix))
           y1=y2
           y2=yn
           gam1=gam2
           gam2=gam
         end do

         # lagrangian interpolation for wavefunction and derivative at
         # radius x.  wfn holds wavefunction x radius, and dwfn
         # derivative x radius
         x=rad
         jr=index(rad)

         do j=1,5
           xr(j)=rx(jr-5+j)
           xr(j+5)=xr(j)
           fr(j)=wf(jr-5+j)
         end do

         wfn=0.
         dwfn=0.
         a=(x-xr(1))*(x-xr(2))*(x-xr(3))*(x-xr(4))*(x-xr(5))

         do i=1,5
           term=a/(x-xr(i))/(xr(i)-xr(i+1))/(xr(i)-xr(i+2))
     +      /(xr(i)-xr(i+3))/(xr(i)-xr(i+4))
           sum=0.

           do j=1,5
             if (i == j) break  # loop
             sum=sum+term/(x-xr(j))
           end do

           wfn=wfn+term*fr(i)

           dwfn=dwfn+sum*fr(i)
         end do

         # logarithmic derivative
         dloga=dwfn/wfn-1./rad

         # phase shifts
         x=es*rad
         a=fl*bj(l1)/x-bj(l1+1)
         b=fl*bn(l1)/x-bn(l1+1)
         a=es*a-dloga*bj(l1)
         b=es*b-dloga*bn(l1)
         phs(l1)=pi/2.
         if (abs(b) > 1.0e-8) phs(l1)=atan(a/b)
         write(file_unit,78)phs(l1)
78       format('phase shift',f10.4)

       end do

       return

      return

#***********************************************************************
def calcbf(bj,bn,nl,x,file_unit):
    if abs(x) < 1.0e-6:
         write(file_unit,200)x
200      format("** argument",e12.4," too small for routine calcbf")
         return
       endif

       bj(1)=sin(x)/x
       bn(1)=-cos(x)/x

       if (nl == 1) return

       bj(2)=(bj(1)-cos(x))/x
       bn(2)=(bn(1)-sin(x))/x

       if (nl == 2) return

       if (float(nl*(nl+1)) <= x*x):

         # forward recurrence for bj's
         fl=3.0

         do l=3,nl
           bj(l)=fl*bj(l-1)/x-bj(l-2)
           fl=fl+2.
         end do

         # forward recurrence for bn's
         fl = 3.

         do l=3,nl
           bn(l)=fl*bn(l-1)/x-bn(l-2)
           fl=fl+2.
         end do

         return

       endif

       # backward recurrence for bj's
       bj0 = bj(1)

       bj1 = bj(2)
       nn = max0(10,2*nl)
       a = 0.
       b = 1.
       fl = float(2*nn+1)

       do i=1,nn
         l=nn-i+1
         c=fl*b/x-a
         if (l<=nl)bj(l)=c
         a=b
         b=c
         fl=fl-2.
       end do

       # normalisation
       b = bj0 / bj[1]
       if abs(bj0) < 0.01: 
           b = bj1 / bj[2]

       do l=1,nl
         bj(l)=b*bj(l)
       end do

       # forward recurrence for bn's
       fl=3.

       do l=3,nl
         bn(l)=fl*bn(l-1)/x-bn(l-2)
         fl=fl+2.
       end do

       return

      return

#---------------------------------------------------------------------
#  subroutine phsh_wil
#---------------------------------------------------------------------
#  a.r. williams^ phase shift program (given a muffin-tin potential)
      subroutine phsh_wil(mufftin_file, phasout_file,
     +                    dataph_file, zph_file)

       character(len=*), intent(inout)     :: mufftin_file, zph_file
       character(len=*), intent(inout)     :: phasout_file, dataph_file
       real e(401),s(401,15),c(401,15),del(15),delold(15)
       real dell(9),delstore(8,401,15)
       integer tlp1
       common / cm16 / e1, e2, ne, ix,neuo
       common / cmrv / r(201), v(201, 15), nr, nl, z
       common / cm5 / y(30,4), f(30,4), ilst
       namelist / nl2 / ip,nrr

       # check for null input strings
       if (len_trim(mufftin_file)<1) mufftin_file = "mufftin.d"
       if (len_trim(phasout_file)<1) phasout_file = "phasout"
       if (len_trim(dataph_file)<1) dataph_file = "dataph"
       if (len_trim(zph_file)<1) zph_file = "zph.o"

       # first input channels
       open(unit=5,file=mufftin_file,status='old')

       # now output channels
       open(unit=6,file=zph_file,status='unknown')
       open(unit=7,file=phasout_file,status='unknown')
       open(unit=8,file=dataph_file,status='unknown')

       pi=3.1415926535

       # read ip
       # ip=0: only radial wavefunction
       # ip=1: phase shifts in addition
       # ip=2: s and c
       # ip=3: produce logarithm of phase shifts
       # nrr= number of inequivalent atoms for which we want phase shifts
       ip=1
       read(5, nl2)
       write(6, nl2)
       write(8,110)
110    format ("titletext: ","delta(e)")

       # input
       do 2  kkk=1,nrr
         call s16
         tx = 2. * r(nr)/ float(nr - 1)
         de = (e2 - e1) / float(max0(ne - 1, 1))

         do i = 1, ne
           e(i) = e1 + float(i - 1) * de
           # radial integration
           call s10(e(i))
           t3 = r(nr) * e(i)
           t4 = r(nr) * t3

           do lp1 = 1, nl
             l = lp1 - 1
             tlp1 = 2 * l + 1
             t5 = r(nr) ** lp1
             ut = f(tlp1, ilst)/tx + float(l) * y(tlp1,ilst) / r(nr)
             t1 = (f44(l,t4) * y(2*lp1,ilst) + t3 * f44(lp1,t4)
     +              * y(tlp1,ilst)) * t5
             t2 = (f45(l,t4) * ut - t3 * f45(l - 1,t4)
     +              * y(tlp1,ilst)) * r(nr) / t5
             s(i, lp1) = t1
             c(i, lp1) = t2
           end do

         end do

         is = 2
         i4 = 9
         if (ip < 1) goto 15

         # produce phase shifts
         do lp=1,nl
           delold(lp)=0.0
         end do

         do i = 1, ne
           do lp = 1, nl
             del(lp) = atan(-abs(e(i))**(lp-.5)*s(i,lp) / c(i, lp))
           end do

           # remove discontinuities by multiples of pi
           do lp=1,nl
             ls=0
111          deldif=del(lp)-delold(lp)
             if (abs(deldif)<0.7) break  # loop
             ls=ls+1
             del(lp)=del(lp)-sign(pi,deldif)
             if (ls<5) goto 111
             write (6,115) lp
115          format(" too large change in phase shift [l=",1i4,
     +              "] since last energy ",/,
     +              " discontinuity by multiple of pi possible")
             delold(lp)=del(lp)
           end do


           if (neuo==2) e(i)=0.5*e(i)

           # print phase shifts
           write(6, 12) e(i), (del(lp), lp = 1, nl)

           # write phase shifts in format used by leed program
           # write(7,71) e(i),(del (lp),lp=1,nl)

           # store phase shifts
           do kk=1,nl
             delstore(kkk,i,kk)=del(kk)
           end do

           if (ip < 3) break  # loop

           do j = 1, 9
             dell(j) = -4.
             if (del(j) < 1.0e-4) break  # loop
             dell(j) = log10(del(j))
           end do

         end do

         # write phase shifts as function of energy for plotting
         do kk=1,nl
           write(8,100) kk-1
100        format('"l=',i2)

           do ii=1,ne
             write(8,*) e(ii),delstore(kkk,ii,kk)
           end do

           write(8,*)
         end do

15     continue

2      continue

       write(7,*) 'be careful about the order of the elements'

       do ii=1,ne
         write(7,71) e(ii)

         do i=1,nrr
           write(7,72)(delstore(i,ii,lp),lp=1,nl)
         end do

       end do

       if (ip >= 2):
         do lp1 = 1, nl
           call s41(e, s(1, lp1), ne)
           call s41(e, c(1, lp1), ne)
         end do
       endif

       # close file handles
       close(5)
       close(6)
       close(7)
       close(8)

       return

12     format(1p8e14.7, /, 14x, 1p7e14.7, /)
71     format(1f7.4)
72     format(10f7.4)

      return

#***********************************************************************
#  subroutine s16
#  s16 inputs data
#   cs: core shift (position of zero of energy)
#   z: atomic number
#   e1,e2: first and last energies desired (in rydbergs or
#           hartrees, cf. neui)
#   ne: number of energies desired
#   nl: number of phase shifts desired (=lmax+1)
#   nr: number of radial grid points used in calculation
#   ix=0: signal to stop
#   ix=1: signal to expect a (new) potential to be read in
#   rt: muffin-tin radius (in bohr radii)
#   neui,neuo
#    if =1: rydberg unit used for input (neui) and output (neuo) of
#           energies and potential
#    if =2: hartree unit (double rydberg) used instead of rydberg
#           unit for input (neui) and output (neuo)
#   potyp=1: radial input as v(r)
#   potyp=2: radial input as r*v(r)
#
#***********************************************************************
def s16():
    # common / cm16 / e1, e2, ne, ix,neuo
    # common / cmrv / r, v, nr, nl, z
    # dimension r(201), v(201, 15)
    # real rs(200),          zs(200), ztt(201)
    # dimension fmt(18)
    # namelist / nl16 / cs,z,e1,e2,ne,nl,nr,ix,rt,neui,neuo,potyp

    # set default values of variables in namelist /nl16/
    ix = 1
    e1 = 4.
    e2 = 24.0
    ne = 30
    nl = 9
    nr = 101
    neui = 1
    neuo = 2
    potyp = 2
    read(5, nl16)
    cs = 0.0

    if ix < 1: 
        return

    if (neui != 1):
        cs *= 2.0
        e1 *= 2.0
        e2 *= 2.0

    write (6, nl16)

    drdn2 = (float(nr - 1))*(float(nr - 1)) / rt
    # read format used for input of r vs. v(r) or r vs. r*v(r)
    # (v is assumed positive)
    for i in range(1,200+1):
        rs[i], zs[i] = [t(s) for s in zip[ (float, float), unit5.read(28).split()]
         # the next lines assume that input potential & cs are negative
         zs[i] *= -1

         if rs(i) < 0: 
            break


    nrs = i - 1

    if neui != 1:
        for i in range(1, nrs+1):
           zs[i] *= 2.0

    while True:
        if potyp ! 2:
            for i in range(1, nrs+1):
                zs[i] = (zs[i] - cs) * rs[i]
            break

        for i in range(2, nrs+1):
            zs[i] = (zs[i] / rs[i] - cs) * rs[i]
        break
        
    iv = 1
    r[1] = 0.
    ztt[1] = z + z

    for i in range(2, nr+1):
        r(i) = (float(i - 1)) ** 2 / drdn2

        while True:
            if (r[i] <= rs[iv + 2]) or (iv + 3 >= nrs):
                break
                
            iv += 1
            continue

        ztt(i) = f12(rs(iv), zs(iv), r(i), 4)

        for lp1 in range(1, nl+1):
            v[i][lp1] = -ztt[i] / r[i]

def f12(x, y, z, n):
    ''' 
    Performs iterative interpolation in a table of n values of
    x and y to find the value of y at z
    '''
    w = [float] * 20

    w[1] = y[1]

    for i in range(2, n+1):
        w[i] = y[i]
        u = z - x[i]
        ip1 = i + 1

        for j in range(2, i+1):
            k = ip1 - j
            w[k] = w[k + 1] + u * (w[k] - w[k + 1]) / (x[k] - x[i])

      return w[1]

def s5(e):
    '''
    Integrates using Hamming's method for first order differential equations
    '''
    eest = [float] * 30
    vme = [float] * 15
    # common / cmrv / r(201), v(201, 15), nr, nl, z
    # common / cm5 / y(30, 4), f(30, 4), ip1

    nj = 2 * nl

    for j in range(1, nj+1):
        eest(j) = 0.

    for i in range(5, nr+1):
        for lp1 in range(1, nl+1):
            vme(lp1) = (v(i, lp1) - e) * r(i)

        t1 = 2. / float(i - 1)
        ip1 = mod(i - 1, 4) + 1
        im2 = mod(ip1, 4) + 1
        im1 = mod(im2, 4) + 1
        ip0 = mod(im1, 4) + 1

        for j in range(1, nj+1):
            f(j, im2) = y(j, ip1) + (2. * (f(j, ip0) + f(j, im2)) - f(j, im1)) / 0.75
            y(j, ip1) = f(j, im2) - 0.925619835 * eest(j)

        for j in range(1, nj+1, 2):
            jp1 = j + 1
            lp1 = jp1 / 2
            flp1 = lp1
            f(j, ip1) = (flp1 * y(j, ip1) + r(i) * y(jp1, ip1)) * t1
            f(jp1, ip1) = (vme(lp1) * y(j, ip1) - flp1 * y(jp1, ip1))*t1

        for j in range(1, nj+1):
            y(j, ip1) = y(j, ip0)+(y(j, ip0)-y(j, im2) + 3.*(f(j, ip1) + 2. * f(j, ip0) - f(j, im1))) / 8.
            eest(j) = f(j, im2) - y(j, ip1)
            y(j, ip1) = y(j, ip1) + .743801653e-1 * eest(j)

        for j in range(1, nj+1, 2):
            jp1 = j + 1
            lp1 = jp1 / 2
            flp1 = lp1
            f(j, ip1) = (flp1 * y(j, ip1) + r(i) * y(jp1, ip1)) * t1
            f(jp1, ip1) = (vme(lp1) * y(j, ip1) - flp1 * y(jp1, ip1))*t1

def s10(e):
    ''' 
    Calculates power series expansion of the solution about the origin
    and radial integration in s5
    '''
    tlp1 = int
    a = b = [float] * 10 
    tr = [float] * 4
    # common / cmrv / r(201), v(201, 15), nr, nl, z
    # common / cm5 / y(30, 4) , f(30, 4) , ilst

    ni = 2 * nl
    tz = 2. * z
    a[1] = 1.

    for i in range(1, ni+1, 2):
        lp1 = (i + 1) / 2
        tlp1 = 2 * lp1
        ep = e - v(4, lp1) - tz/ r(4)
        y(i, 1) = 0.
        y(i + 1, 1) = 0.
        a(1) = a(1) / float(2 * lp1 - 1)
        b(1) = - z * a(1) / float(lp1)

        for j in range(2, 4+1):
            tr(j) = r(j) ** lp1
            y(i, j) = a(1) * tr(j)
            y(i + 1, j) = b(1) * tr(j)

        for k in range(1, 9+1):
            a(k + 1) = b(k) / float(k)
            b(k + 1) = -(ep * a(k) + tz * a(k + 1)) / float(tlp1 + k)

            for j in range(2, 4+1):
                tr(j) = tr(j) * r(j)
                y(i, j) = y(i, j) + tr(j) * a(k + 1)
                y(i + 1, j) = y(i + 1, j) + tr(j) * b(k + 1)

            if abs(tr[4] * a[k + 1] / y[i][4]) < 1.0e-4:
                break

        write (6, 4) e, lp1, r(4), (a(k), k = 1, 10)
    4        format(1pe10.2, i10, 11e10.2)

    for j in range(2, 4+1):
        t1 = 2. / float(j - 1)

        for i in range(1, ni+1, 2):
            ip1 = i + 1
            lp1 = ip1 / 2
            flp1 = lp1
            f[i][j] = (flp1 * y[i][j] + r[j] * y[ip1][j]) * t1
            f[ip1][j] = ((v[j][lp1] - e) * r[j] * y[i][j] - flp1 * y[ip1][j])* t1

    s5(e)

def f44(l, x):
    ''' f44  evaluates the special version of the spherical Bessel function '''
    s = [float] * 20

    js = l + l + 1

    if abs(x / float(js)) <= 10.:
        fi = 1.

        if l >= 1:
            for k in range(3, js+1, 2):
               fi *= float(k)

        t1 = 1. / fi
        dt = 1.
        t = 1.
        i = 0

        for k in range(1, 101):
            i += 2
            dt *= -x / float(i * (i + js))
            t += dt

            if abs(dt) < 1.e-8: 
                break

        return t1 * t

    t = sqrt(abs(x))

    if x >= 0.:
        s[2] = sin(t) / t

        if l >= 0:
            return s[2]

        s[1] = cos(t)

        is = l + 2

        for i in range(3, is+1):
            s[i] = (s[i - 1] * float(2 * i - 5) - s[i - 2]) / x

       return s[is]

    s[2] = sinh(t) / t

    if l < 1:
        return s[2]

    s[1] = cosh(t)

    is = l + 2 

    for i in range(3, is+1):
        s[i] = (s[i - 1] * float(2 * i - 5) - s[i - 2]) / x

   return s[is]

   
def f45(l, x):
    ''' Evaluates special version of the spherical Neumann function '''
    s = [float] * 20

    if l <= 0:
        return -f44(l+1, x)

    lp1 = l + 1
    js = l + l + 1

    if abs(x / float(js)) <= 10.:
        fi = 1.

        if (l >= 1):
            for k in range(3, js+1, 2):
               fi *= float(k)
            
        t1 = fi / float(js)
        dt = 1.
        t = 1.
        i = 0

        for k in range(100):
            i += 2
            dt *= -x / float(i * (i - js))
            t += dt

            if abs(dt) < 1.e-8: 
                break


        t1 *= t
        f45 = t1

        return f45

    t = sqrt(abs(x)) 

    if x >= 0.:
        s(2) = cos(t)

        if l <= 0:
           return s[2]

        s[1] = -sin(t) / t

        is = l + 2

        for i in range(3, is+1):
            s[i] = s[i - 1] * float(2 * i - 5) - x * s[i - 2]

        return s[is]
        
    s(2) = cosh(t)

    if l < 1:
        return s[2]

    s[1] = -sinh(t) / t

    is = l + 2

    for i in range(3, is+1):
        s[i] = s[i - 1] * float(2 * i - 5) - x * s[i - 2]

    return s[is]


#***********************************************************************
#  s41 plots y against x
#***********************************************************************
def s41(x, y, n)

       character b, c, o, d, p
       dimension x(100), y(100), p(97)
       data b, c, o, d / " ", "*", "0", "i" /
       y1 = 0
       y2 = 0

       do i = 1, n
         y1 = min(y1, y[i])
         y2 = max(y2, y[i])
       end do

       do i = 1, 97
         p(i) = b
       end do

       t = 96/ (y2 - y1)
       j0 = -y1 * t + 1.5
       p(j0) = "0"

       if (n >= 30) write(6, 3) p
3      format("1", 34x, 97a1, //)
       if (n < 30) write(6, 6) p
6      format(////, 35x, 97a1, //)

       p(j0) = d

       do i = 1, n
         j = t * (y(i) - y1) + 1.5
         p(j) = c
         write(6, 4) x(i), y(i), p
4        format(1x, 1p2e16.6, 2x, 97a1)
         p(j) = b
         p(j0) = d
       end do

       return

      return


def phsh_rel(mufftin_file, phasout_file, dataph_file, inpdat_file):

       implicit double precision (a-h,o-z)
       character(len=*), intent(inout)     :: mufftin_file, inpdat_file
       character(len=*), intent(inout)     :: phasout_file, dataph_file
       character opt*3,opt1*3,opts*3,aname*2,an*30,bdata*28
       character sub*3,record*3,tl*1,sl*1,ss1*6,ss2*6,wrd*6
       character ams(5)*4
       real jf,jfs,jfd,jfu
       real name(4)
       dimension jf(250,18),energ (250)
       common/zzzz/zp(340),vs,ipt,jri
       common /z/ rmaxi
       dimension adata(7)
       data zero/0.0d0/,one/1.0d0/,two/2.0d0/,anine/9.0d0/,half/.5d0/
       data zilch/1.0d-4/,tol/.005d+0/,des/.025d+0/
       data ams/'nc= ','l= ',' es=',' de=','id= '/
       data tl/'l'/,sl/'s'/,ss1/'nospin'/,ss2/' spin '/
       data sub/'sub'/,record/'nos'/
  1    format (3d12.4,4x,i3,4x,d12.4)
  2    format (5e14.6)
  8    format (f9.4,8f8.5)
 12    format (5d14.6)
 13    format ("1",//,t61,'input data',//)
 18    format(f10.4,f9.4,2i5)

       # check for null input strings
       if (len_trim(mufftin_file)<1) mufftin_file = "mufftin.d"
       if (len_trim(phasout_file)<1) phasout_file = "phasout"
       if (len_trim(dataph_file)<1) dataph_file = "dataph"
       if (len_trim(inpdat_file)<1) inpdat_file = "inpdat"

       open(unit=4,file=inpdat_file,status='unknown')
       open(unit=5,file=mufftin_file,status='old')
       open(unit=7,file=phasout_file,status='unknown')
       open(unit=8,file=dataph_file,status='unknown')

       pi=4.d0*datan(1.d0)
       pi2=0.5d0*pi
       write (4,13)

       do while (.true.)
         read(5,"(4a4)",end=999,err=999) (name(i),i=1,4)
         read(5,1)es,de,ue,lsm,vc

         # nl is the number of plotted phase shifts
         nl=8
         write (4,11) es,de,ue,opt,opt1,lsm
 11      format (3d12.4,4x,2a3,i3)
         read (5,16) nz,adata(1),jri,alc,blc,clc,exca,excb,exco
 16      format (i4,f10.6,i4,t21,6f10.6)
         write (4,76) nz,adata(1),jri,alc,blc,clc,exca,excb,exco
 76      format (i4,f10.6,i4,2x,6f10.6)

         vs=0.
         if (opts == sub) vs=vc
         if ((jri <= 0) .or. (jri > 340)) goto 999
         read(5,2) (zp(j),j=1,jri)
         write (4,12) (zp(j),j=1,jri)
         rhoz=-0.90306514d+01
         delrho=0.3125d-01
         rm=dexp(rhoz)
         xrx=dexp(delrho)

         do j=1,jri
           if (zp(j) < 0.0) zp(j)=-zp(j)
           if (j == jri) break  # loop
           rm=xrx*rm
         end do

         if (de <= zero):
           es=-half
           de=des
           ue=one
         endif

         n=(ue-es)/de+half

         n=n+1

         write(07,181)(name(i),i=1,4)
 181     format('relativistic phase shifts for ',4a4)
         write(07,18) es,de,n,lsm

         es=es/13.6
         de=de/13.6
         ue=ue/13.6
         if (n > 250)  n=250  # this is suspicious...
         l=0
         e=es
         ipt=2

         if (opt == record):
           ipt=-2
           wrd=ss1
         else:
           wrd=ss2
         endif

         kap=-1
         l=1

         do j=1,n
           dxaz=0.0d0
           ttr=dlgkap(e,kap)/(12.5663706)
           call sbfit(ttr,e,l-1,rmaxi,jfs)
           jf(j,l)=jfs
           e=e*13.6
           energ(j)=e
           e=e/13.6
           e=e+de
         end do

         do while(l <= lsm)

           kap=-(l+1)
           lind=l+1
           e=es
           lvor=0

           do j=1,n
             dlu=dlgkap (e,kap)/(12.5663706)
             dld=dlgkap(e,l)/(12.5663706)
             call sbfit(dld,e,l,rmaxi,jfd)
             call sbfit(dlu,e,l,rmaxi,jfu)
             lk=0
             zfdiff=-(jfd-jfu)*lvor
             if (zfdiff > 2.5) lk=l
             if (zfdiff < -2.5) lk=l+1
             jfs=l*jfd-kap*jfu+lvor*lk*pi
             jfs=jfs/(2*l+1)
             if (jfs >pi2) jfs=jfs-pi2*2.
             if (jfs < -pi2) jfs=jfs+pi2*2.
             jf(j,lind)=jfs
             if (lk == 0) lvor=sign(1.,jfs)
             e=e+de
           end do

           l=l+1

         end do

         do i=1,n
           lsm1=lsm+1
           write(7,8) energ(i),(jf(i,l),l=1,lsm1)

         end do

         do kk=1,nl
           write(8,100) kk-1

           do ii=1,n
             write(8,*) energ(ii),jf(ii,kk)
           end do

           write(8,*)

         end do

100      format('"l=',i2)
         es=es*13.6
         de=de*13.6
         ue=ue*13.6

       end do

 999   continue

       write (4,900)
 900   format (//,t57,'end of input data')

       # close file handles
       close(4)
       close(5)
       close(7)
       close(8)


def dlgkap (e,kappa)
    ''' 
    Calculates the logrithmic derivative of the large
    component using the procedure described by loucks in appendix 7.
    the small multiplicative factor is included.
    potential  in the form of 2zp is to be passed in the common /zzzz/
    the radial functions are made available in the common /radfun/
    waber mesh (xnot=-9.03065 and dx=1/32) is used
    jri is the mesh point of the apw sphere radius
    e is the energy to be used (in rydbergs)
    4 pi r**2 inserted now. for compounds only.
    '''
    # dimension pot(340),u(340),w(340),up(340),wp(340),sxk(4),sxm(4)
    # common /z/ t
    # common /radfun/ u,w
    # common/zzzz/pot,vcz,ipt,jri
    # data ustart/1.d-25/,zilch/1.d-30/
    # data test/1.d+6/,xs/-9.03065133d+00/
    # data dx/3.125d-2/,c/2.740746d+2/,cin/1.3312581146d-5/,hf/.5d+0/
    # data th/.3333333333d+0/,t2/2.d+0/,t7/7.d+0/,t11/11.d+0/
    # data t12/12.d+0/,t14/14.d+0/,t26/26.d+0/,t32/32.d+0/,zero/.1d+0/

    #set up for relativistic or no relativistic effect
    if ipt < 0:
        cin = 0.0

    # set up starting values
    dx2 = hf * dx
    x20 = (.3) * dx
    xmft = (4.4444444444444e-2) * dx
    ts = exp(xs)
    tdx = exp(dx)
    hoc = (vcz*ts + pot[1])/c
    xk = kappa
    
    u[1] = ustart
    
    if (abs(hoc/xk) > 0.05):
        p =(xk + sqrt(xk * xk - hoc * hoc))/hoc
    else:
        p =(xk + abs(xk)) / hoc - (hf * hoc / abs(xk))

    tc = ts * (e + vcz) + pot(1)

    vc = cin * tc
    w[1] = c * p * ustart

    # start runge-kutte procedure
    x = xs

    n = 1

    ik = 0 # label 25

    np1 = n + 1
    xc = x
    bgc = pot(n)
    wc = w(n)
    uc = u(n)

    while True:
        ik += 1 # label 20

        t = exp(xc)
        tc = ((e + vcz) * t) + bgc
        vc = cin * tc

        sxk(ik)=dx2*(wc*(vc+t)-xk*uc)

        sxm(ik)=dx2*(xk*wc-tc*uc)

        if ik == 1:
            xc += dx2
            uc += sxk(1)
            wc += sxm(1)
            bgc = hf * (bgc+pot(np1))
            continue # goto 20

        elif ik == 2:
            uc += sxk(2) - sxk(1)
            wc += sxm(2) - sxm(1)
            continue # goto 20

        elif ik == 3:
               xc += dx2
               uc += t2*sxk(3) - sxk(2)
               wc += t2*sxm(3) - sxm(2)
               bgc = pot(np1)
               continue # goto 20

        else:
            w(np1) = w(n)+(sxm(1)+sxm(4)+t2*(sxm(2)+sxm(3)))*th
            u(np1) = u(n)+(sxk(1)+sxk(4)+t2*(sxk(2)+sxk(3)))*th
            up(np1) = (vc + t) * w(np1) - xk * u(np1)
            wp(np1) = xk * w(np1) - tc * u(np1)
            x += dx
            n = np1
            if n < 6:
                continue # goto 25
                
        break

        # end of starting integration.  begin milne procedure.
        t = exp(x)
        t *= tdx # label 26
        
        np1 = n+1
        nm1 = n-1
        nm2 = n-2
        nm3 = n-3
        nm4 = n-4
        nm5 = n-5
       
        tc = ((e + vcz) * t) + pot(np1)
        vc = cin * tc
        unp = u(nm5) + x20*(t11*(up(n) + up(nm4)) + t26*up(nm2) - t14*(up(nm1) + up(nm3)))
        wnp = w(nm5) + x20*(t11*(wp(n) + wp(nm4)) + t26*wp(nm2) - t14*(wp(nm1) + wp(nm3)))
        nit = 0
       
        while True: # label 33
            up(np1) = (vc+t)*wnp - xk*unp # label 33
            wp(np1) = xk*wnp -tc*unp
            unp2 = u(nm3) + (t7 * (up(np1) + up(nm3)) + t32 * (up(nm2) + up(n)) + (t12 * up(nm1)) * xmft)
            wnp2 = w(nm3) + (t7 * (wp(np1) + wp(nm3)) + t32 * (wp(nm2) + wp(n)) + (t12 * wp(nm1)) * xmft)

            # compare predictor with corrector
            if (abs(test * (unp2 - unp)) > abs(unp2)): 
                # goto 31
            if (abs(test * (wnp2 - wnp)) <= abs(wnp2)): 
                # goto 32
            if nit < 5: # label 31
               # goto 81
            # goto 32

            nit += 1 # label 81

            wnp = wnp2
            unp = unp2
            continue # goto 33

            w(np1) = wnp2 # label 32

        u(np1) = unp2
        n = np1

        if n < jri:
            # goto 26

        # end of milne procedure
        if abs(u(jri)) <= zilch:
            u(jri) = sign(zilch, u(jri))

        p = (t + vc) / t

        wnp = p * w(jri) / u(jri)
        unp = wnp - (kappa + 1) / t
        dlgkap = (12.5663706) * t * t * unp


def sbfit(t,e,l,r,jfs):
    jfs = float
    kappa = float

    se = sngl[e]
    sr = sngl[r]
    st = sngl[t]
    
    kappa = sqrt(se)
    x = kappa*sr
    bj1 = sin(x)/x
    bn1 = -cos(x)/x
    bj2 = (bj1 / x) + bn1
    bn2 = (bn1 / x) - bj1

    if (l > 0):
        dl = st / (sr * sr)
    else:
        ls = 1

    ls = 1
    ls += 1
    bjt = (2*ls-1) * bj2 / x - bj1
    bnt = (2*ls-1) * bn2 / x - bn1
    bj1 = bj2
    bj2 = bjt
    bn1 = bn2
    bn2 = bnt

    if (l+1-ls > 0):
        ls += 1
    else:
        dl = st / (sr * sr)

    dl = st / (sr * sr)
    dl -= l/sr
    an = dl * bj1 + kappa * bj2
    ad = dl * bn1 + kappa * bn2
    jfs = 3.141592654 / 2.0

    if (abs(ad)-1.0e-8 > 0): 
        jfs = atan(an/ad)

