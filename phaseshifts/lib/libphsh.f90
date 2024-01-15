!=======================================================================
! libphsh.f90
! 
!!! WARNING: This is an experimental port of libphsh.f to FORTRAN 90 !!! 
!
! Modified code from the Barbieri/Van Hove phase shift (phshift)
! package available at:
!
! http://www.icts.hkbu.edu.hk/vanhove/VanHove_files/leed/leedpack.html
!
!=======================================================================

!-----------------------------------------------------------------------
! hartfock subroutine:
!
!  there are nr grid points, and distances are in bohr radii...
!
!  r(i)=rmin*(rmax/rmin)**(dble(i)/dble(nr)) , i=1,2,3,...nr-1,nr
!
!
!
!  the orbitals are store in phe(), first index goes 1...nr, the
!  second index is the orbital index (i...nel)
!
!  look at the atomic files after printing this out to see everything...
!
!  suffice it to say, that the charge density at radius r(i)
!  in units of electrons per cubic bohr radius is given by
!
!  sum of j=1...nel,
!  occ(j)*phe(i,j)*phe(i,j)/(4.d0*3.14159265....*r(i)*r(i))...
!
!  think of the phe functions as plotting the radial wave-functions
!  as a function of radius...on our logarithmic mesh...
!
!  final note:
!
!  the Dirac equation is solved for the orbitals, whereas their density
!  is treated by setting phe to the square root of Dirac's F*F+G*G
!  times the sign of G...
!
!  so we are doing Dirac-Fock, except that we are not treating exchange
!  exactly, in terms of working with major and minor components of the
!  orbitals, and the phe's give the CORRECT CHARGE DENSITY...
!
!  the above approximation ought to be very small for valence states,
!  so you need not worry about it...
!
!  the Breit interaction has been neglected altogether...it should not
!  have a huge effect on the charge density you are concerned with...
!
! author: Eric Shirley
!
!-----------------------------------------------------------------------
!***********************************************************************     
module WK
      
       real, dimension(3)    :: ra, ga, x, rb, rj
       real, dimension(5)    :: vmm, fr
       real, dimension(20)   :: ex, fac, fnt, nt
       real, dimension(250)  :: wk1, wk2
       real, dimension(3,3)  :: g  
       real, dimension(2000) :: rr, rs
       
       save
       
      end module WK
      
!***********************************************************************      
      module WF
       
       real, dimension(14)      :: wc, lc
       real, dimension(250, 14) :: wfc, wf2
       
       save
       
      end module WF
      
!***********************************************************************
      module CM16
       
       integer NE, ix, neuo
       real e1, e2
       
       save
       
     end module CM16

!***********************************************************************     
     module CMRV
       
       real, dimension(201)     :: r
       real, dimension(201, 15) :: v
       real Z
       integer nr, nl
        
       save
       
     end module CMRV
     
!***********************************************************************
     module CM5
     
      integer ilst, ip1
      real, dimension(30,4) :: y, f 

      save
      
     end module CM5
     
!***********************************************************************

    module ZZZZ
    
       real, dimension(340) :: zp
       real vs, pot, vcz
       integer ipt, jri
    
     save
     
    end module ZZZZ
    
!***********************************************************************
    module Z
    
     real rmaxi, T
    
     save
     
    end module Z

!***********************************************************************    
    module RADFUN
    
     real U, W
    
     save 
    
    end module RADFUN
  
      subroutine hartfock(input_file)
       
       implicit double precision (a-h,o-z)
       character(len=255), intent(in) :: input_file
       integer iorbs, iside, io2, ijive
       integer lmax, ihmax, nrmax, ntmax, npmax
       parameter (iorbs=33,iside=600)
       parameter (io2=iorbs*(iorbs+1)/2)
       parameter (ijive=io2*(io2+1)/2)
       parameter (lmax=4,ihmax=20,nrmax=4000,ntmax=10,npmax=60)
       double precision, dimension(nrmax) :: r, dr, r2, rho, v, q0, xm1, xm2
       integer, dimension(iorbs)          :: no, n
       double precision, dimension(iorbs)        :: ev, occ, ek, xnj
       double precision, dimension(nrmax, iorbs) :: phe, orb
       integer, dimension(4)                     :: njrc(4)
       integer, dimension(iorbs)                 :: is
       double precision, dimension(nrmax,7)      :: vi
       double precision, dimension(33)           :: rhs, co
       double precision, dimension(33,33)        :: w, wi
       double precision, dimension(0:12)         :: xint, rint
       double precision, dimension(11)           :: vav
       double precision, dimension(0:11)         :: pin, sig
       double precision, dimension(nrmax,0:3)    :: vctab
       character ichar
       character(len=11) jive
       character(len=60) jive2
       character(len=70) jive3
    
       open(unit=5, file=trim(input_file), status='old')
       rel = 0.d0

       do while(.true.)
       read (5,"(1a1)") ichar

       if (ichar == 'd') then
         ! RELATIVITY FACTOR.  (0=NR, 1=REL.)
         read (5,*) rel
       endif

       if (ichar == 'x') then
         ! ENTER EXCHANGE CORRELATION (ALPHA).
         ! 0=HARTREE-FOCK, POSITIVE=LDA, NEGATIVE=XALPHA.
         read  (5,*) alfa
       endif

       if (ichar == 'a') then
         call abinitio(etot, nst, rel, alfa, nr, r,                       &
     &                 dr, r2, dl, phe, njrc, vi, zorig, xntot, nel,      &
     &                 no, nl, xnj, ev, occ, is, ek, orb, iuflag)
       endif

       if (ichar == 'i') call initiali(zorig,nr,rmin,rmax,                &
     &    r,dr,r2,dl,njrc,xntot,nel)

       if (ichar == 'q') return

       if (ichar == 'w') then
         ixflag = 1
         iu = -1
         ir = 0
         call hfdisk(iu,ir,etot,nst,rel,nr,rmin,rmax,r,rho,               &
     &      zorig,xntot,ixflag,nel,                                       &
     &      no,nl,xnj,is,ev,ek,occ,njrc,vi,phe,orb)
       endif

       if (ichar =='r') then
         iu = -1
         ir = 1
         call hfdisk(iu,ir,etot,nst,rel,nr,rmin,rmax,r,rho,               &
     &      zorig,xntot,ixflag,nel,                                       &
     &      no,nl,xnj,is,ev,ek,occ,njrc,vi,phe,orb)
         call setgrid(nr,rmin,rmax,r,dr,r2,dl)
       endif

       if (ichar == 'u') then
         ! ENTER IUFLAG. (0=U, 1=SU, 2=R)
         read (5,*) iuflag
       endif

       if (ichar == 'c') then
         ! ENTER ALPHA,RS,RP,RD
         read (5,*) corpol, rs, rp, rd

         do k=1,nr
           fs = (1.d0 - exp(-(r(k) / rs)**2.d0))**2.d0
           fp = (1.d0 - exp(-(r(k) / rp)**2.d0))**2.d0
           fd = (1.d0 - exp(-(r(k) / rd)**2.d0))**2.d0
           vctab(k,0) = -corpol / 2.d0 * fs * fs / r(k)**4.d0
           vctab(k,1) = -corpol / 2.d0 * fp * fp / r(k)**4.d0
           vctab(k,2) = -corpol / 2.d0 * fd * fd / r(k)**4.d0
         end do

       endif

       if (ichar == 'f') then
         ! IUNIT, CORPOL
         read  (5,*) iunit, corpol
         ! ILEV,INUM,EOLD
         read  (5,*) ilev, inum, eold
         xl=nl(ilev)

         if (inum == 1) then
           read (5,*) eav
         else
           read (5,*) e1, e2
           eav = ((e1 * xl) + e2 * (xl + 1.d0)) / (xl + xl + 1.d0 )
         endif

         if (eav < 0.d0) eav = eold + eav
         if (iunit == 2) eav = eav / 2.d0
         if (iunit == 3) eav = eav / 27.2116d0
         if (iunit == 4) eav = eav * 0.000123985d0 / 27.2116d0

         sd = abs(abs(eav) - abs(ev(ilev)))
         rl = 0.d0
         rh = 10.d0
         sl = 0.d0
         sh = 0.d0

         if (sl * sh <= 0.00000001d0) rc = rl + (rh - rl) / 2.d0  ! label 300
         if (sl * sh > 0.00000001d0) rc = rl + (rh - rl) * (sd - sl) / (sh - sl)

         sc = 0.d0

         do i=1,nr
           f = (1.d0 - exp(-(r(i) / rc)**2.d0))**2.d0
           vcpp = corpol / (2.d0 * r(i)**4.d0) * f * f
           sc = sc + dr(i) * phe(i,ilev) * phe(i,ilev) * vcpp
         end do

         if (sc > sd) rl = rc
         if (sc > sd) sl = sc
         if (sc < sd) rh = rc
         if (sc < sd) sh = sc

         write (6,*) rc, sc

         do while (abs(sc - sd) > 0.000001d0)  ! replace goto 300
           if (sl * sh <= 0.00000001d0) rc = rl + (rh - rl) / 2.d0
           if (sl * sh > 0.00000001d0) rc = rl + (rh - rl) * (sd - sl) / (sh - sl)

           sc = 0.d0

           do i=1,nr
             f = (1.d0-exp(-(r(i)/rc)**2.d0))**2.d0
             vcpp = corpol/(2.d0*r(i)**4.d0)*f*f
             sc = sc+dr(i)*phe(i,ilev)*phe(i,ilev)*vcpp
           end do

           if (sc > sd) rl=rc
           if (sc > sd) sl=sc
           if (sc < sd) rh=rc
           if (sc < sd) sh=sc

           write (6,*) rc,sc

         end do

       endif

       if (ichar == 'p') then
         call pseudo(etot,nst,rel,alfa,nr,rmin,rmax,r,dr,r2,dl,           &
     &                phe,orb,njrc,vi,zorig,xntot,nel,                    &
     &                no,nl,xnj,ev,occ,is,ek,iuflag,vctab)
       endif

       if (ichar == 'g') then
         read(5,*) iu
         read(5,"(1x,1a11)") jive
         read(5,"(1x,1a60)") jive2
         read(5,"(1x,1a70)") jive3

         zizv=abs(r(nr-1)*vi(nr-1,1))

         write (iu,"(1x,1a11)") jive
         write (iu,"(1x,1a60)") jive2
         write (iu,"(1x,1a70)") jive3
         write (iu,*) 3,nr,zizv
         write (iu,*) (r(i),i=1,nr)
         write (iu,*) 0,(vi(k,1),k=1,nr)
         write (iu,*) 1,(vi(k,3),k=1,nr)
         write (iu,*) 2,(vi(k,5),k=1,nr)
         write (iu,*) (0.d0,k=1,nr)

         do j=1,nr
           rh=0.d0

           do k=1,nel
             rh=rh+phe(j,k)*phe(j,k)*occ(k)
           end do

           write (iu,*) rh

         end do

       endif

       if (ichar == 'v') then
         rold=0.d0

         ! open file streams
         open(unit=10)
         open(unit=11)
         open(unit=12)

         ! write to file streams
         do k=1,nr
           write (10,*) r(k), vi(k,1) * r(k)
           write (11,*) r(k), vi(k,3) * r(k)
           write (12,*) r(k), vi(k,5) * r(k)
           rold=r(k)
         end do

         ! close file streams
         close(unit=10)
         close(unit=11)
         close(unit=12)

       endif

       if (ichar == 'V') call fourier(nr,r,dr,r2,vi)

       end do

       return ! end

      end subroutine

!-----------------------------------------------------------------------
      subroutine abinitio(etot,nst,rel,alfa,nr,r,dr,r2,dl,                &
     &    phe,njrc,vi,zorig,xntot,nel,no,nl,xnj,                          &
     &    ev,occ,is,ek,orb,iuflag)
     
       implicit double precision (a-h,o-z)
       double precision dl, zorig, xntot
       integer nst, nel, iuflag
       double precision etot, rel, alfa
       integer nr
       integer iorbs, iside, io2, ijive
       integer lmax, ihmax, nrmax, ntmax, npmax
       parameter (iorbs=33,iside=600)
       parameter (io2=iorbs*(iorbs+1)/2)
       parameter (ijive=io2*(io2+1)/2)
       parameter (lmax=4,ihmax=20,nrmax=4000,ntmax=10,npmax=60)
       integer, dimension(4) :: njrc
       integer, dimension(iorbs) :: no, nl, nm, is
       double precision, dimension(iorbs) :: xnj
       double precision, dimension(nrmax) :: r, dr, r2, v, ev, occ, ek
       double precision, dimension(nrmax,iorbs) :: phe, orb
       double precision, dimension(nrmax, 7) :: vi
       double precision, dimension(nrmax,0:15) :: rpower
       
       ! note: this will be good for going up to and including l=3...

       do i=0,7
         xi=i
         do k=1,nr
           rpower(k,i)=r(k)**xi
         end do
       end do

       ! read in nfc, nel.  refer to the documentation for their meanings.

       ! ENTER NFC, NEL, RATIO, ETOL, XNUM
       read (5,*) nfc, nel, ratio, etol, xnum

       ! for all of the electrons, read in the quantum numbers.
       ! get the total Hartree-active charge & initialize eigenvalues.

       xntot = 0.d0

       write (6,*) 'N L M J S OCC ENERGY'

       do i=nfc+1,nel
         read (5,*) no(i), nl(i), nm(i), xnj(i), is(i), occ(i)
         ev(i) = 0.d0
         xntot = xntot + occ(i)

         ! initialise phe and orb
         do j=1,nr
           phe(j,i) = 0.d0
           orb(j,i) = 0.d0

         end do

       end do

       ! initialize the parameters for self-consistency loop.
       ! ratio is the mixture of old and new field mixing.
 110   call atsolve(etot,nst,rel,alfa,eerror,nfc,nr,r,dr,r2,dl,phe,       &
     &              njrc,vi,zorig,xntot,nel,no,nl,nm,xnj,ev,occ,is,ek,    &
     &              ratio,orb,rpower,xnum,etot2,iuflag)                   

       eerror = eerror * (1.d0 - ratio) / ratio

       do while(eerror > etol)
         call atsolve(etot,nst,rel,alfa,eerror,nfc,nr,r,dr,r2,dl,phe,     &
     &                njrc,vi,zorig,xntot,nel,no,nl,nm,xnj,ev,occ,is,ek,  &
     &                ratio,orb,rpower,xnum,etot2,iuflag)

         eerror = eerror * (1.d0 - ratio) / ratio
       end do

       ! write out information about the atom.

       do i=1,nel
         nj=xnj(i)+xnj(i)
         write (6,122) no(i),nl(i),nm(i),nj,'/2',is(i),occ(i),ev(i)
 122     format(1x,2i4,i2,i4,a2,i4,f10.4,f18.6)
       end do

       write (6,132) 'TOTAL ENERGY =  ', etot, etot * 27.2116d0
 132   format (1x,a16,2f14.6)

       return

      end subroutine

!-----------------------------------------------------------------------
      subroutine atsolve(etot,nst,rel,alfa,eerror,nfc,nr,r,dr,r2,dl,      &
     &    phe,njrc,vi,zorig,xntot,nel,no,nl,nm,xnj,ev,occ,is,ek,          &
     &    ratio,orb,rpower,xnum,etot2,iuflag)
     
       implicit double precision (a-h,o-z)
       double precision etot, rel, alfa, eerror, dl, zorig, xntot
       double precision ratio, xnum, etot2
       integer nst, nfc, nr, nel, iuflag
       integer iorbs, iside, io2, ijive
       integer lmax, ihmax, nrmax, ntmax, npmax
       parameter (iorbs=33,iside=600)
       parameter (io2=iorbs*(iorbs+1)/2)
       parameter (ijive=io2*(io2+1)/2)
       parameter (lmax=4,ihmax=20,nrmax=4000,ntmax=10,npmax=60) 
       double precision, dimension(nrmax) :: r, dr, r2, v
       integer, dimension(iorbs)          :: no, nl, nm, is
       double precision, dimension(iorbs) :: ek, ev, occ, xnj
       double precision, dimension(nrmax, iorbs) :: phe, orb
       integer, dimension(4) :: njrc
       double precision, dimension(nrmax, 7) :: vi
       double precision, dimension(nrmax) :: q0, xm1, xm2
       double precision, dimension(nrmax,0:15) :: rpower 

       ! initialize eerror, the biggest change in an eigenvalue, and etot.

       eerror = 0.d0
       etot = 0.d0

       ! run through all the orbitals.  calculate those not in the core.

       do i=1,nel

         if (i > nfc) then

           idoflag = 1
           call setqmm(i,orb,nl(i),is(i),idoflag,v,zeff,zorig,rel,        &
     &                 nr,r,r2,dl,q0,xm1,xm2,njrc,vi)

           xkappa = -1.d0
           if (abs(xnj(i)) > dble(nl(i)) + 0.25d0) xkappa = -nl(i) - 1
           if (abs(xnj(i)) < dble(nl(i)) - 0.25d0) xkappa = nl(i)

           call elsolve(i,occ(i),no(i),nl(i),xkappa,xnj(i),zorig,zeff,    &
     &                  evi,phe(1,i),v,q0,xm1,xm2,nr,r,dr,r2,dl,rel)

           if (abs(ev(i) - evi) > eerror) eerror = abs(ev(i) - evi)
           ev(i) = evi

           ekk = 0.d0
           ll = 2

           do j=nr,1,-1
             dq = phe(j,i) * phe(j,i)
             ekk = ekk + (evi - orb(j,i)) * dr(j) * dq * dble(ll) / 3.d0
             ll = 6 - ll
           end do

           ek(i)=ekk

         endif

         ! add kinetic to total, including frozen core kinetic energy
         etot = etot + (ek(i) * occ(i))

       end do

       call getpot(etot,nst,rel,alfa,dl,nr,dr,r,r2,xntot,                 &
     &             phe,ratio,orb,occ,is,nel,nl,nm,no,xnj,rpower,xnum,     &
     &             etot2,iuflag)
       return

      end subroutine

!----------------------------------------------------------------------
      subroutine getpot(etot,nst,rel,alfa,dl,nr,dr,r,r2,                  &
     &    xntot,phe,ratio,orb,occ,is,                                     &
     &    nel, nl, nm, no, xnj, rpower, xnum, etot2, iuflag)

       implicit double precision (a-h,o-z)
       integer iorbs, iside, io2, ijive
       integer lmax, ihmax, nrmax, ntmax, npmax
       integer nst, nr, nel, iuflag
       double precision etot, etot2, rel, alfa, dl, xntot, ratio, xnum
       parameter (iorbs=33,iside=600)
       parameter (io2=iorbs*(iorbs+1)/2)
       parameter (ijive=io2*(io2+1)/2)
       parameter (lmax=4,ihmax=20,nrmax=4000,ntmax=10,npmax=60)
       double precision, dimension(nrmax)        :: dr, r, r2
       double precision, dimension(nrmax,iorbs)  :: phe, orb
       double precision, dimension(iorbs)        :: occ, xnj
       integer, dimension(iorbs)                 :: is, nl, nm, no
       double precision, dimension(nrmax,0:15)   :: rpower
       double precision, dimension(nrmax)        :: xq1, xq2, xq0
       double precision, dimension(0:8,0:8,0:16) :: pin
       double precision, dimension(2)            :: rsp
       double precision, dimension(nrmax) :: xqj0, xqj1, xqj2, xqi0, xqi1, xqi2
       double precision, dimension(0:6,0:6,0:12,-6:6,-6:6) :: cg

       call clebschgordan(nel,nl,cg)
       call getillls(pin)

       ratio1 = 1.d0 - ratio
       do i=1,nel
         do k=1,nr
           orb(k,i) = ratio1 * orb(k,i)
         end do
       end do

       do i=1,nel

         li = nl(i)
         mi = nm(i)

         jstart = i + 1
         if ((xnj(i) < 0.d0) .or. (occ(i) > 1.d0) .or. (abs(alfa) > 0.001d0)) then
           jstart = i
         endif

         do j=jstart,nel

           if ((occ(i) == 0.d0) .and. (occ(j) == 0.d0)) goto 2990

           lj = nl(j)
           mj = nm(j)

           ! direct coulomb

           lmx = 2 * li
           if (li > lj) lmx = 2 * lj

           ! l=0 is monopole or spherical term for direct coulomb. Therefore,
           ! when we have occ(i) or occ(j) greater than one, set lmx=0.

           if ((occ(i) > 1.d0) .or. (occ(j) > 1.d0) .or.                  &
     &         (xnj(i) < 0.d0) .or. (xnj(j) < 0.d0)) lmx=0

           do la=lmx,0,-2
             lap = la + 1
             coeff = dble((li + li + 1) * (lj + lj + 1)) /dble((la+la+1))**2.d0* &
     &             cg(li,li,la,mi,-mi)*cg(lj,lj,la,mj,-mj)*cg(li,li,la,0 , 0 )*cg(lj,lj,la,0 , 0 )
             if (mi + mj .ne. 2 * ((mi + mj) / 2)) coeff = -coeff
             if (i == j) coeff = coeff / 2.d0
             coeffi = occ(i) * coeff
             coeffj = occ(j) * coeff
             ri = ratio * coeffi
             rj = ratio * coeffj
             rc = coeff * occ(i) * occ(j)

             xouti = 0.d0
             xoutj = 0.d0

             do k=1,nr
               xqi0(k) = dr(k) * phe(k,i) * phe(k,i) / 2.d0
               xqi1(k) = xqi0(k) * rpower(k,la)

               if (rpower(k,lap).ne.0.d0) then
                 xqi2(k) = xqi0(k) / rpower(k,lap)
               else
                 xqi2(k) = 0.d0
               endif

               xouti = xouti + xqi2(k)
               xqj0(k) = dr(k) * phe(k,j) * phe(k,j) / 2.d0
               xqj1(k) = xqj0(k) * rpower(k,la)

               if (rpower(k,lap) .ne. 0.d0) then
                 xqj2(k) = xqj0(k) / rpower(k,lap)
               else
                 xqj2(k) = 0.d0
               endif

               xoutj = xoutj + xqj2(k)
             end do

             xinti = xqi1(1)
             xintj = xqj1(1)
             xouti = (2.d0 * xouti) - xqi2(1)
             xoutj = (2.d0 * xoutj) - xqj2(1)

             do k=2,nr

               xinti = xinti + xqi1(k) + xqi1(k-1)
               xouti = xouti - xqi2(k) - xqi2(k-1)
               vali = xouti * rpower(k,la)
               if (rpower(k,lap) .ne. 0.d0) then
                 vali = vali + (xinti / rpower(k,lap))
               endif
               orb(k,j) = orb(k,j) + (ri * vali)

               xintj = xintj + xqj1(k) + xqj1(k-1)
               xoutj = xoutj - xqj2(k) - xqj2(k-1)
               valj = xoutj * rpower(k,la)
               if (rpower(k,lap) .ne. 0.d0) then
                 valj = valj + xintj / rpower(k,lap)
               endif
               orb(k,i) = orb(k,i) + rj * valj

               etot = etot + (rc * ((xqi0(k) * valj) + (xqj0(k) * vali)))

             end do

           end do

           if ((is(i) .ne. is(j)) .and.                                   &
     &         (occ(i) <= 1.d0) .and.                                     &
     &         (occ(j) <= 1.d0) .and.                                     &
     &         (xnj(i) >= 0.d0) .and.                                     &
     &         (xnj(j) >= 0.d0)) goto 2990
           if (abs(alfa) >= 0.001d0) goto 2990

           ! exchange interaction
           lmx = li + lj
           lmin = abs(mi - mj)
           if ((occ(i) > 1.d0) .or. (occ(j) > 1.d0) .or.                  &
     &         (xnj(i) < 0.d0) .or. (xnj(j) < 0.d0)) lmin=0

           do la=lmx,lmin,-2
             lap = la + 1

             coeff = dble((li+li+1)*(lj+lj+1))/dble((la+la+1))**2.d0* &
     &             (cg(li,lj,la,-mi,mj)*cg(li,lj,la,0,0))**2.d0
             if ((occ(i) > 1.d0) .or. (occ(j) > 1.d0) .or.                &
     &           (xnj(i) < 0.d0) .or. (xnj(j) < 0.d0))                    &
     &         coeff = pin(li,lj,la) / 4.d0
             if (i == j) coeff = coeff/2.d0
             coeffi = occ(i) * coeff
             coeffj = occ(j) * coeff
             ri = ratio * coeffi
             rj = ratio * coeffj
             rc = coeff * occ(i) * occ(j)
             xnum2 = xnum * xnum

             xout = 0.d0

             do k=1,nr
               xq0(k) = dr(k) * phe(k,i) * phe(k,j) / 2.d0
               xq1(k) = xq0(k) * rpower(k,la)

               if (rpower(k,lap).ne.0.d0) then
                 xq2(k) = xq0(k) / rpower(k,lap)
               else
                 xq2(k) = 0.d0
               endif

               xout = xout + xq2(k)
             end do

             xint = xq1(1)
             xout = (2.d0 * xout) - xq2(1)

             do k=2,nr
               xint = xint + xq1(k) + xq1(k-1)
               xout = xout - xq2(k) - xq2(k-1)

               if (xq0(k).ne.0.d0) then
                 val = xout * rpower(k,la)
                 if (rpower(k,lap).ne.0.d0) then
                   val = val + xint / rpower(k,lap)
                 endif
                 etot = etot - (2.d0 * xq0(k) * rc * val)
                 xx = phe(k,j) / phe(k,i)

                 if (abs(xx) / xnum > 1.d0) then
                   orb(k,i) = orb(k,i) - (rj * xnum2 / xx * val)
                 else
                   orb(k,i) = orb(k,i) - (rj * xx * val)
                 endif

                 xx = phe(k,i) / phe(k,j)

                 if (abs(xx)/xnum > 1.d0) then
                   orb(k,j) = orb(k,j) - ((ri * xnum2) / (xx * val))
                 else
                   orb(k,j) = orb(k,j) - (ri * xx * val)
                 endif

               endif

             end do  ! k

           end do  ! la
 
         end do

2990   end do

       ! here we compute the charge density, if needed, for treating
       ! exchange/correlation in a local fashion...

       if (abs(alfa) >= 0.001d0) then
         if (alfa > 0.d0) then
           fx = 1.0d0
           fc = 1.0d0
         else
           fx = 1.5d0 * abs(alfa)
           fc = 0.0d0
         endif

       ! note: we don't deal with spin-polarization in local exchange
       ! picture, since local exchange is totally wrong for such
       ! effects, anyway.  local exchange pretends charge density
       ! is paramagnetic.  also, local exchange treats everything
       ! as spherical.

         fourpi = 16.d0 * atan(1.d0)

         do i=1,nr
           xn = 0.d0

           do j=1,nel
             xn = xn + (phe(i,j) * phe(i,j) * occ(j))
           end do

           xn1 = xn / 2.d0
           xn2 = xn / 2.d0
           nst = 2
           call exchcorr(nst,rel,r2(i),xn1,xn2,ex,ec,ux1,ux2,uc1,uc2)
           exc = (fx * ex) + (fc * ec)
           uxc = (fx * ux1) + (fc * uc1)
           etot = etot + (dr(i) * xn * exc)

           do j=1,nel
             orb(i,j) = orb(i,j) + (uxc * ratio)
           end do

         end do

       endif

       do i=1,nr
         if (iuflag.ne.0) then
           jj = 1
 8960      ii = jj
 8965      if (ii == nel) goto 8970
           icond = 0
           if ((no(jj) == no(ii+1)) .and. (nl(jj) == nl(ii+1)) .and. (iuflag == 2)) icond=1
           if ((no(jj) == no(ii+1)) .and. (nl(jj) == nl(ii+1))            &
     &          .and. (is(jj) == is(ii+1)) .and. (iuflag == 1)) icond=1

           if (icond == 1) then
             ii = ii + 1
             goto 8965
           endif

 8970      orba = 0.d0
           div = 0.d0

           do k=jj,ii
             div = div + occ(k)
             orba = orba + (orb(i,k) * occ(k))
           end do

           if (div.ne.0.d0) then
             orba = orba / div

             do k=jj,ii
               orb(i,k) = orba
             end do

           endif

           if (ii.ne.nel) then
             jj = ii + 1
             goto 8960
           endif

         endif

       end do

       return

      end subroutine

!-----------------------------------------------------------------------
      subroutine elsolve(i,occ,n,l,xkappa,xj,zorig,zeff,e,phi,v,          &
     &                   q0,xm1,xm2,nr,r,dr,r2,dl,rel)

       implicit double precision (a-h,o-z)
       integer iorbs, iside, io2, ijive, lmax, ihmax, nrmax, ntmax, npmax
       integer i, n, l, nr
       double precision xkappa, occ, xj, zorig, zeff, e, rel, dl
       parameter (iorbs=33,iside=600)
       parameter (io2=iorbs*(iorbs+1)/2)
       parameter (ijive=io2*(io2+1)/2)
       parameter (lmax=4,ihmax=20,nrmax=4000,ntmax=10,npmax=60) 
       double precision, dimension(nrmax) :: phi,phi2,v,q0,xm1,xm2,r,dr,r2

       el = -zorig * zorig / dble(n * n)
       eh = 0.d0
       etol = 0.0000000001d0
 155   e = (el + eh) / 2.d0
       istop = 0

       call integ(e,l,xkappa,n,nn,istop,ief,x0,phi,zeff,v,q0,xm1,         &
     &            xm2,nr,r,dr,r2,dl,rel)

       if (nn < n-l-1) ief = -1

       if (ief .ne. 1) then
         el = e
         if (el > -0.001d0) then
           write (6,*) 'MIXING TOO STRONG FOR LEVEL : ', i
           return
         endif
       endif

       if (ief .ne. -1) eh = e
       if (eh-el > etol) goto 155
       if (abs(abs(xj)-abs(dble(l))) > 0.25d0) then
         call augment(e,l,xj,phi,v,nr,r,dl)
       endif
       aa = 0.d0

       do j=1,nr
         aa = aa + (phi(j) * phi(j) * dr(j))
       end do

       xnorm = sqrt(aa)

       do j=1,nr
         phi(j) = phi(j) / xnorm
       end do

       return

      end subroutine

!--------------------------------------------------------------------------
      subroutine augment(e,l,xj,phi,v,nr,r,dl)

       implicit double precision (a-h,o-z)
       integer iorbs, iside, io2, ijive, lmax, ihmax, nrmax, ntmax, npmax
       integer i, n, l, idoflag, nr
       double precision xkappa, occ, xj, aorig, zeff, e, rel, dl
       parameter (iorbs=33,iside=600)
       parameter (io2=iorbs*(iorbs+1)/2)
       parameter (ijive=io2*(io2+1)/2)
       parameter (lmax=4,ihmax=20,nrmax=4000,ntmax=10,npmax=60) 
       double precision, dimension(nrmax) :: phi,phi2,v,r

       c = 137.038d0
       cc = c * c
       c2 = cc + cc
       xkappa = -1

       if (abs(xj) > dble(l) + 0.25d0) xkappa = -l - 1
       if (abs(xj) < dble(l) - 0.25d0) xkappa = l

       do j=4,nr-3
         if (phi(j).ne.0.d0) then
           g0 = phi(j)
           ga = (phi(j+1) - phi(j-1))
           gb = (phi(j+2) - phi(j-2)) / 2.d0
           gc = (phi(j+3) - phi(j-3)) / 3.d0
           gg = (((1.5d0 * ga) - (0.6d0 * gb) + (0.1d0 * gc)) / (2.d0 * dl) + xkappa * g0) / r(j)
           f0 = c * gg / (e - v(j) + c2)
           phi2(j) = sqrt((g0 * g0) + (f0 * f0))
           if (g0 < 0.d0) phi2(j) = -phi2(j)
         else
           phi2(j) = phi(j)
         endif

       end do

       do j=1,3
         phi2(j) = phi(j) * phi(4) / phi2(4)
       end do

       do j=1,nr
         phi(j) = phi2(j)
       end do

       return

      end subroutine

!-----------------------------------------------------------------------
      subroutine setqmm(i, orb, l, ns, idoflag, v, zeff, zorig, rel,      &
     &                  nr, r, r2, dl, q0, xm1, xm2, njrc, vi)

       implicit double precision (a-h,o-z)
       integer iorbs, iside, io2, ijive
       integer lmax, ihmax, nrmax, ntmax, npmax
       integer i, n, l, idoflag, nr, ns
       double precision zeff, zorig, rel, dl
       parameter (iorbs=33,iside=600)
       parameter (io2=iorbs*(iorbs+1)/2)
       parameter (ijive=io2*(io2+1)/2)
       parameter (lmax=4,ihmax=20,nrmax=4000,ntmax=10,npmax=60) 
       double precision, dimension(nrmax) :: v, r, r2, q0, xm1, xm2
       double precision, dimension(nrmax,iorbs) :: orb
       double precision, dimension(nrmax,7), intent(inout) :: vi
       integer, dimension(4) :: njrc

       c = 137.038d0
       alpha = rel / c
       aa = alpha * alpha
       a2 = aa / 2.d0

       lp = l + 1
       lpx = lp
       if (lp > 4) lpx = 4
       lp2 = l + l + 1
       if (lp2 > 7) lp2 = 7
       zeff = zorig
       if (njrc(lpx) > 0) zeff = 0.d0
       zaa = zeff * aa
       za2 = zeff * a2

       if (idoflag.ne.0) then
         if (njrc(lpx) == 0) then

           if (idoflag == 1) then
             do j=1,nr
               v(j) = -zeff / r(j) + orb(j,i)
             end do
           endif

           do j=2,nr-1
             dvdl = (orb(j+1,i) - orb(j-1,i)) / (2.d0 * dl)
             ddvdrr = ((orb(j+1,i)+orb(j-1,i)-2.d0*orb(j,i) )/(dl*dl)-dvdl)/r2(j)
             xm1(j) = -(a2 * dvdl / r(j)) - (za2 / r2(j))
             xm2(j) = -(a2 * ddvdrr) + (zaa / r2(j) / r(j))
           end do

           xm1(nr) = xm1(nr-1)
           xm2(nr) = xm2(nr-1)
           xm1(1) = xm1(2) + za2 / r2(2) - za2 / r2(1)
           xm2(1) = xm2(2) - zaa / r2(2) / r(2) + zaa / r2(1) / r(1)

         else
           if (idoflag == 1) then
             do j=1,nr
               v(j) = vi(j,lp2) + orb(j,i)
             end do
           endif

           do j=2,nr-1
             dvdl = (v(j+1) - v(j-1)) / (2.d0*dl)
             ddvdrr = ((v(j+1) + v(j-1) - 2.d0 * v(j)) / (dl * dl) - dvdl) / r2(j)
             xm1(j) = -a2 * dvdl / r(j)
             xm2(j) = -a2 * ddvdrr
           end do

           xm1(nr) = xm1(nr-1)
           xm2(nr) = xm2(nr-1)
           xm1(1) = xm1(2)
           xm2(1) = xm2(2)

         endif

       endif

       ! figure the (Desclaux-Numerov) effective potential.

       xlb = (dble(l) + 0.5d0)**2.d0 / 2.d0

       do j=1,nr
         vj = v(j)
         q0(j) = vj * (1.d0 - (a2 * vj)) + (xlb / r2(j))
       end do

       return

      end subroutine

!----------------------------------------------------------------------
      subroutine initiali(zorig,nr,rmin,rmax,r,dr,r2,dl,njrc,xntot,nel)

       implicit double precision (a-h,o-z)
       integer iorbs, iside, io2, ijive
       integer lmax, ihmax, nrmax, ntmax, npmax
       integer i, n, l, idoflag, nr, ns, nel
       double precision zorig, rmin, rmax, dl, xntot
       parameter (iorbs=33,iside=600)
       parameter (io2=iorbs*(iorbs+1)/2)
       parameter (ijive=io2*(io2+1)/2)
       parameter (lmax=4,ihmax=20,nrmax=4000,ntmax=10,npmax=60) 
       double precision, dimension(nrmax) :: r, dr, r2
       integer, dimension(4) :: njrc

       ! enter Z, NR
       read (5,*) zorig,nr
       rmin = 0.0001d0 / zorig
       rmax = 800.d0 / sqrt(zorig)
       call setgrid(nr, rmin, rmax, r, dr, r2, dl)

       do j=1,4
         njrc(j) = 0
       end do

       xntot = 0.d0
       nel = 0

       return

      end subroutine

!---------------------------------------------------------------------------
      subroutine setgrid(nr, rmin, rmax, r, dr, r2, dl)

       implicit double precision (a-h,o-z)
       integer iorbs, iside, io2, ijive
       integer lmax, ihmax, nrmax, ntmax, npmax
       integer, intent(in) :: nr
       double precision rmin, rmax, dl
       parameter (iorbs=33,iside=600)
       parameter (io2=iorbs*(iorbs+1)/2)
       parameter (ijive=io2*(io2+1)/2)
       parameter (lmax=4,ihmax=20,nrmax=4000,ntmax=10,npmax=60) 
       double precision, dimension(nrmax) :: r,dr,r2

       ratio = rmax / rmin
       dl = log(ratio) / dble(nr)
       xratio = exp(dl)
       xr1 = sqrt(xratio) - sqrt(1.d0 / xratio)

       do i=1,nr
         r(i) = rmin * xratio**dble(i)
         dr(i) = r(i) * xr1
         r2(i) = r(i) * r(i)
       end do

       return

      end subroutine

!-----------------------------------------------------------------------------
      subroutine integ(e, l, xkappa, n, nn, istop, ief, x0, phi, z, v,     &
     &                 q0, xm1, xm2, nr, r, dr, r2, dl, rel)

       implicit double precision (a-h,o-z)
       integer iorbs, iside, io2, ijive
       integer lmax, ihmax, nrmax, ntmax, npmax
       integer l, n, nn, istop, ief, nr
       double precision e, xkappa, x0, z, rel, dl
       parameter (iorbs=33,iside=600)
       parameter (io2=iorbs*(iorbs+1)/2)
       parameter (ijive=io2*(io2+1)/2)
       parameter (lmax=4,ihmax=20,nrmax=4000,ntmax=10,npmax=60) 
       double precision, dimension(nrmax) :: phi, v
       double precision, dimension(nrmax) :: q0, xm1, xm2, r, dr, r2

       dl2 = dl * dl / 12.d0
       dl5 = 10.d0 * dl2
       c = 137.038d0
       alpha = rel / c
       za2 = z * z * alpha * alpha
       a2 = alpha * alpha / 2.d0
       xl = l
       xlp = l + 1
       xl2 = 0.5d0 + xl
       xl4 = xl2 * xl2

       ! set up the leading power.
       ! adjust for Desclaux's implementation of Numerov.

       if (rel == 0.d0) then
         ss = xlp
       else
         rtest = 1.d0 - za2

         if (rtest < 0.d0) then
           write (6,*) 'Z>137 IS TOO BIG.'
           return
         endif  

         ss = sqrt(rtest)
       endif

       ss2 = ss - 0.5d0

       ! set ief to -1 if energy is too low, +1 if too high.

       ief = 0

       ! see Desclaux and documentation for origin of the equations below.
       ! here, we set up the first two points.

       t = e - v(1)
       xm0 = 1.d0 + (a2 * t)
       tm = xm0 + xm0
       xmx = xm1(1) / xm0
       xk0=r2(1)*(tm*t-xmx*(xkappa/r(1)+0.75d0*xmx)+xm2(1)/tm)-xl4
       dk0=1.d0+dl2*xk0
       p0=dk0
       phi(1)=p0*sqrt(xm0*r(1))/dk0

       t=e-v(2)
       xm=1.d0+a2*t
       tm=xm+xm
       xmx=xm1(2)/xm
       xk2=r2(2)*(tm*t-xmx*(xkappa/r(2)+0.75d0*xmx)+xm2(2)/tm)-xl4
       dk2=1.d0+dl2*xk2
       p1=dk2*((r(2)/r(1))**ss2-(r(2)-r(1))*z/xlp)*sqrt(xm0/xm)
       phi(2)=p1*sqrt(xm*r(2))/dk2

       ! if istop is set, the we know to stop there.  If it is zero, it
       ! shall then be set to the classical turning point.

       is0 = istop
       if (istop == 0) then

         do j=nr-1,2,-1
           if (e > v(j)) goto 15
         end do

         ief = -1
         return

 15      istop = j
       endif

       ! initialize number of nodes, and determine the ideal number.
       nn = 0
       nnideal = n - l - 1

       ! integrate out.
       ! count nodes, and stop along the way if there are too many.
       do i=3,istop+2
         t = e - v(i)
         xm = 1.d0 + (a2 * t)
         tm = xm + xm
         xmx = xm1(i) / xm
         p2 = (2.d0 - (dl5 * xk2)) * p1 / dk2 - p0
         xk2 = r2(i)*((tm * t) - xmx * ((xkappa / r(i)) + 0.75d0 * xmx)+xm2(i)/tm)-xl4
         dk2 = 1.d0 + (dl2 * xk2)
         phi(i) = p2 * sqrt(xm * r(i)) / dk2
         if (abs(p2) > 10000000000.d0) then

           do j=1,i
             phi(j) = phi(j) / p2
           end do

           p0 = p0 / p2
           p1 = p1 / p2
           p2 = p2 / p2
         endif

         if (p2 * p1 < 0.d0) then
           nn = nn + 1
           if (nn > nnideal) then
             ief = 1
             return
           endif
         endif

         p0 = p1
         p1 = p2

       end do

       if (istop > 0) then
         psip2 = (phi(istop+2) - phi(istop-2))
         psip1 = (phi(istop+1) - phi(istop-1))
         psip = ((8.d0 * psip1) - psip2) / (12.d0 * dl * r(istop))
         x0 = psip / phi(istop)
       endif

       if (is0.ne.0) return

       do i=istop+3,nr-1
         t = e - v(i)
         xm = 1.d0 + (a2 * t)
         tm = xm + xm
         xmx = xm1(i) / xm
         p2 = (2.d0 - (dl5 * xk2)) * p1 / dk2 - p0

         if (p2 / p1 > 1.d0) then
           ief = -1
           return
         endif

         xk2 = r2(i)*(tm*t-xmx*(xkappa/r(i)+0.75d0*xmx)+xm2(i)/tm)-xl4
         dk2 = 1.d0 + (dl2 * xk2)
         phi(i) = p2 * sqrt(xm * r(i)) / dk2

         if (abs(p2) > 10000000000.d0) then
           do j=1,i
             phi(j) = phi(j) / p2
           end do

           p0 = p0 / p2
           p1 = p1 / p2
           p2 = p2 / p2
         endif

         if (p2 * p1 < 0.d0) then
           nn = nn + 1
           if (nn > nnideal) then
             ief = 1
             return
           endif
         endif

         p0 = p1
         p1 = p2

       end do

       return

      end subroutine

!-------------------------------------------------------------------------
!  routine to generate Clebsch-Gordan coefficients, in the form of
!  cg(l1,l2,L,m1,m2) = <l1,m1;l2,m2|L,m1+m2>, according to Rose's
!  'Elementary Theory of Angular Momentum', p. 39, Wigner's formula.
!  those coefficients listed are only those for which l1>=l2.
!  coefficients known to be zero because of either the L or M
!  selection rules are not computed, and should not be sought.
      subroutine clebschgordan(nel, nl, cg)

       implicit double precision (a-h, o-z)
       integer iorbs, iside, io2, ijive
       integer lmax, ihmax, nrmax, ntmax, npmax
       integer lmx, nel
       parameter (iorbs=33,iside=600)
       parameter (io2=iorbs*(iorbs+1)/2)
       parameter (ijive=io2*(io2+1)/2)
       parameter (lmax=4,ihmax=20,nrmax=4000,ntmax=10,npmax=60) 
       integer, dimension(iorbs) :: nl
       double precision, dimension(0:6,0:6,0:12,-6:6,-6:6) :: cg
       double precision, dimension(0:32) :: si,fa

       lmx = 0
       do i=1,nel
         if (nl(i) > lmx) lmx = nl(i)
       end do

       si(0) = 1.d0
       fa(0) = 1.d0

       do i=1,32
         si(i) = -si(i-1)
         fa(i) = dble(i) * fa(i-1)
       end do

       do l1=0,lmx

         do l2=0,l1

           do m1=-l1,l1

             do m2=-l2,l2
               m3 = m1 + m2
               lmin = iabs(l1-l2)
               if (lmin < iabs(m3)) then
                 lmin = iabs(m3)
               endif

               do l3=lmin,l1+l2
                 prefactor = dble((2 * l3) + 1)
                 prefactor = prefactor * fa(l3+l1-l2) / fa(l1+l2+l3+1)
                 prefactor = prefactor * fa(l3-l1+l2) / fa(l1-m1)
                 prefactor = prefactor * fa(l1+l2-l3) / fa(l1+m1)
                 prefactor = prefactor * fa(l3+m3) / fa(l2-m2)
                 prefactor = prefactor * fa(l3-m3) / fa(l2+m2)
                 prefactor = sqrt(prefactor)
                 sum = 0.d0
                 numax = l3 - l1 + l2
                
                 if ((l3+m3) < numax) then
                   numax = l3+m3
                 endif
                
                 numin = 0
                
                 if (l1-l2-m3 < numin) then
                   numin = -(l1 - l2 - m3)
                 endif

                 do nu=numin,numax
                   sum = sum+(si(nu+l2+m2)/fa(nu))                        &
     &                  *fa(l2+l3+m1-nu)*fa(l1-m1+nu)                     &
     &                  /fa(l3-l1+l2-nu)/fa(l3+m3-nu)/fa(nu+l1-l2-m3)
                 end do  ! nu

                 cg(l1,l2,l3,m1,m2) = prefactor*sum
                 cg(l2,l1,l3,m2,m1) = si(l1+l2+l3)*prefactor*sum

               end do  ! l3

             end do  ! m2

           end do  ! m1

         end do  ! l2

       end do  ! l1

       return

 !52    format (1x,i3,a3,i3)

      end subroutine

!-----------------------------------------------------------------------
      subroutine pseudo(etot,nst,rel,alfa,nr,rmin,rmax,r,dr,r2,dl,        &
     &  phe,orb,njrc,vi,zorig,xntot,nel,                                  &
     &  no,nl,xnj,ev,occ,is,ek,iuflag,vctab)

       implicit double precision (a-h,o-z)
       integer iorbs, iside, io2, ijive
       integer lmax, ihmax, nrmax, ntmax, npmax
       integer nst, nr, nel, iuflag
       double precision etot, rel, alfa, rmin, rmax, dl, zorig, xntot
       parameter (iorbs=33,iside=600)
       parameter (io2=iorbs*(iorbs+1)/2)
       parameter (ijive=io2*(io2+1)/2)
       parameter (lmax=4,ihmax=20,nrmax=4000,ntmax=10,npmax=60) 
       double precision, dimension(nrmax) :: r, dr, r2, q0, xm1, xm2, phi, v
       double precision, dimension(nrmax,7) :: vi
       double precision, dimension(nrmax,iorbs) :: phe
       integer, dimension(4) :: njrc, njrcdummy
       integer, dimension(iorbs) :: no, nl, nm, is
       double precision, dimension(iorbs) :: xnj, ek, ev, occ
       double precision, dimension(nrmax,0:7) :: rpower
       double precision, dimension(nrmax,iorbs) :: orb
       double precision, dimension(nrmax,0:3) :: vctab

       do i=1,nel
         nm(i) = 0
       end do

       njrcdummy(1) = njrc(1)
       njrcdummy(2) = njrc(2)
       njrcdummy(3) = njrc(3)
       njrcdummy(4) = njrc(4)

       ! ENTER NP,CORPOL,RNORM
       read (5,*) np,corpol,rnorm
       xntot = 0.d0

       do i=np,nel
         write (6,42) 'l=', nl(i), ' ...'
 42      format(1x,1a2,1i1,1a4)
         lp2 = nl(i) + nl(i) + 1
         e = ev(i)

         do j=1,nr
           orb(j,i) = orb(j,i) + vctab(j,nl(i))
         end do

         idoflag = 1
         ns = 1
         call setqmm(i,orb,nl(i),ns,idoflag,vi(1,lp2),zeff,zorig,rel,     &
     &               nr, r, r2, dl, q0, xm1, xm2, njrcdummy, vi)

         ! set j orbitals to zero
         do j=1,nr
           orb(j,i) = 0.d0
         end do

       ! you can replace subroutine pseudize with any type of PP
       ! generation you want... however, kleinman-bylanderization
       ! would take more coding...

         call pseudize(i,orb,e,nl(i),xnj(i),no(i),njrc,zeff,vi(1,lp2),    &
     &                 q0,xm1,xm2,nr,rmin,rmax,r,dr,r2,dl,rel)

         ! WE HAVE GOT THUS FAR...
         no(i) = nl(i) + 1
         ruse = 0.d0
         xkappa = -1.d0
         call elsolve(i, occ(i), no(i), nl(i), xkappa, xnj(i),            &
     &                zorig, zeff, ev(i), phe(1,i), vi(1,lp2),            &
     &                q0, xm1, xm2, nr, r, dr, r2, dl, ruse)
         write (6,*) nl(i),ev(i)
         xntot = xntot + occ(i)

         if (lp2 == 1) exit  ! loop

         do j=1,nr
           vi(j,lp2-1) = vi(j,lp2)
         end do

       end do

       ! everything is pseudized
       do i=np,nel
         inew = 1 + i - np
         no (inew) = no (i)
         nl (inew) = nl (i)
         nm (inew) = nm (i)
         xnj(inew) = xnj(i)
         is (inew) = 1
         ev (inew) = ev (i)
         occ(inew) = occ(i)

         do j=1,nr
           phe(j,inew) = phe(j,i)
         end do

       end do

       nel = 1 + nel - np

       do i=0,7
         xi = i
         do k=1,nr
           rpower(k,i) = r(k)**xi
         end do
       end do

       ! everything is scaled down...ready for unscreening
       xnum = 100.d0
       ratio = 1.d0
       call getpot(etot,nst,rel,alfa,dl,nr,dr,r,r2,                       &
     &             xntot,phe,ratio,orb,occ,is,                            &
     &             nel,nl,nm,no,xnj,rpower,xnum,etot2,iuflag)

       ! screening effects in pseudo atom computed...
       do k=1,nel
         lp2 = nl(k)+nl(k)+1

         do j=1,nr
           vi(j,lp2) = vi(j,lp2)  -orb(j,k)
           if (lp2 > 1) then
             vi(j,lp2-1) = vi(j,lp2-1)-orb(j,k)
           endif
         end do

       end do

       ! we got past the unscreening...

       do j=1,nr
         vl = (vi(j,2)+2.d0*vi(j,3))/3.d0
         vso = 2.d0*(vi(j,3)-vi(j,2))/3.d0
         vi(j,2) = vso
         vi(j,3) = vl
         vl = (2.d0*vi(j,4)+3.d0*vi(j,5))/5.d0
         vso = 2.d0*(vi(j,5)-vi(j,4))/5.d0
         vi(j,4) = vso
         vi(j,5) = vl
         vl = (3.d0*vi(j,6)+4.d0*vi(j,7))/7.d0
         vso = 2.d0*(vi(j,7)-vi(j,6))/7.d0
         vi(j,6) = vso
         vi(j,7) = vl
       end do

       rel = 0.d0

       ! got past the spin-orbit jazz
       izuse = abs(vi(nr-2,1)*r(nr-2))+0.5d0
       zuse = izuse

       do k=1,nr
         if (r(k) > rnorm) then
           videal = -zuse / r(k) - corpol / (2.d0*r(k)**4.d0)
           vi(k,1) = videal
           vi(k,3) = videal
           vi(k,5) = videal
           vi(k,7) = videal
           vi(k,2) = 0.d0
           vi(k,4) = 0.d0
           vi(k,6) = 0.d0
         endif

       end do

       ! we got to the end

       return

      end subroutine

!----------------------------------------------------------------------
      subroutine parabreg(f, fp, fpp, rf, vf)

       implicit double precision (a-h,o-z)
       double precision, dimension(3) :: rf, vf
       double precision f, fp, fpp

       f = vf(2)
       r21 = rf(2) - rf(1)
       r32 = rf(3) - rf(2)
       v21 = vf(2) - vf(1)
       v32 = vf(3) - vf(2)
       fp = (v21 + v32) / (r21 + r32)
       fpp = ((v32 / r32) - (v21 / r21)) / ((r21 + r32) / 2.d0)

       return

      end subroutine

!----------------------------------------------------------------------
      double precision function hb(x, factor)

       implicit double precision (a-h, o-z)
       double precision, intent(in) :: x, factor

       if (x > 3.d0) hb = 0.d0
       if (x <= 3.d0) hb = 0.01d0**((dsinh(x / factor) / 1.1752d0)**2.d0)

       return

      end function
!----------------------------------------------------------------------
      subroutine fitx0(i,orb,rcut,njrc,e,l,xj,n,jrt,xideal,phi,           &
     &                   zeff,v,q0,xm1,xm2,nr,r,dr,r2,dl,rel,factor)

       implicit double precision (a-h,o-z)
       integer, intent(in) :: i, n, jrt, nr, l
       double precision, intent(in) :: xj, xideal, factor
       double precision, intent(in) :: e, rcut, zeff, rel, dl
       integer iorbs, iside, io2, ijive
       integer lmax, ihmax, nrmax, ntmax, npmax
       integer idoflag, ns, ii
       double precision vl, vh, xkappa, xla, xerror, xactual
       double precision dxdla, vmaybe
       parameter (iorbs=33, iside=600)
       parameter (io2=iorbs*(iorbs+1)/2)
       parameter (ijive=io2*(io2+1)/2)
       parameter (lmax=4, ihmax=20, nrmax=4000, ntmax=10, npmax=60) 
       double precision, dimension(nrmax) :: r, dr, r2
       integer, dimension(4), intent(in) :: njrc
       double precision, dimension(nrmax,iorbs) :: orb
       double precision, dimension(nrmax) :: q0, xm1, xm2, phi, v
       !dimension dummy(nrmax,7)

       vl = -1000000.d0
       vh = +1000000.d0

       do while(.True.)
         idoflag = 2
         ns = 1
         xkappa = -1.d0

         call setqmm(i, orb, l, ns, idoflag, v, zeff, dummy, rel,         &
     &               nr, r, r2, dl, q0, xm1, xm2, njrc, dummy)

         call integ(e, l, xkappa, n, nn, jrt, ief, xactual, phi, zeff, v, &
     &              q0, xm1, xm2, nr, r, dr, r2, dl, rel)

         if (nn.ne.0) then
           vl = v(1)
           xla = 1.d0
         else

           if (xactual > xideal) then
             vh = v(1)
           else
             vl = v(1)
           endif

           xerror = xideal - xactual
           if (abs(xerror) < 0.000000001d0) return
           dxdla = 0.d0

           do ii = 1,jrt
             dxdla = dxdla + (dr(ii) * phi(ii) * phi(ii) * hb(r(ii) / rcut,factor))
           end do

           dxdla = 2.d0 * dxdla / (phi(jrt) * phi(jrt))
           xla = xerror / dxdla
         endif

         vmaybe = v(1)+xla
         if ((vmaybe > vh).or.(vmaybe < vl)) then
           xla = (vl + vh) / 2.d0 - v(1)
         endif

         do ii=1,jrt-1
           v(ii) = v(ii) + (xla * hb(r(ii) / rcut,factor))
         end do

       end do

      end subroutine

!----------------------------------------------------------------------
      subroutine pseudize(i, orb, ev, l, xj, n, njrc, zeff, v,            &
     &                    q0, xm1, xm2, nr, rmin, rmax, r, dr, r2, dl, rel)

       implicit double precision (a-h,o-z)
       integer, intent(inout) :: i, l, n, nr
       double precision, intent(inout) :: ev, xj, rel, dl
       double precision, intent(inout) :: rmax, rmin, zeff 
       integer iorbs, iside, io2, ijive
       integer lmax, ihmax, nrmax, ntmax, npmax
       parameter (iorbs=33,iside=600)
       parameter (io2=iorbs*(iorbs+1)/2)
       parameter (ijive=io2*(io2+1)/2)
       parameter (lmax=4,ihmax=20,nrmax=4000,ntmax=10,npmax=60) 
       integer, dimension(4), intent(inout) :: njrc
       double precision, dimension(3) :: rf, vf
       double precision, dimension(nrmax) :: r, dr, r2, phi, v
       double precision, dimension(nrmax) :: phi0, yl, vraw, q0, xm1, xm2
       double precision, dimension(nrmax,iorbs) :: orb

       lp = l+1
       xkappa = -1.d0
       istop = nr
 40    istop = istop - 1

       if (ev <= q0(istop)) goto 40

       call integ(ev, l, xkappa, n, nn, istop, ief, xdummy, phi, zeff, v, &
     &            q0, xm1, xm2, nr, r, dr, r2, dl, rel)

       ! ENTER THE CUTOFF RADIUS, AND FACTOR.
       read (5,*) rcut,factor
       if (rcut < 0.d0) then
         xnodefrac = -rcut
         j = istop
 55      j = j-1
         if (phi(j-1) / phi(j) > 1.d0) goto 55

         if (n > l+1) then
           k = j
 60        k = k - 1
           if (phi(k-1) / phi(k) > 0.d0) goto 60
         else
           k = 1
         endif

         rcut = r(k) + (xnodefrac * (r(j) - r(k)))
       endif

       jrc = 1.d0 + (dble(nr - 1) * log(rcut / rmin) / log(rmax / rmin))
       rcut = r(jrc)
       rtest = 2.d0 * rcut
       jrt = 1.d0 + (dble(nr - 1) * log(rtest / rmin) / log(rmax / rmin))
       njrc(lp) = jrt
       rtest = r(jrt)
       switch = phi(jrt) / abs(phi(jrt))

       write (6,92) 'RCUTOFF = ', rcut, '  JRC = ', jrc
       write (6,92) 'RTEST   = ', rtest, '  JRT = ', jrt
 92    format (1x,1a10,1f8.4,1a8,1i5)
 94    format (1x,2d15.8)

       call integ(ev, l, xkappa, n, nn, jrt, ief, x00, phi, zeff, v,      &
     &            q0, xm1, xm2, nr, r, dr, r2, dl, rel)

       do ii = 1,jrt
         phi(ii) = phi(ii) / phi(jrt)
       end do

       xn00 = 0.d0

       do ii = 1,jrt-1
         xn00 = xn00 + (dr(ii) * phi(ii) * phi(ii))
       end do

       xn00 = xn00 + (dr(jrt) * phi(jrt) * phi(jrt) / 2.d0)
       de = 0.0001d0
       ee = ev + (de / 2.d0)
       call integ(ee, l, xkappa, n, nn, jrt, ief, xp, phi, zeff, v,       &
     &            q0, xm1, xm2, nr, r, dr, r2, dl, rel)

       ee = ev-de/2.d0
       call integ(ee,l,xkappa,n,nn,jrt,ief,xm,phi,zeff,v,                 &
     &            q0,xm1,xm2,nr,r,dr,r2,dl,rel)

       c00 = (xm - xp) / (2.d0 * de)

       write (6,94) c00, x00
       write (6,94) xn00

       ruse = 0.d0
       v0 = v(jrc)
       dvdl   = (8.d0*(v(jrc+1)-v(jrc-1))-(v(jrc+2)-v(jrc-2)))/(12.d0*dl)
       ddvdll = (16.d0*(v(jrc+1)+v(jrc-1))                                &
     &         -30.d0*v(jrc)-v(jrc+2)-v(jrc-2))/(12.d0*dl*dl)
       dldr = 1.d0/r(jrc)
       ddldrr = -1.d0/r2(jrc)
       v1 = dvdl*dldr
       v2 = dvdl*ddldrr+ddvdll*dldr*dldr
       b4 = (v2*rcut-v1)/(8.d0*rcut**3.d0)
       b2 = (v1-4.d0*b4*rcut**3.d0)/(2.d0*rcut)
       b0 = v0-b4*rcut**4.d0-b2*rcut**2.d0

       do ii = 1,jrc
         rr = r(ii)
         v(ii) = b0+b2*rr**2.d0+b4*rr**4.d0
       end do

       call fitx0(i,orb,rcut,njrc,ev,l,xj,lp,jrt,x00,phi,zeff,v,          &
     &            q0,xm1,xm2,nr,r,dr,r2,dl,ruse,factor)

       do ii = 1,jrt
         phi0(ii) = phi(ii)
         vraw(ii) = v(ii)
       end do

       xi0 = 0.d0
       xi1 = 0.d0
       xi2 = 0.d0

       do ii = 1,jrt
         f = hb(r(ii) / rcut,factor)
         ph2 = dr(ii) * phi0(ii) * phi0(ii)
         xi0 = xi0 + ph2

         if (ii <= jrt) then
           xi1 = xi1 + (ph2 * f)
           xi2 = xi2 + (ph2 * f * f)
         endif

       end do

       ph2 = phi0(jrt) * phi0(jrt)
       xi0 = xi0 / ph2
       xi1 = xi1 / ph2
       xi2 = xi2 / ph2
       quant = (xi1 * xi1) + (xi2 * (c00 - xi0))

       if (quant > 0.d0) then
         deltal = (sqrt((xi1 * xi1) + (xi2 * (c00 - xi0))) - xi1) / xi2
       else
         deltal = (c00 - xi0) / (2.d0 * xi1)
       endif

       write (6,222) 'DELTAL = ', deltal
 222   format (1x,1a9,1f11.8)

 225   do ii = 1,jrt
         yl (ii) = phi0(ii) * hb(r(ii) / rcut, factor)
         phi(ii) = phi0(ii) + (deltal * yl(ii))

         if (phi(ii) < 0.d0) then
           write (6,*) 'BIG TROUBLE!!! CROSS AXIS!!!'
           return
         endif

       end do

       do ii = 1,jrt-1
         if ((phi(ii) == 0.).or.(yl(ii) == 0.)) exit  ! loop
         jj = ii
         if (ii == 1) jj = 2

         do j = jj-1,jj+1
           rf(2+j-jj) = r(j)
           vf(2+j-jj) = hb(r(j) / rcut, factor)
         end do

         call parabreg(f, fp, fpp, rf, vf)

         do j = jj-1,jj+1
           vf(2+j-jj) = phi0(j)
         end do

         call parabreg(psi, psip, psipp, rf, vf)
         v(ii) = vraw(ii)+(1.d0-phi0(ii)/phi(ii))*(2.d0*psip/psi*fp/f+fpp/f)/2.d0    
       end do

       call fitx0(i,orb,rcut,njrc,ev,l,xj,lp,jrt,x00,phi,zeff,v,          &
     &            q0,xm1,xm2,nr,r,dr,r2,dl,ruse,factor)
       call integ(ev,l,xkappa,n,nn,jrt,ief,x0,phi,zeff,v,                 &
     &            q0,xm1,xm2,nr,r,dr,r2,dl,ruse)

       do ii = 1,jrt
         phi(ii) = phi(ii)/phi(jrt)
       end do

       xn0 = 0.d0

       do ii = 1,jrt-1
         xn0 = xn0 + (dr(ii) * phi(ii) * phi(ii))
       end do

       xn0 = xn0+dr(jrt)*phi(jrt)*phi(jrt)/2.d0
       de = 0.0001d0
       ee = ev+de/2.d0

       call integ(ee,l,xkappa,n,nn,jrt,ief,xp,phi,zeff,v,                 &
     &            q0,xm1,xm2,nr,r,dr,r2,dl,ruse)
       ee = ev-de/2.d0
       call integ(ee,l,xkappa,n,nn,jrt,ief,xm,phi,zeff,v,                 &
     &            q0,xm1,xm2,nr,r,dr,r2,dl,ruse)
       c0 = (xm-xp)/(2.d0*de)

       write (6,94) c0,x0
       write (6,94) xn0

       if (abs(c0-c00) >= 0.000000001d0) then
         dqddel = 2.*(xi1+deltal*xi2)
         deltal = deltal+(c00-c0)/dqddel
         goto 225
       endif

       write (6,94)  c0,x0
       write (6,*) 'NCPP ACHIEVED !!!'

       return
      end subroutine

!-----------------------------------------------------------------------
      subroutine fourier(nr, r, dr, r2, vi)

       implicit double precision (a-h,o-z)
       integer nr
       integer iorbs, iside, io2, ijive
       integer lmax, ihmax, nrmax, ntmax, npmax
       parameter (iorbs=33,iside=600)
       parameter (io2=iorbs*(iorbs+1)/2)
       parameter (ijive=io2*(io2+1)/2)
       parameter (lmax=4,ihmax=20,nrmax=4000,ntmax=10,npmax=60) 
       double precision, dimension(nrmax) :: r, dr, r2, a, v1, v2
       double precision, dimension(nrmax,7) :: vi

       do l=0,2
         lp2 = l + l + 1
         dl = log(r(2) / r(1))
         dl1 = 12.d0 * dl
         dl2 = 12.d0 * dl * dl

         do i = 1,nr
           a(i) = r(i) * vi(i,lp2)
         end do

         do i = 3,nr-2
           al =(-(a(i+2)-a(i-2))+ 8.d0*(a(i+1)-a(i-1)))/dl1

           ar = al / r(i)
           v1(i) = ar
         end do

         open (unit=20+l,status='unknown')

         do ii = 1,200
           q = dble(ii) / 10.d0
           vq = 0.d0

           do i = 3,nr-2
             vq = vq + (dr(i) * dcos(q * r(i)) * v1(i))
           end do

           write (20+l,*) q,vq
         end do

         close(unit=1)

       end do

       return

      end subroutine

!-----------------------------------------------------------------------
      subroutine getillls(pin)

       implicit double precision (A-H,O-Z)
       integer l, m, n, ia, ib, ic, nn, mm, ll
       double precision af, bf, xcf, xi
       double precision, dimension(0:40) :: fa, si
       double precision, dimension(0:8,0:8,0:16), intent(inout) :: pin

       fa(0) = 1.d0
       si(0) = 1.d0

       do I = 1,32
         fa(I) = dble(I)*fa(I-1)
         si(I) = -si(I-1)
       end do

       do L = 0,8

         do M = 0,8

           do N = M+L,0,-2
             XI = 0.d0
             XF = 2.d0/2.d0**dble(N+L+M)
             nn = (N+1)/2
             MM = (M+1)/2
             LL = (L+1)/2

             do ia = nn,N
               AF = si(ia)*fa(ia+ia)/fa(ia)/fa(N-ia)/fa(ia+ia-N)

               do ib = LL,L
                 BF = si(ib)*fa(ib+ib)/fa(ib)/fa(L-ib)/fa(ib+ib-L)

                 do IC = MM,M
                   XCF = si(IC)*fa(IC+IC)/fa(IC)/fa(M-IC)/fa(IC+IC-M)
                   XI = XI+XF*AF*BF*XCF/dble(ia*2+ib*2+IC*2-N-L-M+1)
                 end do

               end do

             end do

             pin(L,M,N) = XI

           end do

         end do

       end do

       return

      end subroutine

!-----------------------------------------------------------------------
      subroutine hfdisk(iu,ir,etot,nst,rel,nr,rmin,rmax,r,rho,            &
     &                    zorig,xntot,ixflag,nel,                         &
     &                    no,nl,xnj,is,ev,ek,occ,njrc,vi,phe,orb)

       implicit double precision (a-h,o-z)
       integer iorbs, iside, io2, ijive
       integer lmax, ihmax, nrmax, ntmax, npmax
       integer iu, ir, nst, nr, nel, ixflag
       integer i, iu1, iprint
       double precision rmin, rmax, xntot, etot, rel, zorig
       double precision pi, pi4, rden
       parameter (iorbs=33,iside=600)
       parameter (io2=iorbs*(iorbs+1)/2)
       parameter (ijive=io2*(io2+1)/2)
       parameter (lmax=4,ihmax=20,nrmax=4000,ntmax=10,npmax=60)
       integer, dimension(4)     :: njrc 
       integer, dimension(iorbs) :: no, is, nl
       double precision, dimension(iorbs)   :: xnj, ev, ek, occ
       double precision, dimension(nrmax)   :: rho, r
       double precision, dimension(nrmax,7) :: vi
       double precision, dimension(nrmax,iorbs) :: phe, orb
       character(len=255) filename

       pi = 4.*atan(1.)
       pi4 = 4.*pi
       rden = 3.0

       if (iu < 0) then
         ! READ FULL FILENAME
         read (5,52) filename
 52      format (a255)
         iu1 = 1
         open (unit=iu1,status='unknown',file=trim(filename))
       else
         iu1 = iu
         open (unit=iu1,status='unknown')
       endif

       ! define the logarithmic grid
       do i = 1,nr
          r(i) = rmin*(rmax/rmin)**(dble(i)/dble(nr))
       end do

       ! obtain the charge density on the logarithmic grid

       do i = 1,nr
         rho(i) = .0
         do ii = 1,nel
             rho(i) = rho(i)+occ(ii)*phe(i,ii)**2
         end do
       end do

       ! write output file

       iprint = 0
       write(iu1,550)iprint
550     format('RELA'/'RELAT. ATOMIC CHARGE DENSITY'/I1)
       write(iu1,54) rmin,rmax,nr,zorig
54      format (d15.8,d15.8,i5,f5.2)
        nden = nr*(log(rden/rmin)/log(rmax/rmin))
        write(iu1,56) (rho(j),j=1,nr)
56      format (f15.11)
       close (unit=iu1)

       return

      end subroutine

!-----------------------------------------------------------------------
!  exchange correlation routine, via Ceperley-Alder, as parametrized by
!  Perdew and Zunger, Phys. Rev. B 23, 5048.  we use their interpolation
!  between the unpolarized and polarized gas for the correlation part.
      subroutine exchcorr(nst,rel,rr,rh1,rh2,ex,ec,ux1,ux2,uc1,uc2)

       implicit double precision(a-h,o-z)
       integer nst
       double precision rel, rr, rh1, rh2, ex, ec, ux1, ux2, uc1, uc2
       double precision trd, ft, rh, pi, fp, xn1, xn2, xn, xl, rs
       double precision zeta, exchfactor, fe1, fu1, beta, eta, b2

       trd = 1.d0/3.d0
       ft = 4.d0/3.d0

       rh = rh1+rh2

       ! if one spin type, average polarization

       if (nst == 1) then
         rh1 = rh/2.d0
         rh2 = rh/2.d0
       endif

       ! get the n's, and the rs.

       pi = 3.14159265358979d0
       fp = 4.d0*pi
       xn1 = rh1/(rr*fp)
       xn2 = rh2/(rr*fp)
       xn = xn1+xn2

       ! effect cutoff, to avoid overflow

       if ((nst == 3) .or. (xn < 0.00000001d0)) then

         ex = 0.d0
         ec = 0.d0
         ux1 = 0.d0
         ux2 = 0.d0
         uc1 = 0.d0      
         uc2 = 0.d0

       else

         rs = (3.d0/(fp*xn))**trd
         zeta = (xn1-xn2)/xn

         exchfactor = -0.930525546d0

         if (xn1 == 0.d0) then
           fe1 = 1.d0
           fu1 = 1.d0
           ex1 = 0.d0
           ux1 = 0.d0
         else
           beta = 0.028433756d0*xn1**trd
           b2 = beta*beta
           eta = sqrt(1.d0+b2)
           xl = log(beta+eta)
           fe1 = 1.d0-1.5d0*((beta*eta-xl)/b2)**2.d0
           fu1 = -0.5d0+1.5d0*xl/beta/eta
           ex1 = exchfactor*xn1**trd
           ux1 = 4.d0*ex1/3.d0
         endif

         if (xn2 == 0.d0) then
           fe2 = 1.d0
           fu2 = 1.d0
           ex2 = 0.d0
           ux2 = 0.d0
         else
           beta = 0.028433756d0 * xn2**trd
           b2 = beta * beta
           eta = sqrt(1.d0 + b2)
           xl = log(beta + eta)
           fe2 = 1.d0 - 1.5d0 * ((beta * eta - xl) / b2)**2.d0
           fu2 = -0.5d0 + 1.5d0 * xl / beta / eta
           ex2 = exchfactor * xn2**trd
           ux2 = 4.d0 * ex2 / 3.d0
         endif

         ! these next lines do the Ceperley-Alder correlation

         if (rs >= 1.d0) then

           rootr = sqrt(rs)

           gamma = -0.1423d0
           beta1 = 1.0529d0
           beta2 = 0.3334d0
           denom = (1.d0 + (beta1 * rootr) + (beta2 * rs))
           ecu = gamma / denom
           ucu = ecu * (1.d0 + (7.d0 / 6.d0 * beta1 * rootr) + (ft * beta2 * rs)) / denom

           gamma = -0.0843d0
           beta1 = 1.3981d0
           beta2 = 0.2611d0
           denom = (1.d0 + (beta1 * rootr) + (beta2 * rs))
           ecp = gamma / denom
           ucp = ecp * (1.d0 + (7.d0 / 6.d0 * beta1 * rootr) + (ft * beta2 * rs)) / denom

         else

           xlr = log(rs)
           rlr = rs * xlr

           au = 0.0311d0
           bu = -0.048d0
           cu = 0.002d0
           du = -0.0116d0
           ecu = (au * xlr) + bu + (cu * rlr) + (du * rs)
           ucu = au * xlr + (bu - au / 3.d0) + 2.d0 / 3.d0 * cu * rlr + (2.d0 * du - cu) * rs / 3.d0

           ap = 0.01555d0
           bp = -0.0269d0
           cp =  0.0007d0
           dp = -0.0048d0
           ecp = ap*xlr+bp+cp*rlr+dp*rs
           ucp = ap*xlr+(bp-ap/3.d0)+2.d0/3.d0*cp*rlr+(2.d0*dp-cp)*rs/3.d0

         endif

         ! if we are nonrelativistic, turn off MacDonald-Vosko correction
         if (rel == 0.d0) then
           fe1 = 1.d0
           fu1 = 1.d0
           fe2 = 1.d0
           fu2 = 1.d0
         endif

         ! interpolate the correlation energies
         denom = 2.d0**ft - 2.d0
         f = ((1.d0 + zeta)**ft + (1.d0 - zeta)**ft - 2.d0) / denom
         dfdz = ft/denom*((1.d0+zeta)**trd-(1.d0-zeta)**trd)
         ec = ecu + (f * (ecp - ecu))
         uc1 = ucu + f * (ucp - ucu) + (ecp - ecu) * (1.d0 - zeta) * dfdz
         uc2 = ucu + f * (ucp - ucu) + (ecp - ecu) * (-1.d0 - zeta) * dfdz

         ! get the final functional and potential
         ex = (xn1*fe1*ex1+xn2*fe2*ex2)/xn
         ux1 = fu1*ux1
         ux2 = fu2*ux2
         uc1 = uc1
         uc2 = uc2

       endif

       return

      end subroutine

!---------------------------------------------------------------------
      function cavpot(mtz_string, slab_flag, atomic_file,                 &
     &   cluster_file, mufftin_file, output_file, info_file)

       use WK
       use WF
       character(len=*), intent(in) :: atomic_file, mufftin_file
       character(len=*), intent(in) :: cluster_file, output_file
       character(len=*), intent(in) :: mtz_string, info_file
       integer, intent(in)          :: slab_flag
       real                         :: cavpot
       integer nieq, NTOT
       real rx(550), pot(550), RC(3,3)
       real title(20), sum
       real, allocatable    :: rk(:,:), zm(:), trk(:,:), Tzm(:)
       real, allocatable    :: sig(:,:), rho(:,:), vh(:,:)
       real, allocatable    :: vs(:,:), vmad(:,:), AD(:,:)
       real, allocatable    :: Z(:), ZC(:), rmt(:)
       integer, allocatable :: jrmt(:), jrmt2(:), nrr(:)
       integer, allocatable :: ncon(:), nx(:), NA(:,:), IA(:,:)
       real, allocatable    :: name(:,:)
       character(len=4) wfn,wfn0,wfn1,wfn2,wfn3
       data ngrid,MC,pi/250,30,3.1415926536/
       data wfn0,wfn1,wfn2,wfn3/"RELA","HERM","CLEM","POTE"/
       integer pos1, pos2

       index(X)=20.0*(log(X)+8.8)+2.0

       ! First input channels
       open(unit=4, file=atomic_file, status='old')
       open(unit=7, file=cluster_file, status='old')

       ! Now output channels
       open(unit=11, file=info_file, status='unknown')
       open(unit=9, file=mufftin_file, status='unknown')
       open(unit=10, file=output_file, status='unknown')

       ! INITIALISATION OF LOUCKS' EXPONENTIAL MESH
       x = -8.8

       do ix = 1,ngrid
         rx(ix) = exp(X)
         x = x + 0.05
       end do

       read(7,100) title
       write(11,200) title

       ! INPUT OF CRYSTALLOGRAPHIC DATA:
       !     SPA = LATTICE CONSTANT IN A.U.
       !     RC(I,J) = I'TH COORDINATE OF THE J'TH AXIS OF UNIT CELL,
       !               IN UNITS OF SPA
       !     RK(I,J) = I'TH COORDINATE OF THE J'TH ATOM IN UNIT CELL,
       !               IN UNITS OF SPA
       !     NIEQ = NUMBER OF INEQUIVALENT ATOMS IN UNIT CELL
       !
       !     FOR AN ATOM OF TYPE IR:
       !       NRR(IR) = NUMBER IN UNIT CELL
       !       Z(IR) = ATOMIC NUMBER
       !       ZC(IR) = VALENCE CHARGE
       !       RMT(IR) = MUFFIN-TIN RADIUS IN BOHR RADII
       read(7,101) spa
       read(7,101)((RC(I,J),I=1,3),J=1,3)

       do I=1,3
         do J=1,3
           RC(I,J) = spa*RC(I,J)
         end do
       end do

       ! read number of inequivalent atoms in unit cell
       read(7,102) nieq

       ! allocate arrays for nieq atoms
       if(.not. allocated(Z)) allocate(Z(nieq))
       if(.not. allocated(ZC)) allocate(ZC(nieq))
       if(.not. allocated(vh)) allocate(vh(250,nieq))
       if(.not. allocated(sig)) allocate(sig(250,nieq))
       if(.not. allocated(rho)) allocate(rho(250,nieq))
       if(.not. allocated(vs)) allocate(vs(250,nieq))
       if(.not. allocated(vmad)) allocate(vmad(550,nieq))
       if(.not. allocated(AD)) allocate(AD(30,nieq))
       if(.not. allocated(rmt)) allocate(rmt(nieq))
       if(.not. allocated(jrmt)) allocate(jrmt(nieq))
       if(.not. allocated(jrmt2)) allocate(jrmt2(nieq))
       if(.not. allocated(nrr)) allocate(nrr(nieq))
       if(.not. allocated(ncon)) allocate(ncon(nieq))
       if(.not. allocated(nx)) allocate(nx(nieq))
       if(.not. allocated(NA)) allocate(NA(30,nieq))
       if(.not. allocated(IA)) allocate(IA(30,nieq))
       if(.not. allocated(name)) allocate(name(8,nieq))

       ! initialise values
       do ir = 1,nieq
         do i = 1,ngrid
           vh(i,ir) = 0.0
           vs(i,ir) = 0.0
           vmad(i,ir) = 0.0
           sig(i,ir) = 0.0
           rho(i,ir) = 0.0
         end do
       end do

       vhar = 0.0
       vex = 0.0

       JJ = 0
       zz = 0.0
       inrr = 0

       ! loop through each INEQ atom
       do ir=1,nieq
         ! read name and split character string
         read(7,"(10A4)")(name(i,ir),i=1,4)

         ! read parameters NUM_ATOMS, Z, VALENCE, RYDBERG_MUFFTIN_RADIUS
         read(7,103) nrr(ir), Z(ir), ZC(ir), rmt(ir)

         ! reallocate enough space for extra atoms
         if(.not. allocated(rk)) then
           allocate(rk(3,nrr(ir)))
         else
           allocate(trk(3,inrr + nrr(ir)))
           trk(1:3,1:inrr) = rk
           deallocate(rk)
           call move_alloc(trk,rk)
         endif
         if(.not. allocated(zm)) then
           allocate(zm(nrr(ir)))
         else
           allocate(TZM(inrr + nrr(ir)))
           TZM(1:size(zm)) = zm
           deallocate(zm)
           call move_alloc(TZM,zm)
         endif
         inrr = inrr + nrr(ir)

         zz = zz + abs(ZC(ir))
         jrmt(ir) = index(rmt(ir))
         N =nrr(ir)

         do J=1,N
           JJ = JJ + 1
           zm(JJ) = ZC(ir)
           read(7,101)(rk(i,JJ),i=1,3)

           do i=1,3
             rk(i,JJ) = spa * rk(i,JJ)
           end do

         end do

       end do

       ! N = TOTAL NUMBER OF ATOMS IN UNIT CELL
       ! AV = TOTAL VOLUME OF UNIT CELL
       ! OMA = ATOMIC VOLUME
       ! RWS = WIGNER-SEITZ RADIUS
       N = JJ
       RCC1 = RC(2,2) * RC(3,3) - RC(3,2) * RC(2,3)
       RCC2 = RC(3,2) * RC(1,3) - RC(1,2) * RC(3,3)
       RCC3 = RC(1,2) * RC(2,3) - RC(2,2) * RC(1,3)
       AV = abs((RC(1,1) * RCC1) + (RC(2,1) * RCC2) + (RC(3,1) * RCC3))
       OMA = AV / float(N)
       rws =(0.75 * OMA / pi)**(1.0 / 3.0)
       jrws = index(rws)
       write(11,201)((RC(i,J),i=1,3),J=1,3)
       write(11,202)AV, OMA, rws
       JJ=0

       do ir=1,nieq
         write(11,203)ir, (name(i,ir),i=1,4),nrr(ir)
         inr = nrr(ir)

         do iir=1,inr
           JJ=JJ+1
           write(11,204)(rk(i,JJ),i=1,3)
         end do

         write(11,205)Z(ir),ZC(ir),rmt(ir)
       end do

       write(11,216)(rx(ix),ix=1,ngrid)

       ! FOR EACH ATOMIC TYPE, READ IN ATOMIC WAVEfUNCTIONS FOR NEUTRAL
       ! ATOM, IN EITHER HERMAN-SKILLMAN OR CLEMENTI FORM, PRODUCING?
       ! RHO = 4*PI*CHARGE DENSITY * RADIUS**2
       mix=0

       do ir=1,nieq
         read(4,100) wfn
         ! OPTION 0)  RELATIVISTIC CHARGE DENSITY INPUT
         if (wfn==wfn0) call RELA(rho(1,ir),rx,nx(ir),ngrid)
         ! OPTION 1)  HERMAN-SKILLMAN INPUT
         if (wfn==wfn1) call HSIN(rho(1,ir),rx,nx(ir),ngrid)
         ! OPTION 2)  CLEMENTI INPUT
         if (wfn==wfn2) call CLEMIN(rho(1,ir),rx,nx(ir),ngrid)
         ! OPTION 3)  POTENTIAL INPUT
         if (wfn==wfn3) goto 14

         ! RHO IS NORMALISED USING TOTAL ELECTRONIC CHARGE ON THE ATOM
         ! CALCULATED BY THE TRAPEZOIDAL RULE
         nix = nx(ir)
         mix = max(nix,mix)
         sum = 0.0D0
         W1 = 0.025 * rho(1,ir) * rx(1)

         do ix=2,nix
           W2 = 0.025 * rho(ix,ir) * rx(ix)
           sum = sum + W1 + W2
           W1 = W2
         end do

         ZE=sum

         ! SOLVE POISSON'S EQUATION:
         ! SIG = COULOMB POTENTIAL
         ! RHO = 4*PI*CHARGE DENSITY*RADIUS SQUARED
         call poisson(rho(1,ir), Z(ir), nix, sig(1,ir))

         X = -8.8

         do ix=1,nix
           CE = exp(-0.5*X)
           sig(ix,ir) = CE  *((-2.0 * Z(ir) * CE) + sig(ix,ir))
           rho(ix,ir) = rho(ix,ir) / (rx(ix)**2)
           X = X + 0.05
         end do

         write(11,206)(name(I,ir), I=1,4), ZE, rx(nix), nix
         write(11,207)(sig(ix,ir),ix=1,nix)

       end do

       ! DETAILS OF NEIGHBOURING SHELLS FOR EACH ATOMIC TYPE IR?
       ! NCON(IR) = NUMBER OF SHELLS INCLUDED
       ! IA(J,IR) = ATOMIC TYPE IN J'TH SHELL
       ! NA(J,IR) = NUMBER OF ATOMS IN J'TH SHELL
       ! AD(J,IR) = DISTANCE TO J'TH SHELL
       rmax = rx(mix)

       call nbr(IA,NA,ad,ncon,nrr,nieq,RC,rk,N,rmax,MC)
       write(11,208)

       do ir=1,nieq
         write(11,209)ir
         NC = ncon(ir)
         ic = (NC - 1) / 12+1
         kc = 0

         do I=1,ic
           jc = kc + 1
           kc = min0(NC, kc + 12)
           write(11,210)(ad(J,ir),J=jc,kc)
           write(11,211)(NA(J,ir),J=jc,kc)
           write(11,212)(IA(J,ir),J=jc,kc)
         end do

       end do

       read(7,102) nform

       ! CALCULATION OF THE MUFFIN-TIN POTENTIAL FOR EACH NEUTRAL
       ! ATOM, FOLLOWING THE MATTHEISS PRESCRIPTION
       ! READ IN ALPHA FOR THE SLATER EXCHANGE TERM
       read(7,101) alpha
       write(11,215) alpha
       pd = 6.0 / (pi*pi)

       do ir=1,nieq
         jrx = max(jrws,jrmt(ir))
         ! SUMMING THE POTENTIALS FROM NEUTRAL ATOMS
         ! VH = HARTREE POTENTIAL
         call sumax(vh(1,ir),sig,rx,nx,ncon(ir),IA(1,ir),NA(1,ir),        &
     &              ad(1,ir),jrx,ngrid,nieq)

         ! SUMMING THE CHARGE DENSITY ABOUT EACH ATOMIC TYPE
         ! VS = TOTAL CHARGE DENSITY, THEN SLATER EXCHANGE TERM
         call sumax(vs(1,ir),rho,rx,nx,ncon(ir),IA(1,ir),NA(1,ir),        &
     &              ad(1,ir),jrx,ngrid,nieq)

         do ix=1,jrx
           vs(ix,ir) = -1.5 * alpha * (pd * vs(ix,ir))**(1.0 / 3.0)
         end do

       end do

       ! CALCULATE THE MUFFIN-TIN ZERO
       vint = 0.
       read(7,102) nh

       if (nh == 0 .and. nieq == 1) then
         call mtzm(vh(1,1),vs(1,1),rx,ngrid,rmt(1),rws,jrmt(1),jrws,vhar,vex)
       endif

       if (nh.ne.0) then
         call mtz(sig,rho,rx,ngrid,rmt,nrr,nx,nieq,RC,rk,N,vhar,vex,alpha,AV,nh)
       endif

       ! Slab or Bulk calculation?
       ! slab_flag = 1 for Slab or 0 for Bulk (default = 0)
       if (slab_flag == 1) then
         ! Input the MTZ value from the bulk substrate calculation
         if (len_trim(mtz_string) >= 1) then
           read(mtz_string,*) esht
         else
           write(*,*) 'Error: mtz input is invalid for slab calculation'
           close(4)
           close(7)
           close(9)
           close(10)
           close(11)
           return
         endif

         esh=esht-(vhar+vex)

       else
         ! check slab_flag is valid
         if (slab_flag.ne.0) then
           write(*,*) 'Warning: SLAG_FLAG defaulted to 0'
         endif
         ! If you are interested in adatoms on this substrate
         ! rerun a slab calculation with the adatoms
         ! and use this MTZ value as input when asked
         write(10,*) (vhar + vex)
         cavpot = vhar + vex  ! return value of function
       endif

       goto 16

       ! OPTION 3)  READ IN POTENTIAL OF NEUTRAL ATOM, vh, ON RADIAL
       !            GRID, RX, FOR CORRECTION BY MADELUNG SUMMATION
14     read(4,104) ngrid, (rx(ix),ix=1,ngrid)

       do ir=1,nieq
         read(4,104) jrx, (vh(ix,ir),ix=1,jrx)
         jrmt(ir) = jrx
       end do

       ! THE MADELUNG CORRECTION FOR IONIC MATERIALS. SUBROUTINE MAD
       ! COMPUTES THE SPHERICALLY AND SPATIALLY AVERAGED FIELDS FOR
       ! THE LATTICE OF POINT CHARGES ABOUT EACH ATOMIC TYPE
16     if (zz .ne. 0) call mad(vmad,rx,ngrid,rmt,nrr,jrmt,nieq,RC,rk,zm,N,AV)

       ! THE TOTAL MUFFIN-TIN POTENTIAL IS ACCUMULATED INTO SIG,
       ! REFERRED TO THE MUFFIN-TIN ZERO
       vint = vhar + vex
       if (nform == 0) write(9,102) nieq

       do ir=1,nieq
         write(11,213) (name(I,ir),I=1,4), vint, rmt(ir)
         jrx = jrmt(ir)

         do ix=1,jrx
           vh(ix,ir) = vh(ix,ir) - vhar
           vs(ix,ir) = vs(ix,ir) - vex
           sig(ix,ir) = vh(ix,ir) + vs(ix,ir) + vmad(ix,ir)
           write(11,214) rx(ix), vh(ix,ir),vs(ix,ir),vmad(ix,ir),sig(ix,ir)
         end do

       end do

       ! write output in a format to be read by WILLIAMS phase shift
       ! program (NFORM=1), by CAVLEED phase shift program (NFORM=0), or
       ! by the relativistic phase shift program (NFORM=2)

       ! Also prepare to shift the potential by an amount of the order
       ! of the bulk muffintin zero.
       ! This is needed only if the cluster.i file corresponds to a
       ! surface adsorbate
       if (nform == 1) write(9,220) nieq

       if (nform == 2) then
         ! define German grid RX and save old grid in RS
         rm = 60.0
         dx = 0.03125
         nmx = 421
         RS(1) = rx(1)
         rx(1)= rm * exp(dx * (1 - nmx))
         J=1
         rm=exp(dx)

         do while (J < nmx)
           K=J+1
           RS(K)=rx(K)
           rx(K)=rm*rx(J)
           J=K
         end do

       endif

       do ir=1,nieq
         jrx = jrmt(ir)
         if (nform == 0) then
           write(9,217) (name(I,ir),I=1,4)
           write(9,218) Z(ir), rmt(ir), vint
         elseif (nform == 1) then
           write(9,221) Z(ir), rmt(ir)
         else
           ! es=Emin for phase shift calculation (ev)
           ! de=delta E for phase shift calculation (ev)
           ! ue=Emax for phase shift calculation (ev)
           ! lsm=maximum number of phase shifts desired
           es = 20.
           de = 5.
           ue = 300.
           lsm = 12
           write(9,217)(name(I,ir),I=1,4)
           write(9,"(3D12.4,4X,I3,4X,D12.4)") ES,DE,ue,lsm,vint

           ! INTERPOLATION TO GRID RX
           do k=1,jrx
             sig(k,ir) = (sig(k,ir) - esh) * rs(k)
           end do

           nmxx = nmx
           call chgrid(sig(1,ir), RS, jrx, pot, rx, nmxx)
           iz = Z(ir)
           write(9,"(I4,F10.6,I4)") iz, rmt(ir), nmxx
           jrx = nmxx
         endif

         if (nform == 0) write(9,102) jrx
         if (nform == 1) then
           do ix=1,jrx
             write(9,219) rx(ix), rx(ix) * (sig(ix,ir) - esh)
           end do
           rneg = -1.
           write(9,219) rneg
         elseif (nform == 0) then
           do ix=1,jrx
             write(9,219) rx(ix), (sig(ix,ir) - esh)
           end do
         else
           write(9,"(5E14.7)") (pot(ix),ix=1,jrx)
         endif

       end do

       ! close file handles
       close(4)
       close(7)
       close(9)
       close(10)
       close(11)

       return

100    format(20A4)
101    format(3F8.4)
102    format(I4)
103    format(I4,3F8.4)
104    format(I4/(5E14.5))
200    format("MUFFIN-TIN POTENTIAL PROGRAM",5X,20A4)
201    format(///" AXES OF UNIT CELL"/(6X,3F8.4))
202    format(" UNIT CELL VOLUME",F15.4/" ATOMIC VOLUME",F18.4/           &
     & " WIGNER-SEITZ RADIUS",F12.4)
203    format(///" TYPE",I2," ATOM",2X,4A4/I4," ATOMS IN UNIT CELL")
204    format(6X,3F8.4)
205    format(" ATOMIC NUMBER",F15.1/"VALENCE",F21.1/                     &
     & "MUFFIN-TIN RADIUS",F14.4)
206    format(///1X,4A4," ELECTRONIC CHARGE",F12.5/                       &
     & "COULOMB POTENTIAL FOR ISOLATED ATOM, OUT TO RADIUS",              &
     & F12.5,10X,"NX",I4/)
207    format(5(10E12.4/))
208    format("1")
209    format(//"NEAREST NEIGHBOUR SHELLS FOR TYPE",I2," ATOM")
210    format(" DISTANCE",1X,15(F8.4))
211    format(" NUMBER",3X,15(I5,3X))
212    format(" TYPE",5X,15(I5,3X))
213    format("1",4A4,5X,                                                 &
     & "POTENTIALS IN RYDBERGS CORRECT TO MUFFIN-TIN ZERO",F8.4,          &
     & "MUFFIN-TIN RADIUS",F8.4//5X,"RADIUS",5X,"HARTREE POTENTIAL",9X,   &
     & "EXCHANGE",4X,"MADELUNG CORRECTION",5X,"TOTAL POTENTIAL")
214    format(F12.5,4E20.6)
215    format(///"STATISTICAL EXCHANGE PARAMETER, ALPHA",F10.4)
216    format(///"LOUCKS' RADIAL MESH"//5(10F11.5/))
217    format(4A4)
218    format(3F8.4)
219    format(2E14.5)
220    format(" &nl2 nrr=",i2," &end")
221    format(" &nl16 Z=",f7.4,",RT=",f7.4," &end")

      end function

!---------------------------------------------------------------------
!     PIECEWISE QUADRATIC INTERPOLATION FROM GRID X TO GRID Y,  BY
!     AITKEN'S DIVIDED DIFFERENCE SCHEME.   NX,NY ARE ARRAY dimensions
!     NOTE THAT NY IS RESET IN CHGRID
      subroutine chgrid(fx,x,nx,fy,y,ny)
      
       integer, intent(in) :: nx
       integer, intent(inout) :: ny
       double precision, dimension(nx), intent(in)    :: fx, x
       double precision, dimension(ny), intent(in)    :: y
       double precision, dimension(ny), intent(inout) :: fy
       integer ix, iy
       double precision yy, a1, a2, a3, a12, a13
       
       iy=1

       do 2 ix=3,nx
         do while(.True.)
           if (iy > ny) goto 3
           yy=y(iy)
           if (yy > x(ix)) goto 2
           a1=x(ix-2)-yy
           a2=x(ix-1)-yy
           a3=x(ix)-yy
           a12=(fx(ix-2)*a2-fx(ix-1)*a1)/(x(ix-1)-x(ix-2))
           a13=(fx(ix-2)*a3-fx(ix)*a1)/(x(ix)-X(ix-2))
           fy(iy)=(a12*a3-a13*a2)/(x(ix)-x(ix-1))
           if (iy > ny) goto 3
           iy=iy+1
         end do
2      continue

3      ny=iy-1

       return

      end subroutine

!---------------------------------------------------------------------
! ROUTINE FOR INPUT OF WAVEFUNCTIONS IN THE CLEMENTI PARAMETRISED
! FORM, AND CALCULATION OF CHARGE DENSITY ON RADIAL MESH RX
!     RHO = 4*PI*SUM OVER STATES OF (MODULUS(WAVE FN)**2) *
!          RADIUS**2
!     NC = NUMBER OF ATOMIC STATES
! FOR EACH ATOMIC STATE I?
!     LC(I) = ANGULAR MOMENTUM
!     FRAC = FRACTIONAL OCCUPATION
!     WC(I) = NUMBER OF ELECTRONS
!     WFC(IX,I) = WAVEFUNCTION X RADIUS AT GRID POINT IX
      subroutine clemin(rho,rx,nx,ngrid)

       use WF
       use WK
       integer, intent(in) :: ngrid
       integer, intent(inout) :: nx
       double precision, dimension(ngrid) :: rho, rx
       real, dimension(4) :: name
       real sum
       integer ic, ig, ns, j, i, k, kd
       double precision a, b, c

       read(4,100) name
       read(4,101) iprint
       read(4,101) nc
 
       do ic=1,nc
         do ig=1,ngrid
           wfc(ig,ic) = 0.0
         end do
       end do

       ! INPUT OF CLEMENTI PARAMETERS
       ic = 1
2      read(4,101) ns
       if (ns <= 0) goto 8
       do i=1,ns
         read(4,102) nt(i),ex(i)
       end do

       do J=1,ns
         a = 1.0
         b = 2.0
         k = nt(j)
         c = float(k)
         kd = k + k

         do I=2,kd
           a = a * b
           b = b + 1.0
         end do

         fnt(j) = exp(-0.5 * log(a) + (c + 0.5) * log(2.0 * ex(j)))
       end do

5      read(4,101) lc(ic)
       if (lc(ic) < 0) goto 2
       read(4,103)(fac(j),j=1,ns)
       read(4,103) frac
       WC(ic) = 2.0 * float(2 * lc(ic) + 1) * frac

       do ix=1,ngrid
         sum = 0.0d0

         do k=1,ns
           exx = ex(k) * rx(ix)
           if (exx > 80.0) exit  ! loop
           sum = sum + fac(k)*fnt(k)*(rx(ix)**(nt(k)))*exp(-exx)
         end do

         wfc(ix,ic) = sum
       end do

       ic=ic+1
       goto 5

       ! CALCULATION OF CHARGE DENSITY
8      do ix=1,ngrid
         sum = 0.0d0
         do ic=1,nc
           sum = sum + (WC(ic) * wfc(ix,ic) * wfc(ix,ic))
         end do
         rho(ix) = sum

         if (sum < 1.0d-9) exit  ! loop
       end do

       nx=ix

       if (iprint == 0) return

       write(11,200) name

       do ic=1,nc
        write(11,201) lc(ic), (wfc(ix,ic),ix=1,ngrid)
       end do

       write(11,202) rx(nx), nx, (rho(ix),ix=1,nx)

       return

100    format(4A4)
101    format(I4)
102    format(I4,F11.5)
103    format(5F11.5)
200    format("1",4A4," ATOMIC WAVEFUNCTIONS (CLEMENTI)", " X RADIUS")
201    format("L",I3//5(10F11.5/))
202    format("CHARGE DENSITY OUT TO RADIUS",F12.5,10X,                   &
     & "NX",I4//5(10E12.4/))

       return

      end subroutine

!---------------------------------------------------------------------
!  ROUTINE FOR INPUT OF ATOMIC WAVEfunctionS FROM HERMAN-SKILLMAN
!  TABLES, AND CALCULATION OF CHARGE DENSITY ON THE RADIAL MESH RX
!    RHO = 4*PI*SUM OVER STATES OF (MODULUS(WAVE FN)**2) *
!          RADIUS**2
!    NM ? H-S GRID INTERVAL doubleS EVERY NM MESH POINTS
!    NC = NUMBER OF ATOMIC STATES
! FOR EACH ATOMIC STATE I?
!    LC(I) = ANGULAR MOMENTUM
!    FRAC = FRACTIONAL OCCUPATION
!    WC(I) = NUMBER OF ELECTRONS
!    WFC(IX,I) = WAVEfunction X RADIUS AT GRID POINT IX
      subroutine hsin(rho,rx,nx,ngrid)

       use WF
       use WK
       integer, intent(in) :: ngrid
       integer, intent(inout) :: nx
       real, dimension(ngrid) :: rho,rx
       real, dimension(4) :: name
       real sum
       integer ig, ic, iprint


       read(4,100) name, z
       read(4,101) iprint
       read(4,101) nm
       read(4,101) nc

       do ig=1,250
         rs(ig) = 0.0

         do ic=1,nc
           wfc(ig,ic) = 0.0
         end do

       end do

       ! INITIALISATION OF HERMAN-SKILLMAN MESH
       dr = 0.005 * 0.88534138 / exp(log(z) / 3.0)
       rr(1) = 0.0

       do i=2,250
         if (mod(i, nm) == 2) dr=dr+dr
         rr(i) = rr(i-1) + dr
       end do

       ns = 0

       do ic=1,nc
         read(4,101) lc(ic), n, frac
         ns = max0(ns,n)
         WC(ic) = 2.0 * float(2 * lc(ic) + 1) * frac
         
         do ix=1,n
           read(4,102) wfc(ix,ic)
         end do
         
       end do

       ! CALCULATION OF CHARGE DENSITY
       do ix=1,ns
         sum = 0.0d0

         do ic=1,nc
           sum = sum + WC(ic) * wfc(ix,ic) * wfc(ix,ic)
         end do

         rs(ix) = sum
       end do

       ! INTERPOLATION TO GRID RX
       nx = ngrid
       call chgrid(rs, rr, ns, rho, rx, nx)
       if (iprint == 0) return
       write(11,200) name,(RR(ix),ix=1,ns)

       do ic=1,nc
         write(11,201) lc(ic),(wfc(ix,ic),ix=1,ns)
       end do

       do ix=1,nx
         if (rho(ix) < 1.0E-9) exit  ! loop
       end do

       nx = ix
       write(11,202) rx(nx), nx, (rho(ix),ix=1,nx)

       return

100    format(4A4/F9.4)
101    format(2I4,F9.4)
102    format(5F9.4)
200    format("1",4A4," ATOMIC WAVEFUNCTIONS (HERMAN-SKILLMAN)",          &
     & " X RADIUS"//" HERMAN-SKILLMAN MESH"//5(10F12.5/))
201    format("L",I3//5(10F11.5/))
202    format("CHARGE DENSITY OUT TO RADIUS",F12.5,10X,                   &
     & "NX",I4//5(10E12.4/))

      end subroutine

!---------------------------------------------------------------------
!  subroutine MAD CALCULATES THE SPHERIcallY AND SPATIALLY AVERAGED
!  FIELDS FROM A LATTICE OF POINT CHARGES, AND TABULATES THEM ON
!  A RADIAL MESH RX, ABOUT EACH ATOMIC TYPE IN THE UNIT CELL
!  ** NB THIS ROUTINE WORKS IN HARTREES, BUT CONVERTS TO RYDBERGS **
!  RC(I,J) = THE I'TH COORDINATE OF THE J'TH AXIS OF THE UNIT CELL
!  RK(I,J) = THE I'TH COORDINATE OF THE J'TH ATOM IN THE UNIT CELL
!  VMAD(J,IR) = THE J'TH TABULATED VALUE OF THE SPHERIcallY AVERAGED
!  POTENTIAL ABOUT A TYPE-IR ATOM
!  ZM(K)=CHARGE ON THE K'TH ATOM
!  RMT(IR) = MUFFIN-TIN RADIUS OF A TYPE-IR ATOM
!  NR = NUMBER OF INEQUIVALENT ATOMS IN THE CELL
!  AV = VOLUME OF UNIT CELL
!  G(I,J) = I'TH COORDINATE OF THE J'TH RECIPROCAL LATTICE VECTOR
!  VMM(IR) = THE INTEGRAL OF THE POTENTIAL ABOUT A TYPE-IR ATOM
!  OUT TO THE MUFFIN-TIN RADIUS
      subroutine mad(vmad,rx,ngrid,rmt,nrr,nx,nr,rc,rk,zm,n,av)
    
       use WK
       integer, intent(in) :: ngrid, nr, n 
       double precision av
       double precision, dimension(ngrid)     :: rx 
       double precision, dimension(ngrid, nr) :: vmad
       double precision, dimension(3,3)       :: rc
       double precision, dimension(3,n)       :: rk
       double precision, dimension(n)         :: zm
       double precision, dimension(nr)        :: rmt
       integer, dimension(nr) :: nrr,nx
       integer j, ir
       data pi,TEST/3.1415926536,1.0E-4/

       rad(a1,a2,a3) = sqrt((a1 * a1) + (a2 * a2) + (a3 * a3))


       do ir=1,nr
         fr(ir) = 0.0
         do j=1,ngrid
           vmad(j,ir) = 0.0
         end do
       end do

       ! THE RECIPROCAL LATTICE IS DEFINED BY THREE VECTORS, G
       atv = 2.0 * pi / av
       g(1,1) = (rc(2,1)*rc(3,2)-rc(3,1)*rc(2,2))*atv
       g(2,1) = (rc(3,1)*rc(1,2)-rc(1,1)*rc(3,2))*atv
       g(3,1) = (rc(1,1)*rc(2,2)-rc(2,1)*rc(1,2))*atv

       g(1,2) = (rc(2,2)*rc(3,3)-rc(3,2)*rc(2,3))*atv
       g(2,2) = (rc(3,2)*rc(1,3)-rc(1,2)*rc(3,3))*atv
       g(3,2) = (rc(1,2)*rc(2,3)-rc(2,2)*rc(1,3))*atv

       g(1,3) = (rc(2,3)*rc(3,1)-rc(3,3)*rc(2,1))*atv
       g(2,3) = (rc(3,3)*rc(1,1)-rc(1,3)*rc(3,1))*atv
       g(3,3) = (rc(1,3)*rc(2,1)-rc(2,3)*rc(1,1))*atv

       ! MAXIMUM VALUE OF RK, AND MINIMUM VALUES OF RC,G - PRIOR TO
       ! CHOOSING THE SEPARATION CONSTANT AL AND LIMITS FOR SUMMATIONS
       rkmax = 0.0

       do j=1,n
         rkmax = max(rkmax, rad(rk(1,j), rk(2,j), rk(3,j)))
       end do

       rcmin = 1.0E6
       gmin = 1.0E6

       do j=1,3
         rcmin = min(rcmin, rad(rc(1,j), rc(2,j), rc(3,j)))
         gmin = min(gmin, rad(g(1,j), g(2,j), g(3,j)))
       end do

       ! AL IS CHOSEN TO GIVE EQUAL NUMBERS OF TERMS IN real AND
       ! RECIPROCAL SPACE SUMMATIONS
       fac1 = TEST * log(TEST) ** 4
       fac2 = (4.0 * pi * rcmin ** 4) / (av * gmin ** 4)
       al = exp(log(fac1 / fac2) / 6.0)
       itr = 1 + ifix(((al * rkmax) - log(TEST)) / (al * rcmin))
       limr = itr + itr + 1
       fac1 = 4.0 * pi * al * al / (av * gmin ** 4)
       itg = 1 + ifix(exp(log(fac1 / TEST) / 4.0))
       limg = itg + itg + 1
       
       do i=1,3
         do j=1,3
           write(11,200) g(i,j)
         end do
       end do
       
       write(11,201) rcmin, gmin, rkmax, TEST, al

       ! REAL SPACE SUMMATION
       write(11,202) itr
       ! THE PREFACTORS FR FROM THE real SPACE SUMMATION ARE CALCULATED
       as = -float(itr) - 1.0
       ax = as

       do 5 jx=1,limr
         ax = ax + 1.0
         ay = as

         do 5 jy=1,limr
           ay = ay + 1.0
           az = as

           do 5 jz=1,limr
             az = az + 1.0

             do i=1,3
               ra(i) = (ax * rc(i,1)) + (ay * rc(i,2)) + (az * rc(i,3))
             end do

             do 5 j=1,n
               k = 1

               do 5 kr=1,nr
                 r = rad(ra(1) + rk(1,j) - rk(1,k),                       &
     &                   ra(2) + rk(2,j) - rk(2,k),                       &
     &                   ra(3) + rk(3,j) - rk(3,k))
                 if (r < 1.0E-4) goto 5
                 fr(kr) = fr(kr) + zm(j) * exp(-al*r) / r
5                k = k + nrr(kr)

       k = 1

       do kr=1,nr
         x = rmt(kr)
         a = exp(-al*x)
         ai1 = ((1.0-a)/al-X*a)/al
         ai2 = (x * 0.5 * (1.0 / a + a)                                   &
    &           - 0.5 * (1.0 / a - a) / al) / al / al
         vmm(kr) = 4.0 * pi * (zm(k) * ai1 + ai2 * fr(kr))
         nix = nx(kr)

         do j=1,nix
           x = rx(j)
           a = exp(al * x)
           vmad(j,kr) = fr(kr) * 0.5 * (a - 1.0 /a)/(al*x) + zm(k)/(a*x)
         end do

         k = k + nrr(kr)

       end do

       write(11,203)(vmm(kr),kr=1,nr)

       ! NEXT COMES THE SUMMATION IN RECIPROCAL SPACE
       write(11,204) itg
       as = -float(itg) - 1.0
       ax = as

       do 13 jx=1,limg
         ax = ax + 1.0
         ay = as

         do 13 jy=1,limg
           ay = ay + 1.0
           az = as

           do 13 jz=1,limg
             az = az + 1.0

             do i=1,3
               ga(i) = (ax * g(i,1)) + (ay * g(i,2)) + (az * g(i,3))
             end do

             gm = rad(ga(1), ga(2), ga(3))
             gs = gm * gm
             fac1 = 0.0
             if (gs < 1.0E-4) goto 13
             fac1 = 4.0 * pi * al * al /(av * gs * (gs + (al * al)))
             k = 1

             do kr=1,nr
               fac2 = 0.0

               do j=1,n
                 gr = 0.0

                 do i=1,3
                   gr = gr + ga(i) * (rk(i,k) - rk(i,j))
                 end do

                 fac2 = fac2 + cos(gr) * zm(j)
               end do

               x = rmt(kr)
               ai3 = (sin(gm*x) / gm - x * cos(gm * x)) / gs
               vmm(kr) = vmm(kr) + 4.0 * pi * ai3 * fac1 * fac2
               nix = nx(kr)

               do i=1,nix
                 x = rx(i)
                 vmad(i,kr) = vmad(i,kr) + fac1 * fac2 * sin(gm * x) / (gm * x)
               end do

               k = k + nrr(kr)
             end do

13     continue

       write(11,203)(vmm(kr),kr=1,nr)

       ! REFER TO MUFFIN-TIN ZERO
       vm = 0.0
       amt = 0.0

       do ir=1,nr
         vm = vm + float(nrr(ir)) * rmt(ir) ** 3
         amt = amt + float(nrr(ir)) * vmm(ir)
       end do

       amt = amt / (av-4.0*pi*vm/3.0)

       ! EXPRESS THE FINAL POTENTIAL IN RYDBERGS
       amt = -2.0 * amt
       write(11,205) amt

       do kr=1,nr
         nix = nx(kr)
         do j=1,nix
           vmad(j,kr) = 2.0 * vmad(j,kr) - amt
         end do
       end do

       return

200    format(///"MADELUNG CORRECTION"//"RECIPROCAL LATTICE"/(6X,3F8.4))
201    format("RCMIN",F10.4,10X,"GMIN",F10.4,10X,"RKMAX",F10.4,           &
     & 10X,"TEST",E12.4/" SEPARATION CONSTANT",E12.4)
202    format("REAL SPACE SUMMATION",11X,"ITR",I3)
203    format(" VMM (HARTREES) ",5E12.4)
204    format("RECIPROCAL SPACE SUMMATION",5X,"ITG",I3)
205    format("MADELUNG MUFFIN-TIN ZERO",5E12.4)

      end subroutine

!---------------------------------------------------------------------
!  SUBROUTINE FOR CALCULATION OF THE MUFFIN-TIN ZERO LEVEL?
!  THE AVERAGE VALUE OF THE POTENTIAL BETWEEN THE MUFFIN-TIN
!  SPHERES IN THE UNIT CELL
      subroutine mtz(sig, rho, rx, ngrid, rmt, nrr, nx, nr,               &
     &               rc, rk, n, vhar, vex, alpha, av, nh)
       
       use WK
       integer, intent(inout) :: ngrid, nr, n, nh
       real, dimension(ngrid) :: rx
       real, dimension(ngrid,nr) :: sig, rho
       real, dimension(nr)  :: rmt, nrr, nx
       real, dimension(20)  :: vg
       real, dimension(3,3) :: rc
       real, dimension(3,n) :: rk
       real vhar, vex, alpha, av
       data pi,ng/3.14159265358,20/

       ! GRID REFERENCE FOR RADIUS ON LOUCKS' MESH
       index(y) = 20.0 * (log(y) + 8.8) + 1.0
       rad(A1,A2,A3) = sqrt((A1 * A1) + (A2 * A2) + (A3 * A3))

       pd = 6.0 / pi / pi

       do ig=1,ng
         vg(ig) = 0.0
       end do

       ig = 0
       vhar = 0.0
       vex = 0.0
       npoint = 0
       nint = 0
       dh = 1.0 / float(nh)

1      AH = dh / 2.0

       AX = -AH

       do 7 ix=1,nh
         AX = AX + dh
         AY = -AH

         do 7 IY=1,nh
           AY = AY + dh
           AZ = -AH

           do 7 iz=1,nh
             AZ = AZ + dh

             do I=1,3
               X(I) = (AX * RC(I,1)) + (AY * RC(I,2)) + (AZ * RC(I,3))
             end do

             npoint = npoint + 1
             ! GIVES SAMPLE POINT X INSIDE THE UNIT CELL - TEST WHETHER
             ! INTERSTITIAL
             bx = -1.0

             do 4 JX=1,2
               bx = bx + 1.0
               BY = -1.0

               do 4 JY=1,2
                 BY = BY + 1.0
                 bz = -1.0

                 do 4 JZ=1,2
                   bz = bz + 1.0

                   do I=1,3
                     RB(I) = X(I) - (bx * RC(I,1)) - (by * RC(I,2)) - (bz * RC(I,3))
                   end do

                   I = 0

                   do 4 ir=1,nr
                     inr = nrr(ir)

                     do 4 iir=1,inr
                       I = I + 1
                       XR = rad(RB(1)-rk(1,I), RB(2) - rk(2,I), RB(3)-rk(3,I))
                       if (XR < rmt(ir)) goto 7

4                      continue

             ! WE HAVE AN INTERSTITIAL POINT
             nint = nint + 1
             ! SUM COULOMB AND EXCHANGE ENERGIES FROM ATOMS WITHIN 2 UNIT
             ! CELLS AROUND THIS POINT
             sumc = 0.0
             sume = 0.0
             bx = -3.0

             do 6 JX=1,5
               bx = bx + 1.0
               by = -3.0

               do 6 JY=1,5
                 by = by + 1.0
                 bz = -3.0

                 do 6 JZ=1,5
                   bz = bz + 1.0

                   do I=1,3
                     RB(I) = bx*RC(I,1)+by*RC(I,2)+bz*RC(I,3)-X(I)
                   end do

                   J = 0

                   do 6 jr=1,nr
                     jnr = nrr(jr)

                     do 6 jjr=1,jnr
                       J = J + 1
                       XR=rad(RB(1)+rk(1,J),RB(2)+rk(2,J),RB(3)+rk(3,J))
                       J2 = index(XR)
                       if (J2 >= nx(jr)) goto 6
                       J1 = J2 - 1
                       J3 = J2 + 1
                       X1 = RX(J1)
                       X2 = RX(J2)
                       X3 = RX(J3)
                       TERMC = (XR-X2)*(XR-X3)/(X1-X2)/(X1-X3)*sig(J1,jr) &
     &                       + (XR-X1)*(XR-X3)/(X2-X1)/(X2-X3)*sig(J2,jr) &
     &                       + (XR-X2)*(XR-X1)/(X3-X2)/(X3-X1)*sig(J3,jr)
                       TERME = (XR-X2)*(XR-X3)/(X1-X2)/(X1-X3)*rho(J1,jr) &
     &                       + (XR-X1)*(XR-X3)/(X2-X1)/(X2-X3)*rho(J2,jr) &
     &                       + (XR-X2)*(XR-X1)/(X3-X2)/(X3-X1)*rho(J3,jr)
                        sumc = sumc + TERMC
                        sume = sume + TERME

6                       continue


           if (sume <= 1.E-8) then
             sume = 0.0
           else
             sume = -1.5 * alpha * (pd * sume)**(1. / 3.)
           endif

           vhar = vhar + sumc
           vex = vex + sume
           jg = mod(IG, 20) + 1
           vg(jg) = vg(jg) + sumc + sume
           ig = ig + 1

7          continue

       dh = dh / 2.0
       nh = nh + nh
       if (nint == 0) goto 1

       ant = float(nint)
       vhar = vhar / ant
       vex = vex / ant
       vint = vhar + vex

       ! ESTIMATE STANDARD DEVIATION
       if (nint < ng) ng = nint
       NAG = nint / ng
       AG = float(NAG)

       do ig=1,ng
         vg(ig) = vg(ig) / AG
       end do

       VAR = 0.0

       do ig=1,ng
         VAR = VAR + (vint - vg(ig))**2
       end do

       VAR = sqrt(VAR / float(ng * (ng - 1)))
       ! THE CURRENT MONTE-CARLO VOLUME FOR THE INTERSTITIAL REGION
       ! IS VOLC
       VOLC = ant / float(npoint) * AV
      ! VOLT IS THE TRUE VOLUME OF THE REGION BETWEEN MUFFIN-TIN
      ! SPHERES IN THE UNIT CELL
       VM = 0.0

       do ir=1,nr
         VM = VM + float(nrr(ir)) * rmt(ir)**3
       end do

       VOLT = AV - 4.0 * pi * VM / 3.0

       write(11,200) nint, npoint, ng, NAG, VOLT, VOLC
       write(11,201) vhar, vex, vint, VAR

       return

200    format(///"MUFFIN-TIN ZERO CALCULATION, SAMPLING WITH",I6,         &
     & " POINTS FROM GRID OF",I6/" VARIANCE ESTIMATED FROM",              &
     & I4," GROUPS OF",I5//" TRUE VOLUME OF INTERSTITIAL REGION",         &
     & F11.4,5X,"MONTE-CARLO VOLUME",11X,F9.4)
201    format(" AVERAGE HARTREE POTENTIAL",6X,F14.5,5X,                   &
     & "AVERAGE EXCHANGE POTENTIAL",F12.5/                                &
     & "MUFFIN-TIN ZERO",F12.5,10X,"STANDARD DEVIATION",F12.5)

      end subroutine

!---------------------------------------------------------------------
!  SUBROUTINE FOR CALCULATION OF THE MUFFIN-TIN ZERO LEVEL FOR
!  MONOATOMIC CRYSTALS, UsinG A SPHERICAL AVERAGE OF THE POTENTIAL
!  BETWEEN MUFFIN-TIN RADIUS AND WIGNER-SEITZ RADIUS, AS IN EQ 3.31
!  OF LOUCKS, TRANSFORMED TO THE EXPONENTIAL GRID RX?
!                       RX(I)=EXP(-8.8+0.05(I-1))
!  INTEGRATION BY TRAPEZIUM RULE.  JRMT,JRWS ARE GRID POINTS OUTSIDE
!  MUFFIN-TIN RADIUS AND WIGNER-SEITZ RADIUS RESPECTIVELY
      subroutine mtzm(vh,vs,rx,ngrid,rmt,rws,jmrt,jrws,vhar,vex)

       integer, intent(in) :: ngrid, jmrt, jrws
       real, dimension(ngrid) :: vh, vs, rx
       real, intent(in) :: rmt, rws
       real, intent(inout) :: vhar, vex
       double precision sumh, sume

       dx = 0.05
       ddx = 0.5 * dx
       dxx = exp(3. * dx)
       X = log(rx(jrmt) / rmt)
       rdx = X / dx
       xx = rx(jrmt - 1)**3
       xxmt = xx * dxx
       sumh = 0.5 * X * (rdx*xx*vh(jrmt-1) + (2.-rdx)*xxmt*vh(jrmt))
       sume = 0.5 * X * (rdx*xx*vs(jrmt-1) + (2.-rdx)*xxmt*vs(jrmt))
       xx = xxmt
       jrw = jrws - 1

       if (jrmt .ne. jrw) then

         vh1 = ddx * xx * vh(jrmt)
         vx1 = ddx * xx * vs(jrmt)
         jrwm = jrmt + 1

         do J=jrwm,jrw
           xx = xx * dxx
           vh2 = ddx * xx * vh(J)
           vx2 = ddx * xx * vs(J)
           sumh = sumh + vh1 + vh2
           sume = sume + vx1 + vx2
           vh1 = vh2
           vx1 = vx2
         end do

       endif

       X = log(rws / rx(jrw))

       rdx = X / dx
       xxws = xx * dxx
       sumh = sumh + 0.5 * X * ((2.-rdx)*xx*vh(jrw)+rdx*xxws*vh(jrws))
       sume = sume + 0.5 * X * ((2.-rdx)*xx*vs(jrw)+rdx*xxws*vs(jrws))
       c = 3. / ((rws * rws * rws) - (rmt * rmt * rmt))
       vhar = c * sumh
       vex = c * sume
       vint = vhar + vex

       write(11,200) vhar, vex, vint

       return

200    format(///"MUFFIN-TIN ZERO BY SPHERICAL AVERAGE",/                 &
     & " AVERAGE HARTREE POTENTIAL",6X,F14.5,5X,                          &
     & "AVERAGE EXCHANGE POTENTIAL",F12.5,/"MUFFIN-TIN ZERO",F12.5)

      end subroutine

!---------------------------------------------------------------------
!  ROUTINE TO SUPPLY NEAREST NEIGHBOUR DATA FOR ATOMS IN
!  A CRYSTAL STRUCTURE, GIVEN?
!  RC(I,J)? THE I'TH COORDINATE OF THE J'TH AXIS OF THE UNIT CELL
!  RK(I,J)? THE I'TH COORDINATE OF THE J'TH ATOM IN THE UNIT CELL
!  NRR(IR)? THE NUMBER OF TYPE-IR ATOMS IN THE UNIT CELL
!  THE INformatION returnED, FOR A TYPE-IR ATOM, IS
!  NCON(IR)? THE NUMBER OF NEAREST NEIGHBOUR SHELLS OF A TYPE-IR
!  ATOM INCLUDED, OUT TO A DISTANCE OF RMAX, BUT <= MC
!  IA(J,IR)? THE TYPE OF ATOMS IN THE J'TH NEIGHBOURING SHELL
!  NA(J,IR)? THE NUMBER OF ATOMS IN THE J'TH SHELL
!  AD(J,IR)? THE RADIUS OF THE J'TH SHELL
      subroutine nbr(ia, na, ad, ncon, nrr, nr, rc, rk, n, rmax, mc)

       use WK
       integer, intent(in)       :: mc, nr, n
       integer, dimension(nr)    :: ncon
       integer, dimension(mc,nr) :: ia, na
       real, dimension(mc,nr)    :: ad
       real, dimension(nr)       :: nrr
       real, dimension(3,3)      :: rc
       real, dimension(3,n)      :: rk
       real rmax

      ! INITIALISATION
       rad(A1,A2,A3) = sqrt((A1 * A1) + (A2 * A2) + (A3 * A3))

       rcmin = 1.0E6

       do I=1,3
         rcmin = min(rcmin, rad(rc(1,I), rc(2,I), rc(3,I)))
       end do

       do ir=1,nr
         do ic=1,MC
           IA(ic,ir) = 0
           NA(ic,ir) = 0
           AD(ic,ir) = 1.0E6
         end do
       end do

       ! SEARCH OVER ADJACENT UNIT CELLS TO INCLUDE MC NEAREST NEIGHBOURS
       ITC = ifix(rmax / rcmin) + 1
       limc = ITC + ITC + 1
       AS = -float(ITC + 1)
       AX = AS

       do 10 JX=1,limc
         AX = AX + 1.0
         AY = AS

         do 10 JY=1,limc
           AY = AY + 1.0
           AZ = AS

           do 10 JZ=1,limc
             AZ = AZ + 1.0

             do J=1,3
               RJ(J) = (AX * rc(J,1)) + (AY * rc(J,2)) + (AZ * rc(J,3))
             end do

             ! RJ IS CURRENT UNIT CELL ORIGIN.
             ! FOR EACH ATOM IN THIS UNIT CELL FIND DISPLACEMENT R
             ! FROM KR-TYPE ATOM IN BASIC UNIT CELL
             J = 0

             do 10 JR=1,nr
               jnr = nrr(JR)

               do 10 jjr=1,jnr
                 J = J + 1 
                 K = 1

                 do KR=1,nr
                   R = rad(RJ(1) + rk(1,J) - rk(1,K),                     &
     &                     RJ(2) + rk(2,J) - rk(2,K),                     &
     &                     RJ(3) + rk(3,J) - rk(3,K))
                   if (R > rmax) exit  ! loop
                   ! COMPARE R WITH NEAREST NEIGHBOUR DISTANCES
                   ! ALREADY FOUND
                   IC = 0
4                  IC = IC + 1
                   if (IC > MC) exit  ! loop
                   DR = R - AD(IC,KR)
                   if (abs(DR) < 1.0E-4) DR = 0.0
                   if (DR) 6,5,4
5                  if (IA(IC,KR) .ne. jr) goto 4
                   NA(IC,KR) = NA(IC,KR) + 1
                   exit  ! loop
6                  if (IC == MC) goto 8
                   IIC = IC + 1

                   do jjc=IIC,MC
                     jc = MC + IIC - jjc
                     IA(jc,KR) = IA(jc-1,KR)
                     NA(jc,KR) = NA(jc-1,KR)
                     AD(jc,KR) = AD(jc-1,KR)
                   end do

8                  IA(IC,KR) = jr

                   NA(IC,KR) = 1
                   AD(IC,KR) = R
                   K = K + nrr(KR)
                 end do
10     continue

       do 12 ir=1,nr
         ncon(ir) = 0

         do IC=1,MC
           if (NA(IC,ir) == 0) goto 12
           ncon(ir) = ncon(ir) + 1
         end do

12    continue

       return

      end subroutine

!---------------------------------------------------------------------
! TAKEN FROM LOUCKS' BOOK, APPENDIX 1
      subroutine poisson(psq, z, j, w)
    
       integer, intent(in) :: j
       real, intent(in)    :: z
       real, dimension(j)  :: psq, w
       double precision E(250),F(250)
       double precision acc,A,B,C,D,C2

       A = 1.0D0 - 0.0025D0 / 48.0D0

       ! EQ. A1.11
       B = -2.0D0 - 0.025D0 / 48.0D0

       ! EQ. A1.12
       C = 0.0025D0 / 6.0D0
       D = exp(0.025D0)
       C2 = -B / A
       E(1) = 0.0D0

       ! EQ. A1.29
       F(1) = D

       ! EQ.A1.30
       X = -8.75
       J1 = J-1

       do I=2,J1
         acc = C * exp(0.5 * X) * (D * PSQ(I + 1) + 10.0 * PSQ(I) + PSQ(I - 1) / D)
         ! EQS. A1.13, A1.6
         F(I) = C2 - 1.0 / F(I - 1)
         ! EQ. A1.20
         E(I) = ((acc / A) + E(I - 1)) / F(I)
         ! EQ. A1.21
         X = X + 0.05
       end do

       W(J) = 2.0 * Z * exp(-0.5 * X)
       acc = W(J)

       ! EQ.A1.15
       do I=1,J1
         jc = J - I
         acc = E(jc) + (acc / F(jc))
         W(jc) = acc
       end do

       return

      end subroutine

!---------------------------------------------------------------------
!  ROUTINE FOR INPUT OF CHARGE DENSITY FROM RELATIVISTIC ORBITALS
!  (ERIC SHIRLEY PROGRAM), AND CALCULATION OF CHARGE DENSITY ON
!  THE RADIAL MESH RX
!    RHO = 4*PI*SUM OVER STATES OF (MODULUS(WAVE FN)**2) *
!          RADIUS**2
!    RMIN= minimum radial coordinate defining the logarithmic mesh used
!          in relativistic calculation
!    RMAX= maximum radial coordinate defining the logarithmic mesh used
!          in relativistic calculation
!    NR  = number of points in the mesh
! the mesh is defined as r(i)=rmin*(rmax/rmin)**(dble(i)/dble(nr))
! FOR EACH ATOMIC STATE I?
      subroutine rela(rho, rx, nx, ngrid)

       use WK
       integer, intent(in) :: ngrid
       integer, intent(inout) :: nx
       real, dimension(ngrid) :: rho, rx
       real name(4)
       real sum

       read(4,100) name, iprint
       read(4,54) rmin, rmax, nr, z
54     format (d15.8,d15.8,i5,f5.2)

       ! initialization of logarithmic grid
       do i=1,nr
         rr(i) = rmin * (rmax / rmin) ** (dble(i) / dble(nr))
       end do

       ns = nr
       ! read in charge density
       read(4,56) (rs(j),j=1,nr)
56     format (f15.10)

       ! INTERPOLATION TO GRID RX
       nx = ngrid
       call chgrid(rs, rr, ns, rho, rx, nx)
       if (iprint == 0) return
       write(11,200) name, (RR(ix),ix=1,NS)

       do ix=1,nx
         if (rho(ix) < 1.0E-9) exit  ! loop
       end do

       nx = ix

       write(11,202) rx(nx), nx, (rho(ix),ix=1,nx)

       return

100    format(4A4/I4)
200    format("1",4A4," RELAT. WAVEFUNCTIONS (ERIC SHIRLEY)",             &
     & " R RADIUS"," LOGARITHMIC MESH",/(10F12.5/))
202    format("CHARGE DENSITY OUT TO RADIUS",F12.5,10X,                   &
     & "NX",I4//5(10E12.4/))

      end subroutine

!---------------------------------------------------------------------
!  ROUTINE TO PERFORM THE SUMMATION OF CONTRIBUTIONS FROM
!  NEIGHBOURING ATOMS (EQ. 3.22,3.26,3.28).  INTEGRATION BY
!  TRAPEZIUM RULE ON RADIAL GRID  RX(I)=EXP(-8.8+0.05(I-1))
      subroutine sumax(acc, chi, rx, nx, ncon, ia, na, ad, imax, ngrid, nr)

       integer, intent(in)       :: ngrid, nr, ncon, imax
       integer, dimension(nr)    :: nx
       real, dimension(ngrid)    :: acc, rx
       real, dimension(ngrid,nr) :: chi
       real, dimension(ncon)     :: ad
       integer, dimension(ncon)  :: ia, na
       double precision sum

       index(X) = (20. * (log(X) + 8.8)) + 2.
       dx = 0.05
       ddx = 0.5 * dx
       dxx = exp(2. * dx)
       IC = IA(1)

       do I=1,imax
         acc(I) = CHI(I,IC)
       end do

       do 4 JA=2,ncon
         IC = IA(JA)
         nix = nx(IC)

         do 4 I=1,imax
           sum = 0.0D0
           X1 = abs(rx(I) - AD(JA))
           ix1 = index(X1)

           if (ix1 > nix) goto 4

           dx1 = log(rx(ix1) / X1)
           rdx1 = dx1 / dx
           X2 = min((rx(I) + AD(JA)), rx(nix))
           ix2 = min0(index(X2), nix)
           dx2 = log(rx(ix2) / X2)
           rdx2 = dx2 / dx
           xx = rx(ix2 - 1)**2
           xx1 = xx * dxx

           if (ix1 == ix2) goto 3

           sum = sum + 0.5 * dx2 * ((2. - rdx2) * xx * CHI(ix2-1,IC) +    &
     &                      rdx2 * xx1 * CHI(ix2,IC))
           xx = rx(ix1 - 1) ** 2
           xx1 = xx * dxx
           sum = sum + 0.5 * dx1 * (rdx1 * xx * CHI(ix1-1,IC) +           &
     &              (2. - rdx1) * xx1 * CHI(ix1,IC))
           ix1 = ix1 + 1
           if (ix1 == ix2) goto 4
           xx = xx1
           T1 = ddx * xx * CHI(ix1,IC)
           ix2 = ix2 - 1

           do ix=ix1,ix2
             xx = xx * dxx
             T2 = ddx * xx * CHI(ix,IC)
             sum = sum + T1 + T2
             T1 = T2
           end do

           goto 4

3          sum = 0.5 * (dx2 - dx1) * ((rdx1 + rdx2) * xx * CHI(ix1-1,IC) + &
     &              (2. - rdx1 - rdx2) * xx1 * CHI(ix1,IC))

4        acc(I) = acc(I) + 0.5 * sum * float(NA(JA)) / (AD(JA) * rx(I))

       return

      end subroutine

!---------------------------------------------------------------------
!  SUBROUTINE PHSH2CAV
!  POTENTIAL-TO-PHASE-SHIFT CALCULATION(CAVLEED PACKAGE)
!  USES LOUCKS GRID (E.G. AS SUPPLIED BY THE MUFFIN-TIN POTENTIAL
!  PROGRAM).  ENERGIES INPUT IN HARTREES.
      subroutine phsh_cav(mufftin_file, phasout_file,                     &
     &                    dataph_file, zph_file)

       character(len=*), intent(inout) :: mufftin_file, zph_file
       character(len=*), intent(inout) :: phasout_file, dataph_file
       dimension V(250), rx(250), phs(20)
       double precision name(2),mtz,delstore(401,15),estore(401)
       integer unit_muff, unit_phas, unit_data, unit_zph

       ! check for null input strings
       if (len_trim(mufftin_file) < 1) mufftin_file = "mufftin.d"
       if (len_trim(phasout_file) < 1) phasout_file = "phasout"
       if (len_trim(dataph_file) < 1) dataph_file = "dataph"
       if (len_trim(zph_file) < 1) zph_file = "zph.o"

       ! get free file units
       unit_muff = 5
       unit_data = 6
       unit_phas = 7
       unit_zph = 8

       ! First input channels
       open(unit=unit_muff,file=mufftin_file,status='old')

       ! Now output channels
       open(unit=unit_zph,file=zph_file,status='unknown')
       open(unit=unit_phas,file=phasout_file,status='unknown')
       open(unit=unit_data,file=dataph_file,status='unknown')

      ! standard values for phase shifts calculation
       write(unit_data,110)
       emin = 1.
       emax = 12.
       estep = .25
       ianz = (emax - emin) / estep + 1.01
       nl = 12
       read(unit_muff,103) nr

       do KKK=1,nr
         read(unit_muff,100) name
         read(unit_muff,101) Z,rmt, mtz
         read(unit_muff,103) NTAB
         mtz=mtz/2.

         do ix=1,NTAB
           read(unit_muff,219) rx(ix), V(ix)
         end do

         write(unit_zph,200)name,Z,rmt,mtz
         write(unit_data,181)(name(I),I=1,2)
 181     format('NON-RELATIVISTIC PHASE SHIFTS FOR ',2A4)
         write(unit_data,1030) emin,estep,IANZ,nl
 1030    format(2F9.4,2(2X,I3))
         E = EMIN
         ncount = 0
1        E = 2. * E
         ncount = ncount + 1
         call PS(V,rx,NTAB,rmt,E,phs,nl,unit_zph)
         E = 0.5 * E
         write(unit_zph,201) E, (phs(L),L=1,nl)
         write(unit_data,1040) E*27.21, (phs(L),L=1,nl)
 1040    format(F9.4,8F8.4)

         ! store phase shifts
         do kk=1,nl
           delstore(ncount,kk) = phs(kk)
         end do

         estore(ncount)=E
         E=E+ESTEP
         if (E <= EMAX) goto 1

         ! write phase shifts as function of energy for plotting
         do kk=1,nl
           write(unit_data,107) kk-1

           do ii=1,ncount
             write(unit_data,*) estore(II),delstore(II,kk)
           end do

           write(unit_data,*)
         end do

       end do

       ! close file handles
       close(unit_muff)
       close(unit_data)
       close(unit_zph)
       close(unit_phas)

       return

100    format(2A8)
101    format(3F8.4)
103    format(I4)
107    format('"L=',i2)
110    format ("TitleText: ","DELTA(E)")
200    format("PHASE SHIFTS FOR ",2A8,2X,"ATOMIC NUMBER,",F6.1/           &
     & "MUFFIN-TIN RADIUS",F8.4,6X,                                       &
     & " MUFFIN-TIN ZERO LEVEL",F8.4," HARTREES")
201    format("ENERGY",F8.4," HARTREES"/(10F12.5))
219    format(2E14.5)

       return

      end subroutine

!***********************************************************************
!  SUBROUTINE TO CALCULATE NL PHASE SHIFTS (L=0,NL-1) FOR AN
!  ATOMIC POTENTIAL TABULATED ON THE LOUCKS RADIAL GRID.
!     V?  ATOMIC POTENTIAL (RYDBERGS)
!     RX?  LOUCKS' EXPONENTIAL GRID  RX(I)=EXP(-8.8+0.05(I-1))
!     NGRID?  NUMBER OF TABULATED POINTS
!     RAD?  LIMIT OF INTEGRATION OF SCHRODINGER EQUATION (A.U.)
!     E?  ENERGY (RYDBERGS)
!     PHS?  PHASE SHIFTS FOR L=0 TO NL-1
!  REFERENCE? LOUCKS T L, (1967), A.P.W. METHOD, BENJAMIN, NY.
!***********************************************************************
      subroutine ps(v, rx, ngrid, rad, e, phs, nl, file_unit)

       integer, intent(in)    :: ngrid, nl
       integer, intent(in)    :: file_unit
       real, dimension(ngrid) :: v, rx
       real, dimension(nl)    :: phs
       real rad, e
       real, dimension(250)   :: WF
       real, dimension(25)    :: bj, bn
       real, dimension(10)    :: xr
       real, dimension(5)     :: fr
       data pi, dx, DAC, DB/3.141592653589, 0.05, 2.083333E-04, 2.083333E-03/

       index(X)=20.*(log(X)+8.8)+2.
       ! TABULATION OF SPHERICAL BESSEL FUNCTIONS IN BJ AND BN
       es = sqrt(E)
       X = es * rad
       Z = X
       ll = nl + 1
       call calcbf(bj, bn, ll, X, FILE_UNIT)

       ! INTEGRATION OF THE RADIAL SCHRODINGER EQUATION BY THE NUMEROV
       ! METHOD (SEE LOUCKS P56 FF). WF CONTAINS WAVE function X RADIUS,
       ! ON THE LOUCKS GRID
       X1=exp(-8.8)

       do L1=1,nl
         fl = L1 - 1
         fl2 = fl + 0.5
         FF = fl2 * fl2
         Y1 = X1 ** fl2
         Y2 = exp(dx * fl2) * Y1

         write(FILE_UNIT,60)fl,Y1,Y2
60       format('0L',F5.1,5X,'Y1,Y2',2E14.5)

         GAM1 = FF + (rx(1) * rx(1) * (V(1) - E))
         GAM2 = FF + (rx(2) * rx(2) * (V(2) - E))
         WF(1) = Y1 * sqrt(rx(1))
         WF(2) = Y2 * sqrt(rx(2))

         do ix=3,ngrid
           GAM = FF + (rx(ix) * rx(ix) * (V(ix) - E))
           A = 1. - (DAC * GAM)
           B = -2. - (DB * GAM2)
           C = 1. - (DAC * GAM1)
           YN = -((B * Y2) + (C * Y1)) / A
           WF(ix) = YN * sqrt(rx(ix))
           Y1 = Y2
           Y2 = YN
           GAM1 = GAM2
           GAM2 = GAM
         end do

         ! LAGRANGIAN INTERPOLATION FOR WAVEFUNCTION AND DERIVATIVE AT
         ! RADIUS X.  WFN HOLDS WAVEfunction X RADIUS, AND DWFN
         ! DERIVATIVE X RADIUS
         X = rad
         jr = index(rad)

         do J=1,5
           xr(J) = rx(jr-5+J)
           xr(J+5) = xr(J)
           fr(J) = WF(jr-5+J)
         end do

         wfn = 0.
         dwfn = 0.
         A = (X - xr(1)) * (X - xr(2)) * (X - xr(3)) * (X - xr(4)) * (X - xr(5))

         do I=1,5
           TERM=A/(X-xr(I))/(xr(I)-xr(I+1))/(xr(I)-xr(I+2))               &
     &      /(xr(I)-xr(I+3))/(xr(I)-xr(I+4))
           sum=0.

           do J=1,5
             if (I == J) exit  ! loop
             sum = sum + (TERM / (X - xr(J)))
           end do

           wfn = wfn + (TERM * fr(I))

           dwfn = dwfn + (sum * fr(I))
         end do

         ! LOGARITHMIC DERIVATIVE
         dloga = dwfn / wfn - 1. / rad

         ! PHASE SHIFTS
         X = es * rad
         A = fl * bj(L1)/X-bj(L1+1)
         B = fl * bn(L1)/X-bn(L1+1)
         A = es * A - dloga * bj(L1)
         B = es * B - dloga * bn(L1)
         phs(L1)=pi/2.
         if (abs(B) > 1.0E-8) phs(L1) = atan(A / B)
         write(FILE_UNIT,78) phs(L1)
78       format('PHASE SHIFT',F10.4)

       end do

       return

      end subroutine

!***********************************************************************
      subroutine calcbf(bj, bn, nl, x, file_unit)

       integer, intent(in)    :: nl
       integer, intent(in)    :: file_unit
       real, dimension(nl)    :: bj, bn
       real x

       if (abs(X) < 1.0E-6) then
         write(FILE_UNIT,200) X
200      format("** ARGUMENT",E12.4," TOO SMALL FOR ROUTINE CALCBF")
         return
       endif

       bj(1) = sin(X) / X
       bn(1) = -cos(X) / X

       if (nl == 1) return

       bj(2) = (bj(1) -cos(X)) / X
       bn(2) = (bn(1) -sin(X)) / X

       if (nl == 2) return

       if (float(nl * (nl + 1)) <= X * X) then

         ! FORWARD RECURRENCE FOR BJ'S
         fl=3.0

         do L=3,nl
           bj(L) = fl * bj(L-1) / X - bj(L-2)
           fl = fl + 2.
         end do

         ! FORWARD RECURRENCE FOR BN'S
         fl = 3.

         do L=3,nl
           bn(L) = fl * bn(L-1) / X - bn(L-2)
           fl = fl + 2.
         end do

         return

       endif

       ! BACKWARD RECURRENCE FOR BJ'S
       bj0 = bj(1)

       bj1 = bj(2)
       nn = max(10, 2 * nl)
       A = 0.
       B = 1.
       fl = float((2 * nn) + 1)

       do I=1,nn
         L = nn - I + 1
         C = fl * B / X - A
         if (L <= nl) bj(L) = C
         A = B
         B = C
         fl = fl - 2.
       end do

       ! NORMALISATION
       B = bj0 / bj(1)
       if (abs(bj0) < 0.01) B = bj1 / bj(2)

       do L=1,nl
         bj(L) = B * bj(L)
       end do

       ! FORWARD RECURRENCE FOR BN'S
       fl = 3.

       do L=3,nl
         bn(L) = fl * bn(L - 1) / X - bn(L-2)
         fl = fl + 2.
       end do

       return

      end subroutine

!---------------------------------------------------------------------
!  subroutine phsh_wil
!---------------------------------------------------------------------
!  A.R. WILLIAMS^ PHASE SHIFT PROGRAM (GIVEN A MUFFIN-TIN POTENTIAL)
      subroutine phsh_wil(mufftin_file, phasout_file,                     &
     &                    dataph_file, zph_file)

       use CM16
       use CMRV
       use CM5
       character(len=*), intent(inout)     :: mufftin_file, zph_file
       character(len=*), intent(inout)     :: phasout_file, dataph_file
       real E(401), S(401,15), C(401,15), del(15), delold(15)
       real dell(9), delstore(8,401,15)
       integer tlp1
       namelist / nl2 / ip, nrr

       ! check for null input strings
       if (len_trim(mufftin_file) < 1) mufftin_file = "mufftin.d"
       if (len_trim(phasout_file) < 1) phasout_file = "phasout"
       if (len_trim(dataph_file) < 1) dataph_file = "dataph"
       if (len_trim(zph_file) < 1) zph_file = "zph.o"

       ! First input channels
       open(unit=5,file=mufftin_file,status='OLD')

       ! Now output channels
       open(unit=6,file=zph_file,status='UNKNOWN')
       open(unit=7,file=phasout_file,status='UNKNOWN')
       open(unit=8,file=dataph_file,status='UNKNOWN')

       pi = 3.1415926535

       ! READ IP
       ! IP=0: ONLY RADIAL WAVEfunction
       ! IP=1: PHASE SHIFTS IN ADDITION
       ! IP=2: S AND C
       ! IP=3: PRODUCE LOGARITHM OF PHASE SHIFTS
       ! NRR= number of inequivalent atoms for which we want phase shifts
       IP = 1
       read(5, nl2)
       write(6, nl2)
       write(8,110)
110    format ("TitleText: ","DELTA(E)")

       ! INPUT
       do 2  KKK=1,nrr
         call S16
         TX = 2. * R(nr)/ float(nr - 1)
         de = (E2 - E1) / float(max(NE - 1, 1))

         do I = 1, NE
           E(I) = E1 + float(I - 1) * de
           ! RADIAL INTEGRATION
           call S10(E(I))
           T3 = R(nr) * E(I)
           T4 = R(nr) * T3

           do lp1 = 1, nl
             L = lp1 - 1
             tlp1 = 2 * L + 1
             T5 = R(nr) ** lp1
             UT = F(tlp1, ilst)/TX + float(L) * Y(tlp1,ilst) / R(nr)
             T1 = (F44(L,T4) * Y(2*lp1,ilst) + T3 * F44(lp1,T4)           &
     &              * Y(tlp1,ilst)) * T5
             T2 = (F45(L,T4) * UT - T3 * F45(L - 1,T4)                    &
     &              * Y(tlp1,ilst)) * R(nr) / T5
             S(I, lp1) = T1
             C(I, lp1) = T2
           end do

         end do

         IS = 2
         I4 = 9
         if (IP < 1) goto 15

         ! PRODUCE PHASE SHIFTS
         do lp=1,nl
           delold(lp)=0.0
         end do

         do I = 1, NE
           do lp = 1, nl
             del(lp) = atan(-abs(E(I))**(lp-.5)*S(I,lp) / C(I, lp))
           end do

           ! REMOVE DISCONTINUITIES BY MULTIPLES OF PI
           do lp=1,nl
             LS = 0
111          deldif = del(lp) - delold(lp)
             if (abs(deldif) < 0.7) exit  ! loop
             LS = LS + 1
             del(lp) = del(lp) - sign(pi, deldif)
             if (ls < 5) goto 111
             write (6,115) lp
115          format(" TOO LARGE CHANGE IN PHASE SHIFT [L=",1I4,           &
     &              "] SINCE LAST ENERGY ",/,                             &
     &              " DISCONTINUITY BY MULTIPLE OF PI POSSIBLE")
             delold(lp) = del(lp)
           end do


           if (neuo == 2) E(I) = 0.5 * E(I)

           ! PRINT PHASE SHIFTS
           write(6, 12) E(I), (del(lp), lp = 1, nl)

           ! WRITE PHASE SHIFTS IN FORMAT USED BY LEED PROGRAM
           ! write(7,71) E(I),(del (lp),lp=1,nl)

           ! store phase shifts
           do kk=1,nl
             delstore(KKK,I,kk) = del(kk)
           end do

           if (IP < 3) exit  ! loop

           do J = 1, 9
             dell(J) = -4.
             if (del(J) < 1.0E-4) exit  ! loop
             dell(J) = log10(del(J))
           end do

         end do

         ! write phase shifts as function of energy for plotting
         do kk=1,nl
           write(8,100) kk-1
100        format('"L=',i2)

           do ii=1,NE
             write(8,*) E(II), delstore(kkk,II,kk)
           end do

           write(8,*)
         end do

15     continue

2      continue

       write(7,*) 'BE CAREFUL ABOUT THE ORDER OF THE ELEMENTS'

       do ii=1,NE
         write(7,71) E(II)

         do i=1,nrr
           write(7,72)(delstore(i,ii,lp),lp=1,nl)
         end do

       end do

       if (IP >= 2) then
         do lp1 = 1, nl
           call S41(E, S(1, lp1), NE)
           call S41(E, C(1, lp1), NE)
         end do
       endif

       ! close file handles
       close(5)
       close(6)
       close(7)
       close(8)

       return

12     format(1P8E14.7, /, 14X, 1P7E14.7, /)
71     format(1F7.4)
72     format(10F7.4)

      end subroutine

!***********************************************************************
!  subroutine S16
!  S16 INPUTS data
!   CS: CORE SHIFT (POSITION OF ZERO OF ENERGY)
!   Z: ATOMIC NUMBER
!   E1,E2: FIRST AND LAST ENERGIES DESIRED (IN RYDBERGS OR
!           HARTREES, CF. NEUI)
!   NE: NUMBER OF ENERGIES DESIRED
!   NL: NUMBER OF PHASE SHIFTS DESIRED (=LMAX+1)
!   NR: NUMBER OF RADIAL GRID POINTS USED IN CALCULATION
!   IX=0: SIGNAL TO STOP
!   IX=1: SIGNAL TO EXPECT A (NEW) POTENTIAL TO BE READ IN
!   RT: MUFFIN-TIN RADIUS (IN BOHR RADII)
!   NEUI,NEUO
!    if =1: RYDBERG UNIT USED FOR INPUT (NEUI) AND OUTPUT (NEUO) OF
!           ENERGIES AND POTENTIAL
!    if =2: HARTREE UNIT (double RYDBERG) USED INSTEAD OF RYDBERG
!           UNIT FOR INPUT (NEUI) AND OUTPUT (NEUO)
!   POTYP=1: RADIAL INPUT AS V(R)
!   POTYP=2: RADIAL INPUT AS R*V(R)
!
!***********************************************************************
      subroutine S16

       use CM16
       use CMRV
       real, dimension(200) :: rs, zs
       real, dimension(201) :: ztt
       dimension fmt(18)
       namelist / nl16 / cs,Z,E1,E2,NE,nl,nr,ix,RT,neui,neuo,potyp

       ! SET DEFAULT VALUES OF VARIABLES IN NAMELIST /NL16/
       ix = 1
       E1 = 4.
       E2 = 24.0
       NE = 30
       nl = 9
       nr = 101
       neui = 1
       neuo = 2
       potyp = 2
       read(5, nl16)
       cs = 0.0

       if (ix < 1) return

       if (neui .ne. 1) then
         cs = 2.0 * cs
         E1 = 2.0 * E1
         E2 = 2.0 * E2
       endif

       write (6, nl16)

       DRDN2 = (float(nr - 1))*(float(nr - 1)) / RT
       ! READ format USED FOR INPUT OF R VS. V(R) OR R VS. R*V(R)
       ! (V IS ASSUMED POSITIVE)
       do I=1,200
         read(5,111) RS(I), ZS(I)
  111    format(2E14.5)
         ! the next lines assume that input potential & cs are negative
         ZS(I) = -ZS(I)

         if ( RS(I) < 0) exit  ! loop

       end do

       nrs = I - 1

       if (neui .ne. 1) then

         do I=1,nrs
           ZS(I) = 2.0 * ZS(I)
         end do

       endif

       if (potyp .ne. 2) then

         do I=1,nrs
           ZS(I) = (ZS(I) - cs) * RS(I)
         end do

         goto 21

       endif

       do I = 2, nrs
         ZS(I) = (ZS(I) / RS(I) - cs) * RS(I)
       end do

21     IV = 1
       R(1) = 0.
       ZTT(1) = Z + Z

       do I = 2, nr
         R(I) = (float(I - 1)) ** 2 / DRDN2

40       if ((R(I) <= RS(IV + 2)).or.(IV + 3 >= nrs)) goto 50

         IV = IV + 1
         goto 40

50       ZTT(I) = F12(RS(IV), ZS(IV), R(I), 4)

         do lp1 = 1, nl
           V(I, lp1) = -ZTT(I) / R(I)
         end do

       end do



       return

      end subroutine

!***********************************************************************
!  F12 PERFORMS ITERATIVE INTERPOLATION IN A TABLE OF N VALUES OF
!  X AND Y TO FIND THE VALUE OF Y AT Z
!***********************************************************************
      function f12(x, y, z, n)

       integer, intent(in) :: n
       real, intent(in)    :: z
       real, intent(in)    :: x(10), y(10)
       real w(20)

       w(1) = y(1)

       do i=2,n
         w(i) = y(i)
         u = z - x(i)
         ip1 = i + 1

         do j=2,i
           k = ip1 - j
           w(k) = w(k + 1) + u * (w(k) - w(k + 1)) / (x(k) - x(i))
         end do

       end do

       f12 = w(1)

      return

      end function f12

!***********************************************************************
!  S5 -- HAMMING^S METHOD FOR THE INTEGRATION OF SYSTEMS OF FIRST
!  ORDER DIFFERENTIAL EQUATIONS
!***********************************************************************
      subroutine S5(e)

       use CMRV
       use CM5
       real, intent(in) :: e
       real eest(30), vme(15)

       nj = 2 * nl

       do J = 1, nj
         eest(J) = 0.
       end do

       do I = 5, nr

         do lp1 = 1, nl
           vme(lp1) = (V(I, lp1) - E) * R(I)
         end do

         T1 = 2. / float(I - 1)
         IP1 = mod(I - 1, 4) + 1
         IM2 = mod(IP1, 4) + 1
         IM1 = mod(IM2, 4) + 1
         IP0 = mod(IM1, 4) + 1

         do J = 1, nj
           F(J, IM2) = Y(J, IP1) + (2. * (F(J, IP0) + F(J, IM2)) - F(J, IM1)) / 0.75
           Y(J, IP1) = F(J, IM2) - 0.925619835 * eest(J)
         end do

         do J = 1, nj, 2
           jp1 = J + 1
           lp1 = jp1 / 2
           flp1 = lp1
           F(J, IP1) = (flp1 * Y(J, IP1) + R(I) * Y(jp1, IP1)) * T1
           F(jp1, IP1) = (vme(lp1) * Y(J, IP1) - flp1 * Y(jp1, ip1))*T1
         end do

        do J = 1, nj
          Y(J, ip1) = Y(J, ip0)+(Y(J, ip0)-Y(J, IM2) + 3.*(F(J, ip1)      &
     & + 2. * F(J, ip0) - F(J, IM1))) / 8.
          eest(J) = F(J, IM2) - Y(J, ip1)
          Y(J, ip1) = Y(J, ip1) + 0.743801653E-1 * eest(J)
        end do

        do J = 1, nj, 2
          jp1 = J + 1
          lp1 = jp1 / 2
          flp1 = lp1
          F(J, ip1) = (flp1 * Y(J, ip1) + R(I) * Y(jp1, ip1)) * T1
          F(jp1, ip1) = (vme(lp1) * Y(J, ip1) - flp1 * Y(jp1, ip1))*T1
        end do

       end do

       return

      end subroutine
!***********************************************************************
!  S10   POWER SERIES EXPANSION OF THE SOLUTION ABOUT THE ORIGIN
!        AND RADIAL INTEGRATION IN S5
!***********************************************************************
      subroutine s10(e)

       use CMRV
       use CM5
       real, intent(in) :: e
       integer tlp1, i, j, k
       real a(10), b(10), tr(4)
       

       ni = 2 * nl
       tz = 2. * z
       a(1) = 1.

       do i = 1, ni, 2
         lp1 = (I + 1) / 2
         tlp1 = 2 * lp1
         ep = E - V(4, lp1) - TZ/ R(4)
         Y(I, 1) = 0.
         Y(I + 1, 1) = 0.
         A(1) = A(1) / float(2 * lp1 - 1)
         B(1) = - Z * A(1) / float(lp1)

         do j=2,4
           TR(J) = R(J) ** lp1
           Y(I, J) = A(1) * TR(J)
           Y(I + 1, J) = B(1) * TR(J)
         end do

         do k=1,9
           A(K + 1) = B(K) / float(K)
           B(K + 1) = -(ep * A(K) + TZ * A(K + 1)) / float(tlp1 + K)

           do j=2,4
             TR(J) = TR(J) * R(J)
             Y(I, J) = Y(I, J) + TR(J) * A(K + 1)
             Y(I + 1, J) = Y(I + 1, J) + TR(J) * B(K + 1)
           end do

           if (abs(TR(4) * A(K + 1) / Y(I, 4)) < 1.0E-4) exit  ! loop

         end do

         write (6, 4) E, lp1, R(4), (A(K), K = 1, 10)
4        format(1PE10.2, I10, 11E10.2)

       end do

       do j = 2, 4
         T1 = 2. / float(J - 1)

         do i = 1, NI, 2
           ip1 = I + 1
           lp1 = ip1 / 2
           flp1 = lp1
           F(I, J) = (flp1 * Y(I, J) + R(J) * Y(ip1, J)) * T1
           F(ip1,J)=((V(J, lp1)-E) * R(J) * Y(I, J) - flp1 * Y(ip1, J))* T1
         end do

       end do

       call s5(e)

       return

      end subroutine

!***********************************************************************
!  F44  EVALUATES THE SPECIAL VERSION OF THE SPHERICAL BESSEL FUNCT.
!***********************************************************************
      function f44(l, x)

       integer, intent(in) :: l
       real, intent(in) :: x
       real, dimension(20) :: s
       real fi, t1, dt, t
       integer i, js, k, is

       js = l + l + 1

       if (abs(x / float(js)) > 10.) goto 5

       fi = 1.

       if (l >= 1) then

         do k = 3, js, 2
           fi = fi * float(k)
         end do

       endif

       t1 = 1. / fi
       dt = 1.
       t = 1.
       i = 0

       do k = 1, 100
         i = i + 2
         dt = -dt * x / float(i * (i + js))
         t = t + dt

         if (abs(dt) < 1.E-8) exit  ! loop

       end do

       t1 = t1 * t
       f44 = t1

       return

5      t = sqrt(abs(x))

       if (x < 0.) goto 9

       s(2) = sin(t) / t

       if (l > 0) goto 6

11     f44 = s(2)

       return

6      s(1) = cos(t)

       goto 10

9      s(2) = sinh(t) / t

       if (l < 1) goto 11

       s(1) = cosh(t)

10     is = l + 2

       do i=3,is
         s(i) = (s(i - 1) * float(2 * i - 5) - s(i - 2)) / x
       end do

       f44 = s(is)

       return

      end function

!***********************************************************************
!  F45  EVALUATES SPECIAL VERSION OF THE SPHERICAL NEUMANN function
!***********************************************************************
      function F45(l, x)

       real, intent(in) :: l, x
       real s(20)

       if (L < 0) then
         F45 = -F44(L+1, X)
         return
       end if

       lp1 = L + 1
       JS = L + L + 1

       if (abs(X / float(JS)) > 10.) goto 5

       FI = 1.

       if (L >= 1) then

         do K = 3, JS, 2
           FI = FI * float(K)
         end do

       end if

       T1 = FI / float(JS)
       DT = 1.
       T = 1.
       I = 0

       do K = 1, 100
         I = I + 2
         DT = -DT * X / float(I * (I - JS))
         T = T + DT

         if (abs(DT) < 1.E-8) exit  ! loop

       end do

       T1 = T1 * T
       F45 = T1

       return

5      T = sqrt(abs(X))

       if (X < 0.) goto 9

       S(2) = cos(T)

       if (L > 0) goto 6

11     F45 = S(2)

       return

6      S(1) = -sin(T) / T

       goto 10

9      S(2) = cosh(T)

       if (L < 1) goto 11

       S(1) = -sinh(T) / T

10     IS = L + 2

       do I = 3, IS
         S(I) = S(I - 1) * float(2 * I - 5) - X * S(I - 2)
       end do

       F45 = S(IS)

       return

      end function

!***********************************************************************
!  S41 PLOTS Y AGAINST X
!***********************************************************************
      subroutine s41(x, y, n)

       integer, intent(in) :: n
       real, dimension(100), intent(in) :: x, y
       character B, C, O, D
       character, dimension(97)  :: p
       real y1, y2
       integer i
       data B, C, O, D / " ", "*", "0", "I" /
       
       Y1 = 0
       Y2 = 0

       do I = 1, N
         Y1 = min(Y1, Y(I))
         Y2 = max(Y2, Y(I))
       end do

       do I = 1, 97
         P(I) = B
       end do

       T = 96/ (Y2 - Y1)
       J0 = -Y1 * T + 1.5
       P(J0) = "0"

       if (N >= 30) write(6, 3) P
3      format("1", 34X, 97A1, //)
       if (N < 30) write(6, 6) P
6      format(////, 35X, 97A1, //)

       P(J0) = D

       do I = 1, N
         J = T * (Y(I) - Y1) + 1.5
         P(J) = C
         write(6, 4) X(I), Y(I), P
4        format(1X, 1P2E16.6, 2X, 97A1)
         P(J) = B
         P(J0) = D
       end do

       return

      end subroutine

!---------------------------------------------------------------------
!  subroutine phsh_rel
!---------------------------------------------------------------------
      subroutine phsh_rel(mufftin_file, phasout_file, dataph_file, inpdat_file)

       use omp_lib
       use ZZZZ
       use Z
       implicit double precision (A-H,O-Z)
       character(len=*), intent(inout)     :: mufftin_file, inpdat_file
       character(len=*), intent(inout)     :: phasout_file, dataph_file
       character(len=6) :: SS1, SS2, WRD
       character, dimension(4,5) :: AMS
       character(len=3) :: opt, opt1, opts, SUB, record
       character(len=2) :: aname, AN*30, BDATA*28
       character(len=1) :: TL,SL
       integer jfs, jfd, jfu
       real name(4)
       integer, dimension(250,18) :: jf
       double precision, dimension(250) :: energ
       double precision, dimension(7) ::  adata
       data ZERO/0.0D0/,ONE/1.0D0/,TWO/2.0D0/,ANINE/9.0D0/,HALF/.5D0/
       data zilch/1.0D-4/,TOL/.005D+0/,DES/.025D+0/
       data AMS/'NC= ','L= ',' ES=',' DE=','ID= '/
       data TL/'L'/,SL/'S'/,SS1/'NOSPIN'/,SS2/' SPIN '/
       data SUB/'SUB'/,record/'NOS'/
  1    format (3D12.4,4X,I3,4X,D12.4)
  2    format (5E14.6)
  8    format (F9.4,8F8.5)
 12    format (5D14.6)
 13    format ("1",//,T61,'INPUT data',//)
 18    format(F10.4,F9.4,2I5)

       ! check for null input strings
       if (len_trim(mufftin_file) < 1) mufftin_file = "mufftin.d"
       if (len_trim(phasout_file) < 1) phasout_file = "phasout"
       if (len_trim(dataph_file) < 1) dataph_file = "dataph"
       if (len_trim(inpdat_file) < 1) inpdat_file = "inpdat"

       open(unit=4,file=inpdat_file,status='unknown')
       open(unit=5,file=mufftin_file,status='old')
       open(unit=7,file=phasout_file,status='unknown')
       open(unit=8,file=dataph_file,status='unknown')

       pi = 4.D0 * dataN(1.D0)
       pi2 = 0.5D0 * pi
       write (4,13)

       do while (.True.)
         read(5,"(4A4)",end=999,err=999) (name(I),I=1,4)
         read(5,1)ES,de,ue,lsm,VC

         ! nl is the number of plotted phase shifts
         nl = 8
         write (4,11) ES, de, ue, opt, opt1, lsm
 11      format (3D12.4,4X,2A3,I3)
         read (5,16) NZ,Adata(1),jri,ALC,BLC,CLC,EXCA,EXCB,EXCO
 16      format (I4,F10.6,I4,T21,6F10.6)
         write (4,76) NZ,Adata(1),jri,ALC,BLC,CLC,EXCA,EXCB,EXCO
 76      format (I4,F10.6,I4,2X,6F10.6)

         vs = 0.
         if (opts == SUB) vs = VC
         if ((jri <= 0) .or. (jri > 340)) goto 999
         read(5,2) (ZP(J),J=1,jri)
         write (4,12) (ZP(J),J=1,jri)
         rhoz = -0.90306514D+01
         delrho = 0.3125D-01
         RM = exp(rhoz)
         xrx = exp(delrho)

         !$OMP PARALLEL DO
         do J=1,jri
           if (ZP(J) < 0.0) ZP(J) = -ZP(J)
           if (J < jri) RM = xrx * rm
         end do
         !$OMP END PARALLEL DO

         if (de <= ZERO) then
           ES = -HALF
           de = des
           ue = ONE
         endif

         N=(ue-ES)/de+HALF

         N=N+1

         write(07,181)(name(I),I=1,4)
 181     format('RELATIVISTIC PHASE SHIFTS FOR ',4A4)
         write(07,18) ES,de,N,lsm

         ES=ES/13.6
         de=de/13.6
         ue=ue/13.6
         if (N > 250)  N=250  ! this is suspicious...
         L=0
         E=ES
         ipt=2

         if (opt == record) then
           ipt=-2
           WRD=SS1
         else
           WRD=SS2
         endif

         KAP=-1
         L=1

         !$OMP PARALLEL DO
         do J=1,N
           dxaz = 0.0D0
           ttr = dlgkap(E,KAP) / (12.5663706)
           call sbfit(ttr, E, L - 1, rmaxi, jfs)
           JF(J,L) = jfs
           E = E*13.6
           energ(J) = E
           E = E / 13.6
           E = E + de
         end do
         !$OMP END PARALLEL DO

         do while(L <= lsm)

           KAP=-(L+1)
           LIND=L+1
           E=ES
           lvor=0

           !$OMP PARALLEL DO
           do J=1,N
             DLU=dlgkap(E,KAP)/(12.5663706)
             DLD=dlgkap(E,L)/(12.5663706)
             call sbfit(DLD,E,L,rmaxi,jfd)
             call sbfit(DLU,E,L,rmaxi,jfu)
             LK = 0
             zfdiff = -(jfd - jfu) * lvor
             if (zfdiff > 2.5) LK = L
             if (zfdiff < -2.5) LK = L + 1
             jfs = L * jfd - KAP * jfu + lvor * LK * pi
             jfs = jfs / (2 * L + 1)
             if (jfs > pi2) jfs=jfs-pi2*2.
             if (jfs < -pi2) jfs=jfs+pi2*2.
             JF(J,LIND)=jfs
             if (LK == 0) lvor = real(sign(1, jfs))
             E=E+de
           end do
           !$OMP END PARALLEL DO

           L=L+1

         end do

         do I=1,N
           lsm1=lsm+1
           write(7,8) energ(I),(JF(I,L),L=1,lsm1)

         end do

         do kk=1,nl
           write(8,100) kk-1

           do ii=1,N
             write(8,*) energ(II),JF(II,kk)
           end do

           write(8,*)

         end do

100      format('"L=',i2)
         ES=ES*13.6
         de=de*13.6
         ue=ue*13.6

       end do

 999   continue

       write (4,900)
 900   format (//,T57,'end OF INPUT data')

       ! close file handles
       close(4)
       close(5)
       close(7)
       close(8)

       return

      end subroutine

!***********************************************************************
! DLGKAP CALCULATES THE LOGRITHMIC DERIVATIVE OF THE LARGE
! COMPONENT UsinG THE PROCEDURE DESCRIBED BY LOUCKS IN APPENDIX 7.
! THE SMALL MULTIPLICATIVE FACTOR IS INCLUDED.
! POTENTIAL  IN THE FORM OF 2ZP IS TO BE PASSED IN THE ZZZZ MODULE
! THE RADIAL functionS ARE MADE AVAILABLE IN THE RADFUN MODULE
! WABER MESH (XNOT=-9.03065 AND DX=1/32) IS USED
! JRI IS THE MESH POINT OF THE APW SPHERE RADIUS
! E IS THE ENERGY TO BE USED (IN RYDBERGS)
! 4 PI R**2 INSERTED NOW. FOR COMPOUNDS ONLY.
!***********************************************************************
      function dlgkap(e, kappa)

       use RADFUN
       use ZZZZ
       use Z
       implicit double precision(A-H,O-Z)
       double precision e, kappa
       double precision, dimension(4)   :: sxk, sxm
       data ustart/1.D-25/,zilch/1.D-30/
       data TEST/1.D+6/,XS/-9.03065133D+00/
       data dx/3.125D-2/,C/2.740746D+2/,cin/1.3312581146D-5/,hf/.5D+0/
       data TH/.3333333333D+0/,T2/2.D+0/,T7/7.D+0/,T11/11.D+0/
       data T12/12.D+0/,T14/14.D+0/,T26/26.D+0/,T32/32.D+0/,ZERO/.1D+0/

       !SET UP FOR RELATIVISTIC OR NO RELATIVISTIC EFFECT
       if (ipt <= 0) then
         cin = 0.0D00
       end if

       ! SET UP STARTING VALUES
       dx2 = hf * dx
       X20 = 0.3 * dx
       XMFT = 4.4444444444444D-2 * dx
       TS = exp(XS)
       TDX = exp(dx)
       hoc = (vcz * TS + pot(1))/C
       XK = kappa
       U(1) = ustart
       if (abs(hoc / XK) > 0.05) then
         P = (XK + sqrt((XK * XK) - (hoc * hoc))) / hoc
       else
         P = (XK + (abs(XK)) / hoc) - (hf * hoc / abs(XK))
       endif

       TC = ((E + vcz) * TS) + pot(1)

       VC = cin * TC
       W(1) = C * P * ustart

       ! START RUNGE-KUTTE PROCEDURE
       X = XS

       N = 1

 25    IK = 0

       np1 = N + 1
       xc = X
       bgc = pot(N)
       WC = W(N)
       UC = U(N)

 20    IK = IK + 1

       T = exp(xc)
       TC = ((E + vcz) * T) + bgc
       VC = cin * TC

       sxk(IK) = dx2 * ((WC * (VC + T)) - (XK * UC))

       sxm(IK) = dx2 * ((XK * WC) - (TC * UC))

       select case(IK)
         case(1)
           xc = xc + dx2
           UC = UC + sxk(1)
           WC = WC + sxm(1)
           bgc = hf * (bgc + pot(np1))
           goto 20

         case(2)
           UC = UC + sxk(2) - sxk(1)
           WC = WC + sxm(2) - sxm(1)
           goto 20

         case(3)
           xc = xc + dx2
           UC = UC + (T2 * sxk(3)) - sxk(2)
           WC = WC + (T2 * sxm(3)) - sxm(2)
           bgc = pot(np1)
           goto 20

         case default
           W(np1) = W(N)+(sxm(1)+sxm(4)+T2*(sxm(2)+sxm(3)))*TH
           U(np1) = U(N)+(sxk(1)+sxk(4)+T2*(sxk(2)+sxk(3)))*TH
           UP(np1) = (VC+T)*W(np1)-XK*U(np1)
           WP(np1) = XK*W(np1)-TC*U(np1)
           X = X + dx
           N = np1
           if (N < 6) goto 25

       end select

       ! END OF STARTING INTEGRATION.  BEGIN MILNE PROCEDURE.
       T = exp(X)
   26  T = T * TDX
       np1 = N + 1
       nm1 = N - 1
       nm2 = N - 2
       nm3 = N - 3
       nm4 = N - 4
       nm5 = N - 5
       TC = ((E + vcz) * T) + pot(np1)
       VC = cin * TC
       unp = U(nm5) + X20 * (T11 * (UP(N) + UP(nm4)) + T26 * UP(nm2) - T14 * (UP(nm1) + UP(nm3)))
       wnp = W(nm5) + X20 * (T11 * (WP(N) + WP(nm4)) + T26 * WP(nm2) - T14 * (WP(nm1) + WP(nm3)))
       nit = 0
  33   UP(np1) = (VC+T)*wnp-XK*unp
       WP(np1) = XK*wnp-TC*unp
       unp2 = U(nm3) + (T7 * (UP(np1) + UP(nm3)) + T32*(UP(nm2) + UP(N)) + T12 * UP(nm1)) * XMFT
       wnp2 = W(nm3) + (T7 * (WP(np1) + WP(nm3)) + T32*(WP(nm2) + WP(N)) + T12 * WP(nm1)) * XMFT

       ! COMPARE PREDICTOR WITH CORRECTOR
       if (abs(TEST*(unp2 - unp)) > abs(unp2)) goto 31
       if (abs(TEST*(wnp2 - wnp)) <= abs(wnp2)) goto 32
  31   if (nit < 5) goto 81
       goto 32

  81   nit = nit + 1

       wnp = wnp2
       unp = unp2
       goto 33

  32   W(np1) = wnp2

       U(np1) = unp2
       N = np1

       if (N < jri) goto 26

       ! END OF MILNE PROCEDURE
       if (abs(U(jri)) <= zilch) U(jri) = sign(zilch, U(jri))

       P = (T + VC) / T

       wnp = P * W(jri) / U(jri)
       unp = wnp - (kappa + 1) / T
       dlgkap = (12.5663706) * T * T * unp

       return

      end function dlgkap

!***********************************************************************
      subroutine sbfit(T,E,L,R,jfs)

       double precision E,R,T
       real L
       real jfs,kappa

       SE = SNGL(E)
       SR = SNGL(R)
       ST = SNGL(T)
       kappa = sqrt(SE)
       X = kappa * SR
       bj1 = sin(X) / X
       bn1 = -cos(X) / X
       bj2 = (bj1 / X) + bn1
       bn2 = (bn1 / X) - bj1

       if (L > 0) then
         DL = ST / (SR * SR)
       else
         ls = 1
       endif

       ls = 1
       ls = ls + 1
       bjT = ((2 * ls) - 1) * (bj2 / X) - bj1
       bnT = ((2 * ls) - 1) * (bn2 / X) - bn1
       bj1 = bj2
       bj2 = bjT
       bn1 = bn2
       bn2 = bnT

       if (L+1-ls > 0) then
         ls = ls + 1
       else
         DL = ST / (SR * SR)
       endif

       DL = ST / (SR * SR)
       DL = DL - (L / SR)
       AN = (DL * bj1) + (kappa * bj2)
       AD = (DL * bn1) + (kappa * bn2)
       jfs = 3.141592654 / 2.0

       if (abs(AD)-1.0E-8 > 0) jfs = atan(AN / AD)

       return

      end subroutine sbfit
