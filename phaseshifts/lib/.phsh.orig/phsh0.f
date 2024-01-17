C  file PSPROG.AB3  Feb. 5, 1995
C  containing programs PHSH0.FOR and PHSH1.FOR
C
C---------------------------------------------------------------------
C  program PHSH0.FOR
C---------------------------------------------------------------------
c
c  there are nr grid points, and distances are in bohr radii...
c
c  r(i)=rmin*(rmax/rmin)**(dfloat(i)/dfloat(nr)) , i=1,2,3,...nr-1,nr
c
c
c
c  the orbitals are store in phe(), first index goes 1...nr, the
c  second index is the orbital index (i...nel)
c
c  look at the atomic files after printing this out to see everything...
c
c  suffice it to say, that the charge density at radius r(i)
c  in units of electrons per cubic bohr radius is given by
c
c  sum of j=1...nel,
c  occ(j)*phe(i,j)*phe(i,j)/(4.d0*3.14159265....*r(i)*r(i))...
c
c  think of the phe functions as plotting the radial wave-functions
c  as a function of radius...on our logarithmic mesh...
c
c  final note:
c
c  the Dirac equation is solved for the orbitals, whereas their density
c  is treated by setting phe to the square root of Dirac's F*F+G*G
c  times the sign of G...
c
c  so we are doing Dirac-Fock, except that we are not treating exchange
c  exactly, in terms of working with major and minor components of the
c  orbitals, and the phe's give the CORRECT CHARGE DENSITY...
c
c  the above approximation ought to be very small for valence states,
c  so you need not worry about it...
c
c  the Breit interaction has been neglected altogether...it should not
c  have a huge effect on the charge density you are concerned with...
C
C author: Eric Shirley
C
C
C  modified: 26/01/2011 LD - added $1 argument instead of atorb
      program phsh0
      implicit real*8 (a-h,o-z)
      parameter (iorbs=33,iside=600)
      parameter (io2=iorbs*(iorbs+1)/2)
      parameter (ijive=io2*(io2+1)/2)
      parameter (lmax=4,ihmax=20,nrmax=4000,ntmax=10,npmax=60)
      dimension r(nrmax),dr(nrmax),r2(nrmax)
      dimension no(iorbs),nl(iorbs),xnj(iorbs)
      dimension ev(iorbs),occ(iorbs),is(iorbs)
      dimension ek(iorbs),phe(nrmax,iorbs),orb(nrmax,iorbs)
      dimension njrc(4),vi(nrmax,7),rho(nrmax)
      dimension v(nrmax),q0(nrmax),xm1(nrmax),xm2(nrmax)
      dimension w(33,33),wi(33,33),rhs(33),co(33)
      dimension xint(0:12),vav(11),rint(0:12)
      dimension pin(0:11),sig(0:11),vctab(nrmax,0:3)
      character*1 ichar
      character*11 jive
      character*60 jive2
      character*70 jive3
      integer I

c       modified by Liam Deacon
      CHARACTER(len=255)     ::      ARG, input, STRING(3)
      input = "atorb"

      I = 1
      DO
        CALL GETARG(I, ARG)
        IF ((ARG.EQ.'-i').OR.(ARG.EQ.'--input')) THEN
          I = I+1
          CALL GETARG(I, input)
        ENDIF
        IF ((ARG.EQ.'-h').OR.(ARG.EQ.'--help')) THEN
          WRITE(*,*) 'phsh0 usage:-'
          WRITE(*,*)
          STRING(1) = 'phsh0 -i <file>'
          WRITE(*,*) TRIM(STRING(1))
          WRITE(*,*) ''
          WRITE(*,*) 'where:-'
          STRING(1) = '-i or --input <file> specifies '
          STRING(2) = 'the input file path (default "atorb")'
          WRITE(*,*) TRIM(STRING(1))//TRIM(STRING(2))
          STRING(1) = '-h or --help print help and exit'
          WRITE(*,*) TRIM(STRING(1))
          STOP
        ENDIF
        IF (I.GE.IARGC()) THEN
          EXIT
        ENDIF
        I = I+1
      END DO

      call hartfock(input)

      end

c-----------------------------------------------------------------------
      subroutine hartfock(input_file)
      implicit real*8 (a-h,o-z)
      character(len=255), intent(in) :: input_file
      parameter (iorbs=33,iside=600)
      parameter (io2=iorbs*(iorbs+1)/2)
      parameter (ijive=io2*(io2+1)/2)
      parameter (lmax=4,ihmax=20,nrmax=4000,ntmax=10,npmax=60)
      dimension r(nrmax),dr(nrmax),r2(nrmax)
      dimension no(iorbs),nl(iorbs),xnj(iorbs)
      dimension ev(iorbs),occ(iorbs),is(iorbs)
      dimension ek(iorbs),phe(nrmax,iorbs),orb(nrmax,iorbs)
      dimension njrc(4),vi(nrmax,7),rho(nrmax)
      dimension v(nrmax),q0(nrmax),xm1(nrmax),xm2(nrmax)
      dimension w(33,33),wi(33,33),rhs(33),co(33)
      dimension xint(0:12),vav(11),rint(0:12)
      dimension pin(0:11),sig(0:11),vctab(nrmax,0:3)
      character*1 ichar
      character*11 jive
      character*60 jive2
      character*70 jive3
       !implicit real*8 (a-h,o-z)

       open(unit=5, file=input_file, status='old')
       rel=0.d0
 10     read (5,20) ichar
 20     format (1a1)
       if (ichar .eq. 'd') then
         write (6,*) 'PLEASE ENTER RELATIVITY FACTOR.  (0=NR, 1=REL.)'
         read (5,*) rel
       endif
       if (ichar .eq. 'x') then
         write (6,*) 'PLEASE ENTER EXCHANGE CORRELATION (ALPHA).'
         write (6,*) '0=HARTREE-FOCK, POSITIVE=LDA, NEGATIVE=XALPHA.'
         read  (5,*) alfa
       endif
       if (ichar .eq. 'a') then
         call abinitio(etot,nst,rel,alfa,nr,r,
     1    dr,r2,dl,phe,njrc,vi,zorig,xntot,nel,
     1    no,nl,xnj,ev,occ,is,ek,orb,iuflag)
       endif
       if (ichar .eq. 'i') call initiali(zorig,nr,rmin,rmax,
     1    r,dr,r2,dl,njrc,xntot,nel)
       if (ichar .eq. 'q') stop
       if (ichar .eq. 'w') then
         ixflag=1
         iu=-1
         ir=0
         call hfdisk(iu,ir,etot,nst,rel,nr,rmin,rmax,r,rho,
     1      zorig,xntot,ixflag,nel,
     1      no,nl,xnj,is,ev,ek,occ,njrc,vi,phe,orb)
       endif
       if (ichar.eq.'r') then
         iu=-1
         ir=1
         call hfdisk(iu,ir,etot,nst,rel,nr,rmin,rmax,r,rho,
     1      zorig,xntot,ixflag,nel,
     1      no,nl,xnj,is,ev,ek,occ,njrc,vi,phe,orb)
         call setgrid(nr,rmin,rmax,r,dr,r2,dl)
       endif
       if (ichar.eq.'u') then
         write (6,*) 'PLEASE ENTER IUFLAG. (0=U, 1=SU, 2=R).'
         read (5,*) iuflag
       endif
       if (ichar.eq.'c') then
         write (6,*) 'PLEASE ENTER ALPHA,RS,RP,RD.'
         read (5,*) corpol,rs,rp,rd
         do 100 k=1,nr
         fs=(1.d0-dexp(-(r(k)/rs)**2.d0))**2.d0
         fp=(1.d0-dexp(-(r(k)/rp)**2.d0))**2.d0
         fd=(1.d0-dexp(-(r(k)/rd)**2.d0))**2.d0
         vctab(k,0)=-corpol/2.d0*fs*fs/r(k)**4.d0
         vctab(k,1)=-corpol/2.d0*fp*fp/r(k)**4.d0
         vctab(k,2)=-corpol/2.d0*fd*fd/r(k)**4.d0
 100      continue
       endif
       if (ichar.eq.'f') then
         write (6,*) 'PLEASE ENTER IUNIT,CORPOL'
         read  (5,*) iunit,corpol
         write (6,*) 'PLEASE ENTER ILEV,INUM,EOLD'
         read  (5,*) ilev,inum,eold
         xl=nl(ilev)
         if (inum.eq.1) then
           read (5,*) eav
         else
           read (5,*) e1,e2
           eav=(e1*xl+e2*(xl+1.d0))
     1         /(   xl+    xl+1.d0 )
         endif
         if (eav.lt.0.d0) eav=eold+eav
         if (iunit.eq.2) eav=eav/2.d0
         if (iunit.eq.3) eav=eav/27.2116d0
         if (iunit.eq.4) eav=eav*0.000123985d0/27.2116d0
         sd=dabs(dabs(eav)-dabs(ev(ilev)))
         rl= 0.d0
         rh=10.d0
         sl= 0.d0
         sh= 0.d0
 300      if (sl*sh.le.0.00000001d0) rc=rl+(rh-rl)/2.d0
         if (sl*sh.gt.0.00000001d0) rc=rl+(rh-rl)*(sd-sl)/(sh-sl)
         sc=0.d0
         do 320 i=1,nr
         f=(1.d0-dexp(-(r(i)/rc)**2.d0))**2.d0
         vcpp=corpol/(2.d0*r(i)**4.d0)*f*f
         sc=sc+dr(i)*phe(i,ilev)*phe(i,ilev)*vcpp
 320      continue
         if (sc.gt.sd) rl=rc
         if (sc.gt.sd) sl=sc
         if (sc.lt.sd) rh=rc
         if (sc.lt.sd) sh=sc
         write (6,*) rc,sc
         if (dabs(sc-sd).gt.0.000001d0) goto 300
       endif
       if (ichar.eq.'p') then
         call pseudo(etot,nst,rel,alfa,nr,rmin,rmax,r,dr,r2,dl,
     1                phe,orb,njrc,vi,zorig,xntot,nel,
     1                no,nl,xnj,ev,occ,is,ek,iuflag,vctab)
       endif
       if (ichar.eq.'g') then
         read(5,*) iu
         read(5,2202) jive
         read(5,2212) jive2
         read(5,2222) jive3
 2202     format(1x,1a11)
 2212     format(1x,1a60)
 2222     format(1x,1a70)
         zizv=dabs(r(nr-1)*vi(nr-1,1))
         write (iu,2202) jive
         write (iu,2212) jive2
         write (iu,2222) jive3
         write (iu,*) 3,nr,zizv
         write (iu,*) (r(i),i=1,nr)
         write (iu,*) 0,(vi(k,1),k=1,nr)
         write (iu,*) 1,(vi(k,3),k=1,nr)
         write (iu,*) 2,(vi(k,5),k=1,nr)
         write (iu,*) (0.d0,k=1,nr)
         do 500 j=1,nr
         rh=0.d0
         do 480 k=1,nel
         rh=rh+phe(j,k)*phe(j,k)*occ(k)
 480      continue
         write (iu,*) rh
 500      continue
       endif
       if (ichar.eq.'v') then
         rold=0.d0
         open(unit=10)
         open(unit=11)
         open(unit=12)
         do 600 k=1,nr
c         if ((r(k).lt.10.d0).and.((r(k)-rold).gt.0.05d0)) then
           write (10,*) r(k),vi(k,1)*r(k)
           write (11,*) r(k),vi(k,3)*r(k)
           write (12,*) r(k),vi(k,5)*r(k)
           rold=r(k)
c         endif
 600      continue
         close(unit=10)
         close(unit=11)
         close(unit=12)
       endif
       if (ichar.eq.'V') call fourier(nr,r,dr,r2,vi)
       goto 10
       return
       end
c-----------------------------------------------------------------------
      subroutine abinitio(etot,nst,rel,alfa,nr,r,dr,r2,dl,
     1    phe,njrc,vi,zorig,xntot,nel,no,nl,xnj,
     1    ev,occ,is,ek,orb,iuflag)
      implicit real*8 (a-h,o-z)
      parameter (iorbs=33,iside=600)
      parameter (io2=iorbs*(iorbs+1)/2)
      parameter (ijive=io2*(io2+1)/2)
      parameter (lmax=4,ihmax=20,nrmax=4000,ntmax=10,npmax=60)
      dimension r(nrmax),dr(nrmax),r2(nrmax),v(nrmax)
      dimension no(iorbs),nl(iorbs),nm(iorbs),xnj(iorbs)
      dimension ev(iorbs),occ(iorbs),is(iorbs),ek(iorbs)
      dimension phe(nrmax,iorbs),njrc(4),vi(nrmax,7)
      dimension orb(nrmax,iorbs),rpower(nrmax,0:15)
c  this will be good for going up to and including l=3...
       do 10 i=0,7
       xi=i
       do 10 k=1,nr
       rpower(k,i)=r(k)**xi
 10     continue

c  read in nfc, nel.  refer to the documentation for their meanings.

 168    write (6,*) 'PLEASE ENTER NFC, NEL, RATIO, ETOL, XNUM'
       read (5,*) nfc,nel,ratio,etol,xnum

c  for all of the electrons, read in the quantum numbers.
c  get the total Hartree-active charge.  initialize eigenvalues.

       xntot=0.d0

       write (6,*) 'PLEASE ENTER N L M J S OCC.'

       do 100 i=nfc+1,nel
       read (5,*) no(i),nl(i),nm(i),xnj(i),is(i),occ(i)
       ev(i)=0.d0
       xntot=xntot+occ(i)
       do 100 j=1,nr
       phe(j,i)=0.d0
       orb(j,i)=0.d0
 100    continue

c  initialize the parameters for self-consistency loop.
c  ratio is the mixture of old and new field mixing.

 110    call atsolve(etot,nst,rel,alfa,
     1    eerror,nfc,nr,r,dr,r2,dl,phe,
     1    njrc,vi,zorig,xntot,nel,no,nl,nm,xnj,ev,occ,is,ek,
     1    ratio,orb,rpower,xnum,etot2,iuflag)

       eerror=eerror*(1.d0-ratio)/ratio
       write (6,112) eerror,etot
112     format (1x,3f14.6)
       if (eerror.gt.etol) goto 110

c  write out information about the atom.

 120    do 130 i=1,nel
       nj=xnj(i)+xnj(i)
       write (6,122) no(i),nl(i),nm(i),nj,'/2',is(i),occ(i),ev(i)
 122    format(1x,2i4,i2,i4,a2,i4,f10.4,f18.6)
 130    continue

       write (6,132) 'TOTAL ENERGY =  ',etot,etot*27.2116d0
 132    format (1x,a16,2f14.6)

       return

       end

c-----------------------------------------------------------------------
       subroutine atsolve(etot,nst,rel,alfa,eerror,nfc,
     1    nr,r,dr,r2,dl,
     1    phe,njrc,vi,zorig,xntot,nel,no,nl,nm,xnj,ev,occ,is,ek,
     1    ratio,orb,rpower,xnum,etot2,iuflag)
       implicit real*8 (a-h,o-z)
       parameter (iorbs=33,iside=600)
       parameter (io2=iorbs*(iorbs+1)/2)
       parameter (ijive=io2*(io2+1)/2)
       parameter (lmax=4,ihmax=20,nrmax=4000,ntmax=10,npmax=60)
       dimension r(nrmax),dr(nrmax),r2(nrmax),v(nrmax)
       dimension no(iorbs),nl(iorbs),nm(iorbs),xnj(iorbs)
       dimension ek(iorbs),ev(iorbs),occ(iorbs),is(iorbs)
       dimension phe(nrmax,iorbs),njrc(4),vi(nrmax,7)
       dimension q0(nrmax),xm1(nrmax),xm2(nrmax),orb(nrmax,iorbs)
       dimension rpower(nrmax,0:15)

c  initialize eerror, the biggest change in an eigenvalue, and etot.

       eerror=0.d0
       etot=0.d0

c  run through all the orbitals.  calculate those not in the core.

       do 102 i=1,nel

       if (i.gt.nfc) then

         idoflag=1
         call setqmm(i,orb,nl(i),is(i),idoflag,v,zeff,zorig,rel,
     1      nr,r,r2,dl,q0,xm1,xm2,njrc,vi)

         xkappa=-1.d0
         if (dabs(xnj(i)).gt.dfloat(nl(i))+0.25d0) xkappa=-nl(i)-1
         if (dabs(xnj(i)).lt.dfloat(nl(i))-0.25d0) xkappa= nl(i)

         call elsolve(i,occ(i),no(i),nl(i),xkappa,xnj(i),zorig,zeff,
     1      evi,phe(1,i),v,q0,xm1,xm2,nr,r,dr,r2,dl,rel)
         if (dabs(ev(i)-evi).gt.eerror) eerror=dabs(ev(i)-evi)
         ev(i)=evi

         ekk=0.d0
         ll=2
         do 100 j=nr,1,-1
         dq=phe(j,i)*phe(j,i)
         ekk=ekk+(evi-orb(j,i))*dr(j)*dq*dfloat(ll)/3.d0
         ll=6-ll
 100      continue
         ek(i)=ekk

       endif

c  add the kinetic to total, including the frozen core kinetic energy.

       etot=etot+ek(i)*occ(i)
 102    continue
       call getpot(etot,nst,rel,alfa,dl,nr,dr,r,r2,xntot,
     1    phe,ratio,orb,occ,is,nel,nl,nm,no,xnj,rpower,xnum,
     1    etot2,iuflag)
       return
       end
c----------------------------------------------------------------------
       subroutine getpot(etot,nst,rel,alfa,dl,nr,dr,r,r2,
     1    xntot,phe,ratio,orb,occ,is,
     1    nel,nl,nm,no,xnj,rpower,xnum,etot2,iuflag)
       implicit real*8 (a-h,o-z)
       parameter (iorbs=33,iside=600)
       parameter (io2=iorbs*(iorbs+1)/2)
       parameter (ijive=io2*(io2+1)/2)
       parameter (lmax=4,ihmax=20,nrmax=4000,ntmax=10,npmax=60)
       dimension dr(nrmax),r(nrmax),r2(nrmax)
       dimension phe(nrmax,iorbs),occ(iorbs)
       dimension is(iorbs),orb(nrmax,iorbs),nl(iorbs)
       dimension nm(iorbs),xnj(iorbs),no(iorbs)
       dimension rpower(nrmax,0:15)
       dimension xq1(nrmax),xq2(nrmax),xq0(nrmax)
       dimension cg(0:6,0:6,0:12,-6:6,-6:6),pin(0:8,0:8,0:16)
       dimension xqj0(nrmax),xqj1(nrmax)
       dimension xqj2(nrmax),xqi0(nrmax)
       dimension xqi1(nrmax),xqi2(nrmax),rsp(2)

       call clebschgordan(nel,nl,cg)
       call getillls(pin)

       ratio1=1.d0-ratio
       do 2100 i=1,nel
       do 2100 k=1,nr
       orb(k,i)=ratio1*orb(k,i)
 2100   continue

       do 2990 i=1,nel

       li=nl(i)
       mi=nm(i)

       jstart=i+1
       if ((xnj(i).lt.0.d0).or.
     1      (occ(i).gt.1.d0).or.
     1      (dabs(alfa).gt.0.001d0)) jstart=i
       do 2990 j=jstart,nel

       if ((occ(i).eq.0.d0).and.(occ(j).eq.0.d0)) goto 2990

       lj=nl(j)
       mj=nm(j)

c  direct coulomb

       lmx=2*li
       if (li.gt.lj) lmx=2*lj

c  l=0 is monopole or spherical term for direct coulomb.  Therefore,
c  when we have occ(i) or occ(j) greater than one, set lmx=0.

       if ((occ(i).gt.1.d0).or.(occ(j).gt.1.d0).or.
     1      (xnj(i).lt.0.d0).or.(xnj(j).lt.0.d0)) lmx=0

       do 2550 la=lmx,0,-2
       lap=la+1
       coeff=dfloat((li+li+1)*(lj+lj+1))/dfloat((la+la+1))**2.d0*
     1    cg(li,li,la,mi,-mi)*cg(lj,lj,la,mj,-mj)*
     1    cg(li,li,la,0 , 0 )*cg(lj,lj,la,0 , 0 )
       if (mi+mj.ne.2*((mi+mj)/2)) coeff=-coeff
       if (i.eq.j) coeff=coeff/2.d0
       coeffi=occ(i)*coeff
       coeffj=occ(j)*coeff
       ri=ratio*coeffi
       rj=ratio*coeffj
       rc=coeff*occ(i)*occ(j)

       xouti=0.d0
       xoutj=0.d0
       do 2500 k=1,nr
       xqi0(k)=dr(k)*phe(k,i)*phe(k,i)/2.d0
       xqi1(k)=xqi0(k)*rpower(k,la)
       if (rpower(k,lap).ne.0.d0) then
         xqi2(k)=xqi0(k)/rpower(k,lap)
       else
         xqi2(k)=0.d0
       endif
       xouti=xouti+xqi2(k)
       xqj0(k)=dr(k)*phe(k,j)*phe(k,j)/2.d0
       xqj1(k)=xqj0(k)*rpower(k,la)
       if (rpower(k,lap).ne.0.d0) then
         xqj2(k)=xqj0(k)/rpower(k,lap)
       else
         xqj2(k)=0.d0
       endif
       xoutj=xoutj+xqj2(k)
 2500   continue

       xinti=xqi1(1)
       xintj=xqj1(1)
       xouti=2.d0*xouti-xqi2(1)
       xoutj=2.d0*xoutj-xqj2(1)

       do 2550 k=2,nr

       xinti=xinti+xqi1(k)+xqi1(k-1)
       xouti=xouti-xqi2(k)-xqi2(k-1)
       vali=xouti*rpower(k,la)
       if (rpower(k,lap).ne.0.d0) vali=vali+xinti/rpower(k,lap)
       orb(k,j)=orb(k,j)+ri*vali

       xintj=xintj+xqj1(k)+xqj1(k-1)
       xoutj=xoutj-xqj2(k)-xqj2(k-1)
       valj=xoutj*rpower(k,la)
       if (rpower(k,lap).ne.0.d0) valj=valj+xintj/rpower(k,lap)
       orb(k,i)=orb(k,i)+rj*valj

       etot=etot+rc*(xqi0(k)*valj+xqj0(k)*vali)

 2550   continue

       if ((is(i).ne.is(j)).and.
     1      (occ(i).le.1.d0).and.
     1      (occ(j).le.1.d0).and.
     1      (xnj(i).ge.0.d0).and.
     1      (xnj(j).ge.0.d0)     ) goto 2990
       if (dabs(alfa).ge.0.001d0) goto 2990

c  exchange interaction

       lmx=li+lj
       lmin=iabs(mi-mj)
       if ((occ(i).gt.1.d0).or.(occ(j).gt.1.d0).or.
     1      (xnj(i).lt.0.d0).or.(xnj(j).lt.0.d0)) lmin=0
       do 2980 la=lmx,lmin,-2
       lap=la+1

       coeff=dfloat((li+li+1)*(lj+lj+1))/dfloat((la+la+1))**2.d0*
     1    (cg(li,lj,la,-mi,mj)*cg(li,lj,la,0,0))**2.d0
       if ((occ(i).gt.1.d0).or.(occ(j).gt.1.d0).or.
     1      (xnj(i).lt.0.d0).or.(xnj(j).lt.0.d0))
     1     coeff=pin(li,lj,la)/4.d0
       if (i.eq.j) coeff=coeff/2.d0
       coeffi=occ(i)*coeff
       coeffj=occ(j)*coeff
       ri=ratio*coeffi
       rj=ratio*coeffj
       rc=coeff*occ(i)*occ(j)
       xnum2=xnum*xnum

       xout=0.d0
       do 2600 k=1,nr
       xq0(k)=dr(k)*phe(k,i)*phe(k,j)/2.d0
       xq1(k)=xq0(k)*rpower(k,la)
       if (rpower(k,lap).ne.0.d0) then
         xq2(k)=xq0(k)/rpower(k,lap)
       else
         xq2(k)=0.d0
       endif
       xout=xout+xq2(k)
 2600   continue
       xint=xq1(1)
       xout=2.d0*xout-xq2(1)
       do 2610 k=2,nr
       xint=xint+xq1(k)+xq1(k-1)
       xout=xout-xq2(k)-xq2(k-1)
       if (xq0(k).ne.0.d0) then
         val=xout*rpower(k,la)
         if (rpower(k,lap).ne.0.d0) val=val+xint/rpower(k,lap)
         etot=etot-2.d0*xq0(k)*rc*val
         xx=phe(k,j)/phe(k,i)
         if (dabs(xx)/xnum.gt.1.d0) then
           orb(k,i)=orb(k,i)-rj*xnum2/xx*val
         else
           orb(k,i)=orb(k,i)-rj*xx*val
         endif
         xx=phe(k,i)/phe(k,j)
         if (dabs(xx)/xnum.gt.1.d0) then
           orb(k,j)=orb(k,j)-ri*xnum2/xx*val
         else
           orb(k,j)=orb(k,j)-ri*xx*val
         endif
       endif
 2610   continue

 2980   continue

 2990   continue
c
c  here we compute the charge density, if needed, for treating
c  exchange/correlation in a local fashion...
c
       if (dabs(alfa).ge.0.001d0) then
         if (alfa.gt.0.d0) then
           fx=1.0d0
           fc=1.0d0
         else
           fx=1.5d0*dabs(alfa)
           fc=0.0d0
         endif
c
c  note: we don't deal with spin-polarization in local exchange
c  picture, since local exchange is totally wrong for such
c  effects, anyway.  local exchange pretends charge density
c  is paramagnetic.  also, local exchange treats everything
c  as spherical.
c
         fourpi=16.d0*datan(1.d0)
         do 5000 i=1,nr
         xn=0.d0
         do 3200 j=1,nel
         xn=xn+phe(i,j)*phe(i,j)*occ(j)
 3200     continue
         xn1=xn/2.d0
         xn2=xn/2.d0
         nst=2
         call exchcorr(nst,rel,r2(i),xn1,xn2,ex,ec,ux1,ux2,uc1,uc2)
         exc=fx*ex +fc*ec
         uxc=fx*ux1+fc*uc1
         etot=etot+dr(i)*xn*exc
         do 3300 j=1,nel
         orb(i,j)=orb(i,j)+uxc*ratio
 3300     continue
 5000     continue
       endif
c
       do 9000 i=1,nr
       if (iuflag.ne.0) then
         jj=1
 8960     ii=jj
 8965     if (ii.eq.nel) goto 8970
         icond=0
         if ((no(jj).eq.no(ii+1)).and.(nl(jj).eq.nl(ii+1))
     1      .and.(iuflag.eq.2)) icond=1
         if ((no(jj).eq.no(ii+1)).and.(nl(jj).eq.nl(ii+1))
     1      .and.(is(jj).eq.is(ii+1)).and.(iuflag.eq.1)) icond=1
         if (icond.eq.1) then
           ii=ii+1
           goto 8965
         endif
 8970     orba=0.d0
         div=0.d0
         do 8980 k=jj,ii
         div=div+occ(k)
         orba=orba+orb(i,k)*occ(k)
 8980     continue
         if (div.ne.0.d0) then
           orba=orba/div
           do 8990 k=jj,ii
           orb(i,k)=orba
 8990       continue
         endif
         if (ii.ne.nel) then
           jj=ii+1
           goto 8960
         endif
       endif
 9000   continue
       return
       end
c-----------------------------------------------------------------------
       subroutine elsolve(i,occ,n,l,xkappa,xj,zorig,zeff,e,phi,v,
     1    q0,xm1,xm2,nr,r,dr,r2,dl,rel)
       implicit real*8 (a-h,o-z)
       parameter (iorbs=33,iside=600)
       parameter (io2=iorbs*(iorbs+1)/2)
       parameter (ijive=io2*(io2+1)/2)
       parameter (lmax=4,ihmax=20,nrmax=4000,ntmax=10,npmax=60)
       dimension phi(nrmax),phi2(nrmax)
       dimension v(nrmax),q0(nrmax),xm1(nrmax),xm2(nrmax)
       dimension r(nrmax),dr(nrmax),r2(nrmax)
       el=-zorig*zorig/dfloat(n*n)
       eh=0.d0
       etol=0.0000000001d0
 155    e=(el+eh)/2.d0
       istop=0
       call integ(e,l,xkappa,n,nn,istop,ief,x0,phi,zeff,v,q0,xm1,
     1    xm2,nr,r,dr,r2,dl,rel)
       if (nn.lt.n-l-1) ief=-1
 200    if (ief.ne.1) then
         el=e
         if (el.gt.-0.001d0) then
           write (6,*) 'MIXING TOO STRONG FOR LEVEL : ',i
           stop
         endif
       endif
       if (ief.ne.-1) eh=e
       if (eh-el.gt.etol) goto 155
       if (dabs(dabs(xj)-dabs(dfloat(l))).gt.0.25d0)
     1    call augment(e,l,xj,phi,v,nr,r,dl)
       aa=0.d0
       do 6130 j=1,nr
       aa=aa+phi(j)*phi(j)*dr(j)
 6130   continue
       xnorm=dsqrt(aa)
       do 6140 j=1,nr
       phi(j)=phi(j)/xnorm
 6140   continue
       return
       end
c--------------------------------------------------------------------------
       subroutine augment(e,l,xj,phi,v,nr,r,dl)
       implicit real*8 (a-h,o-z)
       parameter (iorbs=33,iside=600)
       parameter (io2=iorbs*(iorbs+1)/2)
       parameter (ijive=io2*(io2+1)/2)
       parameter (lmax=4,ihmax=20,nrmax=4000,ntmax=10,npmax=60)
       dimension phi(nrmax),phi2(nrmax),v(nrmax),r(nrmax)
       c=137.038d0
       cc=c*c
       c2=cc+cc
       xkappa=-1
       if (dabs(xj).gt.dfloat(l)+0.25d0) xkappa=-l-1
       if (dabs(xj).lt.dfloat(l)-0.25d0) xkappa= l
       do 6110 j=4,nr-3
       if (phi(j).ne.0.d0) then
         g0=phi(j)
         ga=(phi(j+1)-phi(j-1))
         gb=(phi(j+2)-phi(j-2))/2.d0
         gc=(phi(j+3)-phi(j-3))/3.d0
         gg=((1.5d0*ga-0.6d0*gb+0.1d0*gc)/(2.d0*dl)+xkappa*g0)/r(j)
         f0=c*gg/(e-v(j)+c2)
         phi2(j)=dsqrt(g0*g0+f0*f0)
         if (g0.lt.0.d0) phi2(j)=-phi2(j)
       else
         phi2(j)=phi(j)
       endif
 6110   continue
       do 6115 j=1,3
       phi2(j)=phi(j)*phi(4)/phi2(4)
 6115   continue
       do 6120 j=1,nr
       phi(j)=phi2(j)
 6120   continue
       return
       end
c-----------------------------------------------------------------------
       subroutine setqmm(i,orb,l,ns,idoflag,v,zeff,zorig,rel,
     1    nr,r,r2,dl,q0,xm1,xm2,njrc,vi)
       implicit real*8 (a-h,o-z)
       parameter (iorbs=33,iside=600)
       parameter (io2=iorbs*(iorbs+1)/2)
       parameter (ijive=io2*(io2+1)/2)
       parameter (lmax=4,ihmax=20,nrmax=4000,ntmax=10,npmax=60)
       dimension v(nrmax),r(nrmax),r2(nrmax),orb(nrmax,iorbs)
       dimension q0(nrmax),xm1(nrmax),xm2(nrmax),njrc(4),vi(nrmax,7)

       c=137.038d0
       alpha=rel/c
       aa=alpha*alpha
       a2=aa/2.d0

       lp=l+1
       lpx=lp
       if (lp.gt.4) lpx=4
       lp2=l+l+1
       if (lp2.gt.7) lp2=7
       zeff=zorig
       if (njrc(lpx).gt.0) zeff=0.d0
       zaa=zeff*aa
       za2=zeff*a2

       if (idoflag.ne.0) then
         if (njrc(lpx).eq.0) then
           if (idoflag.eq.1) then
             do 55 j=1,nr
             v(j)=-zeff/r(j)+orb(j,i)
 55           continue
           endif
           do 65 j=2,nr-1
           dvdl=(orb(j+1,i)-orb(j-1,i))/(2.d0*dl)
           ddvdrr=((orb(j+1,i)
     1         +orb(j-1,i)-2.d0*orb(j,i) )/(dl*dl)-dvdl)/r2(j)
           xm1(j)=-a2*dvdl/r(j)-za2/r2(j)
           xm2(j)=-a2*ddvdrr+zaa/r2(j)/r(j)
 65         continue
           xm1(nr)=xm1(nr-1)
           xm2(nr)=xm2(nr-1)
           xm1(1)=xm1(2)+za2/r2(2)-za2/r2(1)
           xm2(1)=xm2(2)-zaa/r2(2)/r(2)+zaa/r2(1)/r(1)
         else
           if (idoflag.eq.1) then
             do 75 j=1,nr
             v(j)=vi(j,lp2)+orb(j,i)
 75           continue
           endif
           do 85 j=2,nr-1
           dvdl=(v(j+1)-v(j-1))/(2.d0*dl)
           ddvdrr=((v(j+1)+v(j-1)-2.d0*v(j))/(dl*dl)-dvdl)/r2(j)
           xm1(j)=-a2*dvdl/r(j)
           xm2(j)=-a2*ddvdrr
 85         continue
           xm1(nr)=xm1(nr-1)
           xm2(nr)=xm2(nr-1)
           xm1(1)=xm1(2)
           xm2(1)=xm2(2)
         endif
       endif

c  figure the (Desclaux-Numerov) effective potential.

       xlb=(dfloat(l)+0.5d0)**2.d0/2.d0
       do 45 j=1,nr
       vj=v(j)
       q0(j)=vj*(1.d0-a2*vj)+xlb/r2(j)
 45     continue

       return

       end
c----------------------------------------------------------------------
       subroutine initiali(zorig,nr,rmin,rmax,r,dr,r2,dl,njrc,
     1    xntot,nel)
       implicit real*8 (a-h,o-z)
       parameter (iorbs=33,iside=600)
       parameter (io2=iorbs*(iorbs+1)/2)
       parameter (ijive=io2*(io2+1)/2)
       parameter (lmax=4,ihmax=20,nrmax=4000,ntmax=10,npmax=60)
       dimension r(nrmax),dr(nrmax),r2(nrmax),njrc(4)
       write (6,*) 'ENTER Z, NR'
       read (5,*) zorig,nr
       rmin=0.0001d0/zorig
       rmax=800.d0/dsqrt(zorig)
       call setgrid(nr,rmin,rmax,r,dr,r2,dl)
       do 5 j=1,4
       njrc(j)=0
 5      continue
       xntot=0.d0
       nel=0
       return
       end
c---------------------------------------------------------------------------
       subroutine setgrid(nr,rmin,rmax,r,dr,r2,dl)
       implicit real*8 (a-h,o-z)
       parameter (iorbs=33,iside=600)
       parameter (io2=iorbs*(iorbs+1)/2)
       parameter (ijive=io2*(io2+1)/2)
       parameter (lmax=4,ihmax=20,nrmax=4000,ntmax=10,npmax=60)
       dimension r(nrmax),dr(nrmax),r2(nrmax)
       ratio=rmax/rmin
       dl=dlog(ratio)/dfloat(nr)
       xratio=dexp(dl)
       xr1=dsqrt(xratio)-dsqrt(1.d0/xratio)
       do 2010 i=1,nr
       r(i)=rmin*xratio**dfloat(i)
       dr(i)=r(i)*xr1
       r2(i)=r(i)*r(i)
 2010   continue
       return
       end
c-----------------------------------------------------------------------------
       subroutine integ(e,l,xkappa,n,nn,istop,ief,x0,phi,z,v,q0,xm1,
     1    xm2,nr,r,dr,r2,dl,rel)
       implicit real*8 (a-h,o-z)
       parameter (iorbs=33,iside=600)
       parameter (io2=iorbs*(iorbs+1)/2)
       parameter (ijive=io2*(io2+1)/2)
       parameter (lmax=4,ihmax=20,nrmax=4000,ntmax=10,npmax=60)
       dimension phi(nrmax),v(nrmax)
       dimension q0(nrmax),xm1(nrmax),xm2(nrmax)
       dimension r(nrmax),dr(nrmax),r2(nrmax)
       dl2=dl*dl/12.d0
       dl5=10.d0*dl2
       c=137.038d0
       alpha=rel/c
       za2=z*z*alpha*alpha
       a2=alpha*alpha/2.d0
       xl=l
       xlp=l+1
       xl2=0.5d0+xl
       xl4=xl2*xl2

c  then, we set up the leading power.
c  adjust for Desclaux's implementation of Numerov.

       if (rel.eq.0.d0) then
         ss=xlp
       else
         rtest=1.d0-za2
         if (rtest.lt.0.d0) then
           write (6,*) 'Z>137 IS TOO BIG.'
           stop
         endif
         ss=dsqrt(rtest)
       endif
       ss2=ss-0.5d0

c  we shall set ief to -1 if energy is too low, +1 if too high.

       ief=0

c  see Desclaux and documentation to see the origin of the below equations.
c  here, we set up the first two points.

       t=e-v(1)
       xm0=1.d0+a2*t
       tm=xm0+xm0
       xmx=xm1(1)/xm0
       xk0=r2(1)*(tm*t-xmx*(xkappa/r(1)+0.75d0*xmx)+xm2(1)/tm)-xl4
       dk0=1.d0+dl2*xk0
       p0=dk0
       phi(1)=p0*dsqrt(xm0*r(1))/dk0

       t=e-v(2)
       xm=1.d0+a2*t
       tm=xm+xm
       xmx=xm1(2)/xm
       xk2=r2(2)*(tm*t-xmx*(xkappa/r(2)+0.75d0*xmx)+xm2(2)/tm)-xl4
       dk2=1.d0+dl2*xk2
       p1=dk2*((r(2)/r(1))**ss2-(r(2)-r(1))*z/xlp)*dsqrt(xm0/xm)
       phi(2)=p1*dsqrt(xm*r(2))/dk2

c  if istop is set, the we know to stop there.  If it is zero, it shall
c  then be set to the classical turning point.

       is0=istop
       if (istop.eq.0) then
         do 10 j=nr-1,2,-1
         if (e.gt.v(j)) goto 15
 10       continue
         ief=-1
         return
 15       istop=j
       endif

c  initialize number of nodes, and determine the ideal number.

       nn=0
       nnideal=n-l-1

c  integrate out.  count nodes, and stop along the way if there are too many.

       do 50 i=3,istop+2
       t=e-v(i)
       xm=1.d0+a2*t
       tm=xm+xm
       xmx=xm1(i)/xm
       p2=(2.d0-dl5*xk2)*p1/dk2-p0
       xk2=r2(i)*(tm*t-xmx*(xkappa/r(i)+0.75d0*xmx)+xm2(i)/tm)-xl4
       dk2=1.d0+dl2*xk2
       phi(i)=p2*dsqrt(xm*r(i))/dk2
       if (dabs(p2).gt.10000000000.d0) then
         do 20 j=1,i
         phi(j)=phi(j)/p2
 20       continue
         p0=p0/p2
         p1=p1/p2
         p2=p2/p2
       endif
       if (p2*p1.lt.0.d0) then
         nn=nn+1
         if (nn.gt.nnideal) then
           ief=1
           return
         endif
       endif
       p0=p1
       p1=p2
 50     continue

       if (istop.gt.0) then
         psip2=(phi(istop+2)-phi(istop-2))
         psip1=(phi(istop+1)-phi(istop-1))
         psip=(8.d0*psip1-psip2)/(12.d0*dl*r(istop))
         x0=psip/phi(istop)
       endif

       if (is0.ne.0) return

       do 150 i=istop+3,nr-1
       t=e-v(i)
       xm=1.d0+a2*t
       tm=xm+xm
       xmx=xm1(i)/xm
       p2=(2.d0-dl5*xk2)*p1/dk2-p0
       if (p2/p1.gt.1.d0) then
         ief=-1
         return
       endif
       xk2=r2(i)*(tm*t-xmx*(xkappa/r(i)+0.75d0*xmx)+xm2(i)/tm)-xl4
       dk2=1.d0+dl2*xk2
       phi(i)=p2*dsqrt(xm*r(i))/dk2
       if (dabs(p2).gt.10000000000.d0) then
         do 120 j=1,i
         phi(j)=phi(j)/p2
 120      continue
         p0=p0/p2
         p1=p1/p2
         p2=p2/p2
       endif
       if (p2*p1.lt.0.d0) then
         nn=nn+1
         if (nn.gt.nnideal) then
           ief=1
           return
         endif
       endif
       p0=p1
       p1=p2
 150    continue
       return
       end
c-------------------------------------------------------------------------
c  routine to generate Clebsch-Gordan coefficients, in the form of
c  cg(l1,l2,L,m1,m2) = <l1,m1;l2,m2|L,m1+m2>, according to Rose's
c  'Elementary Theory of Angular Momentum', p. 39, Wigner's formula.
c  those coefficients listed are only those for which l1.ge.l2.
c  coefficients known to be zero because of either the L or M
c  selection rules are not computed, and should not be sought.

       subroutine clebschgordan(nel,nl,cg)

       implicit real*8 (a-h,o-z)
       parameter (iorbs=33,iside=600)
       parameter (io2=iorbs*(iorbs+1)/2)
       parameter (ijive=io2*(io2+1)/2)
       parameter (lmax=4,ihmax=20,nrmax=4000,ntmax=10,npmax=60)
       dimension nl(iorbs)
       dimension cg(0:6,0:6,0:12,-6:6,-6:6),si(0:32),fa(0:32)

       lmx=0
       do 14 i=1,nel
       if (nl(i).gt.lmx) lmx=nl(i)
 14     continue

       si(0)=1.d0
       fa(0)=1.d0
       do 50 i=1,32
       si(i)=-si(i-1)
       fa(i)=dfloat(i)*fa(i-1)
 50     continue

       do 100 l1=0,lmx
       do 100 l2=0,l1
 52     format (1x,i3,a3,i3)
       do 100 m1=-l1,l1
       do 100 m2=-l2,l2
       m3=m1+m2
       lmin=iabs(l1-l2)
       if (lmin.lt.iabs(m3)) lmin=iabs(m3)
       do 100 l3=lmin,l1+l2
       prefactor=dfloat(2*l3+1)
       prefactor=prefactor*fa(l3+l1-l2)/fa(l1+l2+l3+1)
       prefactor=prefactor*fa(l3-l1+l2)/fa(l1-m1)
       prefactor=prefactor*fa(l1+l2-l3)/fa(l1+m1)
       prefactor=prefactor*fa(l3+m3)/fa(l2-m2)
       prefactor=prefactor*fa(l3-m3)/fa(l2+m2)
       prefactor=dsqrt(prefactor)
       sum=0.d0
       numax=l3-l1+l2
       if ((l3+m3).lt.numax) numax=l3+m3
       numin=0
       if (l1-l2-m3.lt.numin) numin=-(l1-l2-m3)
       do 90 nu=numin,numax
       sum=sum+(si(nu+l2+m2)/fa(nu))*fa(l2+l3+m1-nu)*fa(l1-m1+nu)
     1    /fa(l3-l1+l2-nu)/fa(l3+m3-nu)/fa(nu+l1-l2-m3)
 90     continue
       cg(l1,l2,l3,m1,m2)=prefactor*sum
       cg(l2,l1,l3,m2,m1)=si(l1+l2+l3)*prefactor*sum
 100    continue

       return

       end
       subroutine pseudo(etot,nst,rel,alfa,
     1  nr,rmin,rmax,r,dr,r2,dl,
     1  phe,orb,njrc,vi,zorig,xntot,nel,
     1  no,nl,xnj,ev,occ,is,ek,iuflag,vctab)
       implicit real*8 (a-h,o-z)
       parameter (iorbs=33,iside=600)
       parameter (io2=iorbs*(iorbs+1)/2)
       parameter (ijive=io2*(io2+1)/2)
       parameter (lmax=4,ihmax=20,nrmax=4000,ntmax=10,npmax=60)
       dimension r(nrmax),dr(nrmax),r2(nrmax)
       dimension q0(nrmax),xm1(nrmax),xm2(nrmax),phi(nrmax),v(nrmax)
       dimension njrc(4),njrcdummy(4),vi(nrmax,7),phe(nrmax,iorbs)
       dimension no(iorbs),nl(iorbs),nm(iorbs),xnj(iorbs)
       dimension ek(iorbs),ev(iorbs),occ(iorbs),is(iorbs)
       dimension rpower(nrmax,0:7),orb(nrmax,iorbs)
       dimension vctab(nrmax,0:3)
       do 10 i=1,nel
       nm(i)=0
 10     continue
       njrcdummy(1)=njrc(1)
       njrcdummy(2)=njrc(2)
       njrcdummy(3)=njrc(3)
       njrcdummy(4)=njrc(4)
       write (6,*) 'PLEASE ENTER NP,CORPOL,RNORM'
       read (5,*) np,corpol,rnorm
       xntot=0.d0
       do 1200 i=np,nel
       write (6,42) 'l=',nl(i),' ...'
 42     format(1x,1a2,1i1,1a4)
       lp2=nl(i)+nl(i)+1
       e=ev(i)
       do 57 j=1,nr
       orb(j,i)=orb(j,i)+vctab(j,nl(i))
 57     continue
       idoflag=1
       ns=1
       call setqmm(i,orb,nl(i),ns,idoflag,vi(1,lp2),zeff,zorig,rel,
     1              nr,r,r2,dl,q0,xm1,xm2,njrcdummy,vi)
       do 60 j=1,nr
       orb(j,i)=0.d0
 60     continue
c
c  you can replacing subroutine pseudize with any type of PP generation
c  you want...however, kleinman-bylanderization would take more coding...
c
       call pseudize(i,orb,e,nl(i),xnj(i),no(i),njrc,zeff,vi(1,lp2),
     1                q0,xm1,xm2,nr,rmin,rmax,r,dr,r2,dl,rel)
       write (6,*) 'WE HAVE GOT THUS FAR...'
       no(i)=nl(i)+1
       ruse=0.d0
       xkappa=-1.d0
       call elsolve(i,occ(i),no(i),nl(i),xkappa,xnj(i),
     1               zorig,zeff,ev(i),phe(1,i),vi(1,lp2),
     1               q0,xm1,xm2,nr,r,dr,r2,dl,ruse)
       write (6,*) nl(i),ev(i)
       xntot=xntot+occ(i)
       if (lp2.eq.1) goto 1200
       do 1170 j=1,nr
       vi(j,lp2-1)=vi(j,lp2)
 1170   continue
 1200   continue
       write (6,*) 'everything is pseudized'
       do 1210 i=np,nel
       inew=1+i-np
       no (inew)=no (i)
       nl (inew)=nl (i)
       nm (inew)=nm (i)
       xnj(inew)=xnj(i)
       is (inew)=1
       ev (inew)=ev (i)
       occ(inew)=occ(i)
       do 1210 j=1,nr
       phe(j,inew)=phe(j,i)
 1210   continue
       nel=1+nel-np
       do 1212 i=0,7
       xi=i
       do 1212 k=1,nr
       rpower(k,i)=r(k)**xi
 1212   continue
       write (6,*) 'everything is scaled down...ready for unscreening'
       xnum=100.d0
       ratio=1.d0
       call getpot(etot,nst,rel,alfa,dl,nr,dr,r,r2,
     1              xntot,phe,ratio,orb,occ,is,
     1              nel,nl,nm,no,xnj,rpower,xnum,etot2,iuflag)
       write (6,*) 'screening effects in pseudo atom computed...'
       do 1250 k=1,nel
       lp2=nl(k)+nl(k)+1
       do 1250 j=1,nr
                    vi(j,lp2  )=vi(j,lp2  )-orb(j,k)
       if (lp2.gt.1) vi(j,lp2-1)=vi(j,lp2-1)-orb(j,k)
 1250   continue
       write (6,*) 'we got past the unscreening...'
       do 1300 j=1,nr
       vl =     (     vi(j,2)+2.d0*vi(j,3))/3.d0
       vso=2.d0*(     vi(j,3)-     vi(j,2))/3.d0
       vi(j,2)=vso
       vi(j,3)=vl
       vl =     (2.d0*vi(j,4)+3.d0*vi(j,5))/5.d0
       vso=2.d0*(     vi(j,5)-     vi(j,4))/5.d0
       vi(j,4)=vso
       vi(j,5)=vl
 2222   format (5f8.4)
       vl =     (3.d0*vi(j,6)+4.d0*vi(j,7))/7.d0
       vso=2.d0*(     vi(j,7)-     vi(j,6))/7.d0
       vi(j,6)=vso
       vi(j,7)=vl
 1300   continue
       rel=0.d0
       write (6,*) 'we got past the spin-orbit jazz'
       izuse=dabs(vi(nr-2,1)*r(nr-2))+0.5d0
        zuse=izuse
       do 2100 k=1,nr
       if (r(k).gt.rnorm) then
         videal=-zuse/r(k)-corpol/(2.d0*r(k)**4.d0)
         vi(k,1)=videal
         vi(k,3)=videal
         vi(k,5)=videal
         vi(k,7)=videal
         vi(k,2)=0.d0
         vi(k,4)=0.d0
         vi(k,6)=0.d0
       endif
 2100   continue
       write (6,*) 'we got to the end'
       return
       end
c----------------------------------------------------------------------
       subroutine parabreg(f,fp,fpp,rf,vf)
       implicit real*8 (a-h,o-z)
       dimension rf(3),vf(3)
       f=vf(2)
       r21=rf(2)-rf(1)
       r32=rf(3)-rf(2)
       v21=vf(2)-vf(1)
       v32=vf(3)-vf(2)
       fp=(v21+v32)/(r21+r32)
       fpp=(v32/r32-v21/r21)/((r21+r32)/2.d0)
       return
       end
c----------------------------------------------------------------------
       real*8 function hb(x,factor)
       implicit real*8 (a-h,o-z)
       if (x.gt.3.d0) hb=0.d0
       if (x.le.3.d0) hb=0.01d0**((dsinh(x/factor)/1.1752d0)**2.d0)
       return
       end
c----------------------------------------------------------------------
       subroutine fitx0(i,orb,rcut,njrc,e,l,xj,n,jrt,xideal,phi,
     1                   zeff,v,q0,xm1,xm2,nr,r,dr,r2,dl,rel,factor)
       implicit real*8 (a-h,o-z)
       parameter (iorbs=33,iside=600)
       parameter (io2=iorbs*(iorbs+1)/2)
       parameter (ijive=io2*(io2+1)/2)
       parameter (lmax=4,ihmax=20,nrmax=4000,ntmax=10,npmax=60)
       dimension r(nrmax),dr(nrmax),r2(nrmax)
       dimension njrc(4),orb(nrmax,iorbs)
       dimension q0(nrmax),xm1(nrmax),xm2(nrmax)
       dimension phi(nrmax),v(nrmax)
       vl=-1000000.d0
       vh=+1000000.d0
 115    idoflag=2
       ns=1
       xkappa=-1.d0
       call setqmm(i,orb,l,ns,idoflag,v,zeff,dummy,rel,
     1              nr,r,r2,dl,q0,xm1,xm2,njrc,dummy)
       call integ(e,l,xkappa,n,nn,jrt,ief,xactual,phi,zeff,v,
     1  q0,xm1,xm2,nr,r,dr,r2,dl,rel)
       if (nn.ne.0) then
         vl=v(1)
         xla=1.d0
       else
         if (xactual.gt.xideal) then
           vh=v(1)
         else
           vl=v(1)
         endif
         xerror=xideal-xactual
         if (dabs(xerror).lt.0.000000001d0) return
         dxdla=0.d0
         do 120 ii=1,jrt
         dxdla=dxdla+dr(ii)*phi(ii)*phi(ii)*hb(r(ii)/rcut,factor)
 120      continue
         dxdla=2.d0*dxdla/(phi(jrt)*phi(jrt))
         xla=xerror/dxdla
       endif
       vmaybe=v(1)+xla
       if ((vmaybe.gt.vh).or.(vmaybe.lt.vl)) xla=(vl+vh)/2.d0-v(1)
       do 130 ii=1,jrt-1
       v(ii)=v(ii)+xla*hb(r(ii)/rcut,factor)
 130    continue
       goto 115
       end
c----------------------------------------------------------------------
       subroutine pseudize(i,orb,ev,l,xj,n,njrc,zeff,v,
     1                      q0,xm1,xm2,nr,rmin,rmax,r,dr,r2,dl,rel)
       implicit real*8 (a-h,o-z)
       parameter (iorbs=33,iside=600)
       parameter (io2=iorbs*(iorbs+1)/2)
       parameter (ijive=io2*(io2+1)/2)
       parameter (lmax=4,ihmax=20,nrmax=4000,ntmax=10,npmax=60)
       dimension r(nrmax),dr(nrmax),r2(nrmax)
       dimension phi(nrmax),v(nrmax),rf(3),vf(3)
       dimension phi0(nrmax),yl(nrmax),vraw(nrmax)
       dimension q0(nrmax),xm1(nrmax),xm2(nrmax)
       dimension njrc(4),orb(nrmax,iorbs)
       lp=l+1
       xkappa=-1.d0
       istop=nr
 40     istop=istop-1
       if (ev.le.q0(istop)) goto 40
       call integ(ev,l,xkappa,n,nn,istop,ief,xdummy,phi,zeff,v,
     1  q0,xm1,xm2,nr,r,dr,r2,dl,rel)
 50     write (6,*) 'PLEASE ENTER THE CUTOFF RADIUS, AND FACTOR.'
       read (5,*) rcut,factor
       if (rcut.lt.0.d0) then
         xnodefrac=-rcut
         j=istop
 55       j=j-1
         if (phi(j-1)/phi(j).gt.1.d0) goto 55
         if (n.gt.l+1) then
           k=j
 60         k=k-1
           if (phi(k-1)/phi(k).gt.0.d0) goto 60
         else
           k=1
         endif
         rcut=r(k)+xnodefrac*(r(j)-r(k))
       endif
       jrc=1.d0+dfloat(nr-1)*dlog(rcut /rmin)/dlog(rmax/rmin)
       rcut=r(jrc)
       rtest=2.d0*rcut
       jrt=1.d0+dfloat(nr-1)*dlog(rtest/rmin)/dlog(rmax/rmin)
       njrc(lp)=jrt
       rtest=r(jrt)
       switch=phi(jrt)/dabs(phi(jrt))
       write (6,92) 'RCUTOFF = ',rcut,'  JRC = ',jrc
       write (6,92) 'RTEST   = ',rtest,  '  JRT = ',jrt
 92     format (1x,1a10,1f8.4,1a8,1i5)
 94     format (1x,2d15.8)
       call integ(ev,l,xkappa,n,nn,jrt,ief,x00,phi,zeff,v,
     1  q0,xm1,xm2,nr,r,dr,r2,dl,rel)
       do 8000 ii=1,jrt
       phi(ii)=phi(ii)/phi(jrt)
 8000   continue
       xn00=0.d0
       do 8010 ii=1,jrt-1
       xn00=xn00+dr( ii)*phi( ii)*phi( ii)
 8010   continue
       xn00=xn00+dr(jrt)*phi(jrt)*phi(jrt)/2.d0
       de=0.0001d0
       ee=ev+de/2.d0
       call integ(ee,l,xkappa,n,nn,jrt,ief,xp,phi,zeff,v,
     1  q0,xm1,xm2,nr,r,dr,r2,dl,rel)
       ee=ev-de/2.d0
       call integ(ee,l,xkappa,n,nn,jrt,ief,xm,phi,zeff,v,
     1  q0,xm1,xm2,nr,r,dr,r2,dl,rel)
       c00=(xm-xp)/(2.d0*de)
       write (6,94)  c00,x00
       write (6,94) xn00
       ruse=0.d0
       v0=v(jrc)
       dvdl  =(8.d0*(v(jrc+1)-v(jrc-1))-(v(jrc+2)-v(jrc-2)))
     1         /(12.d0*dl)
       ddvdll=(16.d0*(v(jrc+1)+v(jrc-1))
     1         -30.d0*v(jrc)-v(jrc+2)-v(jrc-2))
     1         /(12.d0*dl*dl)
       dldr=1.d0/r(jrc)
       ddldrr=-1.d0/r2(jrc)
       v1=dvdl*dldr
       v2=dvdl*ddldrr+ddvdll*dldr*dldr
       b4=(v2*rcut-v1)/(8.d0*rcut**3.d0)
       b2=(v1-4.d0*b4*rcut**3.d0)/(2.d0*rcut)
       b0=v0-b4*rcut**4.d0-b2*rcut**2.d0
       do 110 ii=1,jrc
       rr=r(ii)
       v(ii)=b0+b2*rr**2.d0+b4*rr**4.d0
 110    continue
       call fitx0(i,orb,rcut,njrc,ev,l,xj,lp,jrt,x00,phi,zeff,v,
     1  q0,xm1,xm2,nr,r,dr,r2,dl,ruse,factor)
 180    do 200 ii=1,jrt
       phi0(ii)=phi(ii)
       vraw(ii)=v(ii)
 200    continue
 210    xi0=0.d0
       xi1=0.d0
       xi2=0.d0
       do 220 ii=1,jrt
       f=hb(r(ii)/rcut,factor)
       ph2=dr(ii)*phi0(ii)*phi0(ii)
       xi0=xi0+ph2
       if (ii.le.jrt) then
         xi1=xi1+ph2*f
         xi2=xi2+ph2*f*f
       endif
 220    continue
       ph2=phi0(jrt)*phi0(jrt)
       xi0=xi0/ph2
       xi1=xi1/ph2
       xi2=xi2/ph2
       quant=xi1*xi1+xi2*(c00-xi0)
       if (quant.gt.0.d0) then
         deltal=(dsqrt(xi1*xi1+xi2*(c00-xi0))-xi1)/xi2
       else
         deltal=(c00-xi0)/(2.d0*xi1)
       endif
       write (6,222) 'DELTAL = ',deltal
 222    format (1x,1a9,1f11.8)
 225    do 230 ii=1,jrt
       yl (ii)=phi0(ii)*hb(r(ii)/rcut,factor)
       phi(ii)=phi0(ii)+deltal*yl(ii)
       if (phi(ii).lt.0.d0) then
         write (6,*) 'BIG TROUBLE!!! CROSS AXIS!!!'
         stop
       endif
 230    continue
       do 300 ii=1,jrt-1
       if ((phi(ii).eq.0.).or.(yl(ii).eq.0.)) goto 1170
       jj=ii
       if (ii.eq.1) jj=2
       do 240 j=jj-1,jj+1
       rf(2+j-jj)=r(j)
       vf(2+j-jj)=hb(r(j)/rcut,factor)
 240    continue
       call parabreg(f,fp,fpp,rf,vf)
       do 242 j=jj-1,jj+1
       vf(2+j-jj)=phi0(j)
 242    continue
       call parabreg(psi,psip,psipp,rf,vf)
       v(ii)=vraw(ii)+
     1        (1.d0-phi0(ii)/phi(ii))*(2.d0*psip/psi*fp/f+fpp/f)/2.d0
 300    continue
 1170   call fitx0(i,orb,rcut,njrc,ev,l,xj,lp,jrt,x00,phi,zeff,v,
     1  q0,xm1,xm2,nr,r,dr,r2,dl,ruse,factor)
       call integ(ev,l,xkappa,n,nn,jrt,ief,x0,phi,zeff,v,
     1  q0,xm1,xm2,nr,r,dr,r2,dl,ruse)
       do 8015 ii=1,jrt
       phi(ii)=phi(ii)/phi(jrt)
 8015   continue
       xn0=0.d0
       do 8020 ii=1,jrt-1
       xn0=xn0+dr( ii)*phi( ii)*phi( ii)
 8020   continue
       xn0=xn0+dr(jrt)*phi(jrt)*phi(jrt)/2.d0
       de=0.0001d0
       ee=ev+de/2.d0
       call integ(ee,l,xkappa,n,nn,jrt,ief,xp,phi,zeff,v,
     1  q0,xm1,xm2,nr,r,dr,r2,dl,ruse)
       ee=ev-de/2.d0
       call integ(ee,l,xkappa,n,nn,jrt,ief,xm,phi,zeff,v,
     1  q0,xm1,xm2,nr,r,dr,r2,dl,ruse)
       c0=(xm-xp)/(2.d0*de)
       write (6,94)  c0,x0
       write (6,94) xn0
       if (dabs(c0-c00).ge.0.000000001d0) then
         dqddel=2.*(xi1+deltal*xi2)
         deltal=deltal+(c00-c0)/dqddel
         goto 225
       endif
       write (6,94)  c0,x0
       write (6,*) 'NCPP ACHIEVED !!!'
       return
       end
       subroutine fourier(nr,r,dr,r2,vi)
       implicit real*8 (a-h,o-z)
       parameter (iorbs=33,iside=600)
       parameter (io2=iorbs*(iorbs+1)/2)
       parameter (ijive=io2*(io2+1)/2)
       parameter (lmax=4,ihmax=20,nrmax=4000,ntmax=10,npmax=60)
       dimension r(nrmax),dr(nrmax),r2(nrmax),vi(nrmax,7)
       dimension a(nrmax),v1(nrmax),v2(nrmax)
       do 350 l=0,2
       lp2=l+l+1
       dl=dlog(r(2)/r(1))
       dl1=12.d0*dl
       dl2=12.d0*dl*dl
       do 220 i=1,nr
       a(i)=r(i)*vi(i,lp2)
 220    continue
       do 225 i=3,nr-2
       al =(-(a(i+2)-a(i-2))+ 8.d0*(a(i+1)-a(i-1))           )/dl1
c       all=(-(a(i+2)+a(i-2))+16.d0*(a(i+1)+a(i-1))-30.d0*a(i))/dl2
       ar =al/r(i)
       v1(i)=ar
 225    continue
       open (unit=20+l,status='unknown')
       do 300 ii=1,200
       q=dfloat(ii)/10.d0
       vq=0.d0
       do 250 i=3,nr-2
       vq=vq+dr(i)*dcos(q*r(i))*v1(i)
 250    continue
       write (20+l,*) q,vq
 300    continue
       close(unit=1)
 350    continue
       return
       end
       SUBROUTINE GETILLLS(PIN)

       IMPLICIT DOUBLE PRECISION (A-H,O-Z)

       DIMENSION FA(0:40),SI(0:40),PIN(0:8,0:8,0:16)

       FA(0)=1.D0
       SI(0)=1.D0
       DO 200 I=1,32
       FA(I)=DFLOAT(I)*FA(I-1)
       SI(I)=-SI(I-1)
 200    CONTINUE

       DO 1000 L=0,8
       DO 1000 M=0,8
       DO 1000 N=M+L,0,-2
       XI=0.D0
       XF=2.D0/2.D0**DFLOAT(N+L+M)
       NN=(N+1)/2
       MM=(M+1)/2
       LL=(L+1)/2
       DO 500 IA=NN,N
       AF=SI(IA)*FA(IA+IA)/FA(IA)/FA(N-IA)/FA(IA+IA-N)
       DO 500 IB=LL,L
       BF=SI(IB)*FA(IB+IB)/FA(IB)/FA(L-IB)/FA(IB+IB-L)
       DO 500 IC=MM,M
       XCF=SI(IC)*FA(IC+IC)/FA(IC)/FA(M-IC)/FA(IC+IC-M)
       XI=XI+XF*AF*BF*XCF/DFLOAT(IA*2+IB*2+IC*2-N-L-M+1)
 500    CONTINUE
       PIN(L,M,N)=XI
 1000   CONTINUE

       RETURN

       END
       subroutine hfdisk(iu,ir,etot,nst,rel,nr,rmin,rmax,r,rho,
     1                    zorig,xntot,ixflag,nel,
     1                    no,nl,xnj,is,ev,ek,occ,njrc,vi,phe,orb)
       implicit real*8 (a-h,o-z)
       parameter (iorbs=33,iside=600)
       parameter (io2=iorbs*(iorbs+1)/2)
       parameter (ijive=io2*(io2+1)/2)
       parameter (lmax=4,ihmax=20,nrmax=4000,ntmax=10,npmax=60)
       dimension no(iorbs),nl(iorbs),xnj(iorbs),is(iorbs)
       dimension ev(iorbs),ek(iorbs),occ(iorbs),r(nrmax)
       dimension phe(nrmax,iorbs),orb(nrmax,iorbs)
       dimension njrc(4),vi(nrmax,7),rho(nrmax)
       character*10 filename
       pi=4.*atan(1.)
       pi4=4.*pi
       rden=3.0
       if (iu.lt.0) then
         write (6,*) 'PLEASE ENTER FULL FILENAME.'
         read (5,52) filename
 52       format (a10)
         iu1=1
         open (unit=iu1,status='unknown',file=filename)
       else
         iu1=iu
         open (unit=iu1,status='unknown')
       endif
c
c define the logarithmic grid
c
       do 5 i=1,nr
          r(i)=rmin*(rmax/rmin)**(dfloat(i)/dfloat(nr))
5       continue
c
c obtain the charge density on the logarithmic grid
c
       do 20 i=1,nr
         rho(i)=.0
         do 10 ii=1,nel
             rho(i)=rho(i)+occ(ii)*phe(i,ii)**2
10        continue
20      continue
c
c write output file
c
       iprint=0
       write(iu1,550)iprint
550     format('RELA'/'RELAT. ATOMIC CHARGE DENSITY'/I1)
       write(iu1,54) rmin,rmax,nr,zorig
54      format (d15.8,d15.8,i5,f5.2)
        nden=nr*(log(rden/rmin)/log(rmax/rmin))
        write(iu1,56) (rho(j),j=1,nr)
56      format (f15.11)
       close (unit=iu1)
       return
       end
c------------------------------------------------------------------------
c  exchange correlation routine, via Ceperley-Alder, as parametrized by
c  Perdew and Zunger, Phys. Rev. B 23, 5048.  we use their interpolation
c  between the unpolarized and polarized gas for the correlation part.

       subroutine exchcorr(nst,rel,rr,rh1,rh2,ex,ec,ux1,ux2,uc1,uc2)

       implicit double precision(a-h,o-z)

       trd=1.d0/3.d0
       ft=4.d0/3.d0

       rh=rh1+rh2

c  if one spin type, average polarization

       if (nst.eq.1) then
         rh1=rh/2.d0
         rh2=rh/2.d0
       endif

c  get the n's, and the rs.

       pi=3.14159265358979d0
       fp=4.d0*pi
       xn1=rh1/(rr*fp)
       xn2=rh2/(rr*fp)
       xn=xn1+xn2

c  effect cutoff, to avoid overflow

       if ((nst.eq.3).or.(xn.lt.0.00000001d0)) then

         ex=0.d0
         ec=0.d0
         ux1=0.d0
         ux2=0.d0
         uc1=0.d0
         uc2=0.d0

       else

         rs=(3.d0/(fp*xn))**trd
         zeta=(xn1-xn2)/xn

         exchfactor=-0.930525546d0

         if (xn1.eq.0.d0) then
           fe1=1.d0
           fu1=1.d0
           ex1=0.d0
           ux1=0.d0
         else
           beta=0.028433756d0*xn1**trd
           b2=beta*beta
           eta=dsqrt(1.d0+b2)
           xl=dlog(beta+eta)
           fe1=1.d0-1.5d0*((beta*eta-xl)/b2)**2.d0
           fu1=-0.5d0+1.5d0*xl/beta/eta
           ex1=exchfactor*xn1**trd
           ux1=4.d0*ex1/3.d0
         endif

         if (xn2.eq.0.d0) then
           fe2=1.d0
           fu2=1.d0
           ex2=0.d0
           ux2=0.d0
         else
           beta=0.028433756d0*xn2**trd
           b2=beta*beta
           eta=dsqrt(1.d0+b2)
           xl=dlog(beta+eta)
           fe2=1.d0-1.5d0*((beta*eta-xl)/b2)**2.d0
           fu2=-0.5d0+1.5d0*xl/beta/eta
           ex2=exchfactor*xn2**trd
           ux2=4.d0*ex2/3.d0
         endif

c  these next lines do the Ceperley-Alder correlation

         if (rs.ge.1.d0) then

           rootr=dsqrt(rs)

           gamma=-0.1423d0
           beta1=1.0529d0
           beta2=0.3334d0
           denom=(1.d0+beta1*rootr+beta2*rs)
           ecu=gamma/denom
           ucu=ecu*(1.d0+7.d0/6.d0*beta1*rootr+ft*beta2*rs)/denom

           gamma=-0.0843d0
           beta1=1.3981d0
           beta2=0.2611d0
           denom=(1.d0+beta1*rootr+beta2*rs)
           ecp=gamma/denom
           ucp=ecp*(1.d0+7.d0/6.d0*beta1*rootr+ft*beta2*rs)/denom

         else

           xlr=dlog(rs)
           rlr=rs*xlr

           au= 0.0311d0
           bu=-0.048d0
           cu= 0.002d0
           du=-0.0116d0
           ecu=au*xlr+bu+cu*rlr+du*rs
           ucu=au*xlr+(bu-au/3.d0)+2.d0/3.d0*cu*rlr+(2.d0*du-cu)*rs/3.d0

           ap= 0.01555d0
           bp=-0.0269d0
           cp= 0.0007d0
           dp=-0.0048d0
           ecp=ap*xlr+bp+cp*rlr+dp*rs
           ucp=ap*xlr+(bp-ap/3.d0)+2.d0/3.d0*cp*rlr+(2.d0*dp-cp)*rs/3.d0

         endif

c  if we are nonrelativistic, turn off the MacDonald-Vosko correction.

         if (rel.eq.0.d0) then
           fe1=1.d0
           fu1=1.d0
           fe2=1.d0
           fu2=1.d0
         endif

c  interpolate the correlation energies.

         denom=2.d0**ft-2.d0
         f=((1.d0+zeta)**ft+(1.d0-zeta)**ft-2.d0)/denom
         dfdz=ft/denom*((1.d0+zeta)**trd-(1.d0-zeta)**trd)
         ec=ecu+f*(ecp-ecu)
         uc1=ucu+f*(ucp-ucu)+(ecp-ecu)*(1.d0-zeta)*dfdz
         uc2=ucu+f*(ucp-ucu)+(ecp-ecu)*(-1.d0-zeta)*dfdz

c  get the final functional and potential.

         ex=(xn1*fe1*ex1+xn2*fe2*ex2)/xn
         ux1=fu1*ux1
         ux2=fu2*ux2
         uc1=uc1
         uc2=uc2

       endif

       return

       end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
