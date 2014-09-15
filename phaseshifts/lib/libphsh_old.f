c=======================================================================
c libphsh.f
c
c Modified code from the Barbieri/Van Hove phase shift (phshift) 
c package available at:
c 
c http://www.icts.hkbu.edu.hk/vanhove/VanHove_files/leed/leedpack.html  
c
c=======================================================================

c-----------------------------------------------------------------------
c hartfock subroutine:
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
    
       open(unit=5, file=input_file, status='old')
       rel=0.d0
 10    read (5,20) ichar
 20    format (1a1)

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
         call abinitio(etot,nst,rel,alfa,nr,r,
     1    dr,r2,dl,phe,njrc,vi,zorig,xntot,nel,
     1    no,nl,xnj,ev,occ,is,ek,orb,iuflag)
       endif

       if (ichar == 'i') call initiali(zorig,nr,rmin,rmax,
     1    r,dr,r2,dl,njrc,xntot,nel)

       if (ichar == 'q') return

       if (ichar == 'w') then
         ixflag=1
         iu=-1
         ir=0
         call hfdisk(iu,ir,etot,nst,rel,nr,rmin,rmax,r,rho,
     1      zorig,xntot,ixflag,nel,
     1      no,nl,xnj,is,ev,ek,occ,njrc,vi,phe,orb)
       endif

       if (ichar=='r') then
         iu=-1
         ir=1
         call hfdisk(iu,ir,etot,nst,rel,nr,rmin,rmax,r,rho,
     1      zorig,xntot,ixflag,nel,
     1      no,nl,xnj,is,ev,ek,occ,njrc,vi,phe,orb)
         call setgrid(nr,rmin,rmax,r,dr,r2,dl)
       endif

       if (ichar=='u') then
         ! ENTER IUFLAG. (0=U, 1=SU, 2=R)
         read (5,*) iuflag
       endif

       if (ichar=='c') then
         ! ENTER ALPHA,RS,RP,RD
         read (5,*) corpol,rs,rp,rd

         do k=1,nr
           fs=(1.d0-exp(-(r(k)/rs)**2.d0))**2.d0
           fp=(1.d0-exp(-(r(k)/rp)**2.d0))**2.d0
           fd=(1.d0-exp(-(r(k)/rd)**2.d0))**2.d0
           vctab(k,0)=-corpol/2.d0*fs*fs/r(k)**4.d0
           vctab(k,1)=-corpol/2.d0*fp*fp/r(k)**4.d0
           vctab(k,2)=-corpol/2.d0*fd*fd/r(k)**4.d0
         end do

       endif

       if (ichar == 'f') then
         ! IUNIT, CORPOL
         read  (5,*) iunit,corpol
         ! ILEV,INUM,EOLD
         read  (5,*) ilev,inum,eold
         xl=nl(ilev)

         if (inum == 1) then
           read (5,*) eav
         else
           read (5,*) e1,e2
           eav=(e1*xl+e2*(xl+1.d0))/(xl+xl+1.d0 )
         endif

         if (eav < 0.d0) eav=eold+eav
         if (iunit==2) eav=eav/2.d0
         if (iunit==3) eav=eav/27.2116d0
         if (iunit==4) eav=eav*0.000123985d0/27.2116d0

         sd=abs(abs(eav)-abs(ev(ilev)))
         rl= 0.d0
         rh=10.d0
         sl= 0.d0
         sh= 0.d0

 300     if (sl*sh<=0.00000001d0) rc=rl+(rh-rl)/2.d0
         if (sl*sh>0.00000001d0) rc=rl+(rh-rl)*(sd-sl)/(sh-sl)

         sc=0.d0

         do i=1,nr
           f=(1.d0-exp(-(r(i)/rc)**2.d0))**2.d0
           vcpp=corpol/(2.d0*r(i)**4.d0)*f*f
           sc=sc+dr(i)*phe(i,ilev)*phe(i,ilev)*vcpp
         end do

         if (sc>sd) rl=rc
         if (sc>sd) sl=sc
         if (sc<sd) rh=rc
         if (sc<sd) sh=sc

         write (6,*) rc,sc

         if (abs(sc-sd) > 0.000001d0) goto 300

       endif

       if (ichar == 'p') then
         call pseudo(etot,nst,rel,alfa,nr,rmin,rmax,r,dr,r2,dl,
     1                phe,orb,njrc,vi,zorig,xntot,nel,
     1                no,nl,xnj,ev,occ,is,ek,iuflag,vctab)
       endif

       if (ichar == 'g') then
         read(5,*) iu
         read(5,2202) jive
         read(5,2212) jive2
         read(5,2222) jive3

 2202    format(1x,1a11)
 2212    format(1x,1a60)
 2222    format(1x,1a70)

         zizv=abs(r(nr-1)*vi(nr-1,1))

         write (iu,2202) jive
         write (iu,2212) jive2
         write (iu,2222) jive3
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
           write (10,*) r(k),vi(k,1)*r(k)
           write (11,*) r(k),vi(k,3)*r(k)
           write (12,*) r(k),vi(k,5)*r(k)
           rold=r(k)
         end do

         ! close file streams
         close(unit=10)
         close(unit=11)
         close(unit=12)

       endif

       if (ichar == 'V') call fourier(nr,r,dr,r2,vi)

       goto 10

       return ! end

      end subroutine

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

       ! note: this will be good for going up to and including l=3...

       do i=0,7
         xi=i
         do k=1,nr
           rpower(k,i)=r(k)**xi
         end do
       end do

       ! read in nfc, nel.  refer to the documentation for their meanings.

       ! ENTER NFC, NEL, RATIO, ETOL, XNUM
 168   read (5,*) nfc,nel,ratio,etol,xnum

c     for all of the electrons, read in the quantum numbers.
c     get the total Hartree-active charge & initialize eigenvalues.

       xntot=0.d0

       write (6,*) 'N L M J S OCC ENERGY'

       do i=nfc+1,nel
         read (5,*) no(i),nl(i),nm(i),xnj(i),is(i),occ(i)
         ev(i)=0.d0
         xntot=xntot+occ(i)

         do j=1,nr
           phe(j,i)=0.d0
           orb(j,i)=0.d0

         end do

       end do

       ! initialize the parameters for self-consistency loop.
       ! ratio is the mixture of old and new field mixing.

 110   call atsolve(etot,nst,rel,alfa,eerror,nfc,nr,r,dr,r2,dl,phe,
     1              njrc,vi,zorig,xntot,nel,no,nl,nm,xnj,ev,occ,is,ek,
     1              ratio,orb,rpower,xnum,etot2,iuflag)

       eerror=eerror*(1.d0-ratio)/ratio
       !write (6,112) eerror,etot
112    format (1x,3f14.6)
       if (eerror > etol) goto 110

       ! write out information about the atom.

       do i=1,nel
         nj=xnj(i)+xnj(i)
         write (6,122) no(i),nl(i),nm(i),nj,'/2',is(i),occ(i),ev(i)
 122     format(1x,2i4,i2,i4,a2,i4,f10.4,f18.6)
       end do

       write (6,132) 'TOTAL ENERGY =  ',etot,etot*27.2116d0
 132   format (1x,a16,2f14.6)

       return

      end subroutine

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

       ! initialize eerror, the biggest change in an eigenvalue, and etot.

       eerror=0.d0
       etot=0.d0

       ! run through all the orbitals.  calculate those not in the core.

       do i=1,nel

         if (i>nfc) then

           idoflag=1
           call setqmm(i,orb,nl(i),is(i),idoflag,v,zeff,zorig,rel,
     1                 nr,r,r2,dl,q0,xm1,xm2,njrc,vi)

           xkappa=-1.d0
           if (abs(xnj(i)) > dfloat(nl(i))+0.25d0) xkappa=-nl(i)-1
           if (abs(xnj(i)) < dfloat(nl(i))-0.25d0) xkappa= nl(i)

           call elsolve(i,occ(i),no(i),nl(i),xkappa,xnj(i),zorig,zeff,
     1                  evi,phe(1,i),v,q0,xm1,xm2,nr,r,dr,r2,dl,rel)

           if (abs(ev(i)-evi) > eerror) eerror=abs(ev(i)-evi)
           ev(i)=evi

           ekk=0.d0
           ll=2

           do j=nr,1,-1
             dq=phe(j,i)*phe(j,i)
             ekk=ekk+(evi-orb(j,i))*dr(j)*dq*dfloat(ll)/3.d0
             ll=6-ll
           end do

           ek(i)=ekk

         endif

         ! add kinetic to total, including frozen core kinetic energy
         etot=etot+ek(i)*occ(i)

       end do

       call getpot(etot,nst,rel,alfa,dl,nr,dr,r,r2,xntot,
     1             phe,ratio,orb,occ,is,nel,nl,nm,no,xnj,rpower,xnum,
     1             etot2,iuflag)
       return

      end subroutine

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
       do i=1,nel
         do k=1,nr
           orb(k,i)=ratio1*orb(k,i)
         end do
       end do

       loop_el: do i=1,nel

         li=nl(i)
         mi=nm(i)

         jstart=i+1
         if ((xnj(i)<0.d0).or.
     1       (occ(i)>1.d0).or.
     1       (abs(alfa)>0.001d0)) jstart=i

         do j=jstart,nel

           if ((occ(i) == 0.d0) .and. (occ(j) == 0.d0)) exit loop_el

           lj=nl(j)
           mj=nm(j)

           ! direct coulomb

           lmx=2*li
           if (li > lj) lmx=2*lj

           ! l=0 is monopole or spherical term for direct coulomb. Therefore,
           ! when we have occ(i) or occ(j) greater than one, set lmx=0.

           if ((occ(i) > 1.d0) .or. (occ(j) > 1.d0) .or.
     1         (xnj(i) < 0.d0) .or. (xnj(j) < 0.d0)) lmx=0

           do la=lmx,0,-2
             lap=la+1
             coeff=dfloat((li+li+1)*(lj+lj+1))/dfloat((la+la+1))**2.d0*
     1              cg(li,li,la,mi,-mi)*cg(lj,lj,la,mj,-mj)*
     1              cg(li,li,la,0 , 0 )*cg(lj,lj,la,0 , 0 )
             if (mi+mj .ne. 2*((mi+mj)/2)) coeff=-coeff
             if (i==j) coeff=coeff/2.d0
             coeffi=occ(i)*coeff
             coeffj=occ(j)*coeff
             ri=ratio*coeffi
             rj=ratio*coeffj
             rc=coeff*occ(i)*occ(j)

             xouti=0.d0
             xoutj=0.d0

             do k=1,nr
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

               if (rpower(k,lap) .ne. 0.d0) then
                 xqj2(k)=xqj0(k)/rpower(k,lap)
               else
                 xqj2(k)=0.d0
               endif

               xoutj=xoutj+xqj2(k)
             end do

             xinti=xqi1(1)
             xintj=xqj1(1)
             xouti=2.d0*xouti-xqi2(1)
             xoutj=2.d0*xoutj-xqj2(1)

             do k=2,nr

               xinti=xinti+xqi1(k)+xqi1(k-1)
               xouti=xouti-xqi2(k)-xqi2(k-1)
               vali=xouti*rpower(k,la)
               if (rpower(k,lap) .ne. 0.d0) then
                 vali=vali+xinti/rpower(k,lap)
               endif
               orb(k,j)=orb(k,j)+ri*vali

               xintj=xintj+xqj1(k)+xqj1(k-1)
               xoutj=xoutj-xqj2(k)-xqj2(k-1)
               valj=xoutj*rpower(k,la)
               if (rpower(k,lap) .ne. 0.d0) then
                 valj=valj+xintj/rpower(k,lap)
               endif
               orb(k,i)=orb(k,i)+rj*valj

               etot=etot+rc*(xqi0(k)*valj+xqj0(k)*vali)

             end do

           end do

           if ((is(i) .ne. is(j)) .and.
     1         (occ(i) <= 1.d0) .and.
     1         (occ(j) <= 1.d0) .and.
     1         (xnj(i) >= 0.d0) .and.
     1         (xnj(j) >= 0.d0)) exit loop_el
           if (abs(alfa) >= 0.001d0) exit loop_el

           ! exchange interaction
           lmx=li+lj
           lmin=iabs(mi-mj)
           if ((occ(i) > 1.d0) .or. (occ(j) > 1.d0) .or.
     1         (xnj(i) < 0.d0) .or. (xnj(j) < 0.d0)) lmin=0

           do a=lmx,lmin,-2
             lap=la+1

             coeff=dfloat((li+li+1)*(lj+lj+1))/dfloat((la+la+1))**2.d0*
     1              (cg(li,lj,la,-mi,mj)*cg(li,lj,la,0,0))**2.d0
             if ((occ(i) > 1.d0) .or. (occ(j) > 1.d0) .or.
     1           (xnj(i) < 0.d0) .or. (xnj(j) < 0.d0))
     1           coeff=pin(li,lj,la)/4.d0
             if (i == j) coeff=coeff/2.d0
             coeffi=occ(i)*coeff
             coeffj=occ(j)*coeff
             ri=ratio*coeffi
             rj=ratio*coeffj
             rc=coeff*occ(i)*occ(j)
             xnum2=xnum*xnum

             xout=0.d0

             do k=1,nr
               xq0(k)=dr(k)*phe(k,i)*phe(k,j)/2.d0
               xq1(k)=xq0(k)*rpower(k,la)

               if (rpower(k,lap).ne.0.d0) then
                 xq2(k)=xq0(k)/rpower(k,lap)
               else
                 xq2(k)=0.d0
               endif

               xout=xout+xq2(k)
             end do

             xint=xq1(1)
             xout=2.d0*xout-xq2(1)

             do k=2,nr
               xint=xint+xq1(k)+xq1(k-1)
               xout=xout-xq2(k)-xq2(k-1)

               if (xq0(k).ne.0.d0) then
                 val=xout*rpower(k,la)
                 if (rpower(k,lap).ne.0.d0) val=val+xint/rpower(k,lap)
                 etot=etot-2.d0*xq0(k)*rc*val
                 xx=phe(k,j)/phe(k,i)

                 if (abs(xx)/xnum>1.d0) then
                   orb(k,i)=orb(k,i)-rj*xnum2/xx*val
                 else
                   orb(k,i)=orb(k,i)-rj*xx*val
                 endif

                 xx=phe(k,i)/phe(k,j)

                 if (abs(xx)/xnum>1.d0) then
                   orb(k,j)=orb(k,j)-ri*xnum2/xx*val
                 else
                   orb(k,j)=orb(k,j)-ri*xx*val
                 endif

               endif

             end do

           end do

         end do

       end do loop_el

       ! here we compute the charge density, if needed, for treating
       ! exchange/correlation in a local fashion...

       if (abs(alfa) >= 0.001d0) then
         if (alfa > 0.d0) then
           fx=1.0d0
           fc=1.0d0
         else
           fx=1.5d0*abs(alfa)
           fc=0.0d0
         endif

       ! note: we don't deal with spin-polarization in local exchange
       ! picture, since local exchange is totally wrong for such
       ! effects, anyway.  local exchange pretends charge density
       ! is paramagnetic.  also, local exchange treats everything
       ! as spherical.

         fourpi=16.d0*atan(1.d0)

         do i=1,nr
           xn=0.d0

           do j=1,nel
             xn=xn+phe(i,j)*phe(i,j)*occ(j)
           end do

           xn1=xn/2.d0
           xn2=xn/2.d0
           nst=2
           call exchcorr(nst,rel,r2(i),xn1,xn2,ex,ec,ux1,ux2,uc1,uc2)
           exc=fx*ex +fc*ec
           uxc=fx*ux1+fc*uc1
           etot=etot+dr(i)*xn*exc

           do j=1,nel
             orb(i,j)=orb(i,j)+uxc*ratio
           end do


         end do
       endif

       do i=1,nr
         if (iuflag.ne.0) then
           jj=1
 8960      ii=jj
 8965      if (ii==nel) goto 8970
           icond=0
           if ((no(jj)==no(ii+1)).and.(nl(jj)==nl(ii+1))
     1          .and.(iuflag==2)) icond=1
           if ((no(jj)==no(ii+1)).and.(nl(jj)==nl(ii+1))
     1          .and.(is(jj)==is(ii+1)).and.(iuflag==1)) icond=1

           if (icond==1) then
             ii=ii+1
             goto 8965
           endif

 8970      orba=0.d0
           div=0.d0

           do k=jj,ii
             div=div+occ(k)
             orba=orba+orb(i,k)*occ(k)
           end do

           if (div.ne.0.d0) then
             orba=orba/div

             do k=jj,ii
               orb(i,k)=orba
             end do

           endif

           if (ii.ne.nel) then
             jj=ii+1
             goto 8960
           endif

         endif

       end do

       return

      end subroutine

c-----------------------------------------------------------------------
      subroutine elsolve(i,occ,n,l,xkappa,xj,zorig,zeff,e,phi,v,
     1                   q0,xm1,xm2,nr,r,dr,r2,dl,rel)

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
 155   e=(el+eh)/2.d0
       istop=0

       call integ(e,l,xkappa,n,nn,istop,ief,x0,phi,zeff,v,q0,xm1,
     1            xm2,nr,r,dr,r2,dl,rel)

       if (nn < n-l-1) ief=-1

 200   if (ief .ne. 1) then
         el=e
         if (el > -0.001d0) then
           write (6,*) 'MIXING TOO STRONG FOR LEVEL : ', i
           return
         endif
       endif

       if (ief .ne. -1) eh=e
       if (eh-el > etol) goto 155
       if (abs(abs(xj)-abs(dfloat(l))) > 0.25d0)
     1    call augment(e,l,xj,phi,v,nr,r,dl)
       aa=0.d0

       do j=1,nr
         aa=aa+phi(j)*phi(j)*dr(j)
       end do

       xnorm=sqrt(aa)

       do j=1,nr
         phi(j)=phi(j)/xnorm
       end do

       return

      end subroutine

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

       if (abs(xj)>dfloat(l)+0.25d0) xkappa=-l-1
       if (abs(xj)<dfloat(l)-0.25d0) xkappa= l

       do j=4,nr-3
         if (phi(j).ne.0.d0) then
           g0=phi(j)
           ga=(phi(j+1)-phi(j-1))
           gb=(phi(j+2)-phi(j-2))/2.d0
           gc=(phi(j+3)-phi(j-3))/3.d0
           gg=((1.5d0*ga-0.6d0*gb+0.1d0*gc)/(2.d0*dl)+xkappa*g0)/r(j)
           f0=c*gg/(e-v(j)+c2)
           phi2(j)=sqrt(g0*g0+f0*f0)
           if (g0<0.d0) phi2(j)=-phi2(j)
         else
           phi2(j)=phi(j)
         endif

       end do

       do j=1,3
         phi2(j)=phi(j)*phi(4)/phi2(4)
       end do

       do j=1,nr
         phi(j)=phi2(j)
       end do

       return

      end subroutine

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
       if (lp>4) lpx=4
       lp2=l+l+1
       if (lp2>7) lp2=7
       zeff=zorig
       if (njrc(lpx)>0) zeff=0.d0
       zaa=zeff*aa
       za2=zeff*a2

       if (idoflag.ne.0) then
         if (njrc(lpx)==0) then

           if (idoflag==1) then
             do j=1,nr
               v(j)=-zeff/r(j)+orb(j,i)
             end do
           endif

           do j=2,nr-1
             dvdl=(orb(j+1,i)-orb(j-1,i))/(2.d0*dl)
             ddvdrr=((orb(j+1,i)
     1         +orb(j-1,i)-2.d0*orb(j,i) )/(dl*dl)-dvdl)/r2(j)
             xm1(j)=-a2*dvdl/r(j)-za2/r2(j)
             xm2(j)=-a2*ddvdrr+zaa/r2(j)/r(j)
           end do

           xm1(nr)=xm1(nr-1)
           xm2(nr)=xm2(nr-1)
           xm1(1)=xm1(2)+za2/r2(2)-za2/r2(1)
           xm2(1)=xm2(2)-zaa/r2(2)/r(2)+zaa/r2(1)/r(1)

         else
           if (idoflag==1) then
             do j=1,nr
               v(j)=vi(j,lp2)+orb(j,i)
             end do
           endif

           do j=2,nr-1
             dvdl=(v(j+1)-v(j-1))/(2.d0*dl)
             ddvdrr=((v(j+1)+v(j-1)-2.d0*v(j))/(dl*dl)-dvdl)/r2(j)
             xm1(j)=-a2*dvdl/r(j)
             xm2(j)=-a2*ddvdrr
           end do

           xm1(nr)=xm1(nr-1)
           xm2(nr)=xm2(nr-1)
           xm1(1)=xm1(2)
           xm2(1)=xm2(2)

         endif

       endif

       ! figure out the (Desclaux-Numerov) effective potential.

       xlb=(dfloat(l)+0.5d0)**2.d0/2.d0

       do j=1,nr
         vj=v(j)
         q0(j)=vj*(1.d0-a2*vj)+xlb/r2(j)
       end do

       return

      end subroutine

c----------------------------------------------------------------------
      subroutine initiali(zorig,nr,rmin,rmax,r,dr,r2,dl,njrc,
     1                    xntot,nel)

       implicit real*8 (a-h,o-z)
       parameter (iorbs=33,iside=600)
       parameter (io2=iorbs*(iorbs+1)/2)
       parameter (ijive=io2*(io2+1)/2)
       parameter (lmax=4,ihmax=20,nrmax=4000,ntmax=10,npmax=60) 
       dimension r(nrmax),dr(nrmax),r2(nrmax),njrc(4)

       ! ENTER Z, NR
       read (5,*) zorig,nr
       rmin=0.0001d0/zorig
       rmax=800.d0/sqrt(zorig)
       call setgrid(nr,rmin,rmax,r,dr,r2,dl)

       do j=1,4
         njrc(j)=0
       end do

       xntot=0.d0
       nel=0

       return

      end subroutine

c---------------------------------------------------------------------------
      subroutine setgrid(nr,rmin,rmax,r,dr,r2,dl)

       implicit real*8 (a-h,o-z)
       parameter (iorbs=33,iside=600)
       parameter (io2=iorbs*(iorbs+1)/2)
       parameter (ijive=io2*(io2+1)/2)
       parameter (lmax=4,ihmax=20,nrmax=4000,ntmax=10,npmax=60) 
       dimension r(nrmax),dr(nrmax),r2(nrmax)

       ratio=rmax/rmin
       dl=log(ratio)/dfloat(nr)
       xratio=exp(dl)
       xr1=sqrt(xratio)-sqrt(1.d0/xratio)

       do i=1,nr
         r(i)=rmin*xratio**dfloat(i)
         dr(i)=r(i)*xr1
         r2(i)=r(i)*r(i)
       end do

       return

      end subroutine

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

       ! set up the leading power.
       ! adjust for Desclaux's implementation of Numerov.

       if (rel==0.d0) then
         ss=xlp
       else
         rtest=1.d0-za2

         if (rtest<0.d0) then
           write (6,*) 'Z>137 IS TOO BIG.'
           return
         endif  

         ss=sqrt(rtest)
       endif

       ss2=ss-0.5d0

       ! set ief to -1 if energy is too low, +1 if too high.

       ief=0

       ! see Desclaux and documentation for origin of the equations below.
       ! here, we set up the first two points.

       t=e-v(1)
       xm0=1.d0+a2*t
       tm=xm0+xm0
       xmx=xm1(1)/xm0
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

       is0=istop
       if (istop==0) then

         do j=nr-1,2,-1
           if (e>v(j)) goto 15
         end do

         ief=-1
         return
 15      istop=j
       endif

       ! initialize number of nodes, and determine the ideal number.
       nn=0
       nnideal=n-l-1

       ! integrate out.
       ! count nodes, and stop along the way if there are too many.
       do i=3,istop+2
         t=e-v(i)
         xm=1.d0+a2*t
         tm=xm+xm
         xmx=xm1(i)/xm
         p2=(2.d0-dl5*xk2)*p1/dk2-p0
         xk2=r2(i)*(tm*t-xmx*(xkappa/r(i)+0.75d0*xmx)+xm2(i)/tm)-xl4
         dk2=1.d0+dl2*xk2
         phi(i)=p2*sqrt(xm*r(i))/dk2
         if (abs(p2)>10000000000.d0) then

           do j=1,i
             phi(j)=phi(j)/p2
           end do

           p0=p0/p2
           p1=p1/p2
           p2=p2/p2
         endif

         if (p2*p1<0.d0) then
           nn=nn+1
           if (nn>nnideal) then
             ief=1
             return
           endif
         endif

         p0=p1
         p1=p2

       end do

       if (istop>0) then
         psip2=(phi(istop+2)-phi(istop-2))
         psip1=(phi(istop+1)-phi(istop-1))
         psip=(8.d0*psip1-psip2)/(12.d0*dl*r(istop))
         x0=psip/phi(istop)
       endif

       if (is0.ne.0) return

       do i=istop+3,nr-1
         t=e-v(i)
         xm=1.d0+a2*t
         tm=xm+xm
         xmx=xm1(i)/xm
         p2=(2.d0-dl5*xk2)*p1/dk2-p0

         if (p2/p1>1.d0) then
           ief=-1
           return
         endif

         xk2=r2(i)*(tm*t-xmx*(xkappa/r(i)+0.75d0*xmx)+xm2(i)/tm)-xl4
         dk2=1.d0+dl2*xk2
         phi(i)=p2*sqrt(xm*r(i))/dk2

         if (abs(p2)>10000000000.d0) then
           do j=1,i
             phi(j)=phi(j)/p2
           end do

           p0=p0/p2
           p1=p1/p2
           p2=p2/p2
         endif

         if (p2*p1<0.d0) then
           nn=nn+1
           if (nn>nnideal) then
             ief=1
             return
           endif
         endif

         p0=p1
         p1=p2

       end do

       return

      end subroutine

c-------------------------------------------------------------------------
c  routine to generate Clebsch-Gordan coefficients, in the form of 
c  cg(l1,l2,L,m1,m2) = <l1,m1;l2,m2|L,m1+m2>, according to Rose's 
c  'Elementary Theory of Angular Momentum', p. 39, Wigner's formula.
c  those coefficients listed are only those for which l1>=l2.
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
       do i=1,nel
         if (nl(i)>lmx) lmx=nl(i)
       end do

       si(0)=1.d0
       fa(0)=1.d0
       do i=1,32
         si(i)=-si(i-1)
         fa(i)=dfloat(i)*fa(i-1)
       end do

       do l1=0,lmx

         do l2=0,l1
 52        format (1x,i3,a3,i3)

           do m1=-l1,l1

             do m2=-l2,l2
               m3=m1+m2
               lmin=iabs(l1-l2)
               if (lmin<iabs(m3)) lmin=iabs(m3)

               do l3=lmin,l1+l2
                 prefactor=dfloat(2*l3+1)
                 prefactor=prefactor*fa(l3+l1-l2)/fa(l1+l2+l3+1)
                 prefactor=prefactor*fa(l3-l1+l2)/fa(l1-m1)
                 prefactor=prefactor*fa(l1+l2-l3)/fa(l1+m1)
                 prefactor=prefactor*fa(l3+m3)/fa(l2-m2)
                 prefactor=prefactor*fa(l3-m3)/fa(l2+m2)
                 prefactor=sqrt(prefactor)
                 sum=0.d0
                 numax=l3-l1+l2
                 if ((l3+m3)<numax) numax=l3+m3
                 numin=0
                 if (l1-l2-m3<numin) numin=-(l1-l2-m3)

                 do nu=numin,numax
                   sum=sum+(si(nu+l2+m2)/fa(nu))*fa(l2+l3+m1-nu)*
     1                 fa(l1-m1+nu)/fa(l3-l1+l2-nu)/fa(l3+m3-nu)/
     1                 fa(nu+l1-l2-m3)
                 end do

                 cg(l1,l2,l3,m1,m2)=prefactor*sum
                 cg(l2,l1,l3,m2,m1)=si(l1+l2+l3)*prefactor*sum

               end do  !l3

             end do  !m2

           end do  !m1

         end do  !l2

       end do  !l1

       return

      end subroutine

c-----------------------------------------------------------------------
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

       do i=1,nel
         nm(i)=0
       end do

       njrcdummy(1)=njrc(1)
       njrcdummy(2)=njrc(2)
       njrcdummy(3)=njrc(3)
       njrcdummy(4)=njrc(4)

       ! ENTER NP,CORPOL,RNORM
       read (5,*) np,corpol,rnorm
       xntot=0.d0

       do i=np,nel
         write (6,42) 'l=',nl(i),' ...'
 42      format(1x,1a2,1i1,1a4)
         lp2=nl(i)+nl(i)+1
         e=ev(i)

         do j=1,nr
           orb(j,i)=orb(j,i)+vctab(j,nl(i))
         end do

         idoflag=1
         ns=1
         call setqmm(i,orb,nl(i),ns,idoflag,vi(1,lp2),zeff,zorig,rel,
     1               nr,r,r2,dl,q0,xm1,xm2,njrcdummy,vi)

         do j=1,nr
           orb(j,i)=0.d0
         end do

       ! you can replace subroutine pseudize with any type of PP
       ! generation you want... however, kleinman-bylanderization
       ! would take more coding...

         call pseudize(i,orb,e,nl(i),xnj(i),no(i),njrc,zeff,vi(1,lp2),
     1                 q0,xm1,xm2,nr,rmin,rmax,r,dr,r2,dl,rel)

         ! WE HAVE GOT THUS FAR...
         no(i)=nl(i)+1
         ruse=0.d0
         xkappa=-1.d0
         call elsolve(i,occ(i),no(i),nl(i),xkappa,xnj(i),
     1                zorig,zeff,ev(i),phe(1,i),vi(1,lp2),
     1                q0,xm1,xm2,nr,r,dr,r2,dl,ruse)
         write (6,*) nl(i),ev(i)
         xntot=xntot+occ(i)

         if (lp2==1) exit  !loop

         do j=1,nr
           vi(j,lp2-1)=vi(j,lp2)
         end do

       end do  ! i

       ! everything is pseudized
       do i=np,nel
         inew=1+i-np
         no (inew)=no (i)
         nl (inew)=nl (i)
         nm (inew)=nm (i)
         xnj(inew)=xnj(i)
         is (inew)=1
         ev (inew)=ev (i)
         occ(inew)=occ(i)

         do j=1,nr
           phe(j,inew)=phe(j,i)
         end do

       end do

         nel=1+nel-np

         do 1212 i=0,7
           xi=i
           do 1212 k=1,nr
             rpower(k,i)=r(k)**xi
 1212        continue

         ! everything is scaled down...ready for unscreening
         xnum=100.d0
         ratio=1.d0
         call getpot(etot,nst,rel,alfa,dl,nr,dr,r,r2,
     1               xntot,phe,ratio,orb,occ,is,
     1               nel,nl,nm,no,xnj,rpower,xnum,etot2,iuflag)

         ! screening effects in pseudo atom computed...
         do k=1,nel
           lp2=nl(k)+nl(k)+1

           do j=1,nr
             vi(j,lp2  )=vi(j,lp2  )-orb(j,k)
             if (lp2>1) vi(j,lp2-1)=vi(j,lp2-1)-orb(j,k)
           end do

         end do

           ! we got past the unscreening...

           do j=1,nr
             vl =     (     vi(j,2)+2.d0*vi(j,3))/3.d0
             vso=2.d0*(     vi(j,3)-     vi(j,2))/3.d0
             vi(j,2)=vso
             vi(j,3)=vl
             vl =     (2.d0*vi(j,4)+3.d0*vi(j,5))/5.d0
             vso=2.d0*(     vi(j,5)-     vi(j,4))/5.d0
             vi(j,4)=vso
             vi(j,5)=vl
 2222        format (5f8.4)
             vl =     (3.d0*vi(j,6)+4.d0*vi(j,7))/7.d0
             vso=2.d0*(     vi(j,7)-     vi(j,6))/7.d0
             vi(j,6)=vso
             vi(j,7)=vl
           end do

           rel=0.d0

           ! got past the spin-orbit jazz
           izuse=abs(vi(nr-2,1)*r(nr-2))+0.5d0
           zuse=izuse

           do k=1,nr
             if (r(k)>rnorm) then
               videal=-zuse/r(k)-corpol/(2.d0*r(k)**4.d0)
               vi(k,1)=videal
               vi(k,3)=videal
               vi(k,5)=videal
               vi(k,7)=videal
               vi(k,2)=0.d0
               vi(k,4)=0.d0
               vi(k,6)=0.d0
             endif

           end do

       ! we got to the end

       return

      end subroutine

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

      end subroutine

c----------------------------------------------------------------------
      real*8 function hb(x,factor)
       implicit real*8 (a-h,o-z)

       if (x>3.d0) hb=0.d0
       if (x<=3.d0) hb=0.01d0**((dsinh(x/factor)/1.1752d0)**2.d0)

       return

      end function
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
 115   idoflag=2
       ns=1
       xkappa=-1.d0
       call setqmm(i,orb,l,ns,idoflag,v,zeff,dummy,rel,
     1              nr,r,r2,dl,q0,xm1,xm2,njrc,dummy)
       call integ(e,l,xkappa,n,nn,jrt,ief,xactual,phi,zeff,v,
     1 q0,xm1,xm2,nr,r,dr,r2,dl,rel)

       if (nn.ne.0) then
         vl=v(1)
         xla=1.d0
       else

         if (xactual>xideal) then
           vh=v(1)
         else
           vl=v(1)
         endif

         xerror=xideal-xactual
         if (abs(xerror)<0.000000001d0) return
         dxdla=0.d0

         do ii=1,jrt
           dxdla=dxdla+dr(ii)*phi(ii)*phi(ii)*hb(r(ii)/rcut,factor)
         end do

         dxdla=2.d0*dxdla/(phi(jrt)*phi(jrt))
         xla=xerror/dxdla
       endif

       vmaybe=v(1)+xla
       if ((vmaybe>vh).or.(vmaybe<vl)) xla=(vl+vh)/2.d0-v(1)

       do ii=1,jrt-1
         v(ii)=v(ii)+xla*hb(r(ii)/rcut,factor)
       end do

       goto 115

      end subroutine

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
       istop=istop-1

       do while(ev<=q0(istop))
         istop=istop-1
       end do

       call integ(ev,l,xkappa,n,nn,istop,ief,xdummy,phi,zeff,v,
     1            q0,xm1,xm2,nr,r,dr,r2,dl,rel)

 50    ! ENTER THE CUTOFF RADIUS, AND FACTOR.
       read (5,*) rcut,factor
       if (rcut<0.d0) then
         xnodefrac=-rcut
         j=istop
 55      j=j-1
         if (phi(j-1)/phi(j)>1.d0) goto 55

         if (n>l+1) then
           k=j
 60        k=k-1
           if (phi(k-1)/phi(k)>0.d0) goto 60
         else
           k=1
         endif

         rcut=r(k)+xnodefrac*(r(j)-r(k))
       endif

       jrc=1.d0+dfloat(nr-1)*log(rcut /rmin)/log(rmax/rmin)
       rcut=r(jrc)
       rtest=2.d0*rcut
       jrt=1.d0+dfloat(nr-1)*log(rtest/rmin)/log(rmax/rmin)
       njrc(lp)=jrt
       rtest=r(jrt)
       switch=phi(jrt)/abs(phi(jrt))

       write (6,92) 'RCUTOFF = ',rcut,'  JRC = ',jrc
       write (6,92) 'RTEST   = ',rtest,  '  JRT = ',jrt
 92    format (1x,1a10,1f8.4,1a8,1i5)
 94    format (1x,2d15.8)

       call integ(ev,l,xkappa,n,nn,jrt,ief,x00,phi,zeff,v,
     1            q0,xm1,xm2,nr,r,dr,r2,dl,rel)

       do ii=1,jrt
         phi(ii)=phi(ii)/phi(jrt)
       end do

       xn00=0.d0

       do ii=1,jrt-1
         xn00=xn00+dr( ii)*phi( ii)*phi( ii)
       end do

       xn00=xn00+dr(jrt)*phi(jrt)*phi(jrt)/2.d0
       de=0.0001d0
       ee=ev+de/2.d0
       call integ(ee,l,xkappa,n,nn,jrt,ief,xp,phi,zeff,v,
     1            q0,xm1,xm2,nr,r,dr,r2,dl,rel)
       ee=ev-de/2.d0
       call integ(ee,l,xkappa,n,nn,jrt,ief,xm,phi,zeff,v,
     1            q0,xm1,xm2,nr,r,dr,r2,dl,rel)
       c00=(xm-xp)/(2.d0*de)
       write (6,94) c00,x00
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

       do ii=1,jrc
         rr=r(ii)
         v(ii)=b0+b2*rr**2.d0+b4*rr**4.d0
       end do

       call fitx0(i,orb,rcut,njrc,ev,l,xj,lp,jrt,x00,phi,zeff,v,
     1  q0,xm1,xm2,nr,r,dr,r2,dl,ruse,factor)

 180   do ii=1,jrt
         phi0(ii)=phi(ii)
         vraw(ii)=v(ii)
       end do

 210   xi0=0.d0
       xi1=0.d0
       xi2=0.d0

       do ii=1,jrt
         f=hb(r(ii)/rcut,factor)
         ph2=dr(ii)*phi0(ii)*phi0(ii)
         xi0=xi0+ph2

         if (ii<=jrt) then
           xi1=xi1+ph2*f
           xi2=xi2+ph2*f*f
         endif

       end do

       ph2=phi0(jrt)*phi0(jrt)
       xi0=xi0/ph2
       xi1=xi1/ph2
       xi2=xi2/ph2
       quant=xi1*xi1+xi2*(c00-xi0)

       if (quant > 0.d0) then
         deltal=(sqrt(xi1*xi1+xi2*(c00-xi0))-xi1)/xi2
       else
         deltal=(c00-xi0)/(2.d0*xi1)
       endif

       write (6,222) 'DELTAL = ',deltal
 222   format (1x,1a9,1f11.8)

 225   do ii=1,jrt
         yl (ii)=phi0(ii)*hb(r(ii)/rcut,factor)
         phi(ii)=phi0(ii)+deltal*yl(ii)

         if (phi(ii)<0.d0) then
           write (6,*) 'BIG TROUBLE!!! CROSS AXIS!!!'
           return
         endif

       end do

       do ii=1,jrt-1
         if ((phi(ii) == 0.).or.(yl(ii) == 0.)) goto 1170
         jj=ii
         if (ii==1) jj=2

         do j=jj-1,jj+1
           rf(2+j-jj)=r(j)
           vf(2+j-jj)=hb(r(j)/rcut,factor)
         end do

         call parabreg(f,fp,fpp,rf,vf)

         do j=jj-1,jj+1
           vf(2+j-jj)=phi0(j)
         end do

         call parabreg(psi,psip,psipp,rf,vf)
         v(ii)=vraw(ii)+
     1        (1.d0-phi0(ii)/phi(ii))*(2.d0*psip/psi*fp/f+fpp/f)/2.d0
       end do

 1170  call fitx0(i,orb,rcut,njrc,ev,l,xj,lp,jrt,x00,phi,zeff,v,
     1            q0,xm1,xm2,nr,r,dr,r2,dl,ruse,factor)
       call integ(ev,l,xkappa,n,nn,jrt,ief,x0,phi,zeff,v,
     1            q0,xm1,xm2,nr,r,dr,r2,dl,ruse)

       do ii=1,jrt
         phi(ii)=phi(ii)/phi(jrt)
       end do

       xn0=0.d0

       do ii=1,jrt-1
         xn0=xn0+dr( ii)*phi( ii)*phi( ii)
       end do

       xn0=xn0+dr(jrt)*phi(jrt)*phi(jrt)/2.d0
       de=0.0001d0
       ee=ev+de/2.d0
       call integ(ee,l,xkappa,n,nn,jrt,ief,xp,phi,zeff,v,
     1            q0,xm1,xm2,nr,r,dr,r2,dl,ruse)
       ee=ev-de/2.d0
       call integ(ee,l,xkappa,n,nn,jrt,ief,xm,phi,zeff,v,
     1            q0,xm1,xm2,nr,r,dr,r2,dl,ruse)
       c0=(xm-xp)/(2.d0*de)
       write (6,94) c0,x0
       write (6,94) xn0

       if (abs(c0-c00)>=0.000000001d0) then
         dqddel=2.*(xi1+deltal*xi2)
         deltal=deltal+(c00-c0)/dqddel
         goto 225
       endif

       write (6,94)  c0,x0
       write (6,*) 'NCPP ACHIEVED !!!'

       return
      end subroutine

c-----------------------------------------------------------------------
      subroutine fourier(nr,r,dr,r2,vi)

       implicit real*8 (a-h,o-z)
       parameter (iorbs=33,iside=600)
       parameter (io2=iorbs*(iorbs+1)/2)
       parameter (ijive=io2*(io2+1)/2)
       parameter (lmax=4,ihmax=20,nrmax=4000,ntmax=10,npmax=60) 
       dimension r(nrmax),dr(nrmax),r2(nrmax),vi(nrmax,7)
       dimension a(nrmax),v1(nrmax),v2(nrmax)

       do l=0,2
         lp2=l+l+1
         dl=log(r(2)/r(1))
         dl1=12.d0*dl
         dl2=12.d0*dl*dl

         do i=1,nr
           a(i)=r(i)*vi(i,lp2)
         end do

         do i=3,nr-2
           al =(-(a(i+2)-a(i-2))+ 8.d0*(a(i+1)-a(i-1)))/dl1

           ar =al/r(i)
           v1(i)=ar
         end do

         open (unit=20+l,status='unknown')

         do ii=1,200
           q=dfloat(ii)/10.d0
           vq=0.d0

           do i=3,nr-2
             vq=vq+dr(i)*dcos(q*r(i))*v1(i)
           end do

           write (20+l,*) q,vq
         end do

         close(unit=1)

       end do

       return

      end subroutine

c-----------------------------------------------------------------------
      subroutine GETILLLS(PIN)

       implicit double precision (A-H,O-Z)
       dimension FA(0:40),SI(0:40),PIN(0:8,0:8,0:16)

       FA(0)=1.D0
       SI(0)=1.D0

       do I=1,32
         FA(I)=Dfloat(I)*FA(I-1)
         SI(I)=-SI(I-1)
       end do

       do L=0,8

         do M=0,8

           do N=M+L,0,-2
             XI=0.D0
             XF=2.D0/2.D0**Dfloat(N+L+M)
             NN=(N+1)/2
             MM=(M+1)/2
             LL=(L+1)/2

             do IA=NN,N
               AF=SI(IA)*FA(IA+IA)/FA(IA)/FA(N-IA)/FA(IA+IA-N)

               do IB=LL,L
                 BF=SI(IB)*FA(IB+IB)/FA(IB)/FA(L-IB)/FA(IB+IB-L)

                 do IC=MM,M
                   XCF=SI(IC)*FA(IC+IC)/FA(IC)/FA(M-IC)/FA(IC+IC-M)
                   XI=XI+XF*AF*BF*XCF/Dfloat(IA*2+IB*2+IC*2-N-L-M+1)
                 end do  !IC

               end do  ! IB

             end do  ! IA

             PIN(L,M,N)=XI

           end do  ! N

         end do  ! M

       end do  ! L

       return

      end subroutine

c-----------------------------------------------------------------------
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

       if (iu<0) then
         ! READ FULL FILENAME
         read (5,52) filename
 52      format (a255)
         iu1=1
         open (unit=iu1,status='unknown',file=filename)
       else
         iu1=iu
         open (unit=iu1,status='unknown')
       endif

       ! define the logarithmic grid
       do i=1,nr
          r(i)=rmin*(rmax/rmin)**(dfloat(i)/dfloat(nr))
       end do

       ! obtain the charge density on the logarithmic grid
       do i=1,nr
         rho(i)=.0

         do ii=1,nel
             rho(i)=rho(i)+occ(ii)*phe(i,ii)**2
         end do

       end do

       ! write output file
       iprint=0
       write(iu1,550)iprint
550    format('RELA'/'RELAT. ATOMIC CHARGE DENSITY'/I1)
       write(iu1,54) rmin,rmax,nr,zorig
54     format (d15.8,d15.8,i5,f5.2)
       nden=nr*(log(rden/rmin)/log(rmax/rmin))
       write(iu1,56) (rho(j),j=1,nr)
56     format (f15.11)
       close (unit=iu1)

       return

      end subroutine

c-----------------------------------------------------------------------
c  exchange correlation routine, via Ceperley-Alder, as parametrized by
c  Perdew and Zunger, Phys. Rev. B 23, 5048.  we use their interpolation
c  between the unpolarized and polarized gas for the correlation part.
      subroutine exchcorr(nst,rel,rr,rh1,rh2,ex,ec,ux1,ux2,uc1,uc2)

       implicit double precision(a-h,o-z)

       trd=1.d0/3.d0
       ft=4.d0/3.d0

       rh=rh1+rh2

       ! if one spin type, average polarization
       if (nst == 1) then
         rh1=rh/2.d0
         rh2=rh/2.d0
       endif

       ! get the n's, and the rs.
       pi=3.14159265358979d0
       fp=4.d0*pi
       xn1=rh1/(rr*fp)
       xn2=rh2/(rr*fp)
       xn=xn1+xn2

       ! effect cutoff, to avoid overflow
       if ((nst == 3) .or. (xn < 0.00000001d0)) then

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

         if (xn1 == 0.d0) then
           fe1=1.d0
           fu1=1.d0
           ex1=0.d0
           ux1=0.d0
         else
           beta=0.028433756d0*xn1**trd
           b2=beta*beta
           eta=sqrt(1.d0+b2)
           xl=log(beta+eta)
           fe1=1.d0-1.5d0*((beta*eta-xl)/b2)**2.d0
           fu1=-0.5d0+1.5d0*xl/beta/eta
           ex1=exchfactor*xn1**trd
           ux1=4.d0*ex1/3.d0
         endif

         if (xn2 == 0.d0) then
           fe2=1.d0
           fu2=1.d0
           ex2=0.d0
           ux2=0.d0
         else
           beta=0.028433756d0*xn2**trd
           b2=beta*beta
           eta=sqrt(1.d0+b2)
           xl=log(beta+eta)
           fe2=1.d0-1.5d0*((beta*eta-xl)/b2)**2.d0
           fu2=-0.5d0+1.5d0*xl/beta/eta
           ex2=exchfactor*xn2**trd
           ux2=4.d0*ex2/3.d0
         endif

         ! these next lines do the Ceperley-Alder correlation

         if (rs >= 1.d0) then

           rootr=sqrt(rs)

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

           xlr=log(rs)
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

        ! if we are nonrelativistic, turn off MacDonald-Vosko correction

         if (rel == 0.d0) then
           fe1=1.d0
           fu1=1.d0
           fe2=1.d0
           fu2=1.d0
         endif

         ! interpolate the correlation energies

         denom=2.d0**ft-2.d0
         f=((1.d0+zeta)**ft+(1.d0-zeta)**ft-2.d0)/denom
         dfdz=ft/denom*((1.d0+zeta)**trd-(1.d0-zeta)**trd)
         ec=ecu+f*(ecp-ecu)
         uc1=ucu+f*(ucp-ucu)+(ecp-ecu)*(1.d0-zeta)*dfdz
         uc2=ucu+f*(ucp-ucu)+(ecp-ecu)*(-1.d0-zeta)*dfdz         

        ! get the final functional and potential

         ex=(xn1*fe1*ex1+xn2*fe2*ex2)/xn
         ux1=fu1*ux1
         ux2=fu2*ux2
         uc1=uc1
         uc2=uc2

       endif

       return

      end subroutine

C---------------------------------------------------------------------       
      function CAVPOT(MTZ_STRING, SLAB_FLAG, ATOMIC_FILE,
     1   CLUSTER_FILE, MUFFTIN_FILE, OUTPUT_FILE)

       character(len=*), intent(IN) :: ATOMIC_FILE, MUFFTIN_FILE
       character(len=*), intent(IN) :: CLUSTER_FILE, OUTPUT_FILE
       character(len=*), intent(IN) :: MTZ_STRING
       integer, intent(IN)          :: SLAB_FLAG
       real                         :: CAVPOT
       PARAMETER (NIEQ=10,NTOT=40)
       dimension SIG(250,NIEQ),RHO(250,NIEQ),VH(250,NIEQ)
       dimension VS(250,NIEQ),VMAD(550,NIEQ),RX(550),RS(550),POT(550)
       dimension RC(3,3),RK(3,NTOT),ZM(NTOT),Z(NIEQ),ZC(NIEQ)
       dimension RMT(NIEQ),JRMT(NIEQ),JRMT2(NIEQ),NRR(NIEQ)
       dimension NCON(NIEQ),NX(NIEQ)
       dimension IA(30,NIEQ),NA(30,NIEQ),AD(30,NIEQ)
       real TITLE(20),NAME(4,NIEQ),WFN,WFN0,WFN1,WFN2,WFN3,SUM
       common /WK/ WK1(250),WK2(250)
       common /WF/ WF2(250,14),WC(14),LC(14)
       data NGRID,MC,PI/250,30,3.1415926536/,WFN0,WFN1,WFN2,WFN3/
     + 4HRELA,4HHERM,4HCLEM,4HPOTE/

       INDEX(X)=20.0*(log(X)+8.8)+2.0

       ! First input channels
       open(unit=4, file=ATOMIC_FILE, status='old')
       open(unit=7, file=CLUSTER_FILE, status='old')

       ! Now output channels
       open(unit=11, file='check.o', status='unknown')
       open(unit=9, file=MUFFTIN_FILE, status='unknown')
       open(unit=10, file=OUTPUT_FILE, status='unknown')

       ! INITIALISATION OF LOUCKS' EXPONENTIAL MESH
       X=-8.8

       do IX=1,NGRID
         RX(IX)=exp(X)
         X=X+0.05
       end do

       read(7,100)TITLE
       write(11,200)TITLE

       ! INPUT OF CRYSTALLOGRAPHIC DATA:
       !     SPA = LATTICE CONSTANT IN A.U.
       !     RC(I,J) = I'TH COORDINATE OF THE J'TH AXIS OF UNIT CELL,
       !               IN UNITS OF SPA
       !     RK(I,J) = I'TH COORDINATE OF THE J'TH ATOM IN UNIT CELL,
       !               IN UNITS OF SPA
       !     NR = NUMBER OF INEQUIVALENT ATOMS IN UNIT CELL
       !          FOR AN ATOM OF TYPE IR?
       !     NRR(IR) = NUMBER IN UNIT CELL
       !     Z(IR) = ATOMIC NUMBER
       !     ZC(IR) = VALENCE CHARGE
       !     RMT(IR) = MUFFIN-TIN RADIUS
       read(7,101)SPA
       read(7,101)((RC(I,J),I=1,3),J=1,3)

       do I=1,3
         do J=1,3
           RC(I,J)=SPA*RC(I,J)
         end do
       end do

       read(7,102)NR

       do IR=1,NR
         do I=1,NGRID
           VH(I,IR)=0.0
           VS(I,IR)=0.0
           VMAD(I,IR)=0.0
           SIG(I,IR)=0.0
           RHO(I,IR)=0.0
         end do
       end do

       VHAR=0.0
       VEX=0.0

       JJ=0
       ZZ=0.0

       do IR=1,NR
         read(7,100)(NAME(I,IR),I=1,4)
         read(7,103)NRR(IR),Z(IR),ZC(IR),RMT(IR)
         ZZ=ZZ+abs(ZC(IR))
         JRMT(IR)=INDEX(RMT(IR))
         N=NRR(IR)

         do J=1,N
           JJ=JJ+1
           ZM(JJ)=ZC(IR)
           read(7,101)(RK(I,JJ),I=1,3)

           do I=1,3
             RK(I,JJ)=SPA*RK(I,JJ)
           end do

         end do

       end do

       ! N = TOTAL NUMBER OF ATOMS IN UNIT CELL
       ! AV = TOTAL VOLUME OF UNIT CELL
       ! OMA = ATOMIC VOLUME
       ! RWS = WIGNER-SEITZ RADIUS
       N=JJ
       RCC1=RC(2,2)*RC(3,3)-RC(3,2)*RC(2,3)
       RCC2=RC(3,2)*RC(1,3)-RC(1,2)*RC(3,3)
       RCC3=RC(1,2)*RC(2,3)-RC(2,2)*RC(1,3)
       AV=abs(RC(1,1)*RCC1+RC(2,1)*RCC2+RC(3,1)*RCC3)
       OMA=AV/float(N)
       RWS=(0.75*OMA/PI)**(1.0/3.0)
       JRWS=INDEX(RWS)
       write(11,201)((RC(I,J),I=1,3),J=1,3)
       write(11,202)AV,OMA,RWS
       JJ=0

       do IR=1,NR
         write(11,203)IR,(NAME(I,IR),I=1,4),NRR(IR)
         INR=NRR(IR)

         do IIR=1,INR
           JJ=JJ+1
           write(11,204)(RK(I,JJ),I=1,3)
         end do

         write(11,205)Z(IR),ZC(IR),RMT(IR)
       end do

       write(11,216)(RX(IX),IX=1,NGRID)

       ! FOR EACH ATOMIC TYPE, READ IN ATOMIC WAVEfunctionS FOR NEUTRAL
       ! ATOM, IN EITHER HERMAN-SKILLMAN OR CLEMENTI FORM, PRODUCING?
       ! RHO = 4*PI*CHARGE DENSITY * RADIUS**2
       MIX=0

       do 11 IR=1,NR
         read(4,100) WFN
         ! OPTION 0)  RELATIVISTIC CHARGE DENSITY INPUT
         if (WFN==WFN0) call RELA(RHO(1,IR),RX,NX(IR),NGRID)
         ! OPTION 1)  HERMAN-SKILLMAN INPUT
         if (WFN==WFN1) call HSIN(RHO(1,IR),RX,NX(IR),NGRID)
         ! OPTION 2)  CLEMENTI INPUT
         if (WFN==WFN2) call CLEMIN(RHO(1,IR),RX,NX(IR),NGRID)
         ! OPTION 3)  POTENTIAL INPUT
         if (WFN==WFN3) goto 14

         ! RHO IS NORMALISED USING TOTAL ELECTRONIC CHARGE ON THE ATOM
         ! CALCULATED BY THE TRAPEZOIDAL RULE
7        NIX=NX(IR)
         MIX=MAX0(NIX,MIX)
         SUM=0.0D0
         W1=0.025*RHO(1,IR)*RX(1)

         do IX=2,NIX
           W2=0.025*RHO(IX,IR)*RX(IX)
           SUM=SUM+W1+W2
           W1=W2
         end do

         ZE=SUM

         ! SOLVE POISSON'S EQUATION:
         ! SIG = COULOMB POTENTIAL
         ! RHO = 4*PI*CHARGE DENSITY*RADIUS SQUARED
         call POISON(RHO(1,IR),Z(IR),NIX,SIG(1,IR))

         X=-8.8

         do IX=1,NIX
           CE=exp(-0.5*X)
           SIG(IX,IR)=CE*(-2.0*Z(IR)*CE+SIG(IX,IR))
           RHO(IX,IR)=RHO(IX,IR)/(RX(IX)**2)
           X=X+0.05
         end do

         write(11,206)(NAME(I,IR),I=1,4),ZE,RX(NIX),NIX
         write(11,207)(SIG(IX,IR),IX=1,NIX)

11     end do

       ! DETAILS OF NEIGHBOURING SHELLS FOR EACH ATOMIC TYPE IR?
       ! NCON(IR) = NUMBER OF SHELLS INCLUDED
       ! IA(J,IR) = ATOMIC TYPE IN J'TH SHELL
       ! NA(J,IR) = NUMBER OF ATOMS IN J'TH SHELL
       ! AD(J,IR) = DISTANCE TO J'TH SHELL
       RMAX=RX(MIX)

       call NBR(IA,NA,AD,NCON,NRR,NR,RC,RK,N,RMAX,MC)
       write(11,208)

       do IR=1,NR
         write(11,209)IR
         NC=NCON(IR)
         IC=(NC-1)/12+1
         KC=0

         do I=1,IC
           JC=KC+1
           KC=MIN0(NC,KC+12)
           write(11,210)(AD(J,IR),J=JC,KC)
           write(11,211)(NA(J,IR),J=JC,KC)
           write(11,212)(IA(J,IR),J=JC,KC)
         end do

       end do

       read(7,102) nform

       ! CALCULATION OF THE MUFFIN-TIN POTENTIAL FOR EACH NEUTRAL
       ! ATOM, FOLLOWING THE MATTHEISS PRESCRIPTION
       ! READ IN ALPHA FOR THE SLATER EXCHANGE TERM
       read(7,101)ALPHA
       write(11,215)ALPHA
       PD=6.0/(PI*PI)

       do IR=1,NR
         JRX=MAX0(JRWS,JRMT(IR))
         ! SUMMING THE POTENTIALS FROM NEUTRAL ATOMS
         ! VH = HARTREE POTENTIAL
         call SUMAX(VH(1,IR),SIG,RX,NX,NCON(IR),IA(1,IR),NA(1,IR),
     +              AD(1,IR),JRX,NGRID,NR)

         ! SUMMING THE CHARGE DENSITY ABOUT EACH ATOMIC TYPE
         ! VS = TOTAL CHARGE DENSITY, THEN SLATER EXCHANGE TERM
         call SUMAX(VS(1,IR),RHO,RX,NX,NCON(IR),IA(1,IR),NA(1,IR),
     +              AD(1,IR),JRX,NGRID,NR)

         do IX=1,JRX
           VS(IX,IR)=-1.5*ALPHA*(PD*VS(IX,IR))**(1.0/3.0)
         end do

       end do

       ! CALCULATE THE MUFFIN-TIN ZERO
       VINT=0.
       read(7,102) NH

       if (NH==0.and.NR==1) then
         call MTZM(VH(1,1),VS(1,1),RX,NGRID,
     +             RMT(1),RWS,JRMT(1),JRWS,VHAR,VEX)
       endif

       if (NH.ne.0) then
         call MTZ(SIG,RHO,RX,NGRID,RMT,NRR,NX,NR,RC,RK,N,
     +            VHAR,VEX,ALPHA,AV,NH)
       endif

       ! Slab or Bulk calculation?
       ! SLAB_FLAG = 1 for Slab or 0 for Bulk (default = 0)
       if (SLAB_FLAG==1) then
         ! Input the MTZ value from the bulk substrate calculation
         if (len_trim(MTZ_STRING)>=1) then
           read(MTZ_STRING,*) esht
         else
           write(*,*) 'Error: MTZ_STRING is invalid slab calculation'
           return
         endif

         esh=esht-(vhar+vex)

       else
         ! check SLAB_FLAG is valid
         if (SLAB_FLAG.ne.0) then
           write(*,*) 'Warning: SLAG_FLAG defaulted to 0'
         endif
         ! If you are interested in adatoms on this substrate
         ! rerun a slab calculation with the adatoms
         ! and use this MTZ value as input when asked
         write(10,*) vhar+vex
         CAVPOT = vhar+vex  ! return value of function
       endif

       goto 16

       ! OPTION 3)  READ IN POTENTIAL OF NEUTRAL ATOM, VH, ON RADIAL
       !            GRID, RX, FOR CORRECTION BY MADELUNG SUMMATION
14     read(4,104) NGRID,(RX(IX),IX=1,NGRID)

       do IR=1,NR
         read(4,104)JRX,(VH(IX,IR),IX=1,JRX)
         JRMT(IR)=JRX
       end do

       ! THE MADELUNG CORRECTION FOR IONIC MATERIALS. SUBROUTINE MAD
       ! COMPUTES THE SPHERICALLY AND SPATIALLY AVERAGED FIELDS FOR
       ! THE LATTICE OF POINT CHARGES ABOUT EACH ATOMIC TYPE
16     if (ZZ .ne. 0) call MAD(VMAD,RX,NGRID,RMT,NRR,JRMT,NR,
     +                         RC,RK,ZM,N,AV)

       ! THE TOTAL MUFFIN-TIN POTENTIAL IS ACCUMULATED INTO SIG,
       ! REFERRED TO THE MUFFIN-TIN ZERO
       VINT=VHAR+VEX
       if (nform == 0) write(9,102) NR

       do IR=1,NR
         write(11,213) (NAME(I,IR),I=1,4), VINT, RMT(IR)
         JRX=JRMT(IR)

         do IX=1,JRX
           VH(IX,IR)=VH(IX,IR)-VHAR
           VS(IX,IR)=VS(IX,IR)-VEX
           SIG(IX,IR)=VH(IX,IR)+VS(IX,IR)+VMAD(IX,IR)
17         write(11,214) RX(IX),VH(IX,IR),VS(IX,IR),
     +                   VMAD(IX,IR),SIG(IX,IR)
         end do

       end do

       ! write output in a format to be read by WILLIAMS phase shift
       ! program (NFORM=1), by CAVLEED phase shift program (NFORM=0), or
       ! by the relativistic phase shift program (NFORM=2)

       ! Also prepare to shift the potential by an amount of the order
       ! of the bulk muffintin zero.
       ! This is needed only if the cluster.i file corresponds to a
       ! surface adsorbate
       if (nform==1) write(9,220) NR

       if (nform==2) then
         ! define German grid RX and save old grid in RS
         RM=60.0
         DX=0.03125
         NMX=421
         RS(1)=RX(1)
         RX(1)=RM*exp(DX*(1-NMX))
         J=1
         RM=exp(DX)

         do while (J<NMX)
           K=J+1
           RS(K)=RX(K)
           RX(K)=RM*RX(J)
           J=K
         end do

       endif

       do 18 IR=1,NR
         JRX=JRMT(IR)
         if (nform==0) then
           write(9,217)(NAME(I,IR),I=1,4)
           write(9,218)Z(IR),RMT(IR),VINT
         elseif (nform==1)then
           write(9,221)Z(IR),RMT(IR)
         else
           ! es=Emin for phase shift calculation (ev)
           ! de=delta E for phase shift calculation (ev)
           ! ue=Emax for phase shift calculation (ev)
           ! lsm=maximum number of phase shifts desired
           es=20.
           de=5.
           ue=300.
           lsm=12
           write(9,217)(NAME(I,IR),I=1,4)
           write(9,"(3D12.4,4X,I3,4X,D12.4)") ES,DE,UE,LSM,VINT

           ! INTERPOLATION TO GRID RX
           do k=1,jrx
             sig(k,IR)=(sig(k,IR)-esh)*rs(k)
           end do

           NMXX=NMX
           call CHGRID(SIG(1,IR),RS,JRX,POT,RX,NMXX)
           IZ=Z(IR)
           write(9,"(I4,F10.6,I4)") IZ,RMT(IR),NMXX
           JRX=NMXX
         endif

         if (nform==0) write(9,102)JRX
         if (nform==1) then
           do IX=1,JRX
             write(9,219)RX(IX),RX(IX)*(SIG(IX,IR)-esh)
           end do
           rneg=-1.
           write(9,219) rneg
         elseif (nform==0) then
           do IX=1,JRX
             write(9,219)RX(IX),(SIG(IX,IR)-esh)
           end do
         else
           write(9,"(5E14.7)")(POT(IX),IX=1,JRX)
         endif

18     continue

       return

100    format(20A4)
101    format(3F8.4)
102    format(I4)
103    format(I4,3F8.4)
104    format(I4/(5E14.5))
200    format(30H1MUFFIN-TIN POTENTIAL PROGRAM?,5X,20A4)
201    format(///18H AXES OF UNIT CELL/(6X,3F8.4))
202    format(18H0UNIT CELL VOLUME?,F15.4/
     + 15H ATOMIC VOLUME?,F18.4/
     + 21H WIGNER-SEITZ RADIUS?,F12.4)
203    format(///5H TYPE,I2,6H ATOM?,2X,4A4/
     + I4,19H ATOMS IN UNIT CELL)
204    format(6X,3F8.4)
205    format(15H0ATOMIC NUMBER?,F15.1/9H VALENCE?,F21.1/
     + 19H MUFFIN-TIN RADIUS?,F14.4)
206    format(///1H ,4A4,19H ELECTRONIC CHARGE?,F12.5/
     + 51H0COULOMB POTENTIAL FOR ISOLATED ATOM, OUT TO RADIUS,
     + F12.5,10X,3HNX?,I4/)
207    format(5(10E12.4/))
208    format(1H1)
209    format(//34H0NEAREST NEIGHBOUR SHELLS FOR TYPE,I2,5H ATOM)
210    format(9H DISTANCE,1X,15(F8.4))
211    format(7H NUMBER,3X,15(I5,3X))
212    format(5H TYPE,5X,15(I5,3X))
213    format(1H1,4A4,5X,33HPOTENTIALS IN RYDBERGS CORRECT TO,
     + 17H MUFFIN-TIN ZERO?,F8.4/19H0MUFFIN-TIN RADIUS?,F8.4//
     + 5X,6HRADIUS,5X,17HHARTREE POTENTIAL,9X,
     + 8HEXCHANGE,4X,19HMADELUNG CORRECTION,5X,
     + 15HTOTAL POTENTIAL)
214    format(F12.5,4E20.6)
215    format(///39H0STATISTICAL EXCHANGE PARAMETER, ALPHA?,F10.4)
216    format(///20H0LOUCKS' RADIAL MESH//5(10F11.5/))
217    format(4A4)
218    format(3F8.4)
219    format(2E14.5)
220    format(10H &NL2 NRR=,i2,5H &end)
221    format(9H &NL16 Z=,f7.4,4H,RT=,f7.4,5H &end)

      end function

C---------------------------------------------------------------------
      subroutine CHGRID(FX,X,NX,FY,Y,NY)
       dimension FX(NX),X(NX),FY(NY),Y(NY)
C     PIECEWISE QUADRATIC INTERPOLATION FROM GRID X TO GRID Y,  BY
C     AITKEN'S DIVIDED DifFERENCE SCHEME.   NX,NY ARE ARRAY dimensionS?
C     NOTE THAT NY IS RESET IN CHGRID

       IY=1

       do 2 IX=3,NX
1        if (IY>NY) goto 3
         YY=Y(IY)
         if (YY>X(IX)) goto 2
         A1=X(IX-2)-YY
         A2=X(IX-1)-YY
         A3=X(IX)-YY
         A12=(FX(IX-2)*A2-FX(IX-1)*A1)/(X(IX-1)-X(IX-2))
         A13=(FX(IX-2)*A3-FX(IX)*A1)/(X(IX)-X(IX-2))
         FY(IY)=(A12*A3-A13*A2)/(X(IX)-X(IX-1))
         if (IY>NY) goto 3
         IY=IY+1
         goto 1
2        continue

3      NY=IY-1

       return

      end subroutine

C---------------------------------------------------------------------         
C ROUTINE FOR INPUT OF WAVEFUNCTIONS IN THE CLEMENTI PARAMETRISED
C FORM, AND CALCULATION OF CHARGE DENSITY ON RADIAL MESH RX
C     RHO = 4*PI*SUM OVER STATES OF (MODULUS(WAVE FN)**2) *
C          RADIUS**2
C     NC = NUMBER OF ATOMIC STATES
C FOR EACH ATOMIC STATE I?
C     LC(I) = ANGULAR MOMENTUM
C     FRAC = FRACTIONAL OCCUPATION
C     WC(I) = NUMBER OF ELECTRONS
C     WFC(IX,I) = WAVEFUNCTION X RADIUS AT GRID POINT IX
      subroutine CLEMIN(RHO,RX,NX,NGRID)

       dimension RHO(NGRID),RX(NGRID)
       real NAME(4),SUM
       common /WK/ EX(20),FAC(20),FNT(20),NT(20)
       common /WF/ WFC(250,14),WC(14),LC(14)

       read(4,100)NAME
       read(4,101)IPRINT
       read(4,101)NC
 
       do 1 IC=1,NC
         do 1 IG=1,NGRID
1          WFC(IG,IC)=0.0

      ! INPUT OF CLEMENTI PARAMETERS
       IC=1
2      read(4,101) NS
       if (NS<=0) goto 8
       read(4,102)(NT(I),EX(I),I=1,NS)

       do 4 J=1,NS
         A=1.0
         B=2.0
         K=NT(J)
         C=float(K)
         KD=K+K

         do 3 I=2,KD
           A=A*B
3          B=B+1.0

4        FNT(J)=exp(-0.5*log(A)+(C+0.5)*log(2.0*EX(J)))

5      read(4,101) LC(IC)
       if (LC(IC)<0) goto 2
       read(4,103)(FAC(J),J=1,NS)
       read(4,103)FRAC
       WC(IC)=2.0*float(2*LC(IC)+1)*FRAC

       do 7 IX=1,NGRID
         SUM=0.0D0

         do 6 K=1,NS
           EXX=EX(K)*RX(IX)
           if (EXX>80.0) goto 6
           SUM=SUM+FAC(K)*FNT(K)*(RX(IX)**(NT(K)))*EXP(-EXX)
6          continue

7        WFC(IX,IC)=SUM

       IC=IC+1
       goto 5

      ! CALCULATION OF CHARGE DENSITY
8      do 10 IX=1,NGRID
         SUM=0.0D0
         do 9 IC=1,NC
9          SUM=SUM+WC(IC)*WFC(IX,IC)*WFC(IX,IC)
         RHO(IX)=SUM

         if (SUM<1.0D-9) goto 11
10     continue

11     NX=IX

       if (IPRINT==0) return

       write(11,200)NAME

       do 12 IC=1,NC
12       write(11,201)LC(IC),(WFC(IX,IC),IX=1,NGRID)

       write(11,202)RX(NX),NX,(RHO(IX),IX=1,NX)
       return

100    format(4A4)
101    format(I4)
102    format(I4,F11.5)
103    format(5F11.5)
200    format(1H1,4A4,32H ATOMIC WAVEfunctionS (CLEMENTI),
     + 9H X RADIUS)
201    format(3H0L?,I3//5(10F11.5/))
202    format(29H0CHARGE DENSITY OUT TO RADIUS,F12.5,10X,
     + 3HNX?,I4//5(10E12.4/))

       return

      end subroutine

C---------------------------------------------------------------------
C  ROUTINE FOR INPUT OF ATOMIC WAVEfunctionS FROM HERMAN-SKILLMAN
C  TABLES, AND CALCULATION OF CHARGE DENSITY ON THE RADIAL MESH RX
C    RHO = 4*PI*SUM OVER STATES OF (MODULUS(WAVE FN)**2) *
C          RADIUS**2
C    NM ? H-S GRID INTERVAL doubleS EVERY NM MESH POINTS
C    NC = NUMBER OF ATOMIC STATES
C FOR EACH ATOMIC STATE I?
C    LC(I) = ANGULAR MOMENTUM
C    FRAC = FRACTIONAL OCCUPATION
C    WC(I) = NUMBER OF ELECTRONS
C    WFC(IX,I) = WAVEfunction X RADIUS AT GRID POINT IX
      subroutine Hsin(RHO,RX,NX,NGRID)

       dimension RHO(NGRID),RX(NGRID)
       real NAME(4),SUM
       common /WK/ RR(250),RS(250)
       common /WF/ WFC(250,14),WC(14),LC(14)

       read(4,100)NAME,Z
       read(4,101)IPRINT
       read(4,101)NM
       read(4,101)NC

       do 1 IG=1,250
         RS(IG)=0.0

         do 1 IC=1,NC
1          WFC(IG,IC)=0.0

      ! INITIALISATION OF HERMAN-SKILLMAN MESH
       DR=0.005*0.88534138/exp(log(Z)/3.0)
       RR(1)=0.0

       do 2 I=2,250
         if (mod(I,NM)==2) DR=DR+DR
2        RR(I)=RR(I-1)+DR

       NS=0

       do 3 IC=1,NC
         read(4,101)LC(IC),N,FRAC
         NS=MAX0(NS,N)
         WC(IC)=2.0*float(2*LC(IC)+1)*FRAC
3        read(4,102)(WFC(IX,IC),IX=1,N)

      ! CALCULATION OF CHARGE DENSITY
       do 5 IX=1,NS
         SUM=0.0D0

         do 4 IC=1,NC
4          SUM=SUM+WC(IC)*WFC(IX,IC)*WFC(IX,IC)

5        RS(IX)=SUM

      ! INTERPOLATION TO GRID RX
       NX=NGRID
       call CHGRID(RS,RR,NS,RHO,RX,NX)
       if (IPRINT==0) return
       write(11,200)NAME,(RR(IX),IX=1,NS)

       do 6 IC=1,NC
6        write(11,201)LC(IC),(WFC(IX,IC),IX=1,NS)

       do 7 IX=1,NX
         if (RHO(IX)<1.0E-9)goto 8
7        continue

8      NX=IX
       write(11,202)RX(NX),NX,(RHO(IX),IX=1,NX)

       return

100    format(4A4/F9.4)
101    format(2I4,F9.4)
102    format(5F9.4)
200    format(1H1,4A4,39H ATOMIC WAVEfunctionS (HERMAN-SKILLMAN),
     + 9H X RADIUS//21H HERMAN-SKILLMAN MESH//5(10F12.5/))
201    format(3H0L?,I3//5(10F11.5/))
202    format(29H0CHARGE DENSITY OUT TO RADIUS,F12.5,10X,
     + 3HNX?,I4//5(10E12.4/))

      end subroutine

C---------------------------------------------------------------------
C  subroutine MAD CALCULATES THE SPHERIcallY AND SPATIALLY AVERAGED
C  FIELDS FROM A LATTICE OF POINT CHARGES, AND TABULATES THEM ON
C  A RADIAL MESH RX, ABOUT EACH ATOMIC TYPE IN THE UNIT CELL
C  ** NB THIS ROUTINE WORKS IN HARTREES, BUT CONVERTS TO RYDBERGS **
C  RC(I,J) = THE I'TH COORDINATE OF THE J'TH AXIS OF THE UNIT CELL
C  RK(I,J) = THE I'TH COORDINATE OF THE J'TH ATOM IN THE UNIT CELL
C  VMAD(J,IR) = THE J'TH TABULATED VALUE OF THE SPHERIcallY AVERAGED
C  POTENTIAL ABOUT A TYPE-IR ATOM
C  ZM(K)=CHARGE ON THE K'TH ATOM
C  RMT(IR) = MUFFIN-TIN RADIUS OF A TYPE-IR ATOM
C  NR = NUMBER OF INEQUIVALENT ATOMS IN THE CELL
C  AV = VOLUME OF UNIT CELL
C  G(I,J) = I'TH COORDINATE OF THE J'TH RECIPROCAL LATTICE VECTOR
C  VMM(IR) = THE INTEGRAL OF THE POTENTIAL ABOUT A TYPE-IR ATOM
C  OUT TO THE MUFFIN-TIN RADIUS
      subroutine MAD(VMAD,RX,NGRID,RMT,NRR,NX,NR,RC,RK,ZM,N,AV)

       dimension VMAD(NGRID,NR),RX(NGRID),RC(3,3),RK(3,N),ZM(N),
     + RMT(NR),NRR(NR),NX(NR)
       common /WK/ G(3,3),VMM(5),FR(5),RA(3),GA(3)
       data PI,TEST/3.1415926536,1.0E-4/

       RAD(A1,A2,A3)=sqrt(A1*A1+A2*A2+A3*A3)


       do 1 IR=1,NR
         FR(IR)=0.0
         do 1 J=1,NGRID
1          VMAD(J,IR)=0.0

      ! THE RECIPROCAL LATTICE IS DEFINED BY THREE VECTORS, G
       ATV=2.0*PI/AV
       G(1,1)=(RC(2,1)*RC(3,2)-RC(3,1)*RC(2,2))*ATV
       G(2,1)=(RC(3,1)*RC(1,2)-RC(1,1)*RC(3,2))*ATV
       G(3,1)=(RC(1,1)*RC(2,2)-RC(2,1)*RC(1,2))*ATV

       G(1,2)=(RC(2,2)*RC(3,3)-RC(3,2)*RC(2,3))*ATV
       G(2,2)=(RC(3,2)*RC(1,3)-RC(1,2)*RC(3,3))*ATV
       G(3,2)=(RC(1,2)*RC(2,3)-RC(2,2)*RC(1,3))*ATV

       G(1,3)=(RC(2,3)*RC(3,1)-RC(3,3)*RC(2,1))*ATV
       G(2,3)=(RC(3,3)*RC(1,1)-RC(1,3)*RC(3,1))*ATV
       G(3,3)=(RC(1,3)*RC(2,1)-RC(2,3)*RC(1,1))*ATV

      ! MAXIMUM VALUE OF RK, AND MINIMUM VALUES OF RC,G - PRIOR TO
      ! CHOOSING THE SEPARATION CONSTANT AL AND LIMITS FOR SUMMATIONS
       RKMAX=0.0

       do 2 J=1,N
2        RKMAX=AMAX1(RKMAX,RAD(RK(1,J),RK(2,J),RK(3,J)))

       RCMIN=1.0E6
       GMIN=1.0E6

       do 3 J=1,3
         RCMIN=AMIN1(RCMIN,RAD(RC(1,J),RC(2,J),RC(3,J)))
3        GMIN=AMIN1(GMIN,RAD(G(1,J),G(2,J),G(3,J)))

      ! AL IS CHOSEN TO GIVE EQUAL NUMBERS OF TERMS IN real AND
      ! RECIPROCAL SPACE SUMMATIONS
       FAC1=TEST*log(TEST)**4
       FAC2=(4.0*PI*RCMIN**4)/(AV*GMIN**4)
       AL=exp(log(FAC1/FAC2)/6.0)
       ITR=1+ifIX((AL*RKMAX-log(TEST))/(AL*RCMIN))
       LIMR=ITR+ITR+1
       FAC1=4.0*PI*AL*AL/(AV*GMIN**4)
       ITG=1+ifIX(EXP(log(FAC1/TEST)/4.0))
       LIMG=ITG+ITG+1
       write(11,200)((G(I,J),I=1,3),J=1,3)
       write(11,201)RCMIN,GMIN,RKMAX,TEST,AL

      ! REAL SPACE SUMMATION
       write(11,202)ITR
      ! THE PREFACTORS FR FROM THE real SPACE SUMMATION ARE CALCULATED
       AS=-float(ITR)-1.0
       AX=AS

       do 5 JX=1,LIMR
         AX=AX+1.0
         AY=AS

         do 5 JY=1,LIMR
           AY=AY+1.0
           AZ=AS

           do 5 JZ=1,LIMR
             AZ=AZ+1.0

             do 4 I=1,3
4              RA(I)=AX*RC(I,1)+AY*RC(I,2)+AZ*RC(I,3)

             do 5 J=1,N
               K=1

               do 5 KR=1,NR
                 R=RAD(RA(1)+RK(1,J)-RK(1,K),RA(2)+RK(2,J)-RK(2,K),
     +             RA(3)+RK(3,J)-RK(3,K))
                 if (R<1.0E-4)goto 5
                 FR(KR)=FR(KR)+ZM(J)*exp(-AL*R)/R
5                K=K+NRR(KR)

       K=1

       do 7 KR=1,NR
         X=RMT(KR)
         A=exp(-AL*X)
         AI1=((1.0-A)/AL-X*A)/AL
         AI2=(X*0.5*(1.0/A+A)-0.5*(1.0/A-A)/AL)/AL/AL
         VMM(KR)=4.0*PI*(ZM(K)*AI1+AI2*FR(KR))
         NIX=NX(KR)

         do 6 J=1,NIX
           X=RX(J)
           A=exp(AL*X)
6          VMAD(J,KR)=FR(KR)*0.5*(A-1.0/A)/(AL*X)+ZM(K)/(A*X)

7        K=K+NRR(KR)

       write(11,203)(VMM(KR),KR=1,NR)

      ! NEXT COMES THE SUMMATION IN RECIPROCAL SPACE
       write(11,204)ITG
       AS=-float(ITG)-1.0
       AX=AS

       do 13 JX=1,LIMG
         AX=AX+1.0
         AY=AS

         do 13 JY=1,LIMG
           AY=AY+1.0
           AZ=AS

           do 13 JZ=1,LIMG
             AZ=AZ+1.0

             do 8 I=1,3
8              GA(I)=AX*G(I,1)+AY*G(I,2)+AZ*G(I,3)

             GM=RAD(GA(1),GA(2),GA(3))
             GS=GM*GM
             FAC1=0.0
             if (GS<1.0E-4) goto 13
             FAC1=4.0*PI*AL*AL/(AV*GS*(GS+AL*AL))
             K=1

             do 12 KR=1,NR
               FAC2=0.0

               do 10 J=1,N
                 GR=0.0

                 do 9 I=1,3
9                  GR=GR+GA(I)*(RK(I,K)-RK(I,J))

10               FAC2=FAC2+cos(GR)*ZM(J)

               X=RMT(KR)
               AI3=(sin(GM*X)/GM-X*cos(GM*X))/GS
               VMM(KR)=VMM(KR)+4.0*PI*AI3*FAC1*FAC2
               NIX=NX(KR)

               do 11 I=1,NIX
                 X=RX(I)
11               VMAD(I,KR)=VMAD(I,KR)+FAC1*FAC2*sin(GM*X)/(GM*X)

12             K=K+NRR(KR)

13           continue

       write(11,203)(VMM(KR),KR=1,NR)

      ! REFER TO MUFFIN-TIN ZERO
       VM=0.0
       AMT=0.0

       do 14 IR=1,NR
         VM=VM+float(NRR(IR))*RMT(IR)**3
14       AMT=AMT+float(NRR(IR))*VMM(IR)

       AMT=AMT/(AV-4.0*PI*VM/3.0)

      ! EXPRESS THE FINAL POTENTIAL IN RYDBERGS
       AMT=-2.0*AMT
       write(11,205)AMT

       do 15 KR=1,NR
         NIX=NX(KR)
         do 15 J=1,NIX
15         VMAD(J,KR)=2.0*VMAD(J,KR)-AMT

       return

200    format(///20H0MADELUNG CORRECTION//
     + 19H0RECIPROCAL LATTICE/(6X,3F8.4))
201    format(7H0RCMIN?,F10.4,10X,5HGMIN?,F10.4,10X,6HRKMAX?,F10.4,
     + 10X,5HTEST?,E12.4/21H SEPARATION CONSTANT?,E12.4)
202    format(21H0real SPACE SUMMATION,11X,4HITR?,I3)
203    format(17H VMM (HARTREES) ?,5E12.4)
204    format(27H0RECIPROCAL SPACE SUMMATION,5X,4HITG?,I3)
205    format(26H0MADELUNG MUFFIN-TIN ZERO?,5E12.4)

      end subroutine

C---------------------------------------------------------------------
C  SUBROUTINE FOR CALCULATION OF THE MUFFIN-TIN ZERO LEVEL?
C  THE AVERAGE VALUE OF THE POTENTIAL BETWEEN THE MUFFIN-TIN
C  SPHERES IN THE UNIT CELL
      subroutine MTZ(SIG,RHO,RX,NGRID,RMT,NRR,NX,NR,
     + RC,RK,N,VHAR,VEX,ALPHA,AV,NH)

       dimension SIG(NGRID,NR),RHO(NGRID,NR),RX(NGRID),RMT(NR),
     + NRR(NR),NX(NR),VG(20),RC(3,3),RK(3,N)
       common /WK/ X(3),RB(3)
       data PI,NG/3.14159265358,20/

      ! GRID REFERENCE FOR RADIUS ON LOUCKS' MESH
       INDEX(y)=20.0*(log(y)+8.8)+1.0
       RAD(A1,A2,A3)=sqrt(A1*A1+A2*A2+A3*A3)

       PD=6.0/PI/PI

       do 12 IG=1,NG
         VG(IG)=0.0
12       continue

       IG=0
       VHAR=0.0
       VEX=0.0
       NPOINT=0
       NINT=0
       DH=1.0/float(NH)

1      AH=DH/2.0

       AX=-AH

       do 7 IX=1,NH
         AX=AX+DH
         AY=-AH

         do 7 IY=1,NH
           AY=AY+DH
           AZ=-AH

           do 7 IZ=1,NH
             AZ=AZ+DH

             do 2 I=1,3
2              X(I)=AX*RC(I,1)+AY*RC(I,2)+AZ*RC(I,3)

             NPOINT=NPOINT+1
             ! GIVES SAMPLE POINT X INSIDE THE UNIT CELL - TEST WHETHER
             ! INTERSTITIAL
             BX=-1.0

             do 4 JX=1,2
               BX=BX+1.0
               BY=-1.0

               do 4 JY=1,2
                 BY=BY+1.0
                 BZ=-1.0

                 do 4 JZ=1,2
                   BZ=BZ+1.0

                   do 3 I=1,3
3                    RB(I)=X(I)-BX*RC(I,1)-BY*RC(I,2)-BZ*RC(I,3)

                   I=0

                   do 4 IR=1,NR
                     INR=NRR(IR)

                     do 4 IIR=1,INR
                       I=I+1
                       XR=RAD(RB(1)-RK(1,I),RB(2)-RK(2,I),RB(3)-RK(3,I))
                       if (XR<RMT(IR)) goto 7

4                      continue

             ! WE HAVE AN INTERSTITIAL POINT
             NINT=NINT+1
             ! SUM COULOMB AND EXCHANGE ENERGIES FROM ATOMS WITHIN 2 UNIT
             ! CELLS AROUND THIS POINT
             SUMC=0.0
             SUME=0.0
             BX=-3.0

             do 6 JX=1,5
               BX=BX+1.0
               BY=-3.0

               do 6 JY=1,5
                 BY=BY+1.0
                 BZ=-3.0

                 do 6 JZ=1,5
                   BZ=BZ+1.0

                   do 5 I=1,3
5                    RB(I)=BX*RC(I,1)+BY*RC(I,2)+BZ*RC(I,3)-X(I)

                   J=0

                   do 6 JR=1,NR
                     JNR=NRR(JR)

                     do 6 JJR=1,JNR
                       J=J+1
                       XR=RAD(RB(1)+RK(1,J),RB(2)+RK(2,J),RB(3)+RK(3,J))
                       J2=INDEX(XR)
                       if (J2>=NX(JR)) goto 6
                       J1=J2-1
                       J3=J2+1
                       X1=RX(J1)
                       X2=RX(J2)
                       X3=RX(J3)
                       TERMC=(XR-X2)*(XR-X3)/(X1-X2)/(X1-X3)*SIG(J1,JR)
     1                      +(XR-X1)*(XR-X3)/(X2-X1)/(X2-X3)*SIG(J2,JR)
     1                      +(XR-X2)*(XR-X1)/(X3-X2)/(X3-X1)*SIG(J3,JR)
                       TERME=(XR-X2)*(XR-X3)/(X1-X2)/(X1-X3)*RHO(J1,JR)
     1                      +(XR-X1)*(XR-X3)/(X2-X1)/(X2-X3)*RHO(J2,JR)
     1                      +(XR-X2)*(XR-X1)/(X3-X2)/(X3-X1)*RHO(J3,JR)
                        SUMC=SUMC+TERMC
                        SUME=SUME+TERME

6                       continue


           if (SUME<=1.E-8) then
             SUME=.0
           else
             SUME=-1.5*ALPHA*(PD*SUME)**(1./3.)
           endif

           VHAR=VHAR+SUMC
           VEX=VEX+SUME
           JG=mod(IG,20)+1
           VG(JG)=VG(JG)+SUMC+SUME
           IG=IG+1

7          continue

       DH=DH/2.0
       NH=NH+NH
       if (NINT==0) goto 1

       ANT=float(NINT)
       VHAR=VHAR/ANT
       VEX=VEX/ANT
       VINT=VHAR+VEX

      ! ESTIMATE STANDARD DEVIATION
       if (NINT < NG) NG=NINT
       NAG=NINT/NG
       AG=float(NAG)

       do 10 IG=1,NG
10       VG(IG)=VG(IG)/AG

       VAR=0.0

       do 11 IG=1,NG
11       VAR=VAR+(VINT-VG(IG))**2

       VAR=sqrt(VAR/float(NG*(NG-1)))
      ! THE CURRENT MONTE-CARLO VOLUME FOR THE INTERSTITIAL REGION
      ! IS VOLC
       VOLC=ANT/float(NPOINT)*AV
      ! VOLT IS THE TRUE VOLUME OF THE REGION BETWEEN MUFFIN-TIN
      ! SPHERES IN THE UNIT CELL
       VM=0.0

       do 8 IR=1,NR
8        VM=VM+float(NRR(IR))*RMT(IR)**3

       VOLT=AV-4.0*PI*VM/3.0

       write(11,200)NINT,NPOINT,NG,NAG,VOLT,VOLC
       write(11,201)VHAR,VEX,VINT,VAR

       return

200    format(///43H0MUFFIN-TIN ZERO CALCULATION, SAMPLING WITH,I6,
     + 20H POINTS FROM GRID OF,I6/24H VARIANCE ESTIMATED FROM,
     + I4,10H GROUPS OF,I5//36H TRUE VOLUME OF INTERSTITIAL REGION?,
     + F11.4,5X,19HMONTE-CARLO VOLUME?,11X,F9.4)
201    format(27H AVERAGE HARTREE POTENTIAL?,6X,F14.5,5X,
     + 27HAVERAGE EXCHANGE POTENTIAL?,F12.5/
     + 17H0MUFFIN-TIN ZERO?,F12.5,10X,19HSTANDARD DEVIATION?,F12.5)

      end subroutine

C---------------------------------------------------------------------
C  SUBROUTINE FOR CALCULATION OF THE MUFFIN-TIN ZERO LEVEL FOR
C  MONOATOMIC CRYSTALS, UsinG A SPHERICAL AVERAGE OF THE POTENTIAL
C  BETWEEN MUFFIN-TIN RADIUS AND WIGNER-SEITZ RADIUS, AS IN EQ 3.31
C  OF LOUCKS, TRANSFORMED TO THE EXPONENTIAL GRID RX?
C                       RX(I)=EXP(-8.8+0.05(I-1))
C  INTEGRATION BY TRAPEZIUM RULE.  JRMT,JRWS ARE GRID POINTS OUTSIDE
C  MUFFIN-TIN RADIUS AND WIGNER-SEITZ RADIUS RESPECTIVELY
      subroutine MTZM(VH,VS,RX,NGRID,RMT,RWS,JRMT,JRWS,VHAR,VEX)

       dimension VH(NGRID),VS(NGRID),RX(NGRID)
       double precision SUMH,SUME

       DX=0.05
       DDX=0.5*DX
       DXX=exp(3.*DX)
       X=log(RX(JRMT)/RMT)
       RDX=X/DX
       XX=RX(JRMT-1)**3
       XXMT=XX*DXX
       SUMH=0.5*X*(RDX*XX*VH(JRMT-1)+(2.-RDX)*XXMT*VH(JRMT))
       SUME=0.5*X*(RDX*XX*VS(JRMT-1)+(2.-RDX)*XXMT*VS(JRMT))
       XX=XXMT
       JRW=JRWS-1

       if (JRMT == JRW) goto 2

       VH1=DDX*XX*VH(JRMT)
       VX1=DDX*XX*VS(JRMT)
       JRM=JRMT+1

       do 1 J=JRM,JRW
         XX=XX*DXX
         VH2=DDX*XX*VH(J)
         VX2=DDX*XX*VS(J)
         SUMH=SUMH+VH1+VH2
         SUME=SUME+VX1+VX2
         VH1=VH2
1        VX1=VX2

2      X=log(RWS/RX(JRW))

       RDX=X/DX
       XXWS=XX*DXX
       SUMH=SUMH+0.5*X*((2.-RDX)*XX*VH(JRW)+RDX*XXWS*VH(JRWS))
       SUME=SUME+0.5*X*((2.-RDX)*XX*VS(JRW)+RDX*XXWS*VS(JRWS))
       C=3./(RWS*RWS*RWS-RMT*RMT*RMT)
       VHAR=C*SUMH
       VEX=C*SUME
       VINT=VHAR+VEX

       write(11,200)VHAR,VEX,VINT

       return

200    format(///37H0MUFFIN-TIN ZERO BY SPHERICAL AVERAGE,/
     + 27H AVERAGE HARTREE POTENTIAL?,6X,F14.5,5X,
     + 27HAVERAGE EXCHANGE POTENTIAL?,F12.5,/
     + 17H0MUFFIN-TIN ZERO?,F12.5)

      end subroutine

C---------------------------------------------------------------------
C  ROUTINE TO SUPPLY NEAREST NEIGHBOUR data FOR ATOMS IN
C  A CRYSTAL STRUCTURE, GIVEN?
C  RC(I,J)? THE I'TH COORDINATE OF THE J'TH AXIS OF THE UNIT CELL
C  RK(I,J)? THE I'TH COORDINATE OF THE J'TH ATOM IN THE UNIT CELL
C  NRR(IR)? THE NUMBER OF TYPE-IR ATOMS IN THE UNIT CELL
C  THE INformatION returnED, FOR A TYPE-IR ATOM, IS
C  NCON(IR)? THE NUMBER OF NEAREST NEIGHBOUR SHELLS OF A TYPE-IR
C  ATOM INCLUDED, OUT TO A DISTANCE OF RMAX, BUT <= MC
C  IA(J,IR)? THE TYPE OF ATOMS IN THE J'TH NEIGHBOURING SHELL
C  NA(J,IR)? THE NUMBER OF ATOMS IN THE J'TH SHELL
C  AD(J,IR)? THE RADIUS OF THE J'TH SHELL
      subroutine NBR(IA,NA,AD,NCON,NRR,NR,RC,RK,N,RMAX,MC)

       dimension IA(MC,NR),NA(MC,NR),AD(MC,NR),NCON(NR),NRR(NR),
     + RC(3,3),RK(3,N)
       common /WK/ RJ(3)

      ! INITIALISATION
       RAD(A1,A2,A3)=sqrt(A1*A1+A2*A2+A3*A3)

       RCMIN=1.0E6

       do 1 I=1,3
1        RCMIN=AMIN1(RCMIN,RAD(RC(1,I),RC(2,I),RC(3,I)))

         do 2 IR=1,NR
           do 2 IC=1,MC
             IA(IC,IR)=0
             NA(IC,IR)=0
2            AD(IC,IR)=1.0E6

      ! SEARCH OVER ADJACENT UNIT CELLS TO INCLUDE MC NEAREST NEIGHBOURS
         ITC=ifIX(RMAX/RCMIN)+1
         LIMC=ITC+ITC+1
         AS=-float(ITC+1)
         AX=AS

         do 10 JX=1,LIMC
           AX=AX+1.0
           AY=AS

           do 10 JY=1,LIMC
             AY=AY+1.0
             AZ=AS

             do 10 JZ=1,LIMC
               AZ=AZ+1.0

               do 3 J=1,3
3                RJ(J)=AX*RC(J,1)+AY*RC(J,2)+AZ*RC(J,3)
                 ! RJ IS CURRENT UNIT CELL ORIGIN.
                 ! FOR EACH ATOM IN THIS UNIT CELL FIND DISPLACEMENT R
                 ! FROM KR-TYPE ATOM IN BASIC UNIT CELL
                 J=0

                 do 10 JR=1,NR
                   JNR=NRR(JR)

                   do 10 JJR=1,JNR
                     J=J+1
                     K=1

                     do 9 KR=1,NR
                    R=RAD(RJ(1)+RK(1,J)-RK(1,K),RJ(2)+RK(2,J)-RK(2,K),
     + RJ(3)+RK(3,J)-RK(3,K))
                       if (R > RMAX) goto 9
                       ! COMPARE R WITH NEAREST NEIGHBOUR DISTANCES
                       ! ALREADY FOUND
                       IC=0
4     IC=IC+1
      if (IC>MC)goto 9
      DR=R-AD(IC,KR)
      if (abs(DR) < 1.0E-4) DR=0.0
      if (DR) 6,5,4
5     if (IA(IC,KR) .ne. JR) goto 4
      NA(IC,KR)=NA(IC,KR)+1
      goto 9
6     if (IC == MC) goto 8
      IIC=IC+1
      do 7 JJC=IIC,MC
        JC=MC+IIC-JJC
        IA(JC,KR)=IA(JC-1,KR)
        NA(JC,KR)=NA(JC-1,KR)
7       AD(JC,KR)=AD(JC-1,KR)
8     IA(IC,KR)=JR
      NA(IC,KR)=1
      AD(IC,KR)=R
9     K=K+NRR(KR)
10     continue
      do 12 IR=1,NR
      NCON(IR)=0
      do 11 IC=1,MC
      if (NA(IC,IR)==0)goto 12
11    NCON(IR)=NCON(IR)+1
12    continue

       return

      end subroutine

C---------------------------------------------------------------------
C TAKEN FROM LOUCKS' BOOK, APPENDIX 1
      subroutine POISON(PSQ,Z,J,W)

       dimension PSQ(J),W(J)
       double precision E(250),F(250),ACC,A,B,C,D,C2

       A=1.0D0-0.0025D0/48.0D0

       ! EQ. A1.11
       B=-2.0D0-0.025D0/48.0D0

       ! EQ. A1.12
       C=0.0025D0/6.0D0
       D=exp(0.025D0)
       C2=-B/A
       E(1)=0.0D0

       ! EQ. A1.29
       F(1)=D

       ! EQ.A1.30
       X=-8.75
       J1=J-1

       do 1 I=2,J1
         ACC=C*exp(0.5*X)*(D*PSQ(I+1)+10.0*PSQ(I)+PSQ(I-1)/D)
         ! EQS. A1.13, A1.6
         F(I)=C2-1.0/F(I-1)
         ! EQ. A1.20
         E(I)=(ACC/A+E(I-1))/F(I)
         ! EQ. A1.21
1        X=X+0.05

       W(J)=2.0*Z*exp(-0.5*X)
       ACC=W(J)

       ! EQ.A1.15
       do 2 I=1,J1
         JC=J-I
         ACC=E(JC)+ACC/F(JC)
2        W(JC)=ACC

       return

      end subroutine

C---------------------------------------------------------------------
C  ROUTINE FOR INPUT OF CHARGE DENSITY FROM RELATIVISTIC ORBITALS
C  (ERIC SHIRLEY PROGRAM), AND CALCULATION OF CHARGE DENSITY ON 
C  THE RADIAL MESH RX
C    RHO = 4*PI*SUM OVER STATES OF (MODULUS(WAVE FN)**2) *
C          RADIUS**2
C    RMIN= minimum radial coordinate defining the logarithmic mesh used
C          in relativistic calculation
C    RMAX= maximum radial coordinate defining the logarithmic mesh used
C          in relativistic calculation
C    NR  = number of points in the mesh
C the mesh is defined as r(i)=rmin*(rmax/rmin)**(dfloat(i)/dfloat(nr))
C FOR EACH ATOMIC STATE I?
      subroutine RELA(RHO,RX,NX,NGRID)

       dimension RHO(NGRID),RX(NGRID)
       real NAME(4),SUM
       common /WK/ RR(2000),RS(2000)

       read(4,100) NAME,IPRINT
       read(4,54) rmin,rmax,nr,z
54     format (d15.8,d15.8,i5,f5.2)

       ! initialization of logarithmic grid
       do 5 i=1,nr
         rr(i)=rmin*(rmax/rmin)**(dfloat(i)/dfloat(nr))
5        continue

       NS=nr
       ! read in charge density
       read(4,56) (rs(j),j=1,nr)
56     format (f15.10)

       ! INTERPOLATION TO GRID RX
       NX=NGRID
       call CHGRID(RS,RR,NS,RHO,RX,NX)
       if (IPRINT==0) return
       write(11,200) NAME,(RR(IX),IX=1,NS)

       do 7 IX=1,NX
         if (RHO(IX)<1.0E-9) goto 8
7        continue

8      NX=IX

       write(11,202)RX(NX),NX,(RHO(IX),IX=1,NX)

       return

100    format(4A4/I4)
102    format(5F9.4)
200    format(1H1,4A4,36H RELAT. WAVEfunctionS (ERIC SHIRLEY),
     + 9H R RADIUS,17H LOGARITHMIC MESH,/(10F12.5/))
201    format(3H0L?,I3//5(10F11.5/))
202    format(29H0CHARGE DENSITY OUT TO RADIUS,F12.5,10X,
     + 3HNX?,I4//5(10E12.4/))

      end subroutine

C---------------------------------------------------------------------
C  ROUTINE TO PERFORM THE SUMMATION OF CONTRIBUTIONS FROM
C  NEIGHBOURING ATOMS (EQ. 3.22,3.26,3.28).  INTEGRATION BY
C  TRAPEZIUM RULE ON RADIAL GRID  RX(I)=EXP(-8.8+0.05(I-1))
      subroutine SUMAX(ACC,CHI,RX,NX,NCON,IA,NA,AD,IMAX,NGRID,NR)

       dimension ACC(NGRID),CHI(NGRID,NR),RX(NGRID),NX(NR),
     + IA(NCON),NA(NCON),AD(NCON)
       double precision SUM

       INDEX(X)=20.*(log(X)+8.8)+2.
       DX=0.05
       DDX=0.5*DX
       DXX=exp(2.*DX)
       IC=IA(1)

       do 1 I=1,IMAX
1        ACC(I)=CHI(I,IC)

       do 4 JA=2,NCON
         IC=IA(JA)
         NIX=NX(IC)

         do 4 I=1,IMAX
           SUM=0.0D0
           X1=abs(RX(I)-AD(JA))
           IX1=INDEX(X1)

           if (IX1>NIX) goto 4

           DX1=log(RX(IX1)/X1)
           RDX1=DX1/DX
           X2=AMIN1((RX(I)+AD(JA)),RX(NIX))
           IX2=MIN0(INDEX(X2),NIX)
           DX2=log(RX(IX2)/X2)
           RDX2=DX2/DX
           XX=RX(IX2-1)**2
           XX1=XX*DXX

           if (IX1==IX2) goto 3

           SUM=SUM+0.5*DX2*((2.-RDX2)*XX*CHI(IX2-1,IC)+
     +                      RDX2*XX1*CHI(IX2,IC))
           XX=RX(IX1-1)**2
           XX1=XX*DXX
           SUM=SUM+0.5*DX1*(RDX1*XX*CHI(IX1-1,IC)+
     +              (2.-RDX1)*XX1*CHI(IX1,IC))
           IX1=IX1+1
           if (IX1==IX2) goto 4
           XX=XX1
           T1=DDX*XX*CHI(IX1,IC)
           IX2=IX2-1

           do 2 IX=IX1,IX2
             XX=XX*DXX
             T2=DDX*XX*CHI(IX,IC)
             SUM=SUM+T1+T2
2            T1=T2

           goto 4

3          SUM=0.5*(DX2-DX1)*((RDX1+RDX2)*XX*CHI(IX1-1,IC)+
     +              (2.-RDX1-RDX2)*XX1*CHI(IX1,IC))

4        ACC(I)=ACC(I)+0.5*SUM*float(NA(JA))/(AD(JA)*RX(I))

       return

      end subroutine

C---------------------------------------------------------------------
C  SUBROUTINE PHSH2CAV
C  POTENTIAL-TO-PHASE-SHIFT CALCULATION(CAVLEED PACKAGE)
C  USES LOUCKS GRID (E.G. AS SUPPLIED BY THE MUFFIN-TIN POTENTIAL
C  PROGRAM).  ENERGIES INPUT IN HARTREES.         
      subroutine PHSH_CAV(MUFFTIN_FILE, PHASOUT_FILE, DATAPH_FILE)

       character(len=*), intent(IN)       :: MUFFTIN_FILE
       character(len=*), intent(IN)       :: PHASOUT_FILE, DATAPH_FILE
       dimension V(250),RX(250),PHS(20)
       real*8 NAME(2),MTZ,delstore(401,15),estore(401)

       ! First input channels
       open(unit=4,file=MUFFTIN_FILE,status='OLD')

       ! Now output channels
       open(unit=6,file='zph.o',status='UNKNOWN')
       open(unit=9,file=PHASOUT_FILE,status='UNKNOWN')
       open(unit=8,file=DATAPH_FILE,status='UNKNOWN')

      ! standard values for phase shifts calculation
       write(8,110)
       emin=1.
       emax=12.
       estep=.25
       ianz=(emax-emin)/estep +1.01
       nl=12
       read(4,103)NR

       do 2  KKK=1,NR
         read(4,100)NAME
         read(4,101)Z,RMT,MTZ
         read(4,103)NTAB
         MTZ=MTZ/2.

         do 19 IX=1,NTAB
19         read(4,219)RX(IX),V(IX)

         write(6,200)NAME,Z,RMT,MTZ
         write(9,181)(NAME(I),I=1,2)
 181     format('NON-RELATIVISTIC PHASE SHIFTS FOR ',2A4)
         write(9,1030) emin,estep,IANZ,nl
 1030    format(2F9.4,2(2X,I3))
         E=EMIN
         ncount=0
1        E=2.*E
         ncount=ncount+1
         call PS(V,RX,NTAB,RMT,E,PHS,NL)
         E=0.5*E
         write(6,201)E,(PHS(L),L=1,NL)
         write(9,1040)E*27.21,(PHS(L),L=1,NL)
 1040    format(F9.4,8F8.4)

         ! store phase shifts
         do 144 kk=1,nl
144        delstore(ncount,kk)=phs(kk)

         estore(ncount)=E
         E=E+ESTEP
         if (E <= EMAX) goto 1

         ! write phase shifts as function of energy for plotting
         do 145 kk=1,nl
           write(8,107) kk-1

           do 146 ii=1,ncount
146          write(8,*) estore(II),delstore(II,kk)

145        write(8,*)

2        continue

       return

71     format(1F7.4,/,10F7.4)
100    format(2A8)
101    format(3F8.4)
102    format(I4/(5E14.5))
103    format(I4)
107    format(3H"L=,i2)  
110    format (11HTitleText: ,8HDELTA(E))
200    format(18H1PHASE SHIFTS FOR ,2A8,2X,15HATOMIC NUMBER?,,F6.1/
     + 19H0MUFFIN-TIN RADIUS?,F8.4,6X,
     + 23H MUFFIN-TIN ZERO LEVEL?,F8.4,9H HARTREES)
201    format(8H0ENERGY?,F8.4,9H HARTREES/(10F12.5))
202    format(2A8/2F8.4)
203    format(10F8.4)
219    format(2E14.5)

       return

      end subroutine

C***********************************************************************
C  SUBROUTINE TO CALCULATE NL PHASE SHIFTS (L=0,NL-1) FOR AN
C  ATOMIC POTENTIAL TABULATED ON THE LOUCKS RADIAL GRID.
C     V?  ATOMIC POTENTIAL (RYDBERGS)
C     RX?  LOUCKS' EXPONENTIAL GRID  RX(I)=EXP(-8.8+0.05(I-1))
C     NGRID?  NUMBER OF TABULATED POINTS
C     RAD?  LIMIT OF INTEGRATION OF SCHRODINGER EQUATION (A.U.)
C     E?  ENERGY (RYDBERGS)
C     PHS?  PHASE SHIFTS FOR L=0 TO NL-1
C  REFERENCE? LOUCKS T L, (1967), A.P.W. METHOD, BENJAMIN, NY.
C***********************************************************************
      subroutine PS(V,RX,NGRID,RAD,E,PHS,NL)

       dimension V(NGRID),RX(NGRID),PHS(NL)
       dimension WF(250),BJ(25),BN(25),XR(10),FR(5)
       data PI,DX,DAC,DB/3.141592653589,0.05,2.083333E-04,
     + 2.083333E-03/

       INDEX(X)=20.*(log(X)+8.8)+2.
       ! TABULATION OF SPHERICAL BESSEL FUNCTIONS IN BJ AND BN
       ES=sqrt(E)
       X=ES*RAD
       Z=X
       LL=NL+1
       call CALCBF(BJ,BN,LL,X)

       ! INTEGRATION OF THE RADIAL SCHRODINGER EQUATION BY THE NUMEROV
       ! METHOD (SEE LOUCKS P56 FF). WF CONTAINS WAVE function X RADIUS,
       ! ON THE LOUCKS GRID
       X1=exp(-8.8)

       do 6 L1=1,NL
         FL=L1-1
         FL2=FL+0.5
         FF=FL2*FL2
         Y1=X1**FL2
         Y2=exp(DX*FL2)*Y1

         write(6,60)FL,Y1,Y2
60       format('0L?',F5.1,5X,'Y1,Y2?',2E14.5)

         GAM1=FF+RX(1)*RX(1)*(V(1)-E)
         GAM2=FF+RX(2)*RX(2)*(V(2)-E)
         WF(1)=Y1*sqrt(RX(1))
         WF(2)=Y2*sqrt(RX(2))

         do 2 IX=3,NGRID
           GAM=FF+RX(IX)*RX(IX)*(V(IX)-E)
           A=1.-DAC*GAM
           B=-2.-DB*GAM2
           C=1.-DAC*GAM1
           YN=-(B*Y2+C*Y1)/A
           WF(IX)=YN*sqrt(RX(IX))
           Y1=Y2
           Y2=YN
           GAM1=GAM2
2          GAM2=GAM

         ! LAGRANGIAN INTERPOLATION FOR WAVEFUNCTION AND DERIVATIVE AT
         ! RADIUS X.  WFN HOLDS WAVEfunction X RADIUS, AND DWFN
         ! DERIVATIVE X RADIUS
         X=RAD
         JR=INDEX(RAD)

         do 3 J=1,5
           XR(J)=RX(JR-5+J)
           XR(J+5)=XR(J)
3          FR(J)=WF(JR-5+J)

         WFN=0.
         DWFN=0.
         A=(X-XR(1))*(X-XR(2))*(X-XR(3))*(X-XR(4))*(X-XR(5))

         do 5 I=1,5
           TERM=A/(X-XR(I))/(XR(I)-XR(I+1))/(XR(I)-XR(I+2))
     +      /(XR(I)-XR(I+3))/(XR(I)-XR(I+4))
           SUM=0.

           do 4 J=1,5
             if (I == J) goto 4
             SUM=SUM+TERM/(X-XR(J))
4            continue

           WFN=WFN+TERM*FR(I)

5          DWFN=DWFN+SUM*FR(I)

         ! LOGARITHMIC DERIVATIVE
         DLOGA=DWFN/WFN-1./RAD

         ! PHASE SHIFTS
         X=ES*RAD
         A=FL*BJ(L1)/X-BJ(L1+1)
         B=FL*BN(L1)/X-BN(L1+1)
         A=ES*A-DLOGA*BJ(L1)
         B=ES*B-DLOGA*BN(L1)
         PHS(L1)=PI/2.
         if (abs(B) > 1.0E-8) PHS(L1)=atan(A/B)
         write(6,78)PHS(L1)
78       format('0PHASE SHIFT?',F10.4)

6        continue

       return

      end subroutine

C***********************************************************************
      subroutine CALCBF(BJ,BN,NL,X)

       dimension BJ(NL),BN(NL)

       if (abs(X) < 1.0E-6) goto 7

       BJ(1)=sin(X)/X
       BN(1)=-cos(X)/X

       if (NL == 1) return

       BJ(2)=(BJ(1)-cos(X))/X
       BN(2)=(BN(1)-sin(X))/X

       if (NL == 2) return

       if (float(NL*(NL+1)) > X*X) goto 2

       ! FORWARD RECURRENCE FOR BJ'S
       FL=3.0

       do 1 L=3,NL
         BJ(L)=FL*BJ(L-1)/X-BJ(L-2)
1        FL=FL+2.

       goto 5

       ! BACKWARD RECURRENCE FOR BJ'S
2      BJ0=BJ(1)

       BJ1=BJ(2)
       NN=MAX0(10,2*NL)
       A=0.
       B=1.
       FL=float(2*NN+1)

       do 3 I=1,NN
         L=NN-I+1
         C=FL*B/X-A
         if (L<=NL)BJ(L)=C
         A=B
         B=C
3        FL=FL-2.

       ! NORMALISATION
       B=BJ0/BJ(1)
       if (abs(BJ0) < 0.01) B=BJ1/BJ(2)

       do 4 L=1,NL
4        BJ(L)=B*BJ(L)

       ! FORWARD RECURRENCE FOR BN'S
5      FL=3.

       do 6 L=3,NL
         BN(L)=FL*BN(L-1)/X-BN(L-2)
6        FL=FL+2.
       return

7      write(6,200)X

       return

200    format(13H0** ARGUMENT?,E12.4,29H TOO SMALL FOR ROUTINE CALCBF)

      end subroutine

C---------------------------------------------------------------------
C  subroutine PHSH_WIL
C---------------------------------------------------------------------
C  A.R. WILLIAMS^ PHASE SHIFT PROGRAM (GIVEN A MUFFIN-TIN POTENTIAL)
      subroutine PHSH_WIL(MUFFTIN_FILE, PHASOUT_FILE, DATAPH_FILE)

       character(len=*), intent(IN)       :: MUFFTIN_FILE
       character(len=*), intent(IN)       :: PHASOUT_FILE, DATAPH_FILE
       real E(401),S(401,15),C(401,15),DEL(15),DELOLD(15)
       real DELL(9),delstore(8,401,15)
       integer TLP1
       common / CM16 / E1, E2, NE, IX,NEUO
       common / CMRV / R(201), V(201, 15), NR, NL, Z
       common / CM5 / Y(30,4), F(30,4), ILST
       namelist / NL2 / IP,NRR

       ! First input channels
       open(unit=5,file=MUFFTIN_FILE,status='OLD')

       ! Now output channels
       open(unit=6,file='zph.o',status='UNKNOWN')
       open(unit=7,file=PHASOUT_FILE,status='UNKNOWN')
       open(unit=8,file=DATAPH_FILE,status='UNKNOWN')

       PI=3.1415926535

       ! READ IP
       ! IP=0: ONLY RADIAL WAVEfunction
       ! IP=1: PHASE SHIFTS IN ADDITION
       ! IP=2: S AND C
       ! IP=3: PRODUCE LOGARITHM OF PHASE SHIFTS
       ! NRR= number of inequivalent atoms for which we want phase shifts
       IP=1
       read(5, NL2)
       write(6, NL2)
       write(8,110)
110    format (11HTitleText: ,8HDELTA(E))

       ! INPUT
       do 2  KKK=1,NRR
         call S16
         TX = 2. * R(NR)/ float(NR - 1)
         DE = (E2 - E1) / float(MAX0(NE - 1, 1))

         do 6 I = 1, NE
           E(I) = E1 + float(I - 1) * DE
           ! RADIAL INTEGRATION
           call S10(E(I))
           T3 = R(NR) * E(I)
           T4 = R(NR) * T3

           do 6 LP1 = 1, NL
             L = LP1 - 1
             TLP1 = 2 * L + 1
             T5 = R(NR) ** LP1
             UT = F(TLP1, ILST)/TX + float(L) * Y(TLP1,ILST) / R(NR)
             T1 = (F44(L,T4) * Y(2*LP1,ILST) + T3 * F44(LP1,T4)
     1              * Y(TLP1,ILST)) * T5
             T2 = (F45(L,T4) * UT - T3 * F45(L - 1,T4)
     1              * Y(TLP1,ILST)) * R(NR) / T5
             S(I, LP1) = T1
6            C(I, LP1) = T2

         IS = 2
         I4 = 9
         if (IP < 1) goto 15

         ! PRODUCE PHASE SHIFTS
         do 8 LP=1,NL
8          DELOLD(LP)=0.0

           do 10 I = 1, NE
             do 11 LP = 1, NL
11             DEL(LP) = atan(-abs(E(I))**(LP-.5)*S(I,LP) / C(I, LP))

             ! REMOVE DISCONTINUITIES BY MULTIPLES OF PI
             do 117 LP=1,NL
               LS=0
111            DELDIF=DEL(LP)-DELOLD(LP)
               if (abs(DELDIF)<0.7) goto 117
               LS=LS+1
               DEL(LP)=DEL(LP)-sign(PI,DELDIF)
               if (LS<5) goto 111
               write (6,115) LP
115   format(36H TOO LARGE CHANGE IN PHASE SHIFT [L=,1I4,
     %20H] sinCE LAST ENERGY ,/,41H DISCONTINUITY BY MULTIPLE OF PI POSS
     %IBLE)
117            DELOLD(LP)=DEL(LP)

             if (NEUO==2) E(I)=0.5*E(I)

             ! PRINT PHASE SHIFTS
             write(6, 12) E(I), (DEL(LP), LP = 1, NL)
12           format(1P8E14.7, /, 14X, 1P7E14.7, /)

             ! WRITE PHASE SHIFTS IN FORMAT USED BY LEED PROGRAM
             ! write(7,71) E(I),(DEL (LP),LP=1,NL)
71           format(1F7.4)
72           format(10F7.4)

             ! store phase shifts
             do 144 kk=1,nl
144            delstore(KKK,I,kk)=del(kk)

             if (IP < 3) goto 10

             do 14 J = 1, 9
               DELL(J) = -4.
               if (DEL(J) < 1.0E-4) goto 14
               DELL(J) = log10(DEL(J))
14             continue

10           continue

           ! write phase shifts as function of energy for plotting
           do 145 kk=1,nl
             write(8,100) kk-1
100          format(3H"L=,i2)

             do 146 ii=1,NE
146            write(8,*) E(II),delstore(kkk,II,kk)

145          write(8,*)


15         continue

2        continue

       write(7,*) 'BE CAREFUL ABOUT THE ORDER OF THE ELEMENTS'

       do 148 ii=1,NE
         write(7,71) E(II)

         do 147 i=1,nrr
           write(7,72)(delstore(i,ii,LP),LP=1,NL)
147        continue

148      continue

       if (IP < 2) goto 22
       do 9 LP1 = 1, NL
         call S41(E, S(1, LP1), NE)
9        call S41(E, C(1, LP1), NE)

22     continue
       return

      end subroutine

C***********************************************************************
C  subroutine S16
C  S16 INPUTS data
C   CS: CORE SHIFT (POSITION OF ZERO OF ENERGY)
C   Z: ATOMIC NUMBER
C   E1,E2: FIRST AND LAST ENERGIES DESIRED (IN RYDBERGS OR
C           HARTREES, CF. NEUI)
C   NE: NUMBER OF ENERGIES DESIRED
C   NL: NUMBER OF PHASE SHIFTS DESIRED (=LMAX+1)
C   NR: NUMBER OF RADIAL GRID POINTS USED IN CALCULATION
C   IX=0: SIGNAL TO STOP
C   IX=1: SIGNAL TO EXPECT A (NEW) POTENTIAL TO BE READ IN
C   RT: MUFFIN-TIN RADIUS (IN BOHR RADII)
C   NEUI,NEUO
C    if =1: RYDBERG UNIT USED FOR INPUT (NEUI) AND OUTPUT (NEUO) OF
C           ENERGIES AND POTENTIAL
C    if =2: HARTREE UNIT (double RYDBERG) USED INSTEAD OF RYDBERG
C           UNIT FOR INPUT (NEUI) AND OUTPUT (NEUO)
C   POTYP=1: RADIAL INPUT AS V(R)
C   POTYP=2: RADIAL INPUT AS R*V(R)
C
C***********************************************************************
      subroutine S16

       common / CM16 / E1, E2, NE, IX,NEUO
       common / CMRV / R, V, NR, NL, Z
       dimension R(201), V(201, 15)
       real RS(200),          ZS(200), ZTT(201)
       dimension FMT(18)
       namelist / NL16 / CS,Z,E1,E2,NE,NL,NR,IX,RT,NEUI,NEUO,POTYP

       ! SET DEFAULT VALUES OF VARIABLES IN NAMELIST /NL16/
       IX=1
       E1 = 4.
       E2 = 24.0
       NE = 30
       NL = 9
       NR = 101
       NEUI=1
       NEUO=2
       POTYP=2
       read(5, NL16)
       CS=0.0

       if (IX < 1) return

       if (NEUI == 1) goto 5

       CS=2.0*CS
       E1=2.0*E1
       E2=2.0*E2

5      continue

       write (6, NL16)

       DRDN2 = (float(NR - 1))*(float(NR - 1)) / RT
       ! READ format USED FOR INPUT OF R VS. V(R) OR R VS. R*V(R)
       ! (V IS ASSUMED POSITIVE)
8      format(18A4)
111    format(2E14.5)

       do 16 I=1,200
         read(5,111) RS(I), ZS(I)
         ! the next lines assume that input potential & cs are negative
         ZS(I)=-ZS(I)

         if ( RS(I) < 0) goto 17

16       continue

17     NRS = I - 1

       if (NEUI == 1) goto 174

       do 172 I=1,NRS
172      ZS(I)=2.0*ZS(I)
174      if (POTYP == 2) goto 178

         do 176 I=1,NRS
176        ZS(I)=(ZS(I)-CS)*RS(I)

         goto 21

178      continue

         do 20 I = 2, NRS
20         ZS(I) = (ZS(I) / RS(I) - CS) * RS(I)
21         IV = 1
           R(1) = 0.
           ZTT(1) = Z + Z

           do 1 I = 2, NR
             R(I) = (float(I - 1)) ** 2 / DRDN2
40           if (R(I) <= RS(IV + 2)) goto 50
             if (IV + 3 >= NRS) goto 50
             IV = IV + 1
             goto 40
50           ZTT(I) = F12(RS(IV), ZS(IV), R(I), 4)

2            do 1 LP1 = 1, NL
1              V(I, LP1) = -ZTT(I) / R(I)

       return

      end subroutine

C***********************************************************************
C  F12 PERFORMS ITERATIVE INTERPOLATION IN A TABLE OF N VALUES OF
C  X AND Y TO FIND THE VALUE OF Y AT Z
C***********************************************************************
      function F12(X, Y, Z, N)

       real X(10), Y(10), W(20)

       W(1) = Y(1)

       do 1 I = 2, N
         W(I) = Y(I)
         U = Z - X(I)
         IP1 = I + 1

         do 1 J = 2, I
           K = IP1 - J
1          W(K) = W(K + 1) + U * (W(K) - W(K + 1)) / (X(K) - X(I))

       F12 = W(1)

      return

      end function

C***********************************************************************
C  S5 -- HAMMING^S METHOD FOR THE INTEGRATION OF SYSTEMS OF FIRST
C  ORDER DifFERENTIAL EQUATIONS
C***********************************************************************
      subroutine S5(E)

       real EEST(30), VME(15)
       common / CMRV / R(201), V(201, 15), NR, NL, Z
       common / CM5 / Y(30, 4), F(30, 4), IP1

       NJ = 2 * NL

       do 5 J = 1, NJ
5        EEST(J) = 0.

       do 4 I = 5, NR
         do 6 LP1 = 1, NL
6          VME(LP1) = (V(I, LP1) - E) * R(I)

         T1 = 2. / float(I - 1)
         IP1 = mod(I - 1, 4) + 1
         IM2 = mod(IP1, 4) + 1
         IM1 = mod(IM2, 4) + 1
         IP0 = mod(IM1, 4) + 1

         do 1 J = 1, NJ
           F(J, IM2) = Y(J, IP1) + (2. * (F(J, IP0) + F(J, IM2)) -
     1                  F(J, IM1)) / 0.75
1          Y(J, IP1) = F(J, IM2) - 0.925619835 * EEST(J)

         do 2 J = 1, NJ, 2
           JP1 = J + 1
           LP1 = JP1 / 2
           FLP1 = LP1
           F(J, IP1) = (FLP1 * Y(J, IP1) + R(I) * Y(JP1, IP1)) * T1
2          F(JP1, IP1) = (VME(LP1) * Y(J, IP1) - FLP1 * Y(JP1, IP1))*T1

        do 3 J = 1, NJ
          Y(J, IP1) = Y(J, IP0)+(Y(J, IP0)-Y(J, IM2) + 3.*(F(J, IP1)
     1 + 2. * F(J, IP0) - F(J, IM1))) / 8.
          EEST(J) = F(J, IM2) - Y(J, IP1)
3         Y(J, IP1) = Y(J, IP1) + .743801653E-1 * EEST(J)

        do 4 J = 1, NJ, 2
          JP1 = J + 1
          LP1 = JP1 / 2
          FLP1 = LP1
          F(J, IP1) = (FLP1 * Y(J, IP1) + R(I) * Y(JP1, IP1)) * T1
4         F(JP1, IP1) = (VME(LP1) * Y(J, IP1) - FLP1 * Y(JP1, IP1))*T1

       return

      end subroutine
C***********************************************************************
C  S10   POWER SERIES EXPANSION OF THE SOLUTION ABOUT THE ORIGIN
C        AND RADIAL INTEGRATION IN S5
C***********************************************************************
      subroutine S10(E)

       integer TLP1
       real A(10), B(10), TR(4)
       common / CMRV / R(201), V(201, 15), NR, NL, Z
       common / CM5 / Y(30, 4) , F(30, 4) , ILST

       NI = 2 * NL
       TZ = 2. * Z
       A(1) = 1.

       do 5 I = 1, NI, 2
         LP1 = (I + 1) / 2
         TLP1 = 2 * LP1
         EP = E - V(4, LP1) - TZ/ R(4)
         Y(I, 1) = 0.
         Y(I + 1, 1) = 0.
         A(1) = A(1) / float(2 * LP1 - 1)
         B(1) = - Z * A(1) / float(LP1)

         do 1 J = 2, 4
           TR(J) = R(J) ** LP1
           Y(I, J) = A(1) * TR(J)
1          Y(I + 1, J) = B(1) * TR(J)

         do 3 K = 1, 9
           A(K + 1) = B(K) / float(K)
           B(K + 1) = -(EP * A(K) + TZ * A(K + 1)) / float(TLP1 + K)

           do 2 J = 2, 4
             TR(J) = TR(J) * R(J)
             Y(I, J) = Y(I, J) + TR(J) * A(K + 1)
2            Y(I + 1, J) = Y(I + 1, J) + TR(J) * B(K + 1)

           if (abs(TR(4) * A(K + 1) / Y(I, 4)) < 1.0E-4) goto 5

3          continue

         write (6, 4) E, LP1, R(4), (A(K), K = 1, 10)
4        format(1PE10.2, I10, 11E10.2)

5        continue

       do 6 J = 2, 4
         T1 = 2. / float(J - 1)

         do 6 I = 1, NI, 2
           IP1 = I + 1
           LP1 = IP1 / 2
           FLP1 = LP1
           F(I, J) = (FLP1 * Y(I, J) + R(J) * Y(IP1, J)) * T1
6          F(IP1,J)=((V(J, LP1)-E) * R(J) * Y(I, J) - FLP1 * Y(IP1, J))
     1              * T1

       call S5(E)

       return

      end subroutine

C***********************************************************************
C  F44  EVALUATES THE SPECIAL VERSION OF THE SPHERICAL BESSEL FUNCT.
C***********************************************************************
      function F44(L, X)

       real S(20)

       JS = L + L + 1

       if (abs(X / float(JS)) > 10.) goto 5

       FI = 1.

       if (L < 1) goto 2

       do 1 K = 3, JS, 2
1        FI = FI * float(K)

2      T1 = 1. / FI
       DT = 1.
       T = 1.
       I = 0

       do 3 K = 1, 100
         I = I + 2
         DT = -DT * X / float(I * (I + JS))
         T = T + DT

         if (abs(DT) < 1.E-8) goto 4

3       continue

4      T1 = T1 * T
       F44 = T1

       return

5      T = sqrt(abs(X))

       if (X < 0.) goto 9

       S(2) = sin(T) / T

       if (L > 0) goto 6

11     F44 = S(2)

       return

6      S(1) = cos(T)

       goto 10

9      S(2) = sinh(T) / T

       if (L < 1) goto 11

       S(1) = cosh(T)

10     IS = L + 2

       do 7 I = 3, IS
7        S(I) = (S(I - 1) * float(2 * I - 5) - S(I - 2)) / X

       F44 = S(IS)

       return

      end function

C***********************************************************************
C  F45  EVALUATES SPECIAL VERSION OF THE SPHERICAL NEUMANN function
C***********************************************************************
      function F45(L, X)

       real S(20)

       if (L < 0) goto 8

       LP1 = L + 1
       JS = L + L + 1

       if (abs(X / float(JS)) > 10.) goto 5

       FI = 1.

       if (L < 1) goto 2

       do 1 K = 3, JS, 2
1        FI = FI * float(K)

2      T1 = FI / float(JS)
       DT = 1.
       T = 1.
       I = 0

       do 3 K = 1, 100
         I = I + 2
         DT = -DT * X / float(I * (I - JS))
         T = T + DT

         if (abs(DT) < 1.E-8) goto 4

3        continue

4      T1 = T1 * T
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

       do 7 I = 3, IS
7        S(I) = S(I - 1) * float(2 * I - 5) - X * S(I - 2)

       F45 = S(IS)

       return

8      F45 = - F44(L+1, X)

       return

      end function

C***********************************************************************
C  S41 PLOTS Y AGAINST X
C***********************************************************************
      subroutine S41(X, Y, N)

       dimension X(100), Y(100), P(97)
       data B, C, O, D / 1H , 1H*, 1H0, 1HI /

       Y1 = 0
       Y2 = 0

       do 1 I = 1, N
         Y1 = AMIN1(Y1, Y(I))
1        Y2 = AMAX1(Y2, Y(I))

       do 2 I = 1, 97
2        P(I) = B

       T = 96/ (Y2 - Y1)
       J0 = -Y1 * T + 1.5
       P(J0) = 0

       if (N >= 30) write(6, 3) P
3      format(1H1, 34X, 97A1, //)
       if (N < 30) write(6, 6) P
6      format(////, 35X, 97A1, //)

       P(J0) = D

       do 5 I = 1, N
         J = T * (Y(I) - Y1) + 1.5
         P(J) = C
         write(6, 4) X(I), Y(I), P
4        format(1H , 1P2E16.6, 2X, 97A1)
         P(J) = B
5        P(J0) = D

       return

      end subroutine

C---------------------------------------------------------------------
C  subroutine PHSH_REL
C---------------------------------------------------------------------
      subroutine PHSH_REL(MUFFTIN_FILE, INPDAT_FILE,
     1                    PHASOUT_FILE, DATAPH_FILE)

       implicit double precision (A-H,O-Z)
       character(len=*), intent(in)       :: MUFFTIN_FILE, INPDAT_FILE
       character(len=*), intent(in)       :: PHASOUT_FILE, DATAPH_FILE
       character OPT*3,OPT1*3,OPTS*3,ANAME*2,AN*30,BDATA*28
       character SUB*3,RECORD*3,TL*1,SL*1,SS1*6,SS2*6,WRD*6
       character AMS(5)*4
       real JF,JFS,JFD,JFU
       real NAME(4)
       dimension JF(250,18),ENERG (250)
       common/ZZZZ/ZP(340),VS,IPT,JRI
       common /Z/ RMAXI
       dimension Adata(7)
       data ZERO/0.0D0/,ONE/1.0D0/,TWO/2.0D0/,ANINE/9.0D0/,HALF/.5D0/
       data ZILCH/1.0D-4/,TOL/.005D+0/,DES/.025D+0/
       data AMS/'NC= ','L= ',' ES=',' DE=','ID= '/
       data TL/'L'/,SL/'S'/,SS1/'NOSPIN'/,SS2/' SPIN '/
       data SUB/'SUB'/,RECORD/'NOS'/
  1    format (3D12.4,4X,I3,4X,D12.4)
  2    format (5E14.6)
  3    format    (///,10X,A6,////,10X,18HMUFFIN-TIN RADIUS=,
     1F10.6,4X,14HCONSTANT POT.=,F10.6,//,10X,'ATOMIC data SET FOR Z=',
     2I3,' AND L=',I3,4X,'     ',A2,//,14X,5HE(EV),12X,7HL - 0.5,13X,
     37HL + 0.5,12X,10HS-AVERAGED,/)
  4    format (10X,F10.6,3F20.8)
  8    format (F9.4,8F8.5)
  9    format (1X,A1,I3,3D15.7,22X,A3)
 12    format (5D14.6)
 13    format (1H1,//,T61,'INPUT data',//)
 15    format (A28,T29,A2,T35,A30,T65,F10.7,T76,A3)
 17    format (1H1,/,10X,'RELATIVISTIC PHASE SHIFTS FOR ',A30,//,10X,'EX
     1CA =',F10.6,4X,'EXCB =',F10.6,4X,'EXCO =',F10.6//,10X,'LATTICE CON
     2STANT =',F10.6,' ,',F10.6,' ,',F10.6)
 18    format(F10.4,F9.4,2I5)

       open(unit=4,file=INPDAT_FILE,status='unknown')
       open(unit=5,file=MUFFTIN_FILE,status='old')
       open(unit=7,file=PHASOUT_FILE,status='unknown')
       open(unit=8,file=DATAPH_FILE,status='unknown')

       PI=4.D0*dataN(1.D0)
       PI2=0.5D0*PI
       write (4,13)
10     read(5,217,end=999,err=999) (NAME(I),I=1,4)
217    format(4A4)
       read(5,1)ES,DE,UE,LSM,VC

       ! nl is the number of plotted phase shifts
       nl=8
       write (4,11) ES,DE,UE,OPT,OPT1,LSM
 11    format (3D12.4,4X,2A3,I3)
 75    format (A28,A2,4X,A30,F10.7,1X,A3)
       read (5,16) NZ,Adata(1),JRI,ALC,BLC,CLC,EXCA,EXCB,EXCO
 16    format (I4,F10.6,I4,T21,6F10.6)
       write (4,76) NZ,Adata(1),JRI,ALC,BLC,CLC,EXCA,EXCB,EXCO
 76    format (I4,F10.6,I4,2X,6F10.6)

       VS=0.
       if (OPTS == SUB) VS=VC
       if ((JRI <= 0) .or. (JRI > 340)) goto 999
       read(5,2) (ZP(J),J=1,JRI)
       write (4,12) (ZP(J),J=1,JRI)
       RHOZ=-0.90306514D+01
       DELRHO=0.3125D-01
       RM=DEXP(RHOZ)
       XRX=DEXP(DELRHO)

       do 19 J=1,JRI
         if (ZP(J) < 0.0) ZP(J)=-ZP(J)
         if (J == JRI) goto 21
 19      RM=XRX*RM

 21    continue

       if (DE > ZERO)  goto 20
       ES=-HALF
       DE=DES
       UE=ONE

  20   N=(UE-ES)/DE+HALF

       N=N+1
       write(07,181)(NAME(I),I=1,4)

 181   format('RELATIVISTIC PHASE SHIFTS FOR ',4A4)
       write(07,18) ES,DE,N,LSM

       ES=ES/13.6
       DE=DE/13.6
       UE=UE/13.6
       if (N > 250)  N=250
       L=0
       E=ES
       IPT=2
       if (OPT == RECORD) goto 23
       WRD=SS2
       goto 24

 23    IPT=-2

       WRD=SS1

 24    continue

       KAP=-1
       L=1

       do 30 J=1,N
         DXAZ=0.0D0
         TTR=DLGKAP(E,KAP)/(12.5663706)
         call SBFIT(TTR,E,L-1,RMAXI,JFS)
         JF(J,L)=JFS
         E=E*13.6
         ENERG(J)=E
         E=E/13.6
         E=E+DE
  30     continue

  40   if (L > LSM) goto 80

       KAP=-(L+1)
       LIND=L+1
       E=ES
       LVOR=0

       do 50 J=1,N
         DLU=DLGKAP (E,KAP)/(12.5663706)
         DLD=DLGKAP(E,L)/(12.5663706)
         call SBFIT(DLD,E,L,RMAXI,JFD)
         call SBFIT(DLU,E,L,RMAXI,JFU)
         LK=0
         ZFDIFF=-(JFD-JFU)*LVOR
         if (ZFDIFF > 2.5) LK=L
         if (ZFDIFF < -2.5) LK=L+1
         JFS=L*JFD-KAP*JFU+LVOR*LK*PI
         JFS=JFS/(2*L+1)
         if (JFS >PI2) JFS=JFS-PI2*2.
         if (JFS < -PI2) JFS=JFS+PI2*2.
         JF(J,LIND)=JFS
         if (LK == 0) LVOR=sign(1.,JFS)
         E=E+DE
  50     continue

       L=L+1

       goto 40

  80   do 90 I=1,N
         LSM1=LSM+1
         write(7,8) ENERG(I),(JF(I,L),L=1,LSM1)

  90    continue

       do 145 kk=1,nl
         write(8,100) kk-1

         do 146 ii=1,N
146        write(8,*) ENERG(II),JF(II,kk)

145      write(8,*)

100    format(3H"L=,i2)  
       ES=ES*13.6
       DE=DE*13.6
       UE=UE*13.6
       goto 10

 999   continue

       write (4,900)

 900   format (//,T57,'end OF INPUT data')

       return

      end subroutine

C***********************************************************************
C DLGKAP CALCULATES THE LOGRITHMIC DERIVATIVE OF THE LARGE
C COMPONENT UsinG THE PROCEDURE DESCRIBED BY LOUCKS IN APPENDIX 7.
C THE SMALL MULTIPLICATIVE FACTOR IS INCLUDED.
C POTENTIAL  IN THE FORM OF 2ZP IS TO BE PASSED IN THE common /ZZZZ/
C THE RADIAL functionS ARE MADE AVAILABLE IN THE common /RADFUN/
C WABER MESH (XNOT=-9.03065 AND DX=1/32) IS USED
C JRI IS THE MESH POINT OF THE APW SPHERE RADIUS
C E IS THE ENERGY TO BE USED (IN RYDBERGS)
C 4 PI R**2 INSERTED NOW. FOR COMPOUNDS ONLY.
C***********************************************************************
      function DLGKAP (E,KAPPA)

       implicit double precision(A-H,O-Z)
       dimension POT(340),U(340),W(340),UP(340),WP(340),SXK(4),SXM(4)
       common /Z/ T
       common /RADFUN/ U,W
       common/ZZZZ/POT,VCZ,IPT,JRI
       data USTART/1.D-25/,ZILCH/1.D-30/,TEST/1.D+6/,XS/-9.03065133D+00/,
     1 DX/3.125D-2/,C/2.740746D+2/,CIN/1.3312581146D-5/,HF/.5D+0/,
     2 TH/.3333333333D+0/,T2/2.D+0/,T7/7.D+0/,T11/11.D+0/,T12/12.D+0/,
     3 T14/14.D+0/,T26/26.D+0/,T32/32.D+0/,ZERO/.1D+0/
 83    format (10H HARD TEST,D20.8, I5, 4D20.8 )

       !SET UP FOR RELATIVISTIC OR NO RELATIVISTIC EFFECT
       if (IPT > 0) goto 88
       CIN=0.0D00
 88    continue

       ! SET UP STARTING VALUES
       DX2=HF*DX
       X20=(.3)*DX
       XMFT=(4.4444444444444D-2)*DX
       TS=DEXP(XS)
       TDX=DEXP(DX)
       HOC=(VCZ*TS+POT(1))/C
       XK=KAPPA
 49    U(1)=USTART
       if (abs(HOC/XK) > 0.05) goto 6
       P=(XK+abs(XK))/HOC - HF*HOC/abs(XK)
       goto 7
 6     P=(XK+Dsqrt(XK*XK-HOC*HOC))/HOC

 7     TC=(E+VCZ)*TS+POT(1)

       VC=CIN*TC
       W(1)=C*P*USTART

       ! START RUNGE-KUTTE PROCEDURE
 11    X = XS

       N=1

 25    IK= 0

       NP1=N+1
       XC=X
       BGC=POT(N)
       WC= W(N)
       UC= U(N)

 20    IK=IK+1

       T=DEXP(XC)
       TC=(E+VCZ)*T+BGC
       VC=CIN*TC

 12    SXK(IK)=DX2*(WC*(VC+T)-XK*UC)

       SXM(IK)=DX2*(XK*WC-TC*UC)

 15    goto(16,17,18,19),IK
 16    XC=XC+DX2
       UC=UC+SXK(1)
       WC=WC+SXM(1)
       BGC=HF*(BGC+POT(NP1))
       goto 20
 17    UC=UC+SXK(2)-SXK(1)
       WC=WC+SXM(2)-SXM(1)
       goto 20
 18    XC=XC+DX2
       UC=UC+T2*SXK(3)-SXK(2)
       WC=WC+T2*SXM(3)-SXM(2)
       BGC=POT(NP1)
       goto 20
 19    W(NP1)=W(N)+(SXM(1)+SXM(4)+T2*(SXM(2)+SXM(3)))*TH
       U(NP1)=U(N)+(SXK(1)+SXK(4)+T2*(SXK(2)+SXK(3)))*TH
       UP(NP1)=(VC+T)*W(NP1)-XK*U(NP1)
       WP(NP1)=XK*W(NP1)-TC*U(NP1)
 24    X=X+DX
       N=NP1
       if (N<6) goto 25

       ! END OF STARTING INTEGRATION.  BEGIN MILNE PROCEDURE.
       T=DEXP(X)
   26  T=T*TDX
       NP1=N+1
       NM1=N-1
       NM2=N-2
       NM3=N-3
       NM4=N-4
       NM5=N-5
       TC=(E+VCZ)*T+POT(NP1)
       VC=CIN*TC
  27   UNP=U(NM5)+X20*(T11*(UP(N)+UP(NM4))+T26*UP(NM2)
     +                                 -T14*(UP(NM1)+UP(NM3)))
       WNP=W(NM5)+X20*(T11*(WP(N)+WP(NM4))+T26*WP(NM2)
     +                                      -T14*(WP(NM1)+WP(NM3)))
       NIT=0
  33   UP(NP1)=(VC+T)*WNP-XK*UNP
       WP(NP1)=XK*WNP-TC*UNP
       UNP2=U(NM3) + (T7*(UP(NP1) + UP(NM3)) + T32*(UP(NM2) + UP(N))
     +             + T12 * UP(NM1)) * XMFT
       WNP2=W(NM3) + (T7*(WP(NP1) + WP(NM3)) + T32*(WP(NM2) + WP(N))
     +             + T12 * WP(NM1)) * XMFT

       ! COMPARE PREDICTOR WITH CORRECTOR
       if (abs(TEST*(UNP2 -UNP)) > abs(UNP2)) goto 31
       if (abs(TEST*(WNP2-WNP)) <= abs(WNP2)) goto 32
  31   if (NIT < 5) goto 81
       goto 32

  81   NIT=NIT+1

       WNP=WNP2
       UNP=UNP2
       goto 33

  32   W(NP1)=WNP2

       U(NP1)=UNP2
       N=NP1

       if (N < JRI) goto 26

       ! END OF MILNE PROCEDURE
       if ( abs(U(JRI)) > ZILCH) goto 37
       U(JRI)=sign(ZILCH,U(JRI))

  37   P=(T+VC)/T

       WNP=P*W(JRI)/U(JRI)
       UNP=WNP-(KAPPA+1)/T
       DLGKAP=(12.5663706)*T*T*UNP

   46  return

      end function

C***********************************************************************
      subroutine SBFIT(T,E,L,R,JFS)

       implicit double precision(E,R,T)
       real JFS,KAPPA

       SE=SNGL(E)
       SR=SNGL(R)
       ST=SNGL(T)
       KAPPA=sqrt(SE)
       X=KAPPA*SR
       BJ1=sin(X)/X
       BN1=-cos(X)/X
       BJ2=BJ1/X + BN1
       BN2=BN1/X - BJ1

       if (L > 0) then
         DL=ST/(SR*SR)
       else
         LS=1
       endif

       LS=1
       LS=LS+1
       BJT=(2*LS-1)*BJ2/X - BJ1
       BNT=(2*LS-1)*BN2/X - BN1
       BJ1=BJ2
       BJ2=BJT
       BN1=BN2
       BN2=BNT

       if (L+1-LS > 0) then
         LS=LS+1
       else
         DL=ST/(SR*SR)
       endif

       DL=ST / (SR*SR)
       DL=DL - L/SR
       AN=DL*BJ1 + KAPPA*BJ2
       AD=DL*BN1 + KAPPA*BN2
       JFS=3.141592654/2.0

       if (abs(AD)-1.0E-8 > 0) JFS=atan(AN/AD)

       return

      end subroutine
