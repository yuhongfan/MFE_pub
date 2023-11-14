program getxi
   implicit none
   real, parameter :: pi   = 3.1415926535898D0
   real, parameter :: tiny = 1.0D-99, huge = 1.0D+99

   character*150 gridfile, xifile
   integer :: nt1,nt2
   namelist /input/ gridfile, xifile, nt1, nt2
   integer :: itsub
   character*4 idout

   integer ::       in, jn, kn
   integer ::       is  , js  , ks  , ie  , je  , ke
   integer ::       ism1, jsm1, ksm1, ism2, jsm2, ksm2
   integer ::       isp1, jsp1, ksp1, isp2, jsp2, ksp2
   integer ::       iem1, jem1, kem1, iem2, jem2, kem2
   integer ::       iep1, jep1, kep1, iep2, jep2, kep2
   integer ::       iep3, jep3, kep3
   integer ::       ione, jone, kone
!
   real, allocatable :: x1a(:), x2a(:), x3a(:)
   real, allocatable :: x1b(:), x2b(:), x3b(:)
   real, allocatable :: xi(:,:)
!
   real, parameter :: lcell=0.0025, tcycle=2., dt=0.1
   real :: th1, th2, ph1, ph2
   integer :: iseedx,iseedy,iseedt,i,j,k,np,icount,in_in,jn_in,kn_in
   real, allocatable :: x(:),y(:),tinit(:)
   integer, allocatable :: icycle(:)
   real :: ran2,ph0,delph,area,xp,yp,rr,fcell,time
   real :: rrmin,rrth1,rrth2,rrph1,rrph2
!---------------------------------------------------------

!  read inparam file
   write(6,*) 'read inparam'
   open(unit=20,file="inparam_getxi", &
     status="old",position="rewind")
     read(unit=20, nml=input)
     close(unit=20)
   write(6,*) 'read in inparam'
   write(6,*) 'gridfile= ', trim(gridfile)
   write(6,*) 'xifile= ', trim(xifile)
   write(6,*) 'nt1,nt2= ', nt1, nt2

   write(6,*) 'read gridfile'
   open(unit=13, file=trim(gridfile), &
     form='unformatted',access='stream',status='old')
   read(13) in,jn,kn
   write(6,*) 'in, jn, kn = ', in, jn, kn

   is=3
   js=3
   ks=3
   ie=in-3
   je=jn-3
   ke=kn-3
   ione=1
   jone=1
   kone=1

   ism1=is-ione
   jsm1=js-jone
   ksm1=ks-kone
   ism2=is-ione-ione
   jsm2=js-jone-jone
   ksm2=ks-kone-kone
   isp1=is+ione
   jsp1=js+jone
   ksp1=ks+kone
   isp2=is+ione+ione
   jsp2=js+jone+jone
   ksp2=ks+kone+kone

   iem1=ie-ione
   jem1=je-jone
   kem1=ke-kone
   iem2=ie-ione-ione
   jem2=je-jone-jone
   kem2=ke-kone-kone
   iep1=ie+ione
   jep1=je+jone
   kep1=ke+kone
   iep2=ie+ione+ione
   jep2=je+jone+jone
   kep2=ke+kone+kone
   iep3=ie+ione+ione+ione
   jep3=je+jone+jone+jone
   kep3=ke+kone+kone+kone

   allocate(x1a(in))
   allocate(x2a(jn))
   allocate(x3a(kn))

   allocate(x1b(in))
   allocate(x2b(jn))
   allocate(x3b(kn))

   allocate(xi(jn,kn))

   read(13) (x1a(i),i=ism2,iep3)
   read(13) (x1b(i),i=ism2,iep2)
   read(13) (x2a(j),j=jsm2,jep3)
   read(13) (x2b(j),j=jsm2,jep2)
   read(13) (x3a(k),k=ksm2,kep3)
   read(13) (x3b(k),k=ksm2,kep2)

   close(13)

   th1=x2a(js)
   th2=x2a(jep1)
   ph1=x3a(ks)
   ph2=x3a(kep1)
   ph0=(ph1+ph2)/2.
   delph=ph2-ph1
   area=delph*sin((th1+th2)/2.)*(th2-th1)
   np=area/lcell/lcell
   write(6,*) 'np= ',np

!  initialized x, y, tinit
   write(6,*) 'initialize x, y, tinit'
   allocate(x(np),y(np),tinit(np))
   allocate(icycle(np))
   iseedx=1
   iseedy=2
   iseedt=3
   do i=1,np
     xp=ran2(iseedx)*(th2-th1)+th1
     yp=ph0-delph/2./(sin(xp)+tiny) &
       +ran2(iseedy)*delph/(sin(xp)+tiny)
     tinit(i)=-ran2(iseedt)*tcycle
     icycle(i)=0
     x(i)=xp
     y(i)=yp
   enddo

!  start getting xi slices
   write(6,*) 'start getting xi slices'
   do itsub=nt1,nt2
     time=dt*dble(itsub-1)
     do i=1,np
       icount=(time-tinit(i))/tcycle
       if(icount .gt. icycle(i)) then
         xp=ran2(iseedx)*(th2-th1)+th1
         yp=ph0-delph/2./(sin(xp)+tiny) &
           +ran2(iseedy)*delph/(sin(xp)+tiny)
         icycle(i)=icount
         x(i)=xp
         y(i)=yp
       endif
     enddo

!    compute xi
     xi=0.
     do i=1,np
       if(y(i) .lt. ph2 .and. y(i) .gt. ph1) then
         do j=js,jep1
           do k=ks,kep1
             rr=x1a(is)*sqrt((cos(x2a(j))-cos(x(i)))**2 &
               +(sin(x2a(j))*cos(x3a(k))-sin(x(i))*cos(y(i)))**2 &
               +(sin(x2a(j))*sin(x3a(k))-sin(x(i))*sin(y(i)))**2)
             if(rr .lt. lcell*3.) then
               fcell=exp(-rr**2/lcell**2)
               xi(j,k)=xi(j,k)+fcell*sin((time-tinit(i))/tcycle*2.*pi)
             endif
           enddo
         enddo
       endif
     enddo
!    modify xi near boundary
     do j=js,jep1
       do k=ks,kep1
         rrth1=x1a(is)*sqrt((cos(x2a(j))-cos(th1))**2 &
           +(sin(x2a(j))*cos(x3a(k))-sin(th1)*cos(x3a(k)))**2 &
           +(sin(x2a(j))*sin(x3a(k))-sin(th1)*sin(x3a(k)))**2)
         rrth2=x1a(is)*sqrt((cos(x2a(j))-cos(th2))**2 &
           +(sin(x2a(j))*cos(x3a(k))-sin(th2)*cos(x3a(k)))**2 &
           +(sin(x2a(j))*sin(x3a(k))-sin(th2)*sin(x3a(k)))**2)
         rrph1=x1a(is)*sqrt( &
           (sin(x2a(j))*cos(x3a(k))-sin(x2a(j))*cos(ph1))**2 &
           +(sin(x2a(j))*sin(x3a(k))-sin(x2a(j))*sin(ph1))**2)
         rrph2=x1a(is)*sqrt( &
           (sin(x2a(j))*cos(x3a(k))-sin(x2a(j))*cos(ph2))**2 &
           +(sin(x2a(j))*sin(x3a(k))-sin(x2a(j))*sin(ph2))**2)
         rrmin=min(rrth1,rrth2,rrph1,rrph2)
         xi(j,k)=xi(j,k)*(1.-exp(-rrmin**2/lcell**2))
       enddo
     enddo

     write(idout,'(i4.4)') itsub
     write(6,*) 'writing ', trim(xifile)//idout//'.dat'
     open(unit=13, file=trim(xifile)//idout//'.dat', &
       form='unformatted',access='stream')
     write(13) time
     write(13) jn,kn
     write(13) ((xi(j,k),j=1,jn),k=1,kn)
     close(13)
   enddo

end program getxi

FUNCTION ran2(idum)
      implicit none
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1, &
      IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791, &
      NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
!  (C) Copr. 1986-92 Numerical Recipes Software '%12'%&WR2.
end FUNCTION ran2
