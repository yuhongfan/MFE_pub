module ModUserSetup

  implicit none

  private

  real :: etasqrt,z1,lambda,a0,p0

  public :: etasqrt
  public :: grid, field_init
  public :: getemf1ib, getemf1jb, getemf1kb
  public :: getdiffemf1ib, getdiffemf1jb, getdiffemf1kb
  public :: getemf2ib, getemf2jb, getemf2kb
  public :: getdiffemf2ib, getdiffemf2jb, getdiffemf2kb
  public :: getemf3ib, getemf3jb, getemf3kb
  public :: getdiffemf3ib, getdiffemf3jb, getdiffemf3kb
  public :: geteintib, geteintjb, geteintkb
  public :: getdib, getdjb, getdkb
  public :: getv1ib, getv1jb, getv1kb
  public :: getv2ib, getv2jb, getv2kb
  public :: getv3ib, getv3jb, getv3kb
  public :: getdelsib, getdelsjb, getdelskb
  public :: getbib, getbjb, getbkb
  public :: userstep

contains

subroutine grid
  use ModPar
  use ModPhysics,       ONLY:pi
  use ModGrid
  implicit none

  include "mpif.h"

  integer :: i, j, k
  real :: x1min,x2min,x3min,x1max,x2max,x3max,delx1,delx2,delx3

!-------------------------------------------------------

      x1min=1.D0
      x1max=5.D0
      x2min=0.D0
      x2max=pi
      x3min=-pi
      x3max=pi
      delx1=(x1max-x1min)/dble(inmax-5)
      delx2=(x2max-x2min)/dble(jnmax-5)
      delx3=(x3max-x3min)/dble(knmax-5)

      do i = 1, inmax
        xxa(i)=dble(i-is)*delx1+x1min
        dxxa(i)=delx1
        xxb(i)=xxa(i)+0.5D0*delx1
        dxxb(i)=delx1
      enddo

      do j = 1, jnmax
        yya(j)=dble(j-js)*delx2+x2min
        dyya(j)=delx2
        yyb(j)=yya(j)+0.5D0*delx2
        dyyb(j)=delx2
      enddo

      do k = 1, knmax
        zza(k)=dble(k-ks)*delx3+x3min
        dzza(k)=delx3
        zzb(k)=zza(k)+0.5D0*delx3
        dzzb(k)=delx3
      enddo

      do i = 1, in
        x1a(i)=xxa(i+myid1*(in-5))
        dx1a(i)=dxxa(i+myid1*(in-5))
        x1b(i)=xxb(i+myid1*(in-5))
        dx1b(i)=dxxb(i+myid1*(in-5))
      enddo

      do j = 1, jn
        x2a(j)=yya(j+myid2*(jn-5))
        dx2a(j)=dyya(j+myid2*(jn-5))
        x2b(j)=yyb(j+myid2*(jn-5))
        dx2b(j)=dyyb(j+myid2*(jn-5))
      enddo

      do k = 1, kn
        x3a(k)=zza(k+myid3*(kn-5))
        dx3a(k)=dzza(k+myid3*(kn-5))
        x3b(k)=zzb(k+myid3*(kn-5))
        dx3b(k)=dzzb(k+myid3*(kn-5))
      enddo

      do i = 1, in
        g2a(i)=x1a(i)
        g31a(i)=x1a(i)
        dg2bd1(i)=1.D0
        dg31bd1(i)=1.D0
        g2b(i)=x1b(i)
        g31b(i)=x1b(i)
        dg2ad1(i)=1.D0
        dg31ad1(i)=1.D0
      enddo

      do j = 1, jn
        g32a   (j) = abs(sin ( x2a(j) ))
        g32b   (j) = abs(sin ( x2b(j) ))
      enddo

      do j = 1, jn
        dg32bd2(j) = cos(x2a(j))
        dg32ad2(j) = cos(x2b(j))
      enddo
      if (myid2 .eq. 0) then
        dg32bd2(js) = 0.D0
        dg32bd2(jsm1) = - dg32bd2(jsp1)
        dg32bd2(jsm2) = - dg32bd2(jsp2)
        dg32ad2(jsm1) = - dg32ad2(js)
        dg32ad2(jsm2) = - dg32ad2(jsp1)
      endif
      if (myid2 .eq. nproc2-1) then
        dg32bd2(jep1) = 0.D0
        dg32bd2(jep2) = - dg32bd2(je)
        dg32bd2(jep3) = - dg32bd2(jem1)
        dg32ad2(jep1) = - dg32ad2(je)
        dg32ad2(jep2) = - dg32ad2(jem1)
      endif

      do i = 1, in
        dvl1a(i)=x1b(i)*x1b(i)*dx1a(i)
      enddo
      do i = 1, in
        dvl1b(i)=x1a(i)*x1a(i)*dx1b(i)
      enddo

      do j = 1, jn
        dvl2a(j) = abs(sin(x2b(j)))*dx2a(j)
      enddo
      do j = 1, jn
        dvl2b(j) = abs(sin(x2a(j)))*dx2b(j)
      enddo

      do k= 1, kn
        dvl3a(k) = dx3a(k)
      enddo
      do k=1, kn
        dvl3b(k) = dx3b(k)
      enddo

  do i = ism2, iep3
    dx1bi (i) = 1.D0 / ( dx1b (i) + tiny )
    g2ai  (i) = 1.D0 / ( g2a  (i) + tiny )
    g31ai (i) = 1.D0 / ( g31a (i) + tiny )
    dvl1bi(i) = 1.D0 / ( dvl1b(i) + tiny )
  enddo
  do i = ism2, iep2
    dx1ai (i) = 1.D0 / ( dx1a (i) + tiny )
    g2bi  (i) = 1.D0 / ( g2b  (i) + tiny )
    g31bi (i) = 1.D0 / ( g31b (i) + tiny )
    dvl1ai(i) = 1.D0 / ( dvl1a(i) + tiny )
  enddo

  do j=jsm2,jep3
    dx2bi (j) = 1.D0 / ( dx2b (j) + tiny )
    g32ai (j) = 1.D0 / ( g32a (j) + tiny )
    dvl2bi(j) = 1.D0 / ( dvl2b(j) + tiny )
  enddo
  do j=jsm2,jep2
    dx2ai (j) = 1.D0 / ( dx2a (j) + tiny )
    g32bi (j) = 1.D0 / ( g32b (j) + tiny )
    dvl2ai(j) = 1.D0 / ( dvl2a(j) + tiny )
  enddo
! SAFEGUARD
  if (myid2 .eq. 0) then
      g32ai(js) = 0.D0
      dvl2bi(js) = 0.D0
  endif
  if (myid2 .eq. nproc2-1) then
      g32ai(jep1) = 0.D0
      dvl2bi(jep1) = 0.D0
  endif

  do k=ksm2,kep3
    dx3bi (k) = 1.D0 / ( dx3b (k) + tiny )
    dvl3bi(k) = 1.D0 / ( dvl3b(k) + tiny )
  enddo
  do k=ksm2,kep2
    dx3ai (k) = 1.D0 / ( dx3a (k) + tiny )
    dvl3ai(k) = 1.D0 / ( dvl3a(k) + tiny )
  enddo

end subroutine grid

subroutine field_init
  use ModPar
  use ModPhysics
  use ModGrid
  use ModSundry
  use ModField
  use ModBoundary
  use ModIO
  use ModControl,        ONLY:irst

  implicit none

  include "mpif.h"

  integer :: i,j,k,kbc,ierr
  real :: msol,r_c,gacc_c,rho_c,temp_c
  real :: rc, rs, esh, xsh, xshup
!
!-----------------------------------------------------------------
! specify boundary condition
!-----------------------------------------------------------------
!
       niib=3
       noib=2
       nijb=3
       nojb=3
       nikb=4
       nokb=4
       isplt=-1
!
!-----------------------------------------------------
!      initialize fields or use restart file
!-----------------------------------------------------
!
!      specify gravitaional acceleration gacc
!      specify temp_c and rho_c
!
       etasqrt=0.626583
       z1=6.11676
       lambda=z1/2.D0**(1.D0/7.D0)
       a0=0.0895694
       p0=-func(z1)

       msol=2.D33
       r_c=xxa(is)*unit_len
       gacc_c=g_const*msol/r_c/r_c
       gacc=gacc_c/(unit_v*unit_v/unit_len)*(r_c/unit_len)**2
       rho_c=unit_rho
       temp_c=unit_temp
       rstar=kboltz/(muconst*mproton)

#ifdef THCONDUCT
       kappa0_dim=1.D-6
       kappa0=kappa0_dim*unit_temp**(7.D0/2.D0) &
         /(unit_v**2*unit_rho*unit_len**3) &
         *unit_time*unit_len
       diffcoef0=kappa0_dim*unit_temp**(7.D0/2.D0) &
         /(rstar*unit_temp*unit_rho*unit_len**3) &
         *unit_time*unit_len*(gamma-1.D0)
       aval0=1.D0
#endif
!
!      output physical constants and units
!
       if(myid .eq. 0) then
       open(unit=18,file='physparams.dat', &
         form='unformatted',access='stream')
       write(18) g_const
       write(18) rgas
       write(18) kboltz
       write(18) mproton
       write(18) pi
       write(18) gamma
       write(18) unit_rho
       write(18) unit_len
       write(18) unit_b
       write(18) unit_temp
       write(18) unit_v
       write(18) unit_time
       write(18) msol
       write(18) r_c
       write(18) temp_c
       write(18) rho_c
       write(18) muconst
       close(18)
       endif
!
!      initialize fields or use restart file
!
       if(irst .ne. 1) then

       time=3.19219D0
       ctime=time
       rc=2.D0**(1.D0/7.D0)*etasqrt*time
       rs=(etasqrt*time)**(7.D0/6.D0)
!
!      initialize b1
!
       do k=ksm2,kep2
         do j=jsm2,jep2
           do i=ism2,iep2
             if(x1a(i) .le. rc) then
               b1(i,j,k)=2.D0*a0/x1a(i)**2*(p0 &
                 +func(lambda*x1a(i)/(etasqrt*time))) &
                 *cos(x2b(j))
             else
               b1(i,j,k)=0.D0
             endif
           enddo
         enddo
       enddo
       if (myid2 .eq. 0) then
       do k=ksm2,kep2
         do i=ism2,iep2
           kbc=ks+mod(k-ks+(knmax-5)/2,knmax-5)
           b1(i,jsm2,k)=b1(i,jsp1,kbc)
           b1(i,jsm1,k)=b1(i,js,kbc)
         enddo
       enddo
       endif
       if (myid2 .eq. nproc2-1) then
       do k=ksm2,kep2
         do i=ism2,iep2
           kbc=ks+mod(k-ks+(knmax-5)/2,knmax-5)
           b1(i,jep1,k)=b1(i,je,kbc)
           b1(i,jep2,k)=b1(i,jem1,kbc)
         enddo
       enddo
       endif
!
!      initialize b2
!
       do k=ksm2,kep2
         do j=jsm2,jep2
           do i=ism2,iep2
             if(x1b(i) .le. rc) then
               b2(i,j,k)=-a0*lambda/x1b(i)/(etasqrt*time) &
                 *funcp(lambda*x1b(i)/(etasqrt*time)) &
                 *sin(x2a(j))
             else
               b2(i,j,k)=0.D0
             endif
           enddo
         enddo
       enddo
       if (myid2 .eq. 0) then
       do k=ksm2,kep2
         do i=ism2,iep2
           kbc=ks+mod(k-ks+(knmax-5)/2,knmax-5)
           b2(i,js,k)=0.D0
           b2(i,jsm1,k)=-b2(i,jsp1,kbc)
           b2(i,jsm2,k)=-b2(i,jsp2,kbc)
         enddo
       enddo
       endif
       if (myid2 .eq. nproc2-1) then
       do k=ksm2,kep2
         do i=ism2,iep2
           kbc=ks+mod(k-ks+(knmax-5)/2,knmax-5)
           b2(i,jep1,k)=0.D0
           b2(i,jep2,k)=-b2(i,je,kbc)
           b2(i,jep3,k)=-b2(i,jem1,kbc)
         enddo
       enddo
       endif
!
!      initialize b3
!
       do k=ksm2,kep2
         do j=jsm2,jep2
           do i=ism2,iep2
             if(x1b(i) .le. rc) then
               b3(i,j,k)=a0*lambda/x1b(i)/(etasqrt*time) &
                 *(p0+func(lambda*x1b(i)/(etasqrt*time))) &
                 *sin(x2b(j))
             else
               b3(i,j,k)=0.D0
             endif
           enddo
         enddo
       enddo
       if (myid2 .eq. 0) then
       do k=ksm2,kep2
         do i=ism2,iep2
           kbc=ks+mod(k-ks+(knmax-5)/2,knmax-5)
           b3(i,jsm2,k)=-b3(i,jsp1,kbc)
           b3(i,jsm1,k)=-b3(i,js,kbc)
         enddo
       enddo
       endif
       if (myid2 .eq. nproc2-1) then
       do k=ksm2,kep2
         do i=ism2,iep2
           kbc=ks+mod(k-ks+(knmax-5)/2,knmax-5)
           b3(i,jep1,k)=-b3(i,je,kbc)
           b3(i,jep2,k)=-b3(i,jem1,kbc)
         enddo
       enddo
       endif
!
!      initialize e
!
       if(myid .eq. 0) then
         esh=0.D0
         do i=isp1,inmax-3
           xshup=xxb(is)+4.D0/dble(inmax-5)*dble(i-is)
           xsh=xxb(is)+4.D0/dble(inmax-5)*dble(i-1-is)
           if(xshup .gt. rs .and. xsh .le. rs) then
             esh=1.D0/(gamma-1.D0) &
              *7.D0/6.D0*exp(-2.D0/3.D0/etasqrt**2) &
              *etasqrt**2*(etasqrt*time)**(-4.D0) &
              *exp(2.D0/3.D0*etasqrt**(7.D0) &
              *(time/xsh)**(9.D0))
             write(6,*) 'ish=',i
           endif
         enddo
       endif
       call MPI_BCAST(esh, 1, MPI_DOUBLE_PRECISION, 0, &
         MPI_COMM_WORLD,ierr)

       do k=ksm2,kep2
         do j=jsm2,jep2
           do i=ism2,iep2
             if(x1b(i) .le. rc) then
               eint(i,j,k)=1.D0/4.D0/x1b(i)**4/(gamma-1.D0) &
                 +a0**2*p0/x1b(i)**4/(gamma-1.D0) &
                 *(2.D0-(lambda*x1b(i)/etasqrt/time)**2) &
                 *(p0+func(lambda*x1b(i)/(etasqrt*time))) &
                 *sin(x2b(j))**2
             else
               if(x1b(i) .le. rs) then
                 eint(i,j,k)=1.D0/(gamma-1.D0) &
                   *7.D0/6.D0*exp(-2.D0/3.D0/etasqrt**2) &
                   *etasqrt**2*(etasqrt*time)**(-4.D0) &
                   *exp(2.D0/3.D0*etasqrt**(7.D0) &
                   *(time/x1b(i))**(9.D0))
               else
                 eint(i,j,k)=esh*1.D-3
               endif
             endif
           enddo
         enddo
       enddo
!
!      initialize d
!
       do k=ksm2,kep2
         do j=jsm2,jep2
           do i=ism2,iep2
             if(x1b(i) .le. rc) then
               d(i,j,k)=1.D0/x1b(i)**3 &
                 +2.D0*a0**2*p0/x1b(i)**3 &
                 *(4.D0-(lambda*x1b(i)/etasqrt/time)**2) &
                 *(p0+func(lambda*x1b(i)/(etasqrt*time))) &
                 *sin(x2b(j))**2
             else
               if(x1b(i) .le. rs) then
                 d(i,j,k)=7.D0 &
                   *exp(-2.D0/3.D0/etasqrt**2) &
                   *(etasqrt*time)**(-3.D0) &
                   *(etasqrt*time/x1b(i))**(8.D0) &
                   *exp(2.D0/3.D0*etasqrt**(7.D0) &
                   *(time/x1b(i))**(9.D0))
               else
                 d(i,j,k)=exp(-2.D0/3.D0/etasqrt**2) &
                   *(x1b(i))**(-26.D0/7.D0) &
                   *exp(2.D0/3.D0/etasqrt**2 &
                   /(x1b(i))**(9.D0/7.D0))
               endif
             endif
           enddo
         enddo
       enddo
!
!      initialize v1
!
       do k=ksm2,kep2
         do j=jsm2,jep2
           do i=ism2,iep2
             if(x1a(i) .le. rs) then
               v1(i,j,k)=1.D0/time*x1a(i)
             else
               v1(i,j,k)=0.D0
             endif
           enddo
         enddo
       enddo
!
!      initialize v2 field
!
       do k=ksm2,kep2
         do j=jsm2,jep2
           do i=ism2,iep2
             v2(i,j,k)=0.D0
           enddo
         enddo
       enddo
!
!      initialize v3 field
!
       do k=ksm2,kep2
         do j=jsm2,jep2
           do i=ism2,iep2
             v3(i,j,k)=0.D0
           enddo
         enddo
       enddo

       else

       call readrst_mpi
       ctime=time

       endif
!
!      parameters in sundry
!
       dfloor=1.D-11
       efloor=1.D-11
       if (myid .eq. 0) then
         dt=dx1a(is)
       endif
       call MPI_BCAST(dt, 1, MPI_DOUBLE_PRECISION, 0, &
         MPI_COMM_WORLD,ierr)
       dtmdi2 = 0.D0
       alpha_d = 1.D0
       alpha_v = 1.D0
       alpha_e = 1.D0
       nlf_d = 4
       nlf_v = 4
       nlf_e = 4
!
       itnow=0
       ntcond=0
       courno=0.7D0

end subroutine field_init
!=========================================================
subroutine getemf1ib( v2b3mv3b2 )
  use ModPar
  use ModGrid
  use ModBoundary

  implicit none

  real, intent(in), dimension(in,jn,kn) :: v2b3mv3b2
  integer :: j, k
!---------------------------------------------------------
      if(myid1 .eq. 0) then

      do k=ks,kep1
         do j=js,jep1
           emf1iib1(j,k)=0.D0
           emf1iib2(j,k)=0.D0
         enddo
      enddo

      endif

      if(myid1 .eq. nproc1-1) then

      do k=ks,kep1
        do j=js,jep1
          emf1oib1(j,k) = v2b3mv3b2(ie  ,j,k)
          emf1oib2(j,k) = v2b3mv3b2(ie  ,j,k)
        enddo
      enddo

      endif

end subroutine getemf1ib
!=========================================================
subroutine getdiffemf1ib( v2b3mv3b2 )
  use ModPar
  use ModGrid
  use ModBoundary

  implicit none
  real, intent(in), dimension(in,jn,kn) :: v2b3mv3b2
  integer :: j, k
!---------------------------------------------------------
      if(myid1 .eq. 0) then

      do k=ks,kep1
         do j=js,jep1
           emf1iib1(j,k)=0.D0
           emf1iib2(j,k)=0.D0
         enddo
      enddo

      endif

      if(myid1 .eq. nproc1-1) then

      do k=ks,kep1
        do j=js,jep1
          emf1oib1(j,k) = v2b3mv3b2(ie  ,j,k)
          emf1oib2(j,k) = v2b3mv3b2(ie  ,j,k)
        enddo
      enddo

      endif

end subroutine getdiffemf1ib
!=========================================================
subroutine getemf1jb( v2b3mv3b2 )
  use ModPar
  use ModGrid
  use ModBoundary

  implicit none

  real, intent(in), dimension(in,jn,kn) :: v2b3mv3b2
  integer :: k, i, kbc
!---------------------------------------------------------
      if(myid2 .eq. 0) then

      do i=ism2,iep2
        do k=ks,kep1
          kbc=ks+mod(k-ks+(knmax-5)/2,knmax-5)
          emf1ijb1(k,i) = 0.D0
          emf1ijb2(k,i) = v2b3mv3b2(i,jsp1,kbc)
          emf1ijb3(k,i) = v2b3mv3b2(i,jsp2,kbc)
        enddo
      enddo

      endif

      if(myid2 .eq. nproc2-1) then

      do i=ism2,iep2
        do k=ks,kep1
          kbc=ks+mod(k-ks+(knmax-5)/2,knmax-5)
          emf1ojb1(k,i) = 0.D0
          emf1ojb2(k,i) = v2b3mv3b2(i,je,kbc)
          emf1ojb3(k,i) = v2b3mv3b2(i,jem1,kbc)
        enddo
      enddo

      endif

end subroutine getemf1jb
!=========================================================
subroutine getdiffemf1jb( v2b3mv3b2 )
  use ModPar
  use ModGrid
  use ModBoundary

  implicit none

  real, intent(in), dimension(in,jn,kn) :: v2b3mv3b2
  integer :: k, i, kbc
!---------------------------------------------------------
      if(myid2 .eq. 0) then

      do i=ism2,iep2
        do k=ks,kep1
          kbc=ks+mod(k-ks+(knmax-5)/2,knmax-5)
          emf1ijb1(k,i) = 0.D0
          emf1ijb2(k,i) = v2b3mv3b2(i,jsp1,kbc)
          emf1ijb3(k,i) = v2b3mv3b2(i,jsp2,kbc)
        enddo
      enddo

      endif

      if(myid2 .eq. nproc2-1) then

      do i=ism2,iep2
        do k=ks,kep1
          kbc=ks+mod(k-ks+(knmax-5)/2,knmax-5)
          emf1ojb1(k,i) = 0.D0
          emf1ojb2(k,i) = v2b3mv3b2(i,je,kbc)
          emf1ojb3(k,i) = v2b3mv3b2(i,jem1,kbc)
        enddo
      enddo

      endif

end subroutine getdiffemf1jb
!=========================================================
subroutine getemf1kb( v2b3mv3b2 )
  use ModPar,      ONLY: in,jn,kn
  implicit none
  real, dimension(in,jn,kn) :: v2b3mv3b2

end subroutine getemf1kb
!=========================================================
subroutine getdiffemf1kb( v2b3mv3b2 )
  use ModPar,      ONLY: in,jn,kn
  implicit none
  real, dimension(in,jn,kn) :: v2b3mv3b2

end subroutine getdiffemf1kb
!=========================================================
subroutine getemf2ib( v3b1mv1b3 )
  use ModPar
  use ModGrid
  use ModBoundary
  use ModSundry,      ONLY: ctime

  implicit none

  real, intent(in), dimension(in,jn,kn) :: v3b1mv1b3
  integer :: j, k, kbc
!---------------------------------------------------------
      if(myid1 .eq. 0) then

      do k=ksm2,kep3
        do j=jsm2,jep2
            emf2iib1(j,k)=-a0*lambda &
            /(etasqrt*ctime**2)*(p0 &
            +func(lambda*x1a(is)/(etasqrt*ctime))) &
            *sin(x2b(j))
            emf2iib2(j,k)=-a0*lambda &
            /(etasqrt*ctime**2)*(p0 &
            +func(lambda*x1a(ism1)/(etasqrt*ctime))) &
            *sin(x2b(j))
            emf2iib3(j,k)=-a0*lambda &
            /(etasqrt*ctime**2)*(p0 &
            +func(lambda*x1a(ism2)/(etasqrt*ctime))) &
            *sin(x2b(j))
        enddo
      enddo
!     need to fix the j boundary part because that got over written above
      if(myid2 .eq. 0) then
        do k=ksm2,kep3
          kbc=ks+mod(k-ks+(knmax-5)/2,knmax-5)
          emf2iib1(jsm1,k)=-emf2iib1(js,kbc)
          emf2iib1(jsm2,k)=-emf2iib1(jsp1,kbc)
          emf2iib2(jsm1,k)=-emf2iib2(js,kbc)
          emf2iib2(jsm2,k)=-emf2iib2(jsp1,kbc)
          emf2iib3(jsm1,k)=-emf2iib3(js,kbc)
          emf2iib3(jsm2,k)=-emf2iib3(jsp1,kbc)
        enddo
      endif
      if(myid2 .eq. nproc2-1) then
        do k=ksm2,kep3
          kbc=ks+mod(k-ks+(knmax-5)/2,knmax-5)
          emf2iib1(jep1,k)=-emf2iib1(je,kbc)
          emf2iib1(jep2,k)=-emf2iib1(jem1,kbc)
          emf2iib2(jep1,k)=-emf2iib2(je,kbc)
          emf2iib2(jep2,k)=-emf2iib2(jem1,kbc)
          emf2iib3(jep1,k)=-emf2iib3(je,kbc)
          emf2iib3(jep2,k)=-emf2iib3(jem1,kbc)
        enddo
      endif

      endif

      if(myid1 .eq. nproc1-1) then

      do k=ksm2,kep3
        do j=jsm2,jep2
          emf2oib1(j,k) = v3b1mv1b3(iep1,j,k)
          emf2oib2(j,k) = emf2oib1(j,k) + emf2oib1(j,k) &
               - v3b1mv1b3(ie,j,k)
          emf2oib3(j,k) = emf2oib2(j,k) + emf2oib2(j,k) &
               - emf2oib1(j,k)
        enddo
      enddo

      endif

end subroutine getemf2ib
!=========================================================
subroutine getdiffemf2ib( v3b1mv1b3 )
  use ModPar
  use ModGrid
  use ModBoundary

  implicit none

  real, intent(in), dimension(in,jn,kn) :: v3b1mv1b3

  integer :: j, k
!---------------------------------------------------------
      if(myid1 .eq. 0) then

      do k=ksm2,kep3
        do j=jsm2,jep2
          emf2iib1(j,k)=0.D0
          emf2iib2(j,k)=0.D0
          emf2iib3(j,k)=0.D0
        enddo
      enddo

      endif

      if(myid1 .eq. nproc1-1) then

      do k=ksm2,kep3
        do j=jsm2,jep2
          emf2oib1(j,k) = v3b1mv1b3(iep1,j,k)
          emf2oib2(j,k) = emf2oib1(j,k) + emf2oib1(j,k) &
               - v3b1mv1b3(ie,j,k)
          emf2oib3(j,k) = emf2oib2(j,k) + emf2oib2(j,k) &
               - emf2oib1(j,k)
        enddo
      enddo

      endif

end subroutine getdiffemf2ib
!=========================================================
subroutine getemf2jb( v3b1mv1b3 )
  use ModPar
  use ModGrid
  use ModBoundary

  implicit none

  real, intent(in), dimension(in,jn,kn) :: v3b1mv1b3
  integer :: k, i, kbc
!---------------------------------------------------------
      if(myid2 .eq. 0) then

      do i=is,iep1
        do k=ks,kep1
          kbc=ks+mod(k-ks+(knmax-5)/2,knmax-5)
          emf2ijb1(k,i) = -v3b1mv1b3(i,js,kbc)
          emf2ijb2(k,i) = -v3b1mv1b3(i,jsp1,kbc)
        enddo
      enddo

      endif

      if(myid2 .eq. nproc2-1) then

      do i=is,iep1
        do k=ks,kep1
          kbc=ks+mod(k-ks+(knmax-5)/2,knmax-5)
          emf2ojb1(k,i) = -v3b1mv1b3(i,je,kbc)
          emf2ojb2(k,i) = -v3b1mv1b3(i,jem1,kbc)
        enddo
      enddo

      endif

end subroutine getemf2jb
!=========================================================
subroutine getdiffemf2jb( v3b1mv1b3 )
  use ModPar
  use ModGrid
  use ModBoundary

  implicit none

  real, intent(in), dimension(in,jn,kn) :: v3b1mv1b3
  integer :: k, i, kbc
!---------------------------------------------------------
      if(myid2 .eq. 0) then

      do i=is,iep1
        do k=ks,kep1
          kbc=ks+mod(k-ks+(knmax-5)/2,knmax-5)
          emf2ijb1(k,i) = -v3b1mv1b3(i,js,kbc)
          emf2ijb2(k,i) = -v3b1mv1b3(i,jsp1,kbc)
        enddo
      enddo

      endif

      if(myid2 .eq. nproc2-1) then

      do i=is,iep1
        do k=ks,kep1
          kbc=ks+mod(k-ks+(knmax-5)/2,knmax-5)
          emf2ojb1(k,i) = -v3b1mv1b3(i,je,kbc)
          emf2ojb2(k,i) = -v3b1mv1b3(i,jem1,kbc)
        enddo
      enddo

      endif

end subroutine getdiffemf2jb
!=========================================================
subroutine getemf2kb( v3b1mv1b3 )
  use ModPar,      ONLY: in,jn,kn
  implicit none
  real, dimension(in,jn,kn) :: v3b1mv1b3

end subroutine getemf2kb
!=========================================================
subroutine getdiffemf2kb( v3b1mv1b3 )
  use ModPar,      ONLY: in,jn,kn
  implicit none
  real,  dimension(in,jn,kn) :: v3b1mv1b3

end subroutine getdiffemf2kb
!=========================================================
subroutine getemf3ib( v1b2mv2b1 )
  use ModPar
  use ModGrid
  use ModBoundary
  use ModSundry,      ONLY: ctime

  implicit none

  real, intent(in), dimension(in,jn,kn) :: v1b2mv2b1
  integer :: j, k
!---------------------------------------------------------
      if(myid1 .eq. 0) then

      do k=ksm2,kep2
        do j=js,jep1
          emf3iib1(j,k)=-a0*lambda &
            /(etasqrt*ctime**2) &
            *funcp(lambda*x1a(is)/(etasqrt*ctime)) &
            *sin(x2a(j))
          emf3iib2(j,k)=-a0*lambda &
            /(etasqrt*ctime**2) &
            *funcp(lambda*x1a(ism1)/(etasqrt*ctime)) &
            *sin(x2a(j))
          emf3iib3(j,k)=-a0*lambda &
            /(etasqrt*ctime**2) &
            *funcp(lambda*x1a(ism2)/(etasqrt*ctime)) &
            *sin(x2a(j))
        enddo
      enddo

      endif

      if(myid1 .eq. nproc1-1) then

      do k=ksm2,kep2
        do j=js,jep1
          emf3oib1(j,k) = v1b2mv2b1(iep1,j,k)
          emf3oib2(j,k) = emf3oib1(j,k) + emf3oib1(j,k) &
               - v1b2mv2b1(ie,j,k)
          emf3oib3(j,k) = emf3oib2(j,k) + emf3oib2(j,k) &
               - emf3oib1(j,k)
        enddo
      enddo

      endif

end subroutine getemf3ib
!=========================================================
subroutine getdiffemf3ib( v1b2mv2b1 )
  use ModPar
  use ModGrid
  use ModBoundary

  implicit none
  real, intent(in), dimension(in,jn,kn) :: v1b2mv2b1
  integer :: j, k
!---------------------------------------------------------
      if(myid1 .eq. 0) then

      do k=ksm2,kep2
        do j=js,jep1
          emf3iib1(j,k)=0.D0
          emf3iib2(j,k)=0.D0
          emf3iib3(j,k)=0.D0
        enddo
      enddo

      endif

      if(myid1 .eq. nproc1-1) then

      do k=ksm2,kep2
        do j=js,jep1
          emf3oib1(j,k) = v1b2mv2b1(iep1,j,k)
          emf3oib2(j,k) = emf3oib1(j,k) + emf3oib1(j,k) &
               - v1b2mv2b1(ie,j,k)
          emf3oib3(j,k) = emf3oib2(j,k) + emf3oib2(j,k) &
               - emf3oib1(j,k)
        enddo
      enddo

      endif

end subroutine getdiffemf3ib
!=========================================================
subroutine getemf3jb( v1b2mv2b1 )
  use ModPar
  use ModGrid
  use ModBoundary

  implicit none

  real, intent(in), dimension(in,jn,kn) :: v1b2mv2b1
  integer :: k, i, kbc
!---------------------------------------------------------
      if(myid2 .eq. 0) then

      do i=ism2,iep3
        do k=ksm2,kep2
          kbc=ks+mod(k-ks+(knmax-5)/2,knmax-5)
          emf3ijb1(k,i) = 0.D0
          emf3ijb2(k,i) = -v1b2mv2b1(i,jsp1,kbc)
          emf3ijb3(k,i) = -v1b2mv2b1(i,jsp2,kbc)
        enddo
      enddo

      endif

      if(myid2 .eq. nproc2-1) then

      do i=ism2,iep3
        do k=ksm2,kep2
          kbc=ks+mod(k-ks+(knmax-5)/2,knmax-5)
          emf3ojb1(k,i) = 0.D0
          emf3ojb2(k,i) = -v1b2mv2b1(i,je,kbc)
          emf3ojb3(k,i) = -v1b2mv2b1(i,jem1,kbc)
        enddo
      enddo

      endif

end subroutine getemf3jb
!=========================================================
subroutine getdiffemf3jb( v1b2mv2b1 )
  use ModPar
  use ModGrid
  use ModBoundary

  implicit none

  real, intent(in), dimension(in,jn,kn) :: v1b2mv2b1
  integer :: k, i, kbc
!---------------------------------------------------------
      if(myid2 .eq. 0) then

      do i=ism2,iep3
        do k=ksm2,kep2
          kbc=ks+mod(k-ks+(knmax-5)/2,knmax-5)
          emf3ijb1(k,i) = 0.D0
          emf3ijb2(k,i) = -v1b2mv2b1(i,jsp1,kbc)
          emf3ijb3(k,i) = -v1b2mv2b1(i,jsp2,kbc)
        enddo
      enddo

      endif

      if(myid2 .eq. nproc2-1) then

      do i=ism2,iep3
        do k=ksm2,kep2
          kbc=ks+mod(k-ks+(knmax-5)/2,knmax-5)
          emf3ojb1(k,i) = 0.D0
          emf3ojb2(k,i) = -v1b2mv2b1(i,je,kbc)
          emf3ojb3(k,i) = -v1b2mv2b1(i,jem1,kbc)
        enddo
      enddo

      endif

end subroutine getdiffemf3jb
!=========================================================
subroutine getemf3kb( v1b2mv2b1 )
  use ModPar,      ONLY: in,jn,kn
  implicit none
  real, dimension(in,jn,kn) :: v1b2mv2b1

end subroutine getemf3kb
!=========================================================
subroutine getdiffemf3kb( v1b2mv2b1 )
  use ModPar,      ONLY: in,jn,kn
  implicit none
  real, dimension(in,jn,kn) :: v1b2mv2b1

end subroutine getdiffemf3kb
!=========================================================
subroutine geteintib
  use ModPar
  use ModPhysics,        ONLY:gamma
  use ModGrid
  use ModField,          ONLY: eint
  use ModBoundary
  use ModSundry,         ONLY: ctime

  implicit none

  integer :: j, k
!---------------------------------------------------------
      if(myid1 .eq. 0) then

      do k=ks,ke
        do j=js,je
           eint_iib1(j,k)=0.25D0/x1b(ism1)**4/(gamma-1.D0) &
             +a0**2*p0/x1b(ism1)**4/(gamma-1.D0) &
             *(2.D0-(lambda*x1b(ism1)/etasqrt/ctime)**2) &
             *(p0+func(lambda*x1b(ism1)/(etasqrt*ctime))) &
             *sin(x2b(j))**2
           eint_iib2(j,k)=0.25D0/x1b(ism2)**4/(gamma-1.D0) &
             +a0**2*p0/x1b(ism2)**4/(gamma-1.D0) &
             *(2.D0-(lambda*x1b(ism2)/etasqrt/ctime)**2) &
             *(p0+func(lambda*x1b(ism2)/(etasqrt*ctime))) &
             *sin(x2b(j))**2
        enddo
      enddo

      endif

      if(myid1 .eq. nproc1-1) then

      do k=ks,ke
        do j=js,je
          eint_oib1(j,k) = eint(ie,j,k)
          eint_oib2(j,k) = eint(ie,j,k)
        enddo
      enddo

      endif

end subroutine geteintib
!=========================================================
subroutine geteintjb
  use ModPar
  use ModGrid
  use ModField,          ONLY: eint
  use ModBoundary

  implicit none

  integer :: k, i, kbc
!---------------------------------------------------------
      if(myid2 .eq. 0) then

      do i=ism2,iep2
        do k=ks,ke
          kbc=ks+mod(k-ks+(knmax-5)/2,knmax-5)
          eint_ijb1(k,i) = eint(i,js,kbc)
          eint_ijb2(k,i) = eint(i,jsp1,kbc)
        enddo
      enddo

      endif

      if(myid2 .eq. nproc2-1) then

      do i=ism2,iep2
        do k=ks,ke
          kbc=ks+mod(k-ks+(knmax-5)/2,knmax-5)
          eint_ojb1(k,i) = eint(i,je,kbc)
          eint_ojb2(k,i) = eint(i,jem1,kbc)
        enddo
      enddo

      endif

end subroutine geteintjb
!=========================================================
subroutine geteintkb

end subroutine geteintkb
!=========================================================
subroutine getdib
  use ModPar
  use ModGrid
  use ModField,          ONLY: d
  use ModBoundary
  use ModSundry,         ONLY: ctime

  implicit none

  integer :: j, k
!---------------------------------------------------------
      if(myid1 .eq. 0) then

      do k=ks,ke
        do j=js,je
          diib1(j,k)=1.D0/x1b(ism1)**3 &
            +2.D0*a0**2*p0/x1b(ism1)**3 &
            *(4.D0-(lambda*x1b(ism1)/etasqrt/ctime)**2) &
            *(p0+func(lambda*x1b(ism1)/(etasqrt*ctime))) &
            *sin(x2b(j))**2
          diib2(j,k)=1.D0/x1b(ism2)**3 &
            +2.D0*a0**2*p0/x1b(ism2)**3 &
            *(4.D0-(lambda*x1b(ism2)/etasqrt/ctime)**2) &
            *(p0+func(lambda*x1b(ism2)/(etasqrt*ctime))) &
            *sin(x2b(j))**2
        enddo
      enddo

      endif

      if(myid1 .eq. nproc1-1) then

      do k=ks,ke
        do j=js,je
          doib1(j,k) = d(ie,j,k)
          doib2(j,k) = d(ie,j,k)
        enddo
      enddo

      endif

end subroutine getdib
!=========================================================
subroutine getdjb
  use ModPar
  use ModGrid
  use ModField,          ONLY: d
  use ModBoundary

  implicit none

  integer :: k, i, kbc
!---------------------------------------------------------
      if(myid2 .eq. 0) then

      do i=ism2,iep2
        do k=ks,ke
          kbc=ks+mod(k-ks+(knmax-5)/2,knmax-5)
          dijb1(k,i) = d(i,js,kbc)
          dijb2(k,i) = d(i,jsp1,kbc)
        enddo
      enddo

      endif

      if(myid2 .eq. nproc2-1) then

      do i=ism2,iep2
        do k=ks,ke
          kbc=ks+mod(k-ks+(knmax-5)/2,knmax-5)
          dojb1(k,i) = d(i,je,kbc)
          dojb2(k,i) = d(i,jem1,kbc)
        enddo
      enddo

      endif

end subroutine getdjb
!=========================================================
subroutine getdkb

end subroutine getdkb
!=========================================================
subroutine getv1ib
  use ModPar
  use ModGrid
  use ModField,          ONLY: v1
  use ModBoundary
  use ModSundry,         ONLY: ctime

  implicit none

  integer j, k
  real q1

!----------------------------------------------------------------
      if(myid1 .eq. 0) then

      do k=ks,ke
        do j=js,je
           v1iib1(j,k)=1.D0/ctime*x1a(is)
           v1iib2(j,k)=1.D0/ctime*x1a(ism1)
           v1iib3(j,k)=1.D0/ctime*x1a(ism2)
        enddo
      enddo

      endif

      if(myid1 .eq. nproc1-1) then

      do k=ks,ke
        do j=js,je
          q1           = sign ( 0.5D0, v1(ie,j,k) )
          v1oib1(j,k) = v1(ie,j,k) * ( 0.5D0 + q1 )
          v1oib2(j,k) = v1oib1(j,k)
          v1oib3(j,k) = v1oib2(j,k)
        enddo
      enddo

      endif

end subroutine getv1ib
!=========================================================
subroutine getv1jb
  use ModPar
  use ModGrid
  use ModField,          ONLY: v1
  use ModBoundary

  implicit none

  integer :: k, i, kbc
!---------------------------------------------------------
      if(myid2 .eq. 0) then

      do i=ism2,iep3
        do k=ks,ke
          kbc=ks+mod(k-ks+(knmax-5)/2,knmax-5)
          v1ijb0(k,i) = v1(i,js,k)
          v1ijb1(k,i) = v1(i,js,kbc)
          v1ijb2(k,i) = v1(i,jsp1,kbc)
        enddo
      enddo

      endif

      if(myid2 .eq. nproc2-1) then

      do i=ism2,iep3
        do k=ks,ke
          kbc=ks+mod(k-ks+(knmax-5)/2,knmax-5)
          v1ojb1(k,i) = v1(i,je,kbc)
          v1ojb2(k,i) = v1(i,jem1,kbc)
        enddo
      enddo

      endif

end subroutine getv1jb
!=========================================================
subroutine getv1kb

end subroutine getv1kb
!=========================================================
subroutine getv2ib
  use ModPar
  use ModGrid
  use ModField,          ONLY: v2
  use ModBoundary

  implicit none

  integer j, k

!---------------------------------------------------------
      if(myid1 .eq. 0) then

      do k=ks,ke
        do j=js,jep1
          v2iib0(j,k)=v2(is,j,k)
          v2iib1(j,k)=0.D0
          v2iib2(j,k)=0.D0
        enddo
      enddo

      endif

      if(myid1 .eq. nproc1-1) then

      do k=ks,ke
        do j=js,jep1
          v2oib1(j,k) = v2(ie,j,k)
          v2oib2(j,k) = v2(ie,j,k)
        enddo
      enddo

      endif

end subroutine getv2ib
!=========================================================
subroutine getv2jb
  use ModPar
  use ModGrid
  use ModField,          ONLY: v2
  use ModBoundary

  implicit none

  integer :: k, i, kbc
!---------------------------------------------------------
      if(myid2 .eq. 0) then

      do i=ism2,iep2
        do k=ks,ke
          kbc=ks+mod(k-ks+(knmax-5)/2,knmax-5)
          v2ijb1(k,i) = 0.D0
          v2ijb2(k,i) = -v2(i,jsp1,kbc)
          v2ijb3(k,i) = -v2(i,jsp2,kbc)
        enddo
      enddo

      endif

      if(myid2 .eq. nproc2-1) then

      do i=ism2,iep2
        do k=ks,ke
          kbc=ks+mod(k-ks+(knmax-5)/2,knmax-5)
          v2ojb1(k,i) = 0.D0
          v2ojb2(k,i) = -v2(i,je,kbc)
          v2ojb3(k,i) = -v2(i,jem1,kbc)
        enddo
      enddo

      endif

end subroutine getv2jb
!=========================================================
subroutine getv2kb

end subroutine getv2kb
!=========================================================
subroutine getv3ib
  use ModPar
  use ModGrid
  use ModField,          ONLY: v3
  use ModBoundary

  implicit none

  integer j, k

!----------------------------------------------------------------
      if(myid1 .eq. 0) then

      do k=ks,kep1
        do j=js,je
          v3iib0(j,k)=v3(is,j,k)
          v3iib1(j,k)=0.D0
          v3iib2(j,k)=0.D0
        enddo
      enddo

      endif

      if(myid1 .eq. nproc1-1) then

      do k=ks,kep1
        do j=js,je
          v3oib1(j,k) = v3(ie,j,k)
          v3oib2(j,k) = v3(ie,j,k)
        enddo
      enddo

      endif

end subroutine getv3ib
!=========================================================
subroutine getv3jb
  use ModPar
  use ModGrid
  use ModField,          ONLY: v3
  use ModBoundary

  implicit none

  integer :: k, i, kbc
!---------------------------------------------------------
      if(myid2 .eq. 0) then

      do i=ism2,iep2
        do k=ks,kep1
          kbc=ks+mod(k-ks+(knmax-5)/2,knmax-5)
          v3ijb0(k,i) = v3(i,js,k)
          v3ijb1(k,i) = -v3(i,js,kbc)
          v3ijb2(k,i) = -v3(i,jsp1,kbc)
        enddo
      enddo

      endif

      if(myid2 .eq. nproc2-1) then

      do i=ism2,iep2
        do k=ks,kep1
          kbc=ks+mod(k-ks+(knmax-5)/2,knmax-5)
          v3ojb1(k,i) = -v3(i,je,kbc)
          v3ojb2(k,i) = -v3(i,jem1,kbc)
        enddo
      enddo

      endif

end subroutine getv3jb
!=========================================================
subroutine getv3kb

end subroutine getv3kb
!=========================================================
subroutine getdelsib

end subroutine getdelsib
!=========================================================
subroutine getdelsjb

end subroutine getdelsjb
!=========================================================
subroutine getdelskb

end subroutine getdelskb
!=========================================================
subroutine getbib

end subroutine getbib
!=========================================================
subroutine getbjb

end subroutine getbjb
!=========================================================
subroutine getbkb

end subroutine getbkb
!=========================================================
subroutine userstep

end subroutine userstep
!=========================================================
subroutine coord_xtoxp(rc,thc,phc,tilt,xn,yn,zn,xp,yp,zp)
  implicit none
  real, intent(in) :: rc,thc,phc,tilt,xn,yn,zn
  real, intent(out) :: xp,yp,zp
  real :: a11,a12,a13,a21,a22,a23,a31,a32,a33

       a11=cos(phc)*sin(thc)
       a12=-sin(tilt)*cos(phc)*cos(thc)-cos(tilt)*sin(phc)
       a13=-cos(tilt)*cos(phc)*cos(thc)+sin(tilt)*sin(phc)
       a21=sin(phc)*sin(thc)
       a22=-sin(tilt)*sin(phc)*cos(thc)+cos(tilt)*cos(phc)
       a23=-cos(tilt)*sin(phc)*cos(thc)-sin(tilt)*cos(phc)
       a31=cos(thc)
       a32=sin(tilt)*sin(thc)
       a33=cos(tilt)*sin(thc)

       xp=xn*a11+yn*a21+zn*a31
       yp=xn*a12+yn*a22+zn*a32
       zp=xn*a13+yn*a23+zn*a33

end subroutine coord_xtoxp
!=========================================================
subroutine vect_xptox(rc,thc,phc,tilt,bxp,byp,bzp,bx,by,bz)
  implicit none
  real, intent(in) :: rc,thc,phc,tilt,bxp,byp,bzp
  real, intent(out) :: bx,by,bz
  real :: a11,a12,a13,a21,a22,a23,a31,a32,a33

       a11=sin(thc)*cos(phc)
       a12=sin(thc)*sin(phc)
       a13=cos(thc)
       a21=-cos(tilt)*sin(phc)-sin(tilt)*cos(thc)*cos(phc)
       a22=cos(tilt)*cos(phc)-sin(tilt)*cos(thc)*sin(phc)
       a23=sin(tilt)*sin(thc)
       a31=sin(tilt)*sin(phc)-cos(tilt)*cos(thc)*cos(phc)
       a32=-sin(tilt)*cos(phc)-cos(tilt)*cos(thc)*sin(phc)
       a33=cos(tilt)*sin(thc)

       bx=bxp*a11+byp*a21+bzp*a31
       by=bxp*a12+byp*a22+bzp*a32
       bz=bxp*a13+byp*a23+bzp*a33

end subroutine vect_xptox
!=========================================================
subroutine message_pass_nodes_usr(idtop, idbot, nI1, nJ1, nK1, nI2, nJ2, nK2, &
  sendtop, sendbot, recvtop, recvbot)
    implicit none
    include "mpif.h"
    
    integer, intent(in)  :: idtop, idbot, nI1, nJ1, nK1, nI2, nJ2, nK2
    real,    intent(in)  :: sendtop(nI1,nJ1,nK1), sendbot(nI2,nJ2,nK2)
    real,    intent(out) :: recvtop(nI2,nJ2,nK2), recvbot(nI1,nJ1,nK1)
    
    integer :: ireq(4), istatus_arr(MPI_STATUS_SIZE,4), ierr
    integer :: ncount_top, ncount_bot
    real, allocatable, dimension(:) :: sendtopbuf, sendbotbuf, &
       recvtopbuf, recvbotbuf 
!-----------------------------------------------------------------------
    ncount_top = nI1*nJ1*nK1
    ncount_bot = nI2*nJ2*nK2
    allocate(sendtopbuf(1:ncount_top), sendbotbuf(1:ncount_bot))
    allocate(recvbotbuf(1:ncount_top), recvtopbuf(1:ncount_bot))
    sendtopbuf = reshape(sendtop,(/ncount_top/))
    sendbotbuf = reshape(sendbot,(/ncount_bot/))

    call MPI_IRECV(recvtopbuf, ncount_bot, MPI_DOUBLE_PRECISION, idtop, 1, &
         MPI_COMM_WORLD, ireq(2), ierr)
    call MPI_IRECV(recvbotbuf, ncount_top, MPI_DOUBLE_PRECISION, idbot, 2, &
         MPI_COMM_WORLD, ireq(4), ierr)
    call MPI_ISEND(sendbotbuf, ncount_bot, MPI_DOUBLE_PRECISION, idbot, 1, &
         MPI_COMM_WORLD, ireq(1), ierr)
    call MPI_ISEND(sendtopbuf, ncount_top, MPI_DOUBLE_PRECISION, idtop, 2, &
         MPI_COMM_WORLD, ireq(3), ierr)
    call MPI_WAITALL(4,ireq,istatus_arr,ierr)

    recvtop = reshape(recvtopbuf, (/nI2, nJ2, nK2/))
    recvbot = reshape(recvbotbuf,(/nI1, nJ1, nK1/))

    deallocate(sendtopbuf, sendbotbuf, recvbotbuf, recvtopbuf)

end subroutine message_pass_nodes_usr

real function func(x)
    implicit none

    real, intent(in) :: x
!------------------------------------------------------------------------------

    func=sin(x)/x-cos(x)

end function func

real function funcp(x)
    implicit none

    real, intent(in) :: x
!------------------------------------------------------------------------------

    funcp=cos(x)/x-sin(x)/x**2+sin(x)

end function funcp

end module ModUserSetup
