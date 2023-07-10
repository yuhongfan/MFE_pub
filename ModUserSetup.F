module ModUserSetup

  implicit none

contains

subroutine grid
  use ModPar
  use ModGrid

  implicit none

  integer :: i, j, k
  real :: x1min,x2min,x3min,x1max,x2max,x3max,delx1,delx2,delx3
  real :: PDMB = 4.0
  real*8, dimension(:,:,:,:), allocatable :: vol
  enum, bind(C)
        enumerator :: dir1=1, dir2, dir3, center
  real*8, parameter, dimension(0:3):: fac4th = [-1., 7., 7., -1. ]/12.
  allocate(vol(ism2:iep3,jsm2:jep3,ksm2:kep3,dir1:center))

!-------------------------------------------------------
!
      x1min=0.D0
      x1max=1.D0
      x2min=0.D0
      x2max=1.D0
      x3min=0.D0
      x3max=1.D0
      delx1=x1max/dble(inmax-5)
      delx2=x2max/dble(jnmax-5)
      delx3=x3max/dble(knmax-5)

      do i = 1, in
        x1a(i)=dble(i-is+myid1*(in-5))*delx1+x1min
        dx1a(i)=delx1
        x1b(i)=x1a(i)+0.5D0*delx1
        dx1b(i)=delx1
      enddo

      do j = 1, jn
        x2a(j)=dble(j-js+myid2*(jn-5))*delx2+x2min
        dx2a(j)=delx2
        x2b(j)=x2a(j)+0.5D0*delx2
        dx2b(j)=delx2
      enddo

      do k = 1, kn
        x3a(k)=dble(k-ks+myid3*(kn-5))*delx3+x3min
        dx3a(k)=delx3
        x3b(k)=x3a(k)+0.5D0*delx3
        dx3b(k)=delx3
      enddo

      do i = 1, in
        g2a(i)=1.D0
        g31a(i)=1.D0
        dg2bd1(i)=0.D0
        dg31bd1(i)=0.D0
        g2b(i)=1.D0
        g31b(i)=1.D0
        dg2ad1(i)=0.D0
        dg31ad1(i)=0.D0
      enddo

      do j = 1, jn
        g32a   (j) = 1.D0
        dg32bd2(j) = 0.D0
        g32b   (j) = 1.D0
        dg32ad2(j) = 0.D0
      enddo

      do i = 1, in
        dvl1a(i)=dx1a(i)
      enddo
      do i = 1, in
        dvl1b(i)=dx1b(i)
      enddo

      do j = 1, jn
        dvl2a(j) = dx2a(j)
      enddo
      do j = 1, jn
        dvl2b(j) = dx2b(j)
      enddo

      do k = 1, kn
        dvl3a(k) = dx3a(k)
      enddo
      do k = 1, kn
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

      do k=ksm2,kep3
        dx3bi (k) = 1.D0 / ( dx3b (k) + tiny )
        dvl3bi(k) = 1.D0 / ( dvl3b(k) + tiny )
      enddo
      do k=ksm2,kep2
        dx3ai (k) = 1.D0 / ( dx3a (k) + tiny )
        dvl3ai(k) = 1.D0 / ( dvl3a(k) + tiny )
      enddo

      o k=1,kn
    do j=1,jn
        do i=i,in
            vol(i,j,k,center) = dx1a(i)*dx2a(j)*dx3a(k)
        enddo
    enddo
enddo

do k=ks, ke
        do j=js, je
            do i=is,iep1
                vol(i,j,k,dir1) = sum(vol(i-2:i+1,j,k,center)*fac4th(0:3))
            enddo
        enddo
    enddo
    do k=ks, ke
        do j=js, jep1
            do i=is, ie
                vol(i,j,k,dir2) = sum(vol(i,j-2:j+1,k,center)*fac4th(0:3))
            enddo
        enddo
    enddo
    do k=ks, kep1
        do j=js, je
            do i=is, ie
                vol(i,j,k,dir3) = sum(vol(i,j,k-2:k+1,center)*fac4th(0:3))
            enddo
        enddo
    enddo
end subroutine grid

subroutine field_init
  use ModPar
  use ModPhysics
  use ModGrid
  use ModControl,    ONLY: irst
  use ModSundry
  use ModField
  use ModBoundary
  use ModIO

  implicit none

  include "mpif.h"

  integer :: i,j,k,ierr
  real :: msol,rsol,gacc_c,rho_c,temp_c
  real :: afunc(in,jn,kn)
  real :: b0
!
!-----------------------------------------------------------------
! specify boundary condition
!-----------------------------------------------------------------
!
       niib=4
       noib=4
       nijb=4
       nojb=4
       nikb=4
       nokb=4
       isplt=1
!
!-----------------------------------------------------
!      initialize fields or use restart file
!-----------------------------------------------------
!
!      specify gravitaional acceleration gacc
!      specify temp_c and rho_c
!      (these are not used for this case)
       msol=1.99D33
       rsol=6.96D10
       gacc_c=g_const*msol/rsol/rsol
       rho_c=2.7D-7
       temp_c=5100.D0
       gacc=gacc_c/(unit_v*unit_v/unit_len)
!
       rstar=kboltz/(muconst*mproton)
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
       write(18) rsol
       write(18) temp_c
       write(18) rho_c
       write(18) muconst
       close(18)
       endif

       if(irst .ne. 1) then

       time=0.D0
       ctime=time
!
!      initialize v1 field
!
       do k=ks,ke
         do j=js,je
           do i=is,iep1
             v1(i,j,k)=-sin(2.D0*pi*x2b(j))
           enddo
         enddo
       enddo
!
!      initialize v3 field
!
       do k=ks,kep1
         do j=js,je
           do i=is,ie
             v3(i,j,k)=0.D0
           enddo
         enddo
       enddo
!
!      initialize v2 field
!
       do k=ks,ke
         do j=js,jep1
           do i=is,ie
             v2(i,j,k)=sin(2.D0*pi*x1b(i))
           enddo
         enddo
       enddo
!
!      afunc
!
       b0=1.D0/sqrt(4.D0*pi)
       do k=ksm2,kep3
         do j=jsm2,jep3
           do i=ism2,iep3
             afunc(i,j,k)=b0*(cos(4.D0*pi*x1a(i))/(4.D0*pi) &
               +cos(2.D0*pi*x2a(j))/(2.D0*pi))
           enddo
         enddo
       enddo
!
!      initialize b1 field
!
       do k=ksm2,kep2
         do j=jsm2,jep2
           do i=ism2,iep2
             b1(i,j,k)=(afunc(i,j+1,k)-afunc(i,j,k))*g2ai(i)*dx2ai(j)
           enddo
         enddo
       enddo
!
!      initialize b2 field
!
       do k=ksm2,kep2
         do j=jsm2,jep2
           do i=ism2,iep2
             b2(i,j,k)=-(afunc(i+1,j,k)-afunc(i,j,k))*dx1ai(i)
           enddo
         enddo
       enddo
!
!      initialize b3 field
!
       do k=ksm2,kep2
         do j=jsm2,jep2
           do i=ism2,iep2
             b3(i,j,k)=0.D0
           enddo
         enddo
       enddo
!
!      initialize d
!
       do k=ksm2,kep2
         do j=jsm2,jep2
           do i=ism2,iep2
             d(i,j,k)=25.D0/(36.D0*pi)
           enddo
         enddo
       enddo
!
!      initialize e
!
       do k=ksm2,kep2
         do j=jsm2,jep2
           do i=ism2,iep2
             eint(i,j,k)=5.D0/(12.D0*pi)/(gamma-1.D0)
           enddo
         enddo
       enddo
!
       else
!
!      use restart file
!
       call readrst_mpi
       ctime=time
!
       endif

       dfloor=1.D-7
       efloor=1.D-7
       if (myid .eq. 0) then
         dt=dx1a(is)
       endif
       call MPI_BCAST(dt, 1, MPI_DOUBLE_PRECISION, 0, &
         MPI_COMM_WORLD,ierr)
       dtmdi2 = 0.D0
       alpha_d = 1.D0
       alpha_v = 1.D0
       nlf_d = 4
       nlf_v = 4
!
       itnow=0
       ntcond=0
       courno=0.7D0

end subroutine field_init

subroutine getemf1ib( v2b3mv3b2 )
  use ModPar,      ONLY: in,jn,kn
  implicit none
  real, dimension(in,jn,kn) :: v2b3mv3b2

end subroutine getemf1ib
!---------------------------------------------------------
subroutine getdiffemf1ib( v2b3mv3b2 )
  use ModPar,      ONLY: in,jn,kn
  implicit none
  real, dimension(in,jn,kn) :: v2b3mv3b2

end subroutine getdiffemf1ib
!---------------------------------------------------------
subroutine getemf1jb( v2b3mv3b2 )
  use ModPar,      ONLY: in,jn,kn
  implicit none
  real, dimension(in,jn,kn) :: v2b3mv3b2

end subroutine getemf1jb
!---------------------------------------------------------
subroutine getdiffemf1jb( v2b3mv3b2 )
  use ModPar,      ONLY: in,jn,kn
  implicit none
  real, dimension(in,jn,kn) :: v2b3mv3b2

end subroutine getdiffemf1jb
!---------------------------------------------------------
subroutine getemf1kb( v2b3mv3b2 )
  use ModPar,      ONLY: in,jn,kn
  implicit none
  real, dimension(in,jn,kn) :: v2b3mv3b2

end subroutine getemf1kb
!---------------------------------------------------------
subroutine getdiffemf1kb( v2b3mv3b2 )
  use ModPar,      ONLY: in,jn,kn
  implicit none
  real, dimension(in,jn,kn) :: v2b3mv3b2

end subroutine getdiffemf1kb
!---------------------------------------------------------
subroutine getemf2ib( v3b1mv1b3 )
  use ModPar,      ONLY: in,jn,kn
  implicit none
  real, dimension(in,jn,kn) :: v3b1mv1b3

end subroutine getemf2ib
!---------------------------------------------------------
subroutine getdiffemf2ib( v3b1mv1b3 )
  use ModPar,      ONLY: in,jn,kn
  implicit none
  real, dimension(in,jn,kn) :: v3b1mv1b3

end subroutine getdiffemf2ib
!---------------------------------------------------------
subroutine getemf2jb( v3b1mv1b3 )
  use ModPar,      ONLY: in,jn,kn
  implicit none
  real, dimension(in,jn,kn) :: v3b1mv1b3

end subroutine getemf2jb
!---------------------------------------------------------
subroutine getdiffemf2jb( v3b1mv1b3 )
  use ModPar,      ONLY: in,jn,kn
  implicit none
  real, dimension(in,jn,kn) :: v3b1mv1b3

end subroutine getdiffemf2jb
!---------------------------------------------------------
subroutine getemf2kb( v3b1mv1b3 )
  use ModPar,      ONLY: in,jn,kn
  implicit none
  real, dimension(in,jn,kn) :: v3b1mv1b3

end subroutine getemf2kb
!---------------------------------------------------------
subroutine getdiffemf2kb( v3b1mv1b3 )
  use ModPar,      ONLY: in,jn,kn
  implicit none
  real, dimension(in,jn,kn) :: v3b1mv1b3

end subroutine getdiffemf2kb
!---------------------------------------------------------
subroutine getemf3ib( v1b2mv2b1 )
  use ModPar,      ONLY: in,jn,kn
  implicit none
  real, dimension(in,jn,kn) :: v1b2mv2b1

end subroutine getemf3ib
!---------------------------------------------------------
subroutine getdiffemf3ib( v1b2mv2b1 )
  use ModPar,      ONLY: in,jn,kn
  implicit none
  real, dimension(in,jn,kn) :: v1b2mv2b1

end subroutine getdiffemf3ib
!---------------------------------------------------------
subroutine getemf3jb( v1b2mv2b1 )
  use ModPar,      ONLY: in,jn,kn
  implicit none
  real, dimension(in,jn,kn) :: v1b2mv2b1

end subroutine getemf3jb
!---------------------------------------------------------
subroutine getdiffemf3jb( v1b2mv2b1 )
  use ModPar,      ONLY: in,jn,kn
  implicit none
  real, dimension(in,jn,kn) :: v1b2mv2b1

end subroutine getdiffemf3jb
!---------------------------------------------------------
subroutine getemf3kb( v1b2mv2b1 )
  use ModPar,      ONLY: in,jn,kn
  implicit none
  real, dimension(in,jn,kn) :: v1b2mv2b1

end subroutine getemf3kb
!---------------------------------------------------------
subroutine getdiffemf3kb( v1b2mv2b1 )
  use ModPar,      ONLY: in,jn,kn
  implicit none
  real, dimension(in,jn,kn) :: v1b2mv2b1

end subroutine getdiffemf3kb
!---------------------------------------------------------
subroutine getbib

end subroutine getbib
!---------------------------------------------------------
subroutine getbjb

end subroutine getbjb
!---------------------------------------------------------
subroutine getbkb

end subroutine getbkb
!---------------------------------------------------------
subroutine userstep

end subroutine userstep
!---------------------------------------------------------
subroutine geteib

end subroutine geteib
!---------------------------------------------------------
subroutine getejb

end subroutine getejb
!---------------------------------------------------------
subroutine getekb

end subroutine getekb
!---------------------------------------------------------
subroutine geteintib

end subroutine geteintib
!---------------------------------------------------------
subroutine geteintjb

end subroutine geteintjb
!---------------------------------------------------------
subroutine geteintkb

end subroutine geteintkb
!---------------------------------------------------------
subroutine getdib

end subroutine getdib
!---------------------------------------------------------
subroutine getdjb

end subroutine getdjb
!---------------------------------------------------------
subroutine getdkb

end subroutine getdkb
!---------------------------------------------------------
subroutine getv1ib

end subroutine getv1ib
!---------------------------------------------------------
subroutine getv1jb

end subroutine getv1jb
!---------------------------------------------------------
subroutine getv1kb

end subroutine getv1kb
!---------------------------------------------------------
subroutine getv2ib

end subroutine getv2ib
!---------------------------------------------------------
subroutine getv2jb

end subroutine getv2jb
!---------------------------------------------------------
subroutine getv2kb

end subroutine getv2kb
!---------------------------------------------------------
subroutine getv3ib

end subroutine getv3ib
!---------------------------------------------------------
subroutine getv3jb

end subroutine getv3jb
!---------------------------------------------------------
subroutine getv3kb

end subroutine getv3kb
!---------------------------------------------------------
subroutine getdelsib

end subroutine getdelsib
!---------------------------------------------------------
subroutine getdelsjb

end subroutine getdelsjb
!---------------------------------------------------------
subroutine getdelskb

end subroutine getdelskb
!---------------------------------------------------------

end module ModUserSetup
