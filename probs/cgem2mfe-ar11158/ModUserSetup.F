module ModUserSetup

  implicit none

  private

  real :: ovhp0, as2, d0, r_c, tis0, tau_tr
  real, allocatable :: ebase0(:,:), ebase1(:,:), ebase2(:,:), &
                       dbase0(:,:), dbase1(:,:), dbase2(:,:)

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
  use ModGrid
  implicit none

  include "mpif.h"

  integer :: i, j, k, in_in, jn_in, kn_in, ierr

!-------------------------------------------------------

  if(myid .eq. 0) then
    open(unit=13, file='grid_in.dat', &
      form='unformatted',access='stream',status='old')
    read(13) in_in,jn_in,kn_in
    if(in_in .ne. inmax .or. jn_in .ne. jnmax &
      .or. kn_in .ne. knmax) then
      write(6,*) 'in_in jn_in kn_in= ', &
        in_in, jn_in, kn_in, 'but inmax jnmax knmax = ', &
        inmax,jnmax,knmax
      call MPI_ABORT(MPI_COMM_WORLD, 1,ierr)
    endif
    read(13) (xxa(i),i=1,inmax)
    read(13) (xxb(i),i=1,inmax-1)
    read(13) (yya(j),j=1,jnmax)
    read(13) (yyb(j),j=1,jnmax-1)
    read(13) (zza(k),k=1,knmax)
    read(13) (zzb(k),k=1,knmax-1)

    read(13) (dxxa(i),i=1,inmax-1)
    read(13) (dxxb(i),i=1,inmax)
    read(13) (dyya(j),j=1,jnmax-1)
    read(13) (dyyb(j),j=1,jnmax)
    read(13) (dzza(k),k=1,knmax-1)
    read(13) (dzzb(k),k=1,knmax)

    dxxa(inmax)=dxxa(inmax-1)
    xxb(inmax)=xxa(inmax)+0.5D0*dxxa(inmax)
    dyya(jnmax)=dyya(jnmax-1)
    yyb(jnmax)=yya(jnmax)+0.5D0*dyya(jnmax)
    dzza(knmax)=dzza(knmax-1)
    zzb(knmax)=zza(knmax)+0.5D0*dzza(knmax)

    close(13)
  endif
  call MPI_BCAST(xxa, inmax, MPI_DOUBLE_PRECISION, 0, &
       MPI_COMM_WORLD,ierr)
  call MPI_BCAST(xxb, inmax, MPI_DOUBLE_PRECISION, 0, &
       MPI_COMM_WORLD,ierr)
  call MPI_BCAST(yya, jnmax, MPI_DOUBLE_PRECISION, 0, &
       MPI_COMM_WORLD,ierr)
  call MPI_BCAST(yyb, jnmax, MPI_DOUBLE_PRECISION, 0, &
       MPI_COMM_WORLD,ierr)
  call MPI_BCAST(zza, knmax, MPI_DOUBLE_PRECISION, 0, &
       MPI_COMM_WORLD,ierr)
  call MPI_BCAST(zzb, knmax, MPI_DOUBLE_PRECISION, 0, &
       MPI_COMM_WORLD,ierr)
  call MPI_BCAST(dxxa, inmax, MPI_DOUBLE_PRECISION, 0, &
       MPI_COMM_WORLD,ierr)
  call MPI_BCAST(dxxb, inmax, MPI_DOUBLE_PRECISION, 0, &
       MPI_COMM_WORLD,ierr)
  call MPI_BCAST(dyya, jnmax, MPI_DOUBLE_PRECISION, 0, &
       MPI_COMM_WORLD,ierr)
  call MPI_BCAST(dyyb, jnmax, MPI_DOUBLE_PRECISION, 0, &
       MPI_COMM_WORLD,ierr)
  call MPI_BCAST(dzza, knmax, MPI_DOUBLE_PRECISION, 0, &
       MPI_COMM_WORLD,ierr)
  call MPI_BCAST(dzzb, knmax, MPI_DOUBLE_PRECISION, 0, &
       MPI_COMM_WORLD,ierr)

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
    dg32bd2(j) = cos(x2a(j))
    dg32ad2(j) = cos(x2b(j))
  enddo

  do i = 1, in
    dvl1a(i)=x1b(i)*x1b(i)*dx1a(i)
  enddo
  do i = 1, in
    dvl1b(i)=x1a(i)*x1a(i)*dx1b(i)
  enddo

  do j = 1, jn
    dvl2a(j) = sin(x2b(j))*dx2a(j)
  enddo
  do j = 1, jn
    dvl2b(j) = sin(x2a(j))*dx2b(j)
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

  integer :: i,j,k,ierr
  real :: msol,gacc_c,rho_c,temp_c,as2_c,hp_c
!
!-----------------------------------------------------------------
! specify boundary condition
!-----------------------------------------------------------------
!
       niib=3
       noib=2
       nijb=2
       nojb=2
       nikb=2
       nokb=2
       isplt=1
!
!-----------------------------------------------------
!      initialize fields or use restart file
!-----------------------------------------------------
!
!      specify gravitaional acceleration gacc
!      specify temp_c and rho_c
!
       msol=1.99D33
       r_c=xxa(is)*unit_len
       gacc_c=g_const*msol/r_c/r_c
       gacc=gacc_c/(unit_v*unit_v/unit_len)*(r_c/unit_len)**2
       d0=2.D0
       rho_c=d0*unit_rho
       temp_c=1.6D6
       rstar=kboltz/(muconst*mproton)
       as2_c=rstar*temp_c
       hp_c=as2_c/gacc_c

       as2=as2_c/unit_v/unit_v
       ovhp0=unit_len/hp_c

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
!      use restart file
!
#ifdef THCONDUCT
       q=0.D0
#endif
       call readrst_mpi
       ctime=time
       do k=ks,ke
         do j=js,je
           do i=is,ie
             tfield(i,j,k)=eint(i,j,k)*(gamma-1.D0) &
                /d(i,j,k)*unit_v**2/(rstar*unit_temp)
           enddo
         enddo
       enddo
!
!      set tau_tr and tis0
!
       tau_tr=0.1D0
       tis0=temp_c/unit_temp &
           *(1.D0-(1.D0-1.D0/gamma)*ovhp0*r_c/unit_len &
           *(1.D0-r_c/(xxb(is+1)*unit_len)))
!
!      allocate space for lower boundary arrays
!
       allocate(ebase0(jn,kn), ebase1(jn,kn),ebase2(jn,kn), &
         dbase0(jn,kn), dbase1(jn,kn),dbase2(jn,kn))
!
!      set current values of ghostzones for e and d
!
       irst=1
       do k=ks,ke
        do j=js,je
          ebase0(j,k)=eint(is,j,k)
          dbase0(j,k)=d(is,j,k)
          ebase1(j,k)=eint(is-1,j,k)
          dbase1(j,k)=d(is-1,j,k)
          ebase2(j,k)=eint(is-2,j,k)
          dbase2(j,k)=d(is-2,j,k)
        enddo
       enddo
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
       alpha_d = 1.25D0
       alpha_v = 1.D0
       alpha_e = 1.25D0
       nlf_d = 2
       nlf_v = 0
       nlf_e = 2
       clight = 8.D0
       ovc2 = 1.D0/(clight**2)
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
  use ModPar,      ONLY: in,jn,kn
  implicit none
  real, dimension(in,jn,kn) :: v2b3mv3b2

end subroutine getemf1jb
!=========================================================
subroutine getdiffemf1jb( v2b3mv3b2 )
  use ModPar,      ONLY: in,jn,kn
  implicit none
  real, dimension(in,jn,kn) :: v2b3mv3b2

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
  use ModPar,      ONLY: in,jn,kn
  implicit none
  real, dimension(in,jn,kn) :: v3b1mv1b3

end subroutine getemf2jb
!=========================================================
subroutine getdiffemf2jb( v3b1mv1b3 )
  use ModPar,      ONLY: in,jn,kn
  implicit none
  real, dimension(in,jn,kn) :: v3b1mv1b3

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
  use ModPar,      ONLY: in,jn,kn
  implicit none
  real, dimension(in,jn,kn) :: v1b2mv2b1

end subroutine getemf3jb
!=========================================================
subroutine getdiffemf3jb( v1b2mv2b1 )
  use ModPar,      ONLY: in,jn,kn
  implicit none
  real, dimension(in,jn,kn) :: v1b2mv2b1

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
  use ModPhysics
  use ModGrid
  use ModField,          ONLY: eint, tfield
  use ModBoundary
  use ModSundry,         ONLY: time, ctime
  use ModControl,        ONLY: irst

  implicit none

  integer :: j, k
  real :: eval,eval0,eval1,eval2
!---------------------------------------------------------
      if(myid1 .eq. 0) then

      if(irst .ne. 1) then
      eval0 = 1.D0/(gamma-1.D0) &
            *as2*d0*(1.D0-(1.D0-1.D0/gamma)*ovhp0*r_c/unit_len &
            *(1.D0-r_c/(x1b(is)*unit_len))) &
            **(gamma/(gamma-1.D0))
      eval1 = 1.D0/(gamma-1.D0) &
            *as2*d0*(1.D0-(1.D0-1.D0/gamma)*ovhp0*r_c/unit_len &
            *(1.D0-r_c/(x1b(ism1)*unit_len))) &
            **(gamma/(gamma-1.D0))
      eval2 = 1.D0/(gamma-1.D0) &
            *as2*d0*(1.D0-(1.D0-1.D0/gamma)*ovhp0*r_c/unit_len &
            *(1.D0-r_c/(x1b(ism2)*unit_len))) &
            **(gamma/(gamma-1.D0))
      do k=ks,ke
        do j=js,je
          eval=eval0*tfield(is+1,j,k)/tis0
          eint(is,j,k)=ebase0(j,k) &
             +min(1.D0,(ctime-time)/tau_tr)*(eval-ebase0(j,k))
          eval=eval1*tfield(is+1,j,k)/tis0
          eint_iib1(j,k)=ebase1(j,k) &
             +min(1.D0,(ctime-time)/tau_tr)*(eval-ebase1(j,k))
          eval=eval2*tfield(is+1,j,k)/tis0
          eint_iib2(j,k)=ebase2(j,k) &
             +min(1.D0,(ctime-time)/tau_tr)*(eval-ebase2(j,k))
        enddo
      enddo
      else
      do k=ks,ke
        do j=js,je
          eint(is,j,k)=ebase0(j,k)
          eint_iib1(j,k)=ebase1(j,k)
          eint_iib2(j,k)=ebase2(j,k)
        enddo
      enddo
      endif

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

end subroutine geteintjb
!=========================================================
subroutine geteintkb

end subroutine geteintkb
!=========================================================
subroutine getdib
  use ModPar
  use ModPhysics
  use ModGrid
  use ModField,          ONLY: d, tfield
  use ModBoundary
  use ModSundry,         ONLY: time, ctime
  use ModControl,        ONLY: irst

  implicit none

  integer :: j, k
  real :: dval,dval0,dval1,dval2
!---------------------------------------------------------
      if(myid1 .eq. 0) then

      if(irst .ne. 1) then
      dval0 = d0*(1.D0-(1.D0-1.D0/gamma)*ovhp0*r_c/unit_len &
            *(1.D0-r_c/(x1b(is)*unit_len))) &
            **(1.D0/(gamma-1.D0))
      dval1 = d0*(1.D0-(1.D0-1.D0/gamma)*ovhp0*r_c/unit_len &
            *(1.D0-r_c/(x1b(ism1)*unit_len))) &
            **(1.D0/(gamma-1.D0))
      dval2 = d0*(1.D0-(1.D0-1.D0/gamma)*ovhp0*r_c/unit_len &
            *(1.D0-r_c/(x1b(ism2)*unit_len))) &
            **(1.D0/(gamma-1.D0))
      do k=ks,ke
        do j=js,je
          dval=dval0*tfield(is+1,j,k)/tis0
          d(is,j,k)=dbase0(j,k) &
             +min(1.D0,(ctime-time)/tau_tr)*(dval-dbase0(j,k))
          dval=dval1*tfield(is+1,j,k)/tis0
          diib1(j,k)=dbase1(j,k) &
             +min(1.D0,(ctime-time)/tau_tr)*(dval-dbase1(j,k))
          dval=dval2*tfield(is+1,j,k)/tis0
          diib2(j,k)=dbase2(j,k) &
             +min(1.D0,(ctime-time)/tau_tr)*(dval-dbase2(j,k))
        enddo
      enddo
      else
      do k=ks,ke
        do j=js,je
          d(is,j,k)=dbase0(j,k)
          diib1(j,k)=dbase1(j,k)
          diib2(j,k)=dbase2(j,k)
        enddo
      enddo
      endif

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

end subroutine getdjb
!=========================================================
subroutine getdkb

end subroutine getdkb
!=========================================================
subroutine getv1ib
  use ModPar
  use ModGrid
  use ModField
  use ModBoundary

  implicit none

  integer j, k
  real q1

!----------------------------------------------------------------
      if(myid1 .eq. 0) then

      do k=ks,ke
        do j=js,je
            v1iib1(j,k)=0.D0
            v1iib2(j,k)=0.D0
            v1iib3(j,k)=0.D0
        enddo
      enddo

      endif

      if(myid1 .eq. nproc1-1) then

      do k=ks,ke
        do j=js,je
          q1           = sign ( 0.5D0, v1(ie,j,k) )
          v1oib1(j,k) = v1(ie,j,k) * ( 0.5D0 + q1 )
          v1oib2(j,k) = v1oib1(j,k)
        enddo
      enddo

      endif

end subroutine getv1ib
!=========================================================
subroutine getv1jb

end subroutine getv1jb
!=========================================================
subroutine getv1kb

end subroutine getv1kb
!=========================================================
subroutine getv2ib
  use ModPar
  use ModGrid
  use ModField
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

end subroutine getv2jb
!=========================================================
subroutine getv2kb

end subroutine getv2kb
!=========================================================
subroutine getv3ib
  use ModPar
  use ModGrid
  use ModField
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

end subroutine getv3jb
!=========================================================
subroutine getv3kb

end subroutine getv3kb
!=========================================================
subroutine getdelsib
  use ModPar
  use ModGrid
  use ModDel

  implicit none

  integer i, j, k
!---------------------------------------------------------
       if(myid1 .eq. 0) then

       do k=ksm2,kep2
         do j=jsm2,jep2
           dels1(is,j,k) = 0.D0
           dels2(is,j,k) = 0.D0
           dels3(is,j,k) = 0.D0
         enddo
       enddo

       endif

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
  use ModPhysics
  use ModGrid
  use ModField
  
  implicit none
  integer :: i, j, k
  real, parameter :: tceiling_dim=2.D8
!---------------------------------------------------------

  do k=ksm2,kep2
    do j=jsm2,jep2
      do i=ism2,iep2
        eint(i,j,k)=eint(i,j,k) &
          *min(1.D0,tceiling_dim/(tfield(i,j,k)*unit_temp))
      enddo
    enddo
  enddo
  
! use this to set current values of ghostzones for e and d
       
       do k=ks,ke
        do j=js,je
          ebase0(j,k)=eint(is,j,k)
          dbase0(j,k)=d(is,j,k)
          ebase1(j,k)=eint(is-1,j,k)
          dbase1(j,k)=d(is-1,j,k)
          ebase2(j,k)=eint(is-2,j,k)
          dbase2(j,k)=d(is-2,j,k)
        enddo
       enddo

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

end module ModUserSetup
