module ModUserSetup

  implicit none

  private

  real :: as2, d0
  real, allocatable :: ebase0(:,:), dbase0(:,:)

  real, allocatable :: allemf2(:,:,:), allemf3(:,:,:)
  real, allocatable :: twistemf2(:,:,:), twistemf3(:,:,:)
  real, parameter :: pdfidt_dim=12.D0*60.D0
  real, parameter :: vtwist_dim=2.5D5
  integer, parameter :: ipdfistart=479
  integer :: ipdfi

  integer, parameter :: ntsub=50
  real, parameter :: ampxi=3.D-6
  real, allocatable :: xislc1(:,:,:), psi(:,:), bsurf(:,:,:)
  integer ifileread

  integer :: nttot=0
  real, allocatable :: tsub_arr(:)
  integer, allocatable :: itsub_arr(:), itcount_arr(:), ifilesub_arr(:)

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
  real :: msol,gacc_s,rho_c,temp_c,as2_c,hp_c
  real :: ovhp0, rsol
  real :: fctr_dim,pbase_dim,dbase

  character(len=80) :: fileallemf2, fileallemf3
  character(len=80) :: filetwistemf2, filetwistemf3
  character(len=4) :: idoutpdfi
!
!-----------------------------------------------------------------
! specify boundary condition
!-----------------------------------------------------------------
!
       niib=3
       noib=2
       nijb=1
       nojb=1
       nikb=1
       nokb=1
       isplt=1
       call get_timetable
!
!-----------------------------------------------------
!      initialize fields or use restart file
!-----------------------------------------------------
!
!      specify gravitaional acceleration gacc
!      specify temp_c and rho_c
!
       msol=1.99D33
       rsol=xxa(is)*unit_len
       gacc_s=g_const*msol/rsol/rsol
       gacc=gacc_s/(unit_v*unit_v/unit_len)*(rsol/unit_len)**2
       rstar=kboltz/(muconst*mproton)
       temp_c=2.D4
       rho_c=1.D12*mproton
       as2_c=rstar*temp_c
       hp_c=as2_c/gacc_s

       as2=as2_c/unit_v/unit_v
       ovhp0=unit_len/hp_c
       d0=rho_c/unit_rho

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
       write(18) rsol
       write(18) temp_c
       write(18) rho_c
       write(18) muconst
       close(18)
       endif
!
!      use restart file
!
       deltat = 12.D0*60.D0/unit_time
       q=0.D0
       call readrst_mpi
       if(irst .ne. 1) then
         time=deltat*ipdfistart
         ctime=time
       else
         ctime=time
       endif
       do k=ksm2,kep2
         do j=jsm2,jep2
           do i=ism2,iep2
             tfield(i,j,k)=eint(i,j,k)*(gamma-1.D0) &
                /d(i,j,k)*unit_v**2/(rstar*unit_temp)
           enddo
         enddo
       enddo
!
!      allocate space for lower boundary arrays
!
       allocate(ebase0(jn,kn), dbase0(jn,kn))
       allocate(allemf2(1,jn,kn),allemf3(1,jn,kn))
       allocate(twistemf2(1,jn,kn),twistemf3(1,jn,kn))
!
!      set current values of ghostzones for e and d
!
       do k=ksm2,kep2
        do j=jsm2,jep2
          ebase0(j,k)=eint(is,j,k)
          dbase0(j,k)=d(is,j,k)
        enddo
       enddo
       irst=1
!
!      reset b1, b2, b3 in ghostzones
!
       call getbib
!
!      parameters in sundry
!
       dfloor=1.D-5
       efloor=1.D-11
       if (myid .eq. 0) then
         dt=dx1a(is)
       endif
       call MPI_BCAST(dt, 1, MPI_DOUBLE_PRECISION, 0, &
         MPI_COMM_WORLD,ierr)
       dtmdi2 = 0.D0
       alpha_d = 1.25D0
       alpha_v = 1.D0
       alpha_e = 1.D0
       nlf_d = 2
       nlf_v = 0
       nlf_e = 0
       clight = 4.D0
       ovc2 = 1.D0/(clight**2)
!
!      read in first pdfi electric field
!
       ipdfi=int(ctime/deltat)
       write(idoutpdfi,'(i4.4)') ipdfi
       fileallemf2='allemf2_'//idoutpdfi//'.dat'
       fileallemf3='allemf3_'//idoutpdfi//'.dat'
       filetwistemf2='twistemf2_'//idoutpdfi//'.dat'
       filetwistemf3='twistemf3_'//idoutpdfi//'.dat'
       call read23emf2_mpi(allemf2,fileallemf2)
       call read23emf3_mpi(allemf3,fileallemf3)
       call read23emf2_mpi(twistemf2,filetwistemf2)
       call read23emf3_mpi(twistemf3,filetwistemf3)
!
!      set values for the e2kbs,e2kbe,e3jbs,e3jbe
!
       do j=jsm2,jep2
         e2kbs(j) = allemf2(1,j,ks)*xxa(is) &
           *pdfidt_dim/(deltat*unit_time)
         e2kbe(j) = allemf2(1,j,kep1)*xxa(is) &
           *pdfidt_dim/(deltat*unit_time)
       enddo
       do k=ksm2,kep2
         e3jbs(k) = allemf3(1,js,k)*xxa(is) &
           *pdfidt_dim/(deltat*unit_time)
         e3jbe(k) = allemf3(1,jep1,k)*xxa(is) &
           *pdfidt_dim/(deltat*unit_time)
       enddo
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
  use ModPhysics,  ONLY: unit_time, unit_v
  use ModGrid
  use ModField,    ONLY: b1
  use ModSundry,   ONLY: ctime, deltat
  use ModInterp,   ONLY: findpos
  use ModBoundary
  use ModIO

  implicit none

  real, intent(in), dimension(in,jn,kn) :: v3b1mv1b3

  character(len=80) :: filexislc1
  character(len=4) :: idoutsub

  integer :: j, k
  integer :: itsub, itcount, ifilesub
  real :: ctimetable

  integer :: idtop, idbot
  integer :: isendtop, isendbot, irecvtop, irecvbot
  real, allocatable :: recvtop(:,:,:), recvbot(:,:,:)
  real, allocatable :: sendtop(:,:,:), sendbot(:,:,:)
!---------------------------------------------------------

      ctimetable=ctime-deltat*ipdfistart
      call findpos(nttot,tsub_arr,itsub_arr,itcount_arr, &
        ifilesub_arr,ctimetable,itsub,itcount,ifilesub)
      if (ifileread .ne. ifilesub) then
        ifileread=ifilesub
        write(idoutsub,'(i4.4)') ifileread
        filexislc1='xislc1_23.'//idoutsub//'.dat'
        call read23qty_mpi(ntsub,xislc1,filexislc1)
      endif

      do k=ksm2,kep2
        do j=jsm2,jep2
          bsurf(is,j,k)=b1(is,j,k)
        enddo
      enddo

      idtop=myid1+(myid2+1)*nproc1+myid3*(nproc1*nproc2)
      idbot=myid1+(myid2-1)*nproc1+myid3*(nproc1*nproc2)
      isendtop=2
      isendbot=1
      irecvtop=1
      irecvbot=2
      if(myid2 .eq. 0) idbot = myid1 + (nproc2 - 1)*nproc1 &
                             + myid3*(nproc1*nproc2)
      if(myid2 .eq. nproc2-1) idtop = myid1 + 0*nproc1 &
                                    + myid3*(nproc1*nproc2)
      allocate(recvtop(is:is,jep3:jep3,ksm2:kep2), &
        recvbot(is:is,jsm2-1:jsm2-1,ksm2:kep2), &
        sendtop(is:is,jem2:jem2,ksm2:kep2), &
        sendbot(is:is,jsp2:jsp2,ksm2:kep2))
      sendtop=bsurf(is:is,jem2:jem2,ksm2:kep2)
      sendbot=bsurf(is:is,jsp2:jsp2,ksm2:kep2)
      call message_pass_nodes(idtop, idbot, &
        1, 1, kep2-ksm2+1, 1, 1, kep2-ksm2+1, &
        sendtop, sendbot, recvtop, recvbot, &
        isendtop, isendbot, irecvtop, irecvbot)
      bsurf(is:is,jep3:jep3,ksm2:kep2)=recvtop
      bsurf(is:is,jsm2-1:jsm2-1,ksm2:kep2)=recvbot
      deallocate(recvtop, recvbot, sendtop, sendbot)

      idtop = myid1 + myid2*nproc1 + (myid3+1)*(nproc1*nproc2)
      idbot = myid1 + myid2*nproc1 + (myid3-1)*(nproc1*nproc2)
      isendtop=2
      isendbot=1
      irecvtop=1
      irecvbot=2
      if(myid3==0) idbot = myid1 + myid2*nproc1 + (nproc3-1)*(nproc1*nproc2)
      if(myid3==nproc3-1) idtop = myid1 + myid2*nproc1 + 0*(nproc1*nproc2)
      allocate(recvtop(is:is,jsm2-1:jep3,kep3:kep3), &
        recvbot(is:is,jsm2-1:jep3,ksm2-1:ksm2-1))
      allocate(sendtop(is:is,jsm2-1:jep3,kem2:kem2), &
        sendbot(is:is,jsm2-1:jep3,ksp2:ksp2))
      sendtop=bsurf(is:is,jsm2-1:jep3,kem2:kem2)
      sendbot=bsurf(is:is,jsm2-1:jep3,ksp2:ksp2)
      call message_pass_nodes(idtop, idbot, &
        1, jep3-jsm2+2, 1, 1, jep3-jsm2+2, 1, &
        sendtop, sendbot, recvtop, recvbot, &
        isendtop, isendbot, irecvtop, irecvbot)
      bsurf(is:is,jsm2-1:jep3,kep3:kep3) = recvtop
      bsurf(is:is,jsm2-1:jep3,ksm2-1:ksm2-1) = recvbot
      deallocate(recvtop, recvbot, sendtop, sendbot)

      do k=ksm2,kep3
        do j=jsm2,jep3
          psi(j,k)=(xislc1(itcount,j,k)*(ctimetable-tsub_arr(itsub-1)) &
            +xislc1(itcount-1,j,k)*(tsub_arr(itsub)-ctimetable)) &
            /(tsub_arr(itsub)-tsub_arr(itsub-1)) &
            *0.25D0*(bsurf(is,j,k)+bsurf(is,j-1,k) &
         +bsurf(is,j,k-1)+bsurf(is,j-1,k-1))*ampxi
        enddo
      enddo

      if(myid1 .eq. 0) then

      do k=ksm2,kep3
        do j=jsm2,jep2
            emf2iib1(j,k)=-(psi(j+1,k)-psi(j,k))*dx2ai(j)*g2ai(is) &
              +allemf2(1,j,k)*pdfidt_dim/(deltat*unit_time) &
              +twistemf2(1,j,k)*vtwist_dim/unit_v
            emf2iib2(j,k)=emf2iib1(j,k)
            emf2iib3(j,k)=emf2iib1(j,k)
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
  use ModPhysics,   ONLY: unit_time, unit_v
  use ModGrid
  use ModSundry,   ONLY: deltat
  use ModBoundary

  implicit none

  real, intent(in), dimension(in,jn,kn) :: v1b2mv2b1
  integer :: j, k
!---------------------------------------------------------
      if(myid1 .eq. 0) then

      do k=ksm2,kep2
        do j=js,jep1
            emf3iib1(j,k)=-(psi(j,k+1)-psi(j,k)) &
              *dx3ai(k)*g31ai(is)*g32ai(j) &
              +allemf3(1,j,k)*pdfidt_dim/(deltat*unit_time) &
              +twistemf3(1,j,k)*vtwist_dim/unit_v
            emf3iib2(j,k)=emf3iib1(j,k)
            emf3iib3(j,k)=emf3iib1(j,k)
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
  use ModField,          ONLY: eint
  use ModBoundary

  implicit none

  integer :: j, k
!---------------------------------------------------------
      if(myid1 .eq. 0) then

      do k=ks,ke
        do j=js,je
          eint(is,j,k)=ebase0(j,k)
          eint_iib1(j,k)=ebase0(j,k)
          eint_iib2(j,k)=ebase0(j,k)
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

end subroutine geteintjb
!=========================================================
subroutine geteintkb

end subroutine geteintkb
!=========================================================
subroutine getdib
  use ModPar
  use ModPhysics
  use ModGrid
  use ModField,          ONLY: d
  use ModBoundary
  implicit none

  integer :: j, k
!---------------------------------------------------------
      if(myid1 .eq. 0) then

      do k=ks,ke
        do j=js,je
          d(is,j,k)=dbase0(j,k)
          diib1(j,k)=dbase0(j,k)
          diib2(j,k)=dbase0(j,k)
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
  use ModPar,      ONLY: myid1
  use ModGrid
  use ModField,    ONLY: b1, b2, b3

  implicit none

  integer :: j, k
!---------------------------------------------------------

  if(myid1 .eq. 0) then
    do k=ksm2,kep2
      do j=jsm2,jep2
        b1(ism1,j,k)=b1(is,j,k)
        b2(ism1,j,k)=b2(is,j,k)
        b3(ism1,j,k)=b3(is,j,k)
        b1(ism2,j,k)=b1(is,j,k)
        b2(ism2,j,k)=b2(is,j,k)
        b3(ism2,j,k)=b3(is,j,k)
      enddo
    enddo
  endif

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
  use ModSundry,        ONLY: time, dt, deltat, clight
  use ModBoundary
  use ModIO
  
  implicit none
  integer :: i, j, k, nnow
  real, parameter :: tceiling_dim=2.D8

  character(len=80) :: fileallemf2, fileallemf3
  character(len=80) :: filetwistemf2, filetwistemf3
  character(len=4) :: idoutpdfi
!---------------------------------------------------------

  do k=ksm2,kep2
    do j=jsm2,jep2
      do i=ism2,iep2
        eint(i,j,k)=eint(i,j,k) &
          *min(1.D0,tceiling_dim/(tfield(i,j,k)*unit_temp))
      enddo
    enddo
  enddo

  do k=ksm2,kep2
    do j=jsm2,jep2
      do i=ism2,iep3
        v1(i,j,k)=sign(1.,v1(i,j,k))*min(0.5*clight,abs(v1(i,j,k)))
      enddo
    enddo
  enddo

  do k=ksm2,kep2
    do j=jsm2,jep3
      do i=ism2,iep2
        v2(i,j,k)=sign(1.,v2(i,j,k))*min(0.5*clight,abs(v2(i,j,k)))
      enddo
    enddo
  enddo

  do k=ksm2,kep3
    do j=jsm2,jep2
      do i=ism2,iep2
        v3(i,j,k)=sign(1.,v3(i,j,k))*min(0.5*clight,abs(v3(i,j,k)))
      enddo
    enddo
  enddo
!
!      read in the pdfi electric field at time+dt
!
       nnow=int((time+dt)/deltat)
       if (nnow .ne. ipdfi) then
         ipdfi=nnow
         write(idoutpdfi,'(i4.4)') ipdfi
         fileallemf2='allemf2_'//idoutpdfi//'.dat'
         fileallemf3='allemf3_'//idoutpdfi//'.dat'
         filetwistemf2='twistemf2_'//idoutpdfi//'.dat'
         filetwistemf3='twistemf3_'//idoutpdfi//'.dat'
         call read23emf2_mpi(allemf2,fileallemf2)
         call read23emf3_mpi(allemf3,fileallemf3)
         call read23emf2_mpi(twistemf2,filetwistemf2)
         call read23emf3_mpi(twistemf3,filetwistemf3)
!
!      set values for the e2kbs,e2kbe,e3jbs,e3jbe
!
         do j=jsm2,jep2
           e2kbs(j) = allemf2(1,j,ks)*xxa(is) &
             *pdfidt_dim/(deltat*unit_time)
           e2kbe(j) = allemf2(1,j,kep1)*xxa(is) &
             *pdfidt_dim/(deltat*unit_time)
         enddo
         do k=ksm2,kep2
           e3jbs(k) = allemf3(1,js,k)*xxa(is) &
             *pdfidt_dim/(deltat*unit_time)
           e3jbe(k) = allemf3(1,jep1,k)*xxa(is) &
             *pdfidt_dim/(deltat*unit_time)
         enddo

       endif

end subroutine userstep
!=========================================================
subroutine get_timetable
  use ModPar
  use ModGrid
  use ModIO
  implicit none

  include "mpif.h"

  integer :: stat=0
  integer :: itsub, itcount, ifilesub, i, ierr
  real :: tsub

  character(len=80) :: filexislc1
  character(len=4) :: idoutsub
!-----------------------------------------------------------
  if(myid .eq. 0) then
    nttot=0
    open(20,file='timetable',status='old')
    do while(stat == 0)
      read(20, '(e22.14,3i8)', iostat=stat) &
        tsub, itsub, itcount, ifilesub
      if(stat == 0) then
        nttot=nttot+1
        write(6,'(e22.14,3i8)') tsub, itsub, itcount, ifilesub
      endif
    end do
    close(20)
    write(6,*) 'there are ', nttot, 'lines in the file'
    call flush(6)
  endif
  call MPI_BCAST(nttot,1,MPI_INTEGER,0, &
    MPI_COMM_WORLD,ierr)
  if(nttot .le. 1) then
    call MPI_ABORT(MPI_COMM_WORLD, 1,ierr)
  endif

  allocate(tsub_arr(nttot), itsub_arr(nttot), &
    itcount_arr(nttot), ifilesub_arr(nttot))
  open(20,file='timetable')
  do i=1,nttot
    read(20, '(e22.14,3i8)') &
      tsub, itsub, itcount, ifilesub
    tsub_arr(i)=tsub
    itsub_arr(i)=itsub
    itcount_arr(i)=itcount
    ifilesub_arr(i)=ifilesub
  enddo
  close(20)

  if(myid .eq. 0) then
    write(6,*) 'tsub_arr,itsub_arr,itcount_arr,ifilesub_arr(1)', &
      tsub_arr(1),itsub_arr(1),itcount_arr(1),ifilesub_arr(1)
    write(6,*) 'tsub_arr,itsub_arr,itcount_arr,ifilesub_arr(nttot)', &
      tsub_arr(nttot),itsub_arr(nttot),itcount_arr(nttot), &
      ifilesub_arr(nttot)
    write(6,*) 'inmax,jnmax,knmax,ntsub= ', inmax, jnmax, knmax, ntsub
    call flush(6)
  endif
  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  allocate(xislc1(ntsub,jn,kn))
  allocate(psi(jn,kn))
  allocate(bsurf(is:is,jsm2-1:jep3,ksm2-1:kep3))

  ifileread=ifilesub_arr(1)
  write(idoutsub,'(i4.4)') ifileread
  filexislc1='xislc1_23.'//idoutsub//'.dat'

  if(myid .eq. 0) then
    write(6,*) 'reading ', filexislc1
    call flush(6)
  endif

  call read23qty_mpi(ntsub,xislc1,filexislc1)

end subroutine get_timetable
!=========================================================
subroutine message_pass_nodes(idtop, idbot, nI1, nJ1, nK1, nI2, nJ2, nK2, &
  sendtop, sendbot, recvtop, recvbot, &
  isendtop, isendbot, irecvtop, irecvbot)
    implicit none
    include "mpif.h"
    
    integer, intent(in)  :: idtop, idbot, nI1, nJ1, nK1, nI2, nJ2, nK2
    integer, intent(in)  :: isendtop, isendbot, irecvtop, irecvbot
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

    call MPI_IRECV(recvtopbuf, ncount_bot, MPI_DOUBLE_PRECISION, idtop, &
         irecvtop, MPI_COMM_WORLD, ireq(2), ierr)
    call MPI_IRECV(recvbotbuf, ncount_top, MPI_DOUBLE_PRECISION, idbot, &
         irecvbot, MPI_COMM_WORLD, ireq(4), ierr)
    call MPI_ISEND(sendbotbuf, ncount_bot, MPI_DOUBLE_PRECISION, idbot, &
         isendbot, MPI_COMM_WORLD, ireq(1), ierr)
    call MPI_ISEND(sendtopbuf, ncount_top, MPI_DOUBLE_PRECISION, idtop, &
         isendtop, MPI_COMM_WORLD, ireq(3), ierr)
    call MPI_WAITALL(4,ireq,istatus_arr,ierr)

    recvtop = reshape(recvtopbuf, (/nI2, nJ2, nK2/))
    recvbot = reshape(recvbotbuf,(/nI1, nJ1, nK1/))

    deallocate(sendtopbuf, sendbotbuf, recvbotbuf, recvtopbuf)

end subroutine message_pass_nodes

end module ModUserSetup
