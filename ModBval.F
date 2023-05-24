module ModBval
  
  implicit none
  
  private
  
  public :: bvalv
  public :: bvald
  public :: bvalei
  public :: bvalemf1
  public :: bvalemf2
  public :: bvalemf3
  public :: bvaldiffemf1
  public :: bvaldiffemf2
  public :: bvaldiffemf3
  public :: bvaleta
  
contains
!-----------------------------------------------------------------------
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

subroutine bvalv
  use ModPar
  use ModGrid
  use ModField,     ONLY: v1, v2, v3
  use ModBoundary
  use ModUserSetup, ONLY: getv1ib,getv1jb,getv1kb,&
                          getv2ib,getv2jb,getv2kb,&
                          getv3ib,getv3jb,getv3kb
  implicit none

  integer idtop, idbot
  integer :: isendtop, isendbot, irecvtop, irecvbot
  real, allocatable :: sendtop(:,:,:), sendbot(:,:,:)
  real, allocatable :: recvtop(:,:,:), recvbot(:,:,:)
  real q1
  integer i, j, k
  !-------------------------------------------------------
  !------------------ i - boundary ----------------------
  idtop = myid1 + 1 + myid2*nproc1 + myid3*(nproc1*nproc2)
  idbot = myid1 - 1 + myid2*nproc1 + myid3*(nproc1*nproc2)
  isendtop=2
  isendbot=1
  irecvtop=1
  irecvbot=2
  if(myid1 .eq. 0) idbot = nproc1 - 1 + nproc1*myid2 &
                         + (nproc1*nproc2)*myid3
  if(myid1 .eq. nproc1-1) idtop = 0 + nproc1*myid2 &
                                + (nproc1*nproc2)*myid3
  ! v1
  allocate(recvtop(iep2:iep3,js:je,ks:ke),recvbot(ism2:ism1,js:je,ks:ke))
  allocate(sendtop(iem1:ie,js:je,ks:ke),sendbot(isp1:isp2,js:je,ks:ke))
  sendtop=v1(iem1:ie,js:je,ks:ke)
  sendbot=v1(isp1:isp2,js:je,ks:ke)
  call message_pass_nodes(idtop, idbot, &
    2, je-js+1, ke-ks+1, 2, je-js+1, ke-ks+1, &
    sendtop, sendbot, recvtop, recvbot, &
    isendtop, isendbot, irecvtop, irecvbot)
  v1(iep2:iep3,js:je,ks:ke) = recvtop
  v1(ism2:ism1,js:je,ks:ke) = recvbot
  deallocate(recvtop, recvbot, sendtop, sendbot)
  ! v2
  allocate(recvtop(iep1:iep2,js:jep1,ks:ke),recvbot(ism2:ism1,js:jep1,ks:ke))
  allocate(sendtop(iem1:ie,js:jep1,ks:ke), sendbot(is:isp1,js:jep1,ks:ke))
  sendtop=v2(iem1:ie,js:jep1,ks:ke)
  sendbot=v2(is:isp1,js:jep1,ks:ke)
  call message_pass_nodes(idtop, idbot, &
    2, jep1-js+1, ke-ks+1, 2, jep1-js+1, ke-ks+1, &
    sendtop, sendbot, recvtop, recvbot, &
    isendtop, isendbot, irecvtop, irecvbot)
  v2(iep1:iep2,js:jep1,ks:ke) = recvtop
  v2(ism2:ism1,js:jep1,ks:ke) = recvbot
  deallocate(recvtop, recvbot, sendtop, sendbot)
  ! v3
  allocate(recvtop(iep1:iep2,js:je,ks:kep1),recvbot(ism2:ism1,js:je,ks:kep1))
  allocate(sendtop(iem1:ie,js:je,ks:kep1),sendbot(is:isp1,js:je,ks:kep1))
  sendtop=v3(iem1:ie,js:je,ks:kep1)
  sendbot=v3(is:isp1,js:je,ks:kep1)
  call message_pass_nodes(idtop, idbot, &
    2, je-js+1, kep1-ks+1, 2, je-js+1, kep1-ks+1, &
    sendtop, sendbot, recvtop, recvbot, &
    isendtop, isendbot, irecvtop, irecvbot)
  v3(iep1:iep2,js:je,ks:kep1) = recvtop
  v3(ism2:ism1,js:je,ks:kep1) = recvbot
  deallocate(recvtop, recvbot, sendtop, sendbot)
  ! Inner i boundary.
  if(niib .eq. 3 .or. noib .eq. 3) then
    call getv1ib
    call getv2ib
    call getv3ib
  endif
  if(myid1 .eq. 0) then
    if ( niib .eq. 1) then
      do k=ks,ke
        do j=js,je
          v1(is  ,j,k) = 0.
          v1(ism1,j,k) = - v1(isp1,j,k)
          v1(ism2,j,k) = - v1(isp2,j,k)
        enddo
      enddo
      do k=ks,ke
        do j=js,jep1
          v2(ism1,j,k) = v2(is  ,j,k)
          v2(ism2,j,k) = v2(isp1,j,k)
        enddo
      enddo
      do k=ks,kep1
        do j=js,je
          v3(ism1,j,k) = v3(is  ,j,k)
          v3(ism2,j,k) = v3(isp1,j,k)
        enddo
      enddo
    endif
    if ( niib .eq. 2) then
      do k=ks,ke
        do j=js,je
          q1           = sign ( 0.5, v1(isp1,j,k) )
          v1(is  ,j,k) = v1(isp1,j,k) * ( 0.5 - q1 )
          v1(ism1,j,k) = v1(is  ,j,k)
          v1(ism2,j,k) = v1(is  ,j,k)
        enddo
      enddo
      do k=ks,ke
        do j=js,jep1
          v2(ism1,j,k) = v2(is  ,j,k)
          v2(ism2,j,k) = v2(ism1,j,k)
        enddo
      enddo
      do k=ks,kep1
        do j=js,je
          v3(ism1,j,k) = v3(is  ,j,k)
          v3(ism2,j,k) = v3(ism1,j,k)
        enddo
      enddo
    endif
    if ( niib .eq. 3) then
      do k=ks,ke
        do j=js,je
          v1(is  ,j,k) = v1iib1(j,k)
          v1(ism1,j,k) = v1iib2(j,k)
          v1(ism2,j,k) = v1iib3(j,k)
        enddo
      enddo
      do k=ks,ke
        do j=js,jep1
          v2(is,j,k) = v2iib0(j,k)
          v2(ism1,j,k) = v2iib1(j,k)
          v2(ism2,j,k) = v2iib2(j,k)
        enddo
      enddo
      do k=ks,kep1
        do j=js,je
          v3(is,j,k) = v3iib0(j,k)
          v3(ism1,j,k) = v3iib1(j,k)
          v3(ism2,j,k) = v3iib2(j,k)
        enddo
      enddo
    endif
    if ( niib .eq. 5) then
      do k=ks,ke
        do j=js,je
          v1(is  ,j,k) = 0.
          v1(ism1,j,k) = - v1(isp1,j,k)
          v1(ism2,j,k) = - v1(isp2,j,k)
        enddo
      enddo
      do k=ks,ke
        do j=js,jep1
          v2(ism1,j,k) = v2(is  ,j,k)
          v2(ism2,j,k) = v2(isp1,j,k)
        enddo
      enddo
      do k=ks,kep1
        do j=js,je
          v3(ism1,j,k) = v3(is  ,j,k)
          v3(ism2,j,k) = v3(isp1,j,k)
        enddo
      enddo
    endif
  endif
  ! Outer i boundary.
  if(myid1 .eq. nproc1-1) then
    if (noib .eq. 1) then
      do k=ks,ke
        do j=js,je
          v1(iep1,j,k) = 0.
          v1(iep2,j,k) = - v1(ie,j,k)
          v1(iep3,j,k) = - v1(iem1,j,k)
        enddo
      enddo
      do k=ks,ke
        do j=js,jep1
          v2(iep1,j,k) = v2(ie  ,j,k)
          v2(iep2,j,k) = v2(iem1,j,k)
        enddo
      enddo
      do k=ks,kep1
        do j=js,je
          v3(iep1,j,k) = v3(ie  ,j,k)
          v3(iep2,j,k) = v3(iem1,j,k)
        enddo
      enddo
    endif
    if (noib .eq. 2) then
      do k=ks,ke
        do j=js,je
          q1           = sign ( 0.5, v1(ie,j,k) )
          v1(iep1,j,k) = v1(ie,j,k) * ( 0.5 + q1 )
          v1(iep2,j,k) = v1(iep1,j,k)
          v1(iep3,j,k) = v1(iep2,j,k)
        enddo
      enddo
      do k=ks,ke
        do j=js,jep1
          v2(iep1,j,k) = v2(ie,j,k)
          v2(iep2,j,k) = v2(iep1,j,k)
        enddo
      enddo
      do k=ks,kep1
        do j=js,je
          v3(iep1,j,k) = v3(ie,j,k)
          v3(iep2,j,k) = v3(iep1,j,k)
        enddo
      enddo
    endif
    if (noib .eq. 3) then
      do k=ks,ke
        do j=js,je
          v1(iep1,j,k) = v1oib1(j,k)
          v1(iep2,j,k) = v1oib2(j,k)
          v1(iep3,j,k) = v1oib3(j,k)
        enddo
      enddo
      do k=ks,ke
        do j=js,jep1
          v2(iep1,j,k) = v2oib1(j,k)
          v2(iep2,j,k) = v2oib2(j,k)
        enddo
      enddo
      do k=ks,kep1
        do j=js,je
          v3(iep1,j,k) = v3oib1(j,k)
          v3(iep2,j,k) = v3oib2(j,k)
        enddo
      enddo
    endif
    if (noib .eq. 5) then
      do k=ks,ke
        do j=js,je
          v1(iep1,j,k) = 0.
          v1(iep2,j,k) = - v1(ie,j,k)
          v1(iep3,j,k) = - v1(iem1,j,k)
        enddo
      enddo
      do k=ks,ke
        do j=js,jep1
          v2(iep1,j,k) = v2(ie  ,j,k)
          v2(iep2,j,k) = v2(iem1,j,k)
        enddo
      enddo
      do k=ks,kep1
        do j=js,je
          v3(iep1,j,k) = v3(ie  ,j,k)
          v3(iep2,j,k) = v3(iem1,j,k)
        enddo
      enddo
    endif
  endif
  !------------------ j - boundary -----------------
  idtop = myid1 + (myid2 + 1)*nproc1 + myid3*(nproc1*nproc2)
  idbot = myid1 + (myid2 - 1)*nproc1 + myid3*(nproc1*nproc2)
  isendtop=2
  isendbot=1
  irecvtop=1
  irecvbot=2
  if(myid2 .eq. 0) idbot = myid1 + (nproc2 - 1)*nproc1 &
                         + myid3*(nproc1*nproc2)
  if(myid2 .eq. nproc2-1) idtop = myid1 + 0*nproc1 &
                                + myid3*(nproc1*nproc2)
  ! v1
  allocate(recvtop(ism2:iep3,jep1:jep2,ks:ke),recvbot(ism2:iep3,jsm2:jsm1,ks:ke))
  allocate(sendtop(ism2:iep3,jem1:je,ks:ke),sendbot(ism2:iep3,js:jsp1,ks:ke))
  sendtop=v1(ism2:iep3,jem1:je,ks:ke)
  sendbot=v1(ism2:iep3,js:jsp1,ks:ke)
  call message_pass_nodes(idtop, idbot, &
    iep3-ism2+1, 2, ke-ks+1, iep3-ism2+1, 2, ke-ks+1, &
    sendtop, sendbot, recvtop, recvbot, &
    isendtop, isendbot, irecvtop, irecvbot)
  v1(ism2:iep3,jep1:jep2,ks:ke) = recvtop
  v1(ism2:iep3,jsm2:jsm1,ks:ke) = recvbot
  deallocate(recvtop,recvbot,sendtop,sendbot)
  ! v2
  allocate(recvtop(ism2:iep2,jep2:jep3,ks:ke),recvbot(ism2:iep2,jsm2:jsm1,ks:ke))
  allocate(sendtop(ism2:iep2,jem1:je,ks:ke),sendbot(ism2:iep2,jsp1:jsp2,ks:ke))
  sendtop=v2(ism2:iep2,jem1:je,ks:ke)
  sendbot=v2(ism2:iep2,jsp1:jsp2,ks:ke)
  call message_pass_nodes(idtop, idbot, &
    iep2-ism2+1, 2, ke-ks+1, iep2-ism2+1, 2, ke-ks+1, &
    sendtop, sendbot, recvtop, recvbot, &
    isendtop, isendbot, irecvtop, irecvbot)
  v2(ism2:iep2,jep2:jep3,ks:ke) = recvtop
  v2(ism2:iep2,jsm2:jsm1,ks:ke) = recvbot
  deallocate(recvtop,recvbot,sendtop,sendbot)
  ! v3
  allocate(recvtop(ism2:iep2,jep1:jep2,ks:kep1),recvbot(ism2:iep2,jsm2:jsm1,ks:kep1))
  allocate(sendtop(ism2:iep2,jem1:je,ks:kep1),sendbot(ism2:iep2,js:jsp1,ks:kep1))
  sendtop=v3(ism2:iep2,jem1:je,ks:kep1)
  sendbot=v3(ism2:iep2,js:jsp1,ks:kep1)
  call message_pass_nodes(idtop, idbot, &
    iep2-ism2+1, 2, kep1-ks+1, iep2-ism2+1, 2, kep1-ks+1, &
    sendtop, sendbot, recvtop, recvbot, &
    isendtop, isendbot, irecvtop, irecvbot)
  v3(ism2:iep2,jep1:jep2,ks:kep1) = recvtop
  v3(ism2:iep2,jsm2:jsm1,ks:kep1) = recvbot
  deallocate(recvtop,recvbot,sendtop,sendbot)
  ! Inner j boundary
  if(nijb .eq. 3 .or. nojb .eq. 3) then
    call getv1jb
    call getv2jb
    call getv3jb
  endif
  if(myid2 .eq. 0) then
    if ( nijb .eq. 1) then
      do i=ism2,iep3
        do k=ks,ke
          v1(i,jsm1,k) = v1(i,js  ,k)
          v1(i,jsm2,k) = v1(i,jsp1,k)
        enddo
      enddo
      do i=ism2,iep2
        do k=ks,ke
          v2(i,js  ,k) = 0.
          v2(i,jsm1,k) = - v2(i,jsp1,k)
          v2(i,jsm2,k) = - v2(i,jsp2,k)
        enddo
      enddo
      do i=ism2,iep2
        do k=ks,kep1
          v3(i,jsm1,k) = v3(i,js  ,k)
          v3(i,jsm2,k) = v3(i,jsp1,k)
        enddo
      enddo
    endif
    if ( nijb .eq. 2) then
      do i=ism2,iep3
        do k=ks,ke
          v1(i,jsm1,k) = v1(i,js  ,k)
          v1(i,jsm2,k) = v1(i,jsm1,k)
        enddo
      enddo
      do i=ism2,iep2
        do k=ks,ke
          q1           = sign ( 0.5, v2(i,jsp1,k) )
          v2(i,js  ,k) = v2(i,jsp1,k) * ( 0.5 - q1 )
          v2(i,jsm1,k) = v2(i,js  ,k)
          v2(i,jsm2,k) = v2(i,js  ,k)
        enddo
      enddo
      do i=ism2,iep2
        do k=ks,kep1
          v3(i,jsm1,k) = v3(i,js  ,k)
          v3(i,jsm2,k) = v3(i,jsm1,k)
        enddo
      enddo
    endif
    if ( nijb .eq. 3) then
      do i=ism2,iep3
        do k=ks,ke
          v1(i,js,k) = v1ijb0(k,i)
          v1(i,jsm1,k) = v1ijb1(k,i)
          v1(i,jsm2,k) = v1ijb2(k,i)
        enddo
      enddo
      do i=ism2,iep2
        do k=ks,ke
          v2(i,js  ,k) = v2ijb1(k,i)
          v2(i,jsm1,k) = v2ijb2(k,i)
          v2(i,jsm2,k) = v2ijb3(k,i)
        enddo
      enddo
      do i=ism2,iep2
        do k=ks,kep1
          v3(i,js,k) = v3ijb0(k,i)
          v3(i,jsm1,k) = v3ijb1(k,i)
          v3(i,jsm2,k) = v3ijb2(k,i)
        enddo
      enddo
    endif
    if ( nijb .eq. 5) then
      do i=ism2,iep3
        do k=ks,ke
          v1(i,jsm1,k) = v1(i,js  ,k)
          v1(i,jsm2,k) = v1(i,jsp1,k)
        enddo
      enddo
      do i=ism2,iep2
        do k=ks,ke
          v2(i,js  ,k) = 0.
          v2(i,jsm1,k) = - v2(i,jsp1,k)
          v2(i,jsm2,k) = - v2(i,jsp2,k)
        enddo
      enddo
      do i=ism2,iep2
        do k=ks,kep1
          v3(i,jsm1,k) = v3(i,js  ,k)
          v3(i,jsm2,k) = v3(i,jsp1,k)
        enddo
      enddo
    endif
  endif
  ! Outer j boundary
  if(myid2 .eq. nproc2-1) then
    if (nojb .eq. 1) then
      do i=ism2,iep3
        do k=ks,ke
          v1(i,jep1,k) = v1(i,je  ,k)
          v1(i,jep2,k) = v1(i,jem1,k)
        enddo
      enddo
      do i=ism2,iep2
        do k=ks,ke
          v2(i,jep1,k) = 0.
          v2(i,jep2,k) = - v2(i,je,k)
          v2(i,jep3,k) = - v2(i,jem1,k)
        enddo
      enddo
      do i=ism2,iep2
        do k=ks,kep1
          v3(i,jep1,k) = v3(i,je  ,k)
          v3(i,jep2,k) = v3(i,jem1,k)
        enddo
      enddo
    endif
    if (nojb .eq. 2) then
      do i=ism2,iep3
        do k=ks,ke
          v1(i,jep1,k) = v1(i,je  ,k)
          v1(i,jep2,k) = v1(i,jep1,k)
        enddo
      enddo
      do i=ism2,iep2
        do k=ks,ke
          q1           = sign ( 0.5, v2(i,je,k) )
          v2(i,jep1,k) = v2(i,je  ,k) * ( 0.5 + q1 )
          v2(i,jep2,k) = v2(i,jep1,k)
          v2(i,jep3,k) = v2(i,jep2,k)
        enddo
      enddo
      do i=ism2,iep2
        do k=ks,kep1
          v3(i,jep1,k) = v3(i,je  ,k)
          v3(i,jep2,k) = v3(i,jep1,k)
        enddo
      enddo
    endif
    if (nojb .eq. 3) then
      do i=ism2,iep3
        do k=ks,ke
          v1(i,jep1,k) = v1ojb1(k,i)
          v1(i,jep2,k) = v1ojb2(k,i)
        enddo
      enddo
      do i=ism2,iep2
        do k=ks,ke
          v2(i,jep1,k) = v2ojb1(k,i)
          v2(i,jep2,k) = v2ojb2(k,i)
          v2(i,jep3,k) = v2ojb3(k,i)
        enddo
      enddo
      do i=ism2,iep2
        do k=ks,kep1
          v3(i,jep1,k) = v3ojb1(k,i)
          v3(i,jep2,k) = v3ojb2(k,i)
        enddo
      enddo
    endif
    if (nojb .eq. 5) then
      do i=ism2,iep3
        do k=ks,ke
          v1(i,jep1,k) = v1(i,je  ,k)
          v1(i,jep2,k) = v1(i,jem1,k)
        enddo
      enddo
      do i=ism2,iep2
        do k=ks,ke
          v2(i,jep1,k) = 0.
          v2(i,jep2,k) = - v2(i,je,k)
          v2(i,jep3,k) = - v2(i,jem1,k)
        enddo
      enddo
      do i=ism2,iep2
        do k=ks,kep1
          v3(i,jep1,k) = v3(i,je  ,k)
          v3(i,jep2,k) = v3(i,jem1,k)
        enddo
      enddo
    endif
  endif
  !------------------ k - boundary -----------------
  idtop = myid1 + myid2*nproc1 + (myid3+1)*(nproc1*nproc2)
  idbot = myid1 + myid2*nproc1 + (myid3-1)*(nproc1*nproc2)
  isendtop=2
  isendbot=1
  irecvtop=1
  irecvbot=2
  if(myid3==0) idbot = myid1 + myid2*nproc1 + (nproc3-1)*(nproc1*nproc2)
  if(myid3==nproc3-1) idtop = myid1 + myid2*nproc1 + 0*(nproc1*nproc2)
  ! v1
  allocate(recvtop(ism2:iep3,jsm2:jep2,kep1:kep2),recvbot(ism2:iep3,jsm2:jep2,ksm2:ksm1))
  allocate(sendtop(ism2:iep3,jsm2:jep2,kem1:ke),sendbot(ism2:iep3,jsm2:jep2,ks:ksp1))
  sendtop=v1(ism2:iep3,jsm2:jep2,kem1:ke)
  sendbot=v1(ism2:iep3,jsm2:jep2,ks:ksp1)
  call message_pass_nodes(idtop, idbot, &
    iep3-ism2+1, jep2-jsm2+1, 2, iep3-ism2+1, jep2-jsm2+1, 2, &
    sendtop, sendbot, recvtop, recvbot, &
    isendtop, isendbot, irecvtop, irecvbot)
  v1(ism2:iep3,jsm2:jep2,kep1:kep2) = recvtop
  v1(ism2:iep3,jsm2:jep2,ksm2:ksm1) = recvbot
  deallocate(recvtop,recvbot,sendtop,sendbot)
  ! v2
  allocate(recvtop(ism2:iep2,jsm2:jep3,kep1:kep2),recvbot(ism2:iep2,jsm2:jep3,ksm2:ksm1))
  allocate(sendtop(ism2:iep2,jsm2:jep3,kem1:ke),sendbot(ism2:iep2,jsm2:jep3,ks:ksp1))
  sendtop=v2(ism2:iep2,jsm2:jep3,kem1:ke)
  sendbot=v2(ism2:iep2,jsm2:jep3,ks:ksp1)
  call message_pass_nodes(idtop, idbot, &
    iep2-ism2+1, jep3-jsm2+1, 2, iep2-ism2+1, jep3-jsm2+1, 2, &
    sendtop, sendbot, recvtop, recvbot, &
    isendtop, isendbot, irecvtop, irecvbot)
  v2(ism2:iep2,jsm2:jep3,kep1:kep2) = recvtop
  v2(ism2:iep2,jsm2:jep3,ksm2:ksm1) = recvbot
  deallocate(recvtop,recvbot,sendtop,sendbot)
  ! v3
  allocate(recvtop(ism2:iep2,jsm2:jep2,kep2:kep3),recvbot(ism2:iep2,jsm2:jep2,ksm2:ksm1))
  allocate(sendtop(ism2:iep2,jsm2:jep2,kem1:ke),sendbot(ism2:iep2,jsm2:jep2,ksp1:ksp2))
  sendtop=v3(ism2:iep2,jsm2:jep2,kem1:ke)
  sendbot=v3(ism2:iep2,jsm2:jep2,ksp1:ksp2)
  call message_pass_nodes(idtop, idbot, &
    iep2-ism2+1, jep2-jsm2+1, 2, iep2-ism2+1, jep2-jsm2+1, 2, &
    sendtop, sendbot, recvtop, recvbot, &
    isendtop, isendbot, irecvtop, irecvbot)
  v3(ism2:iep2,jsm2:jep2,kep2:kep3) = recvtop
  v3(ism2:iep2,jsm2:jep2,ksm2:ksm1) = recvbot
  deallocate(recvtop,recvbot,sendtop,sendbot)
  ! Inner k boundary
  if(nikb .eq. 3 .or. nokb .eq. 3) then
    call getv1kb
    call getv2kb
    call getv3kb
  endif
  if(myid3 .eq. 0) then
    if (nikb .eq. 1) then
      do j=jsm2,jep2
        do i=ism2,iep3
          v1(i,j,ksm1) = v1(i,j,ks  )
          v1(i,j,ksm2) = v1(i,j,ksp1)
        enddo
      enddo
      do j=jsm2,jep3
        do i=ism2,iep2
          v2(i,j,ksm1) = v2(i,j,ks  )
          v2(i,j,ksm2) = v2(i,j,ksp1)
        enddo
      enddo
      do j=jsm2,jep2
        do i=ism2,iep2
          v3(i,j,ks  ) = 0.
          v3(i,j,ksm1) = - v3(i,j,ksp1)
          v3(i,j,ksm2) = - v3(i,j,ksp2)
        enddo
      enddo
    endif
    if (nikb .eq. 2) then
      do j=jsm2,jep2
        do i=ism2,iep3
          v1(i,j,ksm1) = v1(i,j,ks  )
          v1(i,j,ksm2) = v1(i,j,ksm1)
        enddo
      enddo
      do j=jsm2,jep3
        do i=ism2,iep2
          v2(i,j,ksm1) = v2(i,j,ks  )
          v2(i,j,ksm2) = v2(i,j,ksm1)
        enddo
      enddo
      do j=jsm2,jep2
        do i=ism2,iep2
          q1           = sign ( 0.5, v3(i,j,ksp1) )
          v3(i,j,ks  ) = v3(i,j,ksp1) * ( 0.5 - q1 )
          v3(i,j,ksm1) = v3(i,j,ks  )
          v3(i,j,ksm2) = v3(i,j,ks  )
        enddo
      enddo
    endif
    if (nikb .eq. 3) then
      do j=jsm2,jep2
        do i=ism2,iep3
          v1(i,j,ks) = v1ikb0(i,j)
          v1(i,j,ksm1) = v1ikb1(i,j)
          v1(i,j,ksm2) = v1ikb2(i,j)
        enddo
      enddo
      do j=jsm2,jep3
        do i=ism2,iep2
          v2(i,j,ks) = v2ikb0(i,j)
          v2(i,j,ksm1) = v2ikb1(i,j)
          v2(i,j,ksm2) = v2ikb2(i,j)
        enddo
      enddo
      do j=jsm2,jep2
        do i=ism2,iep2
          v3(i,j,ks  ) = v3ikb1(i,j)
          v3(i,j,ksm1) = v3ikb2(i,j)
          v3(i,j,ksm2) = v3ikb3(i,j)
        enddo
      enddo
    endif
    if (nikb .eq. 5) then
      do j=jsm2,jep2
        do i=ism2,iep3
          v1(i,j,ksm1) = v1(i,j,ks  )
          v1(i,j,ksm2) = v1(i,j,ksp1)
        enddo
      enddo
      do j=jsm2,jep3
        do i=ism2,iep2
          v2(i,j,ksm1) = v2(i,j,ks  )
          v2(i,j,ksm2) = v2(i,j,ksp1)
        enddo
      enddo
      do j=jsm2,jep2
        do i=ism2,iep2
          v3(i,j,ks  ) = 0.
          v3(i,j,ksm1) = - v3(i,j,ksp1)
          v3(i,j,ksm2) = - v3(i,j,ksp2)
        enddo
      enddo
    endif
  endif
  ! Outer k boundary
  if(myid3 .eq. nproc3-1) then
    if (nokb .eq. 1) then
      do j=jsm2,jep2
        do i=ism2,iep3
          v1(i,j,kep1) = v1(i,j,ke  )
          v1(i,j,kep2) = v1(i,j,kem1)
        enddo
      enddo
      do j=jsm2,jep3
        do i=ism2,iep2
          v2(i,j,kep1) = v2(i,j,ke  )
          v2(i,j,kep2) = v2(i,j,kem1)
        enddo
      enddo
      do j=jsm2,jep2
        do i=ism2,iep2
          v3(i,j,kep1) = 0.
          v3(i,j,kep2) = - v3(i,j,ke)
          v3(i,j,kep3) = - v3(i,j,kem1)
        enddo
      enddo
    endif
    if (nokb .eq. 2) then
      do j=jsm2,jep2
        do i=ism2,iep3
          v1(i,j,kep1) = v1(i,j,ke  )
          v1(i,j,kep2) = v1(i,j,kep1)
        enddo
      enddo
      do j=jsm2,jep3
        do i=ism2,iep2
          v2(i,j,kep1) = v2(i,j,ke  )
          v2(i,j,kep2) = v2(i,j,kep1)
        enddo
      enddo
      do j=jsm2,jep2
        do i=ism2,iep2
          q1           = sign ( 0.5, v3(i,j,ke) )
          v3(i,j,kep1) = v3(i,j,ke  ) * ( 0.5 + q1 )
          v3(i,j,kep2) = v3(i,j,kep1)
          v3(i,j,kep3) = v3(i,j,kep2)
        enddo
      enddo
    endif
    if (nokb .eq. 3) then
      do j=jsm2,jep2
        do i=ism2,iep3
          v1(i,j,kep1) = v1okb1(i,j)
          v1(i,j,kep2) = v1okb2(i,j)
        enddo
      enddo
      do j=jsm2,jep3
        do i=ism2,iep2
          v2(i,j,kep1) = v2okb1(i,j)
          v2(i,j,kep2) = v2okb2(i,j)
        enddo
      enddo
      do j=jsm2,jep2
        do i=ism2,iep2
          v3(i,j,kep1) = v3okb1(i,j)
          v3(i,j,kep2) = v3okb2(i,j)
          v3(i,j,kep3) = v3okb3(i,j)
        enddo
      enddo
    endif
    if (nokb .eq. 5) then
      do j=jsm2,jep2
        do i=ism2,iep3
          v1(i,j,kep1) = v1(i,j,ke  )
          v1(i,j,kep2) = v1(i,j,kem1)
        enddo
      enddo
      do j=jsm2,jep3
        do i=ism2,iep2
          v2(i,j,kep1) = v2(i,j,ke  )
          v2(i,j,kep2) = v2(i,j,kem1)
        enddo
      enddo
      do j=jsm2,jep2
        do i=ism2,iep2
          v3(i,j,kep1) = 0.
          v3(i,j,kep2) = - v3(i,j,ke)
          v3(i,j,kep3) = - v3(i,j,kem1)
        enddo
      enddo
    endif
  endif

end subroutine bvalv

subroutine bvald
  use ModPar
  use ModGrid
  use ModField,        ONLY: d
  use ModBoundary
  use ModUserSetup,    ONLY: getdib, getdjb, getdkb
  implicit none

  integer :: idtop, idbot
  integer :: isendtop, isendbot, irecvtop, irecvbot
  real, allocatable :: sendtop(:,:,:), sendbot(:,:,:)
  real, allocatable :: recvtop(:,:,:), recvbot(:,:,:)
  integer i, j, k
  !---------------------------------------------------------
  !------------------ i - boundary -----------------
  idtop = myid1 + 1 + myid2*nproc1 + myid3*(nproc1*nproc2)
  idbot = myid1 - 1 + myid2*nproc1 + myid3*(nproc1*nproc2)
  isendtop=2
  isendbot=1
  irecvtop=1
  irecvbot=2
  if(myid1 .eq. 0) idbot = nproc1 - 1 + nproc1*myid2 &
                         + (nproc1*nproc2)*myid3
  if(myid1 .eq. nproc1-1) idtop = 0 + nproc1*myid2 &
                                + (nproc1*nproc2)*myid3
  allocate(recvtop(iep1:iep2,js:je,ks:ke), recvbot(ism2:ism1,js:je,ks:ke))
  allocate(sendtop(iem1:ie,js:je,ks:ke), sendbot(is:isp1,js:je,ks:ke))
  sendtop=d(iem1:ie,js:je,ks:ke)
  sendbot=d(is:isp1,js:je,ks:ke)
  call message_pass_nodes(idtop, idbot, &
    2, je-js+1, ke-ks+1, 2, je-js+1, ke-ks+1, &
    sendtop, sendbot, recvtop, recvbot, &
    isendtop, isendbot, irecvtop, irecvbot)
  d(iep1:iep2,js:je,ks:ke) = recvtop
  d(ism2:ism1,js:je,ks:ke) = recvbot
  deallocate(recvtop, recvbot, sendtop, sendbot)
  ! Inner i boundary
  if(niib .eq. 3 .or. noib .eq. 3) then
    call getdib
  endif
  if(myid1 .eq. 0) then
    if(niib .eq. 1) then
      do k=ks,ke
        do j=js,je
          d(ism1,j,k) = d(is  ,j,k)
          d(ism2,j,k) = d(isp1,j,k)
        enddo
      enddo
    endif
    if(niib .eq. 2) then
      do k=ks,ke
        do j=js,je
          d(ism1,j,k) = d(is  ,j,k)
          d(ism2,j,k) = d(ism1,j,k)
        enddo
      enddo
    endif
    if(niib .eq. 3) then
      do k=ks,ke
        do j=js,je
          d(ism1,j,k) = diib1(j,k)
          d(ism2,j,k) = diib2(j,k)
        enddo
      enddo
    endif
    if(niib .eq. 5) then
      do k=ks,ke
        do j=js,je
          d(ism1,j,k) = d(is  ,j,k)
          d(ism2,j,k) = d(isp1,j,k)
        enddo
      enddo
    endif
  endif
  ! Outer i boundary
  if(myid1 .eq. nproc1-1) then
    if(noib .eq. 1) then
      do k=ks,ke
        do j=js,je
          d(iep1,j,k) = d(ie,j,k)
          d(iep2,j,k) = d(iem1,j,k)
        enddo
      enddo
    endif
    if(noib .eq. 2) then
      do k=ks,ke
        do j=js,je
          d(iep1,j,k) = d(ie,j,k)
          d(iep2,j,k) = d(iep1,j,k)
        enddo
      enddo
    endif
    if(noib .eq. 3) then
      do k=ks,ke
        do j=js,je
          d(iep1,j,k) = doib1(j,k)
          d(iep2,j,k) = doib2(j,k)
        enddo
      enddo
    endif
    if(noib .eq. 5) then
      do k=ks,ke
        do j=js,je
          d(iep1,j,k) = d(ie  ,j,k)
          d(iep2,j,k) = d(iem1,j,k)
        enddo
      enddo
    endif
  endif

  !------------------ j - boundary -----------------
  idtop = myid1 + (myid2 + 1)*nproc1 + myid3*(nproc1*nproc2)
  idbot = myid1 + (myid2 - 1)*nproc1 + myid3*(nproc1*nproc2)
  isendtop=2
  isendbot=1
  irecvtop=1
  irecvbot=2
  if(myid2 .eq. 0) idbot = myid1 + (nproc2 - 1)*nproc1 &
                         + myid3*(nproc1*nproc2)
  if(myid2 .eq. nproc2-1) idtop = myid1 + 0*nproc1 &
                                + myid3*(nproc1*nproc2)
  allocate(recvtop(ism2:iep2,jep1:jep2,ks:ke),recvbot(ism2:iep2,jsm2:jsm1,ks:ke))
  allocate(sendtop(ism2:iep2,jem1:je,ks:ke), sendbot(ism2:iep2,js:jsp1,ks:ke))
  sendtop=d(ism2:iep2,jem1:je,ks:ke)
  sendbot=d(ism2:iep2,js:jsp1,ks:ke)
  call message_pass_nodes(idtop, idbot, &
    iep2-ism2+1, 2, ke-ks+1, iep2-ism2+1, 2, ke-ks+1, &
    sendtop, sendbot, recvtop, recvbot, &
    isendtop, isendbot, irecvtop, irecvbot)
  d(ism2:iep2,jep1:jep2,ks:ke) = recvtop
  d(ism2:iep2,jsm2:jsm1,ks:ke) = recvbot
  deallocate(recvtop, recvbot,sendtop,sendbot)
  ! Inner j boundary
  if(nijb .eq. 3 .or. nojb .eq. 3) then
    call getdjb
  endif
  if(myid2 .eq. 0) then
    if ( nijb .eq. 1) then
      do i=ism2,iep2
        do k=ks,ke
          d(i,jsm1,k) = d(i,js  ,k)
          d(i,jsm2,k) = d(i,jsp1,k)
        enddo
      enddo
    endif
    if ( nijb .eq. 2) then
      do i=ism2,iep2
        do k=ks,ke
          d(i,jsm1,k) = d(i,js  ,k)
          d(i,jsm2,k) = d(i,jsm1,k)
        enddo
      enddo
    endif
    if ( nijb .eq. 3) then
      do i=ism2,iep2
        do k=ks,ke
          d(i,jsm1,k) = dijb1(k,i)
          d(i,jsm2,k) = dijb2(k,i)
        enddo
      enddo
    endif
    if ( nijb .eq. 5) then
      do i=ism2,iep2
        do k=ks,ke
          d(i,jsm1,k) = d(i,js  ,k)
          d(i,jsm2,k) = d(i,jsp1,k)
        enddo
      enddo
    endif
  endif
  ! Outer j boundary
  if(myid2 .eq. nproc2-1) then
    if ( nojb .eq. 1) then
      do i=ism2,iep2
        do k=ks,ke
          d(i,jep1,k) = d(i,je  ,k)
          d(i,jep2,k) = d(i,jem1,k)
        enddo
      enddo
    endif
    if ( nojb .eq. 2) then
      do i=ism2,iep2
        do k=ks,ke
          d(i,jep1,k) = d(i,je  ,k)
          d(i,jep2,k) = d(i,jep1,k)
        enddo
      enddo
    endif
    if ( nojb .eq. 3) then
      do i=ism2,iep2
        do k=ks,ke
          d(i,jep1,k) = dojb1(k,i)
          d(i,jep2,k) = dojb2(k,i)
        enddo
      enddo
    endif
    if ( nojb .eq. 5) then
      do i=ism2,iep2
        do k=ks,ke
          d(i,jep1,k) = d(i,je  ,k)
          d(i,jep2,k) = d(i,jem1,k)
        enddo
      enddo
    endif
  endif
  !------------------ k - boundary -----------------
  idtop = myid1 + myid2*nproc1 + (myid3+1)*(nproc1*nproc2)
  idbot = myid1 + myid2*nproc1 + (myid3-1)*(nproc1*nproc2)
  isendtop=2
  isendbot=1
  irecvtop=1
  irecvbot=2
  if(myid3==0) idbot = myid1 + myid2*nproc1 + (nproc3-1)*(nproc1*nproc2)
  if(myid3==nproc3-1) idtop = myid1 + myid2*nproc1 + 0*(nproc1*nproc2)
  allocate(recvtop(ism2:iep2,jsm2:jep2,kep1:kep2),recvbot(ism2:iep2,jsm2:jep2,ksm2:ksm1))
  allocate(sendtop(ism2:iep2,jsm2:jep2,kem1:ke), sendbot(ism2:iep2,jsm2:jep2,ks:ksp1))
  sendtop=d(ism2:iep2,jsm2:jep2,kem1:ke)
  sendbot=d(ism2:iep2,jsm2:jep2,ks:ksp1)
  call message_pass_nodes(idtop, idbot, &
    iep2-ism2+1, jep2-jsm2+1, 2, iep2-ism2+1, jep2-jsm2+1, 2, &
    sendtop, sendbot, recvtop, recvbot, &
    isendtop, isendbot, irecvtop, irecvbot)
  d(ism2:iep2,jsm2:jep2,kep1:kep2) = recvtop
  d(ism2:iep2,jsm2:jep2,ksm2:ksm1) = recvbot
  deallocate(recvtop, recvbot, sendtop, sendbot)
  ! Inner k boundary
  if(nikb .eq. 3 .or. nokb .eq. 3) then
    call getdkb
  endif
  if(myid3 .eq. 0) then
    if ( nikb .eq. 1) then
      do j=jsm2,jep2
        do i=ism2,iep2
          d(i,j,ksm1) = d(i,j,ks  )
          d(i,j,ksm2) = d(i,j,ksp1)
        enddo
      enddo
    endif
    if ( nikb .eq. 2) then
      do j=jsm2,jep2
        do i=ism2,iep2
          d(i,j,ksm1) = d(i,j,ks  )
          d(i,j,ksm2) = d(i,j,ksm1)
        enddo
      enddo
    endif
    if ( nikb .eq. 3) then
      do j=jsm2,jep2
        do i=ism2,iep2
          d(i,j,ksm1) = dikb1(i,j)
          d(i,j,ksm2) = dikb2(i,j)
        enddo
      enddo
    endif
    if ( nikb .eq. 5) then
      do j=jsm2,jep2
        do i=ism2,iep2
          d(i,j,ksm1) = d(i,j,ks  )
          d(i,j,ksm2) = d(i,j,ksp1)
        enddo
      enddo
    endif
  endif
  ! Outer k boundary
  if(myid3 .eq. nproc3-1) then
    if ( nokb .eq. 1) then
      do j=jsm2,jep2
        do i=ism2,iep2
          d(i,j,kep1) = d(i,j,ke  )
          d(i,j,kep2) = d(i,j,kem1)
        enddo
      enddo
    endif
    if ( nokb .eq. 2) then
      do j=jsm2,jep2
        do i=ism2,iep2
          d(i,j,kep1) = d(i,j,ke  )
          d(i,j,kep2) = d(i,j,kep1)
        enddo
      enddo
    endif
    if ( nokb .eq. 3) then
      do j=jsm2,jep2
        do i=ism2,iep2
          d(i,j,kep1) = dokb1(i,j)
          d(i,j,kep2) = dokb2(i,j)
        enddo
      enddo
    endif
    if ( nokb .eq. 5) then
      do j=jsm2,jep2
        do i=ism2,iep2
          d(i,j,kep1) = d(i,j,ke  )
          d(i,j,kep2) = d(i,j,kem1)
        enddo
      enddo
    endif
  endif

end subroutine bvald

subroutine bvalei
#ifndef ISOMHD
  use ModPar
  use ModGrid
  use ModField,        ONLY: eint
  use ModBoundary
  use ModUserSetup,    ONLY: geteintib, geteintjb, geteintkb
  implicit none

  integer :: idtop, idbot
  integer :: isendtop, isendbot, irecvtop, irecvbot
  real, allocatable :: sendtop(:,:,:), sendbot(:,:,:)
  real, allocatable :: recvtop(:,:,:), recvbot(:,:,:)
  integer i, j, k
  !---------------------------------------------------------
  !------------------ i - boundary -----------------
  idtop = myid1 + 1 + myid2*nproc1 + myid3*(nproc1*nproc2)
  idbot = myid1 - 1 + myid2*nproc1 + myid3*(nproc1*nproc2)
  isendtop=2
  isendbot=1
  irecvtop=1
  irecvbot=2
  if(myid1 .eq. 0) idbot = nproc1 - 1 + nproc1*myid2 &
                         + (nproc1*nproc2)*myid3
  if(myid1 .eq. nproc1-1) idtop = 0 + nproc1*myid2 &
                                + (nproc1*nproc2)*myid3
  allocate(recvtop(iep1:iep2,js:je,ks:ke), recvbot(ism2:ism1,js:je,ks:ke))
  allocate(sendtop(iem1:ie,js:je,ks:ke), sendbot(is:isp1,js:je,ks:ke))
  sendtop=eint(iem1:ie,js:je,ks:ke)
  sendbot=eint(is:isp1,js:je,ks:ke)
  call message_pass_nodes(idtop, idbot, &
    2, je-js+1, ke-ks+1, 2, je-js+1, ke-ks+1, &
    sendtop, sendbot, recvtop, recvbot, &
    isendtop, isendbot, irecvtop, irecvbot)
  eint(iep1:iep2,js:je,ks:ke) = recvtop
  eint(ism2:ism1,js:je,ks:ke) = recvbot
  deallocate(recvtop, recvbot, sendtop, sendbot)
  ! Inner i boundary
  if(niib .eq. 3 .or. noib .eq. 3) then
    call geteintib
  endif
  if(myid1 .eq. 0) then
    if(niib .eq. 1) then
      do k=ks,ke
        do j=js,je
          eint(ism1,j,k) = eint(is  ,j,k)
          eint(ism2,j,k) = eint(isp1,j,k)
        enddo
      enddo
    endif
    if(niib .eq. 2) then
      do k=ks,ke
        do j=js,je
          eint(ism1,j,k) = eint(is  ,j,k)
          eint(ism2,j,k) = eint(ism1,j,k)
        enddo
      enddo
    endif
    if(niib .eq. 3) then
      do k=ks,ke
        do j=js,je
          eint(ism1,j,k) = eint_iib1(j,k)
          eint(ism2,j,k) = eint_iib2(j,k)
        enddo
      enddo
    endif
    if(niib .eq. 5) then
      do k=ks,ke
        do j=js,je
          eint(ism1,j,k) = eint(is  ,j,k)
          eint(ism2,j,k) = eint(isp1,j,k)
        enddo
      enddo
    endif
  endif
  ! Outer i boundary
  if(myid1 .eq. nproc1-1) then
    if(noib .eq. 1) then
      do k=ks,ke
        do j=js,je
          eint(iep1,j,k) = eint(ie,j,k)
          eint(iep2,j,k) = eint(iem1,j,k)
        enddo
      enddo
    endif
    if(noib .eq. 2) then
      do k=ks,ke
        do j=js,je
          eint(iep1,j,k) = eint(ie,j,k)
          eint(iep2,j,k) = eint(iep1,j,k)
        enddo
      enddo
    endif
    if(noib .eq. 3) then
      do k=ks,ke
        do j=js,je
          eint(iep1,j,k) = eint_oib1(j,k)
          eint(iep2,j,k) = eint_oib2(j,k)
        enddo
      enddo
    endif
    if(noib .eq. 5) then
      do k=ks,ke
        do j=js,je
          eint(iep1,j,k) = eint(ie  ,j,k)
          eint(iep2,j,k) = eint(iem1,j,k)
        enddo
      enddo
    endif
  endif

  !------------------ j - boundary -----------------
  idtop = myid1 + (myid2 + 1)*nproc1 + myid3*(nproc1*nproc2)
  idbot = myid1 + (myid2 - 1)*nproc1 + myid3*(nproc1*nproc2)
  isendtop=2
  isendbot=1
  irecvtop=1
  irecvbot=2
  if(myid2 .eq. 0) idbot = myid1 + (nproc2 - 1)*nproc1 &
                         + myid3*(nproc1*nproc2)
  if(myid2 .eq. nproc2-1) idtop = myid1 + 0*nproc1 &
                                + myid3*(nproc1*nproc2)
  allocate(recvtop(ism2:iep2,jep1:jep2,ks:ke),recvbot(ism2:iep2,jsm2:jsm1,ks:ke))
  allocate(sendtop(ism2:iep2,jem1:je,ks:ke), sendbot(ism2:iep2,js:jsp1,ks:ke))
  sendtop=eint(ism2:iep2,jem1:je,ks:ke)
  sendbot=eint(ism2:iep2,js:jsp1,ks:ke)
  call message_pass_nodes(idtop, idbot, &
    iep2-ism2+1, 2, ke-ks+1, iep2-ism2+1, 2, ke-ks+1, &
    sendtop, sendbot, recvtop, recvbot, &
    isendtop, isendbot, irecvtop, irecvbot)
  eint(ism2:iep2,jep1:jep2,ks:ke) = recvtop
  eint(ism2:iep2,jsm2:jsm1,ks:ke) = recvbot
  deallocate(recvtop, recvbot, sendtop, sendbot)
  ! Inner j boundary
  if(nijb .eq. 3 .or. nojb .eq. 3) then
    call geteintjb
  endif
  if(myid2 .eq. 0) then
    if ( nijb .eq. 1) then
      do i=ism2,iep2
        do k=ks,ke
          eint(i,jsm1,k) = eint(i,js  ,k)
          eint(i,jsm2,k) = eint(i,jsp1,k)
        enddo
      enddo
    endif
    if ( nijb .eq. 2) then
      do i=ism2,iep2
        do k=ks,ke
          eint(i,jsm1,k) = eint(i,js  ,k)
          eint(i,jsm2,k) = eint(i,jsm1,k)
        enddo
      enddo
    endif
    if ( nijb .eq. 3) then
      do i=ism2,iep2
        do k=ks,ke
          eint(i,jsm1,k) = eint_ijb1(k,i)
          eint(i,jsm2,k) = eint_ijb2(k,i)
        enddo
      enddo
    endif
    if ( nijb .eq. 5) then
      do i=ism2,iep2
        do k=ks,ke
          eint(i,jsm1,k) = eint(i,js  ,k)
          eint(i,jsm2,k) = eint(i,jsp1,k)
        enddo
      enddo
    endif
  endif
  ! Outer j boundary
  if(myid2 .eq. nproc2-1) then
    if ( nojb .eq. 1) then
      do i=ism2,iep2
        do k=ks,ke
          eint(i,jep1,k) = eint(i,je  ,k)
          eint(i,jep2,k) = eint(i,jem1,k)
        enddo
      enddo
    endif
    if ( nojb .eq. 2) then
      do i=ism2,iep2
        do k=ks,ke
          eint(i,jep1,k) = eint(i,je  ,k)
          eint(i,jep2,k) = eint(i,jep1,k)
        enddo
      enddo
    endif
    if ( nojb .eq. 3) then
      do i=ism2,iep2
        do k=ks,ke
          eint(i,jep1,k) = eint_ojb1(k,i)
          eint(i,jep2,k) = eint_ojb2(k,i)
        enddo
      enddo
    endif
    if ( nojb .eq. 5) then
      do i=ism2,iep2
        do k=ks,ke
          eint(i,jep1,k) = eint(i,je  ,k)
          eint(i,jep2,k) = eint(i,jem1,k)
        enddo
      enddo
    endif
  endif
  !------------------ k - boundary -----------------
  idtop = myid1 + myid2*nproc1 + (myid3+1)*(nproc1*nproc2)
  idbot = myid1 + myid2*nproc1 + (myid3-1)*(nproc1*nproc2)
  isendtop=2
  isendbot=1
  irecvtop=1
  irecvbot=2
  if(myid3==0) idbot = myid1 + myid2*nproc1 + (nproc3-1)*(nproc1*nproc2)
  if(myid3==nproc3-1) idtop = myid1 + myid2*nproc1 + 0*(nproc1*nproc2)
  allocate(recvtop(ism2:iep2,jsm2:jep2,kep1:kep2),recvbot(ism2:iep2,jsm2:jep2,ksm2:ksm1))
  allocate(sendtop(ism2:iep2,jsm2:jep2,kem1:ke),sendbot(ism2:iep2,jsm2:jep2,ks:ksp1))
  sendtop=eint(ism2:iep2,jsm2:jep2,kem1:ke)
  sendbot=eint(ism2:iep2,jsm2:jep2,ks:ksp1)
  call message_pass_nodes(idtop, idbot, &
    iep2-ism2+1, jep2-jsm2+1, 2, iep2-ism2+1, jep2-jsm2+1, 2, &
    sendtop, sendbot, recvtop, recvbot, &
    isendtop, isendbot, irecvtop, irecvbot)
  eint(ism2:iep2,jsm2:jep2,kep1:kep2) = recvtop
  eint(ism2:iep2,jsm2:jep2,ksm2:ksm1) = recvbot
  deallocate(recvtop, recvbot, sendtop, sendbot)
  ! Inner k boundary
  if(nikb .eq. 3 .or. nokb .eq. 3) then
    call geteintkb
  endif
  if(myid3 .eq. 0) then
    if ( nikb .eq. 1) then
      do j=jsm2,jep2
        do i=ism2,iep2
          eint(i,j,ksm1) = eint(i,j,ks  )
          eint(i,j,ksm2) = eint(i,j,ksp1)
        enddo
      enddo
    endif
    if ( nikb .eq. 2) then
      do j=jsm2,jep2
        do i=ism2,iep2
          eint(i,j,ksm1) = eint(i,j,ks  )
          eint(i,j,ksm2) = eint(i,j,ksm1)
        enddo
      enddo
    endif
    if ( nikb .eq. 3) then
      do j=jsm2,jep2
        do i=ism2,iep2
          eint(i,j,ksm1) = eint_ikb1(i,j)
          eint(i,j,ksm2) = eint_ikb2(i,j)
        enddo
      enddo
    endif
    if ( nikb .eq. 5) then
      do j=jsm2,jep2
        do i=ism2,iep2
          eint(i,j,ksm1) = eint(i,j,ks  )
          eint(i,j,ksm2) = eint(i,j,ksp1)
        enddo
      enddo
    endif
  endif
  ! Outer k boundary
  if(myid3 .eq. nproc3-1) then
    if ( nokb .eq. 1) then
      do j=jsm2,jep2
        do i=ism2,iep2
          eint(i,j,kep1) = eint(i,j,ke  )
          eint(i,j,kep2) = eint(i,j,kem1)
        enddo
      enddo
    endif
    if ( nokb .eq. 2) then
      do j=jsm2,jep2
        do i=ism2,iep2
          eint(i,j,kep1) = eint(i,j,ke  )
          eint(i,j,kep2) = eint(i,j,kep1)
        enddo
      enddo
    endif
    if ( nokb .eq. 3) then
      do j=jsm2,jep2
        do i=ism2,iep2
          eint(i,j,kep1) = eint_okb1(i,j)
          eint(i,j,kep2) = eint_okb2(i,j)
        enddo
      enddo
    endif
    if ( nokb .eq. 5) then
      do j=jsm2,jep2
        do i=ism2,iep2
          eint(i,j,kep1) = eint(i,j,ke  )
          eint(i,j,kep2) = eint(i,j,kem1)
        enddo
      enddo
    endif
  endif

#endif
end subroutine bvalei

subroutine bvalemf1(v2b3mv3b2)
  use ModPar
  use ModGrid
  use ModBoundary
  use ModUserSetup, ONLY: getemf1ib, getemf1jb, getemf1kb
  implicit none

  real, intent(inout), dimension(in,jn,kn) :: v2b3mv3b2

  integer idtop, idbot
  integer :: isendtop, isendbot, irecvtop, irecvbot
  real, allocatable :: sendtop(:,:,:), sendbot(:,:,:)
  real, allocatable :: recvtop(:,:,:), recvbot(:,:,:)
  integer i, j, k
  !-------------------------------------------------------
  !------------------ i - boundary ----------------------
  idtop = myid1 + 1 + myid2*nproc1 + myid3*(nproc1*nproc2)
  idbot = myid1 - 1 + myid2*nproc1 + myid3*(nproc1*nproc2)
  isendtop=2
  isendbot=1
  irecvtop=1
  irecvbot=2
  if(myid1 .eq. 0) idbot = nproc1 - 1 + nproc1*myid2 &
                         + (nproc1*nproc2)*myid3
  if(myid1 .eq. nproc1-1) idtop = 0 + nproc1*myid2 &
                                + (nproc1*nproc2)*myid3
  allocate(recvtop(iep1:iep2,js:jep1,ks:kep1), &
           recvbot(ism2:ism1,js:jep1,ks:kep1))
  allocate(sendtop(iem1:ie,js:jep1,ks:kep1), &
           sendbot(is:isp1,js:jep1,ks:kep1))
  sendtop=v2b3mv3b2(iem1:ie,js:jep1,ks:kep1)
  sendbot=v2b3mv3b2(is:isp1,js:jep1,ks:kep1)
  call message_pass_nodes(idtop, idbot, &
    2, jep1-js+1, kep1-ks+1, 2, jep1-js+1, kep1-ks+1, &
    sendtop, sendbot, recvtop, recvbot, &
    isendtop, isendbot, irecvtop, irecvbot)
  v2b3mv3b2(iep1:iep2,js:jep1,ks:kep1) = recvtop
  v2b3mv3b2(ism2:ism1,js:jep1,ks:kep1) = recvbot
  deallocate(recvtop, recvbot, sendtop, sendbot)
  ! Inner i boundary.
  if(niib .eq. 3 .or. noib .eq. 3) then
    call getemf1ib(v2b3mv3b2)
  endif
  if(myid1 .eq. 0) then
    if ( niib .eq. 1) then
      do k=ks,kep1
        do j=js,jep1
          v2b3mv3b2(ism1,j,k) = v2b3mv3b2(is  ,j,k)
          v2b3mv3b2(ism2,j,k) = v2b3mv3b2(isp1,j,k)
        enddo
      enddo
    endif
    if ( niib .eq. 2) then
      do k=ks,kep1
        do j=js,jep1
          v2b3mv3b2(ism1,j,k) = v2b3mv3b2(is  ,j,k)*dx1a(is)*dx1ai(ism1)
          v2b3mv3b2(ism2,j,k) = v2b3mv3b2(is  ,j,k)*dx1a(is)*dx1ai(ism2)
        enddo
      enddo
    endif
    if ( niib .eq. 3) then
      do k=ks,kep1
        do j=js,jep1
          v2b3mv3b2(ism1,j,k) = emf1iib1(j,k)
          v2b3mv3b2(ism2,j,k) = emf1iib2(j,k)
        enddo
      enddo
    endif
    if ( niib .eq. 5) then
      do k=ks,kep1
        do j=js,jep1
          v2b3mv3b2(ism1,j,k) =-v2b3mv3b2(is  ,j,k)
          v2b3mv3b2(ism2,j,k) =-v2b3mv3b2(isp1,j,k)
        enddo
      enddo
    endif
  endif
  ! Outer i boundary.
  if(myid1 .eq. nproc1-1) then
    if ( noib .eq. 1) then
      do k=ks,kep1
        do j=js,jep1
          v2b3mv3b2(iep1,j,k) = v2b3mv3b2(ie  ,j,k)
          v2b3mv3b2(iep2,j,k) = v2b3mv3b2(iem1,j,k)
        enddo
      enddo
    endif
    if ( noib .eq. 2) then
      do k=ks,kep1
        do j=js,jep1
          v2b3mv3b2(iep1,j,k) = v2b3mv3b2(ie  ,j,k)*dx1a(ie)*dx1ai(iep1)
          v2b3mv3b2(iep2,j,k) = v2b3mv3b2(ie  ,j,k)*dx1a(ie)*dx1ai(iep2)
        enddo
      enddo
    endif
    if ( noib .eq. 3) then
      do k=ks,kep1
        do j=js,jep1
          v2b3mv3b2(iep1,j,k) = emf1oib1(j,k)
          v2b3mv3b2(iep2,j,k) = emf1oib2(j,k)
        enddo
      enddo
    endif
    if ( noib .eq. 5) then
      do k=ks,kep1
        do j=js,jep1
          v2b3mv3b2(iep1,j,k) =-v2b3mv3b2(ie  ,j,k)
          v2b3mv3b2(iep2,j,k) =-v2b3mv3b2(iem1,j,k)
        enddo
      enddo
    endif
  endif
  !------------------ j - boundary -----------------
  idtop = myid1 + (myid2 + 1)*nproc1 + myid3*(nproc1*nproc2)
  idbot = myid1 + (myid2 - 1)*nproc1 + myid3*(nproc1*nproc2)
  isendtop=2
  isendbot=1
  irecvtop=1
  irecvbot=2
  if(myid2 .eq. 0) idbot = myid1 + (nproc2 - 1)*nproc1 &
                         + myid3*(nproc1*nproc2)
  if(myid2 .eq. nproc2-1) idtop = myid1 + 0*nproc1 &
                                + myid3*(nproc1*nproc2)
  allocate(recvtop(ism2:iep2,jep2:jep3,ks:kep1), &
           recvbot(ism2:iep2,jsm2:jsm1,ks:kep1))
  allocate(sendtop(ism2:iep2,jem1:je,ks:kep1), &
           sendbot(ism2:iep2,jsp1:jsp2,ks:kep1))
  sendtop=v2b3mv3b2(ism2:iep2,jem1:je,ks:kep1)
  sendbot=v2b3mv3b2(ism2:iep2,jsp1:jsp2,ks:kep1)
  call message_pass_nodes(idtop, idbot, &
    iep2-ism2+1, 2, kep1-ks+1, iep2-ism2+1, 2, kep1-ks+1, &
    sendtop, sendbot, recvtop, recvbot, &
    isendtop, isendbot, irecvtop, irecvbot)
  v2b3mv3b2(ism2:iep2,jep2:jep3,ks:kep1) = recvtop
  v2b3mv3b2(ism2:iep2,jsm2:jsm1,ks:kep1) = recvbot
  deallocate(recvtop,recvbot,sendtop,sendbot)
  ! Inner j boundary.
  if(nijb .eq. 3 .or. nojb .eq. 3) then
    call getemf1jb(v2b3mv3b2)
  endif
  if(myid2 .eq. 0) then
    if (nijb .eq. 1) then
      do i=ism2,iep2
        do k=ks,kep1
          v2b3mv3b2(i,js  ,k) = e1jbs(i)
          v2b3mv3b2(i,jsm1,k) = 2.*e1jbs(i)-v2b3mv3b2(i,jsp1,k)
          v2b3mv3b2(i,jsm2,k) = 2.*e1jbs(i)-v2b3mv3b2(i,jsp2,k)
        enddo
      enddo
    endif
    if (nijb .eq. 2) then
      do i=ism2,iep2
        do k=ks,kep1
          v2b3mv3b2(i,jsm1,k) = v2b3mv3b2(i,js  ,k) &
            + v2b3mv3b2(i,js,k) - v2b3mv3b2(i,jsp1,k)
          v2b3mv3b2(i,jsm2,k) = v2b3mv3b2(i,jsm1  ,k) &
            + v2b3mv3b2(i,jsm1,k) - v2b3mv3b2(i,js,k)
        enddo
      enddo
    endif
    if (nijb .eq. 3) then
      do i=ism2,iep2
        do k=ks,kep1
          v2b3mv3b2(i,js  ,k) = emf1ijb1(k,i)
          v2b3mv3b2(i,jsm1,k) = emf1ijb2(k,i)
          v2b3mv3b2(i,jsm2,k) = emf1ijb3(k,i)
        enddo
      enddo
    endif
    if (nijb .eq. 5) then
      do i=ism2,iep2
        do k=ks,kep1
          v2b3mv3b2(i,js  ,k) = 0.
          v2b3mv3b2(i,jsm1,k) = v2b3mv3b2(i,jsp1,k)
          v2b3mv3b2(i,jsm2,k) = v2b3mv3b2(i,jsp2,k)
        enddo
      enddo
    endif
  endif
  ! Outer j boundary.
  if(myid2 .eq. nproc2-1) then
    if (nojb .eq. 1) then
      do i=ism2,iep2
        do k=ks,kep1
          v2b3mv3b2(i,jep1,k) = e1jbe(i)
          v2b3mv3b2(i,jep2,k) = 2.*e1jbe(i)-v2b3mv3b2(i,je  ,k)
          v2b3mv3b2(i,jep3,k) = 2.*e1jbe(i)-v2b3mv3b2(i,jem1,k)
        enddo
      enddo
    endif
    if (nojb .eq. 2) then
      do i=ism2,iep2
        do k=ks,kep1
          v2b3mv3b2(i,jep2,k) = v2b3mv3b2(i,jep1,k) &
            + v2b3mv3b2(i,jep1,k) - v2b3mv3b2(i,je,k)
          v2b3mv3b2(i,jep3,k) = v2b3mv3b2(i,jep2,k) &
            + v2b3mv3b2(i,jep2,k) - v2b3mv3b2(i,jep1,k)
        enddo
      enddo
    endif
    if (nojb .eq. 3) then
      do i=ism2,iep2
        do k=ks,kep1
          v2b3mv3b2(i,jep1,k) = emf1ojb1(k,i)
          v2b3mv3b2(i,jep2,k) = emf1ojb2(k,i)
          v2b3mv3b2(i,jep3,k) = emf1ojb3(k,i)
        enddo
      enddo
    endif
    if (nojb .eq. 5) then
      do i=ism2,iep2
        do k=ks,kep1
          v2b3mv3b2(i,jep1,k) = 0.
          v2b3mv3b2(i,jep2,k) = v2b3mv3b2(i,je  ,k)
          v2b3mv3b2(i,jep3,k) = v2b3mv3b2(i,jem1,k)
        enddo
      enddo
    endif
  endif
  !------------------ k - boundary -----------------
  idtop = myid1 + myid2*nproc1 + (myid3+1)*(nproc1*nproc2)
  idbot = myid1 + myid2*nproc1 + (myid3-1)*(nproc1*nproc2)
  isendtop=2
  isendbot=1
  irecvtop=1
  irecvbot=2
  if(myid3==0) idbot = myid1 + myid2*nproc1 + (nproc3-1)*(nproc1*nproc2)
  if(myid3==nproc3-1) idtop = myid1 + myid2*nproc1 + 0*(nproc1*nproc2)
  allocate(recvtop(ism2:iep2,jsm2:jep3,kep2:kep3), &
           recvbot(ism2:iep2,jsm2:jep3,ksm2:ksm1))
  allocate(sendtop(ism2:iep2,jsm2:jep3,kem1:ke), &
           sendbot(ism2:iep2,jsm2:jep3,ksp1:ksp2))
  sendtop=v2b3mv3b2(ism2:iep2,jsm2:jep3,kem1:ke)
  sendbot=v2b3mv3b2(ism2:iep2,jsm2:jep3,ksp1:ksp2)
  call message_pass_nodes(idtop, idbot, &
    iep2-ism2+1, jep3-jsm2+1, 2, iep2-ism2+1, jep3-jsm2+1, 2, &
    sendtop, sendbot, recvtop, recvbot, &
    isendtop, isendbot, irecvtop, irecvbot)
  v2b3mv3b2(ism2:iep2,jsm2:jep3,kep2:kep3) = recvtop
  v2b3mv3b2(ism2:iep2,jsm2:jep3,ksm2:ksm1) = recvbot
  deallocate(recvtop,recvbot,sendtop,sendbot)
  ! Inner k boundary
  if(nikb .eq. 3 .or. nokb .eq. 3) then
    call getemf1kb(v2b3mv3b2)
  endif
  if(myid3 .eq. 0) then
    if ( nikb .eq. 1) then
      do j=jsm2,jep3
        do i=ism2,iep2
          v2b3mv3b2(i,j,ks  ) = e1kbs(i)
          v2b3mv3b2(i,j,ksm1) = 2.*e1kbs(i)-v2b3mv3b2(i,j,ksp1)
          v2b3mv3b2(i,j,ksm2) = 2.*e1kbs(i)-v2b3mv3b2(i,j,ksp2)
        enddo
      enddo
    endif
    if ( nikb .eq. 2) then
      do j=jsm2,jep3
        do i=ism2,iep2
          v2b3mv3b2(i,j,ksm1) = v2b3mv3b2(i,j,ks  ) &
            + v2b3mv3b2(i,j,ks) - v2b3mv3b2(i,j,ksp1)
          v2b3mv3b2(i,j,ksm2) = v2b3mv3b2(i,j,ksm1  ) &
            + v2b3mv3b2(i,j,ksm1) - v2b3mv3b2(i,j,ks)
        enddo
      enddo
    endif
    if ( nikb .eq. 3) then
      do j=jsm2,jep3
        do i=ism2,iep2
          v2b3mv3b2(i,j,ks  ) = emf1ikb1(i,j)
          v2b3mv3b2(i,j,ksm1) = emf1ikb2(i,j)
          v2b3mv3b2(i,j,ksm2) = emf1ikb3(i,j)
        enddo
      enddo
    endif
    if ( nikb .eq. 5) then
      do j=jsm2,jep3
        do i=ism2,iep2
          v2b3mv3b2(i,j,ks  ) = 0.
          v2b3mv3b2(i,j,ksm1) = v2b3mv3b2(i,j,ksp1)
          v2b3mv3b2(i,j,ksm2) = v2b3mv3b2(i,j,ksp2)
        enddo
      enddo
    endif
  endif
  ! Outer k boundary.
  if(myid3 .eq. nproc3-1) then
    if ( nokb .eq. 1) then
      do j=jsm2,jep3
        do i=ism2,iep2
          v2b3mv3b2(i,j,kep1) = e1kbe(i)
          v2b3mv3b2(i,j,kep2) = 2.*e1kbe(i)-v2b3mv3b2(i,j,ke  )
          v2b3mv3b2(i,j,kep3) = 2.*e1kbe(i)-v2b3mv3b2(i,j,kem1)
        enddo
      enddo
    endif
    if ( nokb .eq. 2) then
      do j=jsm2,jep3
        do i=ism2,iep2
          v2b3mv3b2(i,j,kep2) = v2b3mv3b2(i,j,kep1) &
            + v2b3mv3b2(i,j,kep1) - v2b3mv3b2(i,j,ke)
          v2b3mv3b2(i,j,kep3) = v2b3mv3b2(i,j,kep2) &
            + v2b3mv3b2(i,j,kep2) - v2b3mv3b2(i,j,kep1)
        enddo
      enddo
    endif
    if ( nokb .eq. 3) then
      do j=jsm2,jep3
        do i=ism2,iep2
          v2b3mv3b2(i,j,kep1) = emf1okb1(i,j)
          v2b3mv3b2(i,j,kep2) = emf1okb2(i,j)
          v2b3mv3b2(i,j,kep3) = emf1okb3(i,j)
        enddo
      enddo
    endif
    if ( nokb .eq. 5) then
      do j=jsm2,jep3
        do i=ism2,iep2
          v2b3mv3b2(i,j,kep1) = 0.
          v2b3mv3b2(i,j,kep2) = v2b3mv3b2(i,j,ke  )
          v2b3mv3b2(i,j,kep3) = v2b3mv3b2(i,j,kem1)
        enddo
      enddo
    endif
  endif

end subroutine bvalemf1

subroutine bvalemf2(v3b1mv1b3)
  use ModPar
  use ModGrid
  use ModBoundary
  use ModUserSetup,   ONLY: getemf2ib, getemf2jb, getemf2kb
  implicit none

  real, intent(inout), dimension(in,jn,kn) :: v3b1mv1b3

  integer idtop, idbot
  integer :: isendtop, isendbot, irecvtop, irecvbot
  real, allocatable :: sendtop(:,:,:), sendbot(:,:,:)
  real, allocatable :: recvtop(:,:,:), recvbot(:,:,:)
  integer i, j, k
  !-----------------------------------------------------------------------------
  !------------------ j - boundary ----------------- 
  idtop = myid1 + (myid2 + 1)*nproc1 + myid3*(nproc1*nproc2)
  idbot = myid1 + (myid2 - 1)*nproc1 + myid3*(nproc1*nproc2)
  isendtop=2
  isendbot=1
  irecvtop=1
  irecvbot=2
  if(myid2 .eq. 0) idbot = myid1 + (nproc2 - 1)*nproc1 &
                         + myid3*(nproc1*nproc2)
  if(myid2 .eq. nproc2-1) idtop = myid1 + 0*nproc1 &
                                + myid3*(nproc1*nproc2)
  allocate(recvtop(is:iep1,jep1:jep2,ks:kep1), &
           recvbot(is:iep1,jsm2:jsm1,ks:kep1))
  allocate(sendtop(is:iep1,jem1:je, ks:kep1), &
           sendbot(is:iep1,js:jsp1, ks:kep1))
  sendtop=v3b1mv1b3(is:iep1,jem1:je, ks:kep1)
  sendbot=v3b1mv1b3(is:iep1,js:jsp1, ks:kep1)
  call message_pass_nodes(idtop, idbot, &
    iep1-is+1, 2, kep1-ks+1, iep1-is+1, 2, kep1-ks+1, &
    sendtop, sendbot, recvtop, recvbot, &
    isendtop, isendbot, irecvtop, irecvbot)
    v3b1mv1b3(is:iep1,jep1:jep2,ks:kep1) = recvtop
    v3b1mv1b3(is:iep1,jsm2:jsm1,ks:kep1) = recvbot
  deallocate(recvtop,recvbot,sendtop,sendbot)
  ! Inner j boundary.
  if(nijb .eq. 3 .or. nojb .eq. 3) then
    call getemf2jb(v3b1mv1b3)
  endif
  if(myid2 .eq. 0) then
    if (nijb .eq. 1) then
      do i=is,iep1
        do k=ks,kep1
          v3b1mv1b3(i,jsm1,k) = v3b1mv1b3(i,js  ,k)
          v3b1mv1b3(i,jsm2,k) = v3b1mv1b3(i,jsp1,k)
        enddo
      enddo
    endif
    if (nijb .eq. 2) then
      do i=is,iep1
        do k=ks,kep1
          v3b1mv1b3(i,jsm1,k) = v3b1mv1b3(i,js  ,k)*dx2a(js)*dx2ai(jsm1)
          v3b1mv1b3(i,jsm2,k) = v3b1mv1b3(i,js  ,k)*dx2a(js)*dx2ai(jsm2)
        enddo
      enddo
    endif
    if (nijb .eq. 3) then
      do i=is,iep1
        do k=ks,kep1
             v3b1mv1b3(i,jsm1,k) = emf2ijb1(k,i)
             v3b1mv1b3(i,jsm2,k) = emf2ijb2(k,i)
        enddo
      enddo
    endif
    if (nijb .eq. 5) then
      do i=is,iep1
        do k=ks,kep1
          v3b1mv1b3(i,jsm1,k) =-v3b1mv1b3(i,js  ,k)
          v3b1mv1b3(i,jsm2,k) =-v3b1mv1b3(i,jsp1,k)
        enddo
      enddo
    endif
  endif
  ! Outer j boundary.
  if(myid2 .eq. nproc2-1) then
    if (nojb .eq. 1) then
      do i=is,iep1
        do k=ks,kep1
          v3b1mv1b3(i,jep1,k) = v3b1mv1b3(i,je  ,k)
          v3b1mv1b3(i,jep2,k) = v3b1mv1b3(i,jem1,k)
        enddo
      enddo
    endif
    if (nojb .eq. 2) then
      do i=is,iep1
        do k=ks,kep1
          v3b1mv1b3(i,jep1,k) = v3b1mv1b3(i,je  ,k)*dx2a(je)*dx2ai(jep1)
          v3b1mv1b3(i,jep2,k) = v3b1mv1b3(i,je  ,k)*dx2a(je)*dx2ai(jep2)
        enddo
      enddo
    endif
    if (nojb .eq. 3) then
      do i=is,iep1
        do k=ks,kep1
          v3b1mv1b3(i,jep1,k) = emf2ojb1(k,i)
          v3b1mv1b3(i,jep2,k) = emf2ojb2(k,i)
        enddo
      enddo
    endif
    if (nojb .eq. 5) then
      do i=is,iep1
        do k=ks,kep1
          v3b1mv1b3(i,jep1,k) =-v3b1mv1b3(i,je  ,k)
          v3b1mv1b3(i,jep2,k) =-v3b1mv1b3(i,jem1,k)
        enddo
      enddo
    endif
  endif
  !------------------ k - boundary -----------------
  idtop = myid1 + myid2*nproc1 + (myid3+1)*(nproc1*nproc2)
  idbot = myid1 + myid2*nproc1 + (myid3-1)*(nproc1*nproc2)
  isendtop=2
  isendbot=1
  irecvtop=1
  irecvbot=2
  if(myid3==0) idbot = myid1 + myid2*nproc1 + (nproc3-1)*(nproc1*nproc2)
  if(myid3==nproc3-1) idtop = myid1 + myid2*nproc1 + 0*(nproc1*nproc2)
  allocate(recvtop(is:iep1,jsm2:jep2,kep2:kep3), &
           recvbot(is:iep1,jsm2:jep2,ksm2:ksm1))
  allocate(sendtop(is:iep1,jsm2:jep2,kem1:ke), &
           sendbot(is:iep1,jsm2:jep2,ksp1:ksp2))
  sendtop=v3b1mv1b3(is:iep1,jsm2:jep2,kem1:ke)
  sendbot=v3b1mv1b3(is:iep1,jsm2:jep2,ksp1:ksp2)
  call message_pass_nodes(idtop, idbot, &
    iep1-is+1, jep2-jsm2+1, 2, iep1-is+1, jep2-jsm2+1, 2, &
    sendtop, sendbot, recvtop, recvbot, &
    isendtop, isendbot, irecvtop, irecvbot)
  v3b1mv1b3(is:iep1,jsm2:jep2,kep2:kep3) = recvtop
  v3b1mv1b3(is:iep1,jsm2:jep2,ksm2:ksm1) = recvbot
  deallocate(recvtop, recvbot, sendtop, sendbot)
  ! Inner k boundary
  if(nikb .eq. 3 .or. nokb .eq. 3) then
    call getemf2kb(v3b1mv1b3)
  endif
  if(myid3 .eq. 0) then
    if (nikb .eq. 1) then
      do j=jsm2,jep2
        do i=is,iep1
          v3b1mv1b3(i,j,ks  ) = e2kbs(j)*g2ai(i)
          v3b1mv1b3(i,j,ksm1) = 2.*e2kbs(j)*g2ai(i)-v3b1mv1b3(i,j,ksp1)
          v3b1mv1b3(i,j,ksm2) = 2.*e2kbs(j)*g2ai(i)-v3b1mv1b3(i,j,ksp2)
        enddo
      enddo
    endif
    if (nikb .eq. 2) then
      do j=jsm2,jep2
        do i=is,iep1
          v3b1mv1b3(i,j,ksm1) = v3b1mv1b3(i,j,ks  ) &
            + v3b1mv1b3(i,j,ks) - v3b1mv1b3(i,j,ksp1)
          v3b1mv1b3(i,j,ksm2) = v3b1mv1b3(i,j,ksm1  ) &
            + v3b1mv1b3(i,j,ksm1) - v3b1mv1b3(i,j,ks)
        enddo
      enddo
    endif
    if (nikb .eq. 3) then
      do j=jsm2,jep2
        do i=is,iep1
          v3b1mv1b3(i,j,ks  ) = emf2ikb1(i,j)
          v3b1mv1b3(i,j,ksm1) = emf2ikb2(i,j)
          v3b1mv1b3(i,j,ksm2) = emf2ikb3(i,j)
        enddo
      enddo
    endif
    if (nikb .eq. 5) then
      do j=jsm2,jep2
        do i=is,iep1
          v3b1mv1b3(i,j,ks  ) = 0.
          v3b1mv1b3(i,j,ksm1) = v3b1mv1b3(i,j,ksp1)
          v3b1mv1b3(i,j,ksm2) = v3b1mv1b3(i,j,ksp2)
        enddo
      enddo
    endif
  endif
  ! Outer k boundary
  if(myid3 .eq. nproc3-1) then
    if (nokb .eq. 1) then
      do j=jsm2,jep2
        do i=is,iep1
          v3b1mv1b3(i,j,kep1) = e2kbe(j)*g2ai(i)
          v3b1mv1b3(i,j,kep2) = 2.*e2kbe(j)*g2ai(i)-v3b1mv1b3(i,j,ke  )
          v3b1mv1b3(i,j,kep3) = 2.*e2kbe(j)*g2ai(i)-v3b1mv1b3(i,j,kem1)
        enddo
      enddo
    endif
    if (nokb .eq. 2) then
      do j=jsm2,jep2
        do i=is,iep1
          v3b1mv1b3(i,j,kep2) = v3b1mv1b3(i,j,kep1) &
            + v3b1mv1b3(i,j,kep1) - v3b1mv1b3(i,j,ke)
          v3b1mv1b3(i,j,kep3) = v3b1mv1b3(i,j,kep2) &
            + v3b1mv1b3(i,j,kep2) - v3b1mv1b3(i,j,kep1)
        enddo
      enddo
    endif
    if (nokb .eq. 3) then
      do j=jsm2,jep2
        do i=is,iep1
          v3b1mv1b3(i,j,kep1) = emf2okb1(i,j)
          v3b1mv1b3(i,j,kep2) = emf2okb2(i,j)
          v3b1mv1b3(i,j,kep3) = emf2okb3(i,j)
        enddo
      enddo
    endif
    if (nokb .eq. 5) then
      do j=jsm2,jep2
        do i=is,iep1
          v3b1mv1b3(i,j,kep1) = 0.
          v3b1mv1b3(i,j,kep2) = v3b1mv1b3(i,j,ke  )
          v3b1mv1b3(i,j,kep3) = v3b1mv1b3(i,j,kem1)
        enddo
      enddo
    endif
  endif
  !-------------------------------------------------------
  !------------------ i - boundary ----------------------
  idtop = myid1 + 1 + myid2*nproc1 + myid3*(nproc1*nproc2)
  idbot = myid1 - 1 + myid2*nproc1 + myid3*(nproc1*nproc2)
  isendtop=2
  isendbot=1
  irecvtop=1
  irecvbot=2
  if(myid1 .eq. 0) idbot = nproc1 - 1 + nproc1*myid2 &
                         + (nproc1*nproc2)*myid3
  if(myid1 .eq. nproc1-1) idtop = 0 + nproc1*myid2 &
                                + (nproc1*nproc2)*myid3
  allocate(recvtop(iep2:iep3,jsm2:jep2,ksm2:kep3), &
           recvbot(ism2:ism1,jsm2:jep2,ksm2:kep3))
  allocate(sendtop(iem1:ie,jsm2:jep2,ksm2:kep3), &
           sendbot(isp1:isp2,jsm2:jep2,ksm2:kep3))
  sendtop=v3b1mv1b3(iem1:ie,jsm2:jep2,ksm2:kep3)
  sendbot=v3b1mv1b3(isp1:isp2,jsm2:jep2,ksm2:kep3)
  call message_pass_nodes(idtop, idbot, &
    2, jep2-jsm2+1, kep3-ksm2+1, 2, jep2-jsm2+1, kep3-ksm2+1, &
    sendtop, sendbot, recvtop, recvbot, &
    isendtop, isendbot, irecvtop, irecvbot)
  v3b1mv1b3(iep2:iep3,jsm2:jep2,ksm2:kep3) = recvtop
  v3b1mv1b3(ism2:ism1,jsm2:jep2,ksm2:kep3) = recvbot
  deallocate(recvtop, recvbot, sendtop, sendbot)
  ! Inner i boundary.
  if(niib .eq. 3 .or. noib .eq. 3) then
    call getemf2ib(v3b1mv1b3)
  endif
  if(myid1 .eq. 0) then
    if (niib .eq. 1) then
      do k=ksm2,kep3
        do j=jsm2,jep2
             v3b1mv1b3(is  ,j,k) = e2ibs(j)
             v3b1mv1b3(ism1,j,k) = 2.*e2ibs(j)-v3b1mv1b3(isp1,j,k)
             v3b1mv1b3(ism2,j,k) = 2.*e2ibs(j)-v3b1mv1b3(isp2,j,k)
        enddo
      enddo
    endif
    if (niib .eq. 2) then
      do k=ksm2,kep3
        do j=jsm2,jep2
          v3b1mv1b3(ism1,j,k) = v3b1mv1b3(is  ,j,k)*g2a(is)*g2ai(ism1) &
            + v3b1mv1b3(is  ,j,k)*g2a(is)*g2ai(ism1) &
            - v3b1mv1b3(isp1,j,k)*g2a(isp1)*g2ai(ism1)
          v3b1mv1b3(ism2,j,k) = v3b1mv1b3(ism1  ,j,k)*g2a(ism1)*g2ai(ism2) &
            + v3b1mv1b3(ism1  ,j,k)*g2a(ism1)*g2ai(ism2) &
            - v3b1mv1b3(is,j,k)*g2a(is)*g2ai(ism2)
        enddo
      enddo
    endif
    if (niib .eq. 3) then
      do k=ksm2,kep3
        do j=jsm2,jep2
          v3b1mv1b3(is  ,j,k) = emf2iib1(j,k)
          v3b1mv1b3(ism1,j,k) = emf2iib2(j,k)
          v3b1mv1b3(ism2,j,k) = emf2iib3(j,k)
        enddo
      enddo
    endif
    if (niib .eq. 5) then
      do k=ksm2,kep3
        do j=jsm2,jep2
          v3b1mv1b3(is  ,j,k) = 0.
          v3b1mv1b3(ism1,j,k) = v3b1mv1b3(isp1,j,k)
          v3b1mv1b3(ism2,j,k) = v3b1mv1b3(isp2,j,k)
        enddo
      enddo
    endif

  endif
  ! Outer i boundary.
  if(myid1 .eq. nproc1-1) then
    if (noib .eq. 1) then
      do k=ksm2,kep3
        do j=jsm2,jep2
          v3b1mv1b3(iep1,j,k) = e2ibe(j)
          v3b1mv1b3(iep2,j,k) = 2.*e2ibe(j)-v3b1mv1b3(ie  ,j,k)
          v3b1mv1b3(iep3,j,k) = 2.*e2ibe(j)-v3b1mv1b3(iem1,j,k)
        enddo
      enddo
    endif
    if (noib .eq. 2) then
      do k=ksm2,kep3
        do j=jsm2,jep2
          v3b1mv1b3(iep2,j,k) = v3b1mv1b3(iep1,j,k)*g2a(iep1)*g2ai(iep2) &
            + v3b1mv1b3(iep1,j,k)*g2a(iep1)*g2ai(iep2) &
            - v3b1mv1b3(ie,j,k)*g2a(ie)*g2ai(iep2)
          v3b1mv1b3(iep3,j,k) = v3b1mv1b3(iep2,j,k)*g2a(iep2)*g2ai(iep3) &
            + v3b1mv1b3(iep2,j,k)*g2a(iep2)*g2ai(iep3) &
            - v3b1mv1b3(iep1,j,k)*g2a(iep1)*g2ai(iep3)
        enddo
      enddo
    endif
    if (noib .eq. 3) then
      do k=ksm2,kep3
        do j=jsm2,jep2
          v3b1mv1b3(iep1,j,k) = emf2oib1(j,k)
          v3b1mv1b3(iep2,j,k) = emf2oib2(j,k)
          v3b1mv1b3(iep3,j,k) = emf2oib3(j,k)
        enddo
      enddo
    endif
    if (noib .eq. 5) then
      do k=ksm2,kep3
        do j=jsm2,jep2
          v3b1mv1b3(iep1,j,k) = 0.
          v3b1mv1b3(iep2,j,k) = v3b1mv1b3(ie  ,j,k)
          v3b1mv1b3(iep3,j,k) = v3b1mv1b3(iem1,j,k)
        enddo
      enddo
    endif
  endif

end subroutine bvalemf2

subroutine bvalemf3(v1b2mv2b1)
  use ModPar
  use ModGrid
  use ModBoundary
  use ModUserSetup,   ONLY: getemf3ib, getemf3jb, getemf3kb
  implicit none

  real, intent(inout), dimension(in,jn,kn) :: v1b2mv2b1

  integer idtop, idbot
  integer :: isendtop, isendbot, irecvtop, irecvbot
  real, allocatable :: sendtop(:,:,:), sendbot(:,:,:)
  real, allocatable :: recvtop(:,:,:), recvbot(:,:,:)
  integer i, j, k
  !-----------------------------------------------------------------------------
  !------------------ k - boundary -----------------
  idtop = myid1 + myid2*nproc1 + (myid3+1)*(nproc1*nproc2)
  idbot = myid1 + myid2*nproc1 + (myid3-1)*(nproc1*nproc2)
  isendtop=2
  isendbot=1
  irecvtop=1
  irecvbot=2
  if(myid3==0) idbot = myid1 + myid2*nproc1 + (nproc3-1)*(nproc1*nproc2)
  if(myid3==nproc3-1) idtop = myid1 + myid2*nproc1 + 0*(nproc1*nproc2)
  allocate(recvtop(is:iep1,js:jep1,kep1:kep2), &
           recvbot(is:iep1,js:jep1,ksm2:ksm1))
  allocate(sendtop(is:iep1,js:jep1,kem1:ke), &
           sendbot(is:iep1,js:jep1,ks:ksp1))
  sendtop=v1b2mv2b1(is:iep1,js:jep1,kem1:ke)
  sendbot=v1b2mv2b1(is:iep1,js:jep1,ks:ksp1)
  call message_pass_nodes(idtop, idbot, &
    iep1-is+1, jep1-js+1, 2, iep1-is+1, jep1-js+1, 2, &
    sendtop, sendbot, recvtop, recvbot, &
    isendtop, isendbot, irecvtop, irecvbot)
  v1b2mv2b1(is:iep1,js:jep1,kep1:kep2) = recvtop
  v1b2mv2b1(is:iep1,js:jep1,ksm2:ksm1) = recvbot
  deallocate(recvtop, recvbot, sendtop, sendbot)
  ! Inner k boundary
  if(nikb .eq. 3 .or. nokb .eq. 3) then
    call getemf3kb( v1b2mv2b1 )
  endif
  if(myid3 .eq. 0) then
    if (nikb .eq. 1) then
      do j=js,jep1
        do i=is,iep1
          v1b2mv2b1(i,j,ksm1) = v1b2mv2b1(i,j,ks  )
          v1b2mv2b1(i,j,ksm2) = v1b2mv2b1(i,j,ksp1)
        enddo
      enddo
    endif
    if (nikb .eq. 2) then
      do j=js,jep1
        do i=is,iep1
          v1b2mv2b1(i,j,ksm1) = v1b2mv2b1(i,j,ks  )*dx3a(ks)*dx3ai(ksm1)
          v1b2mv2b1(i,j,ksm2) = v1b2mv2b1(i,j,ks  )*dx3a(ks)*dx3ai(ksm2)
        enddo
      enddo
    endif
    if (nikb .eq. 3) then
      do j=js,jep1
        do i=is,iep1
          v1b2mv2b1(i,j,ksm1) = emf3ikb1(i,j)
          v1b2mv2b1(i,j,ksm2) = emf3ikb2(i,j)
        enddo
      enddo
    endif
    if (nikb .eq. 5) then
      do j=js,jep1
        do i=is,iep1
          v1b2mv2b1(i,j,ksm1) =-v1b2mv2b1(i,j,ks  )
          v1b2mv2b1(i,j,ksm2) =-v1b2mv2b1(i,j,ksp1)
        enddo
      enddo
    endif
  endif
  ! Outer k boundary
  if(myid3 .eq. nproc3-1) then
    if (nokb .eq. 1) then
      do j=js,jep1
        do i=is,iep1
          v1b2mv2b1(i,j,kep1) = v1b2mv2b1(i,j,ke  )
          v1b2mv2b1(i,j,kep2) = v1b2mv2b1(i,j,kem1)
        enddo
      enddo
    endif
    if (nokb .eq. 2) then
      do j=js,jep1
        do i=is,iep1
          v1b2mv2b1(i,j,kep1) = v1b2mv2b1(i,j,ke  )*dx3a(ke)*dx3ai(kep1)
          v1b2mv2b1(i,j,kep2) = v1b2mv2b1(i,j,ke  )*dx3a(ke)*dx3ai(kep2)
        enddo
      enddo
    endif
    if (nokb .eq. 3) then
      do j=js,jep1
        do i=is,iep1
          v1b2mv2b1(i,j,kep1) = emf3okb1(i,j)
          v1b2mv2b1(i,j,kep2) = emf3okb2(i,j)
        enddo
      enddo
    endif
    if (nokb .eq. 5) then
      do j=js,jep1
        do i=is,iep1
          v1b2mv2b1(i,j,kep1) =-v1b2mv2b1(i,j,ke  )
          v1b2mv2b1(i,j,kep2) =-v1b2mv2b1(i,j,kem1)
        enddo
      enddo
    endif
  endif
  !------------------ i - boundary -----------------
  idtop = myid1 + 1 + myid2*nproc1 + myid3*(nproc1*nproc2)
  idbot = myid1 - 1 + myid2*nproc1 + myid3*(nproc1*nproc2)
  isendtop=2
  isendbot=1
  irecvtop=1
  irecvbot=2
  if(myid1 .eq. 0) idbot = nproc1 - 1 + nproc1*myid2 &
                         + (nproc1*nproc2)*myid3
  if(myid1 .eq. nproc1-1) idtop = 0 + nproc1*myid2 &
                                + (nproc1*nproc2)*myid3
  allocate(recvtop(iep2:iep3,js:jep1,ksm2:kep2), &
           recvbot(ism2:ism1,js:jep1,ksm2:kep2))
  allocate(sendtop(iem1:ie,js:jep1,ksm2:kep2), &
           sendbot(isp1:isp2,js:jep1,ksm2:kep2))
  sendtop=v1b2mv2b1(iem1:ie,js:jep1,ksm2:kep2)
  sendbot=v1b2mv2b1(isp1:isp2,js:jep1,ksm2:kep2)
  call message_pass_nodes(idtop, idbot, &
    2, jep1-js+1, kep2-ksm2+1, 2, jep1-js+1, kep2-ksm2+1, &
    sendtop, sendbot, recvtop, recvbot, &
    isendtop, isendbot, irecvtop, irecvbot)
  v1b2mv2b1(iep2:iep3,js:jep1,ksm2:kep2) = recvtop
  v1b2mv2b1(ism2:ism1,js:jep1,ksm2:kep2) = recvbot
  deallocate(recvtop, recvbot, sendtop, sendbot)
  ! Inner i boundary
  if(niib .eq. 3 .or. noib .eq. 3) then
    call getemf3ib( v1b2mv2b1 )
  endif
  if(myid1 .eq. 0) then
    if (niib .eq. 1) then
      do k=ksm2,kep2
        do j=js,jep1
          v1b2mv2b1(is  ,j,k) = e3ibs(k)*g32ai(j)
          v1b2mv2b1(ism1,j,k) = 2.*e3ibs(k)*g32ai(j)-v1b2mv2b1(isp1,j,k)
          v1b2mv2b1(ism2,j,k) = 2.*e3ibs(k)*g32ai(j)-v1b2mv2b1(isp2,j,k)
        enddo
      enddo
    endif
    if (niib .eq. 2) then
      do k=ksm2,kep2
        do j=js,jep1
          v1b2mv2b1(ism1,j,k) = v1b2mv2b1(is  ,j,k)*g31a(is)*g31ai(ism1) &
            + v1b2mv2b1(is  ,j,k)*g31a(is)*g31ai(ism1) &
            - v1b2mv2b1(isp1,j,k)*g31a(isp1)*g31ai(ism1)
          v1b2mv2b1(ism2,j,k) = v1b2mv2b1(ism1  ,j,k)*g31a(ism1)*g31ai(ism2) &
            + v1b2mv2b1(ism1  ,j,k)*g31a(ism1)*g31ai(ism2) &
            - v1b2mv2b1(is,j,k)*g31a(is)*g31ai(ism2)
        enddo
      enddo
    endif
    if (niib .eq. 3) then
      do k=ksm2,kep2
        do j=js,jep1
          v1b2mv2b1(is  ,j,k) = emf3iib1(j,k)
          v1b2mv2b1(ism1,j,k) = emf3iib2(j,k)
          v1b2mv2b1(ism2,j,k) = emf3iib3(j,k)
        enddo
      enddo
    endif
    if (niib .eq. 5) then
      do k=ksm2,kep2
        do j=js,jep1
          v1b2mv2b1(is  ,j,k) = 0.
          v1b2mv2b1(ism1,j,k) = v1b2mv2b1(isp1,j,k)
          v1b2mv2b1(ism2,j,k) = v1b2mv2b1(isp2,j,k)
        enddo
      enddo
    endif
  endif
  ! Outer i boundary
  if(myid1 .eq. nproc1-1) then
    if (noib .eq. 1) then
      do k=ksm2,kep2
        do j=js,jep1
          v1b2mv2b1(iep1,j,k) = e3ibe(k)*g32ai(j)
          v1b2mv2b1(iep2,j,k) = 2.*e3ibe(k)*g32ai(j)-v1b2mv2b1(ie  ,j,k)
          v1b2mv2b1(iep3,j,k) = 2.*e3ibe(k)*g32ai(j)-v1b2mv2b1(iem1,j,k)
        enddo
      enddo
    endif
    if (noib .eq. 2) then
      do k=ksm2,kep2
        do j=js,jep1
             v1b2mv2b1(iep2,j,k) = v1b2mv2b1(iep1,j,k)*g31a(iep1)*g31ai(iep2) &
               + v1b2mv2b1(iep1,j,k)*g31a(iep1)*g31ai(iep2) &
               - v1b2mv2b1(ie,j,k)*g31a(ie)*g31ai(iep2)
             v1b2mv2b1(iep3,j,k) = v1b2mv2b1(iep2,j,k)*g31a(iep2)*g31ai(iep3) &
               + v1b2mv2b1(iep2,j,k)*g31a(iep2)*g31ai(iep3) &
               - v1b2mv2b1(iep1,j,k)*g31a(iep1)*g31ai(iep3)
        enddo
      enddo
    endif
    if (noib .eq. 3) then
      do k=ksm2,kep2
        do j=js,jep1
          v1b2mv2b1(iep1,j,k) = emf3oib1(j,k)
          v1b2mv2b1(iep2,j,k) = emf3oib2(j,k)
          v1b2mv2b1(iep3,j,k) = emf3oib3(j,k)
        enddo
      enddo
    endif
    if (noib .eq. 5) then
      do k=ksm2,kep2
        do j=js,jep1
          v1b2mv2b1(iep1,j,k) = 0.
          v1b2mv2b1(iep2,j,k) = v1b2mv2b1(ie  ,j,k)
          v1b2mv2b1(iep3,j,k) = v1b2mv2b1(iem1,j,k)
        enddo
      enddo
    endif
  endif
  !------------------ j - boundary -----------------
  idtop = myid1 + (myid2 + 1)*nproc1 + myid3*(nproc1*nproc2)
  idbot = myid1 + (myid2 - 1)*nproc1 + myid3*(nproc1*nproc2)
  isendtop=2
  isendbot=1
  irecvtop=1
  irecvbot=2
  if(myid2 .eq. 0) idbot = myid1 + (nproc2 - 1)*nproc1 &
                         + myid3*(nproc1*nproc2)
  if(myid2 .eq. nproc2-1) idtop = myid1 + 0*nproc1 &
                                + myid3*(nproc1*nproc2)
  allocate(recvtop(ism2:iep3,jep2:jep3,ksm2:kep2), &
           recvbot(ism2:iep3,jsm2:jsm1,ksm2:kep2))
  allocate(sendtop(ism2:iep3,jem1:je,ksm2:kep2), &
           sendbot(ism2:iep3,jsp1:jsp2,ksm2:kep2))
  sendtop=v1b2mv2b1(ism2:iep3,jem1:je,ksm2:kep2)
  sendbot=v1b2mv2b1(ism2:iep3,jsp1:jsp2,ksm2:kep2)
  call message_pass_nodes(idtop, idbot, &
    iep3-ism2+1, 2, kep2-ksm2+1, iep3-ism2+1, 2, kep2-ksm2+1, &
    sendtop, sendbot, recvtop, recvbot, &
    isendtop, isendbot, irecvtop, irecvbot)
  v1b2mv2b1(ism2:iep3,jep2:jep3,ksm2:kep2) = recvtop
  v1b2mv2b1(ism2:iep3,jsm2:jsm1,ksm2:kep2) = recvbot
  deallocate(recvtop,recvbot,sendtop,sendbot)
  ! Inner j boundary.   
  if(nijb .eq. 3 .or. nojb .eq. 3) then
    call getemf3jb( v1b2mv2b1 )
  endif
  if(myid2 .eq. 0) then
    if ( nijb .eq. 1) then
      do i=ism2,iep3
        do k=ksm2,kep2
          v1b2mv2b1(i,js  ,k) = e3jbs(k)*g31ai(i)
          v1b2mv2b1(i,jsm1,k) = 2.*e3jbs(k)*g31ai(i)-v1b2mv2b1(i,jsp1,k)
          v1b2mv2b1(i,jsm2,k) = 2.*e3jbs(k)*g31ai(i)-v1b2mv2b1(i,jsp2,k)
        enddo
      enddo
    endif
    if ( nijb .eq. 2) then
      do i=ism2,iep3
        do k=ksm2,kep2
             v1b2mv2b1(i,jsm1,k) = v1b2mv2b1(i,js  ,k)*g32a(js)*g32ai(jsm1) &
               + v1b2mv2b1(i,js,k)*g32a(js)*g32ai(jsm1) &
               - v1b2mv2b1(i,jsp1,k)*g32a(jsp1)*g32ai(jsm1)
             v1b2mv2b1(i,jsm2,k) = v1b2mv2b1(i,jsm1  ,k)*g32a(jsm1)*g32ai(jsm2) &
               + v1b2mv2b1(i,jsm1,k)*g32a(jsm1)*g32ai(jsm2) &
               - v1b2mv2b1(i,js,k)*g32a(js)*g32ai(jsm2)
        enddo
      enddo
    endif
    if ( nijb .eq. 3) then
      do i=ism2,iep3
        do k=ksm2,kep2
          v1b2mv2b1(i,js  ,k) = emf3ijb1(k,i)
          v1b2mv2b1(i,jsm1,k) = emf3ijb2(k,i)
          v1b2mv2b1(i,jsm2,k) = emf3ijb3(k,i)
        enddo
      enddo
    endif
    if ( nijb .eq. 5) then
      do i=ism2,iep3
        do k=ksm2,kep2
          v1b2mv2b1(i,js  ,k) = 0.
          v1b2mv2b1(i,jsm1,k) = v1b2mv2b1(i,jsp1,k)
          v1b2mv2b1(i,jsm2,k) = v1b2mv2b1(i,jsp2,k)
        enddo
      enddo
    endif
  endif
  ! Outer j boundary.                      
  if(myid2 .eq. nproc2-1) then
    if (nojb .eq. 1) then
      do i=ism2,iep3
        do k=ksm2,kep2
             v1b2mv2b1(i,jep1,k) = e3jbe(k)*g31ai(i)
             v1b2mv2b1(i,jep2,k) = 2.*e3jbe(k)*g31ai(i)-v1b2mv2b1(i,je  ,k)
             v1b2mv2b1(i,jep3,k) = 2.*e3jbe(k)*g31ai(i)-v1b2mv2b1(i,jem1,k)
        enddo
      enddo
    endif
    if (nojb .eq. 2) then
      do i=ism2,iep3
        do k=ksm2,kep2
             v1b2mv2b1(i,jep2,k) = v1b2mv2b1(i,jep1,k)*g32a(jep1)*g32ai(jep2) &
               + v1b2mv2b1(i,jep1,k)*g32a(jep1)*g32ai(jep2) &
               - v1b2mv2b1(i,je,k)*g32a(je)*g32ai(jep2)
             v1b2mv2b1(i,jep3,k) = v1b2mv2b1(i,jep2,k)*g32a(jep2)*g32ai(jep3) &
               + v1b2mv2b1(i,jep2,k)*g32a(jep2)*g32ai(jep3) &
               - v1b2mv2b1(i,jep1,k)*g32a(jep1)*g32ai(jep3)
        enddo
      enddo
    endif
    if (nojb .eq. 3) then
      do i=ism2,iep3
        do k=ksm2,kep2
          v1b2mv2b1(i,jep1,k) = emf3ojb1(k,i)
          v1b2mv2b1(i,jep2,k) = emf3ojb2(k,i)
          v1b2mv2b1(i,jep3,k) = emf3ojb3(k,i)
        enddo
      enddo
    endif
    if (nojb .eq. 5) then
      do i=ism2,iep3
        do k=ksm2,kep2
          v1b2mv2b1(i,jep1,k) = 0.
          v1b2mv2b1(i,jep2,k) = v1b2mv2b1(i,je  ,k)
          v1b2mv2b1(i,jep3,k) = v1b2mv2b1(i,jem1,k)
        enddo
      enddo
    endif
  endif

end subroutine bvalemf3

subroutine bvaldiffemf1(v2b3mv3b2)
  use ModPar
  use ModGrid
  use ModBoundary
  use ModUserSetup,   ONLY: getdiffemf1ib, getdiffemf1jb, getdiffemf1kb
  implicit none

  real, intent(inout), dimension(in,jn,kn) :: v2b3mv3b2

  integer idtop, idbot
  integer :: isendtop, isendbot, irecvtop, irecvbot
  real, allocatable :: sendtop(:,:,:), sendbot(:,:,:)
  real, allocatable :: recvtop(:,:,:), recvbot(:,:,:)
  integer i, j, k
  !-------------------------------------------------------
  !------------------ i - boundary ----------------------
  idtop = myid1 + 1 + myid2*nproc1 + myid3*(nproc1*nproc2)
  idbot = myid1 - 1 + myid2*nproc1 + myid3*(nproc1*nproc2)
  isendtop=2
  isendbot=1
  irecvtop=1
  irecvbot=2
  if(myid1 .eq. 0) idbot = nproc1 - 1 + nproc1*myid2 &
                         + (nproc1*nproc2)*myid3
  if(myid1 .eq. nproc1-1) idtop = 0 + nproc1*myid2 &
                                + (nproc1*nproc2)*myid3
  allocate(recvtop(iep1:iep2,js:jep1,ks:kep1), &
           recvbot(ism2:ism1,js:jep1,ks:kep1))
  allocate(sendtop(iem1:ie,js:jep1,ks:kep1), &
           sendbot(is:isp1,js:jep1,ks:kep1))
  sendtop=v2b3mv3b2(iem1:ie,js:jep1,ks:kep1)
  sendbot=v2b3mv3b2(is:isp1,js:jep1,ks:kep1)
  call message_pass_nodes(idtop, idbot, &
    2, jep1-js+1, kep1-ks+1, 2, jep1-js+1, kep1-ks+1, &
    sendtop, sendbot, recvtop, recvbot, &
    isendtop, isendbot, irecvtop, irecvbot)
  v2b3mv3b2(iep1:iep2,js:jep1,ks:kep1) = recvtop
  v2b3mv3b2(ism2:ism1,js:jep1,ks:kep1) = recvbot
  deallocate(recvtop, recvbot, sendtop, sendbot)
  ! Inner i boundary.
  if(niib .eq. 3 .or. noib .eq. 3) then
    call getdiffemf1ib(v2b3mv3b2)
  endif
  if(myid1 .eq. 0) then
    if ( niib .eq. 1) then
      do k=ks,kep1
        do j=js,jep1
          v2b3mv3b2(ism1,j,k) = v2b3mv3b2(is  ,j,k)
          v2b3mv3b2(ism2,j,k) = v2b3mv3b2(isp1,j,k)
        enddo
      enddo
    endif
    if ( niib .eq. 2) then
      do k=ks,kep1
        do j=js,jep1
          v2b3mv3b2(ism1,j,k) = v2b3mv3b2(is  ,j,k)*dx1a(is)*dx1ai(ism1)
          v2b3mv3b2(ism2,j,k) = v2b3mv3b2(is  ,j,k)*dx1a(is)*dx1ai(ism2)
        enddo
      enddo
    endif
    if ( niib .eq. 3) then
      do k=ks,kep1
        do j=js,jep1
          v2b3mv3b2(ism1,j,k) = emf1iib1(j,k)
          v2b3mv3b2(ism2,j,k) = emf1iib2(j,k)
        enddo
      enddo
    endif
    if ( niib .eq. 5) then
      do k=ks,kep1
        do j=js,jep1
          v2b3mv3b2(ism1,j,k) =-v2b3mv3b2(is  ,j,k)
          v2b3mv3b2(ism2,j,k) =-v2b3mv3b2(isp1,j,k)
        enddo
      enddo
    endif
  endif
  ! Outer i boundary.
  if(myid1 .eq. nproc1-1) then
    if ( noib .eq. 1) then
      do k=ks,kep1
        do j=js,jep1
          v2b3mv3b2(iep1,j,k) = v2b3mv3b2(ie  ,j,k)
          v2b3mv3b2(iep2,j,k) = v2b3mv3b2(iem1,j,k)
        enddo
      enddo
    endif
    if ( noib .eq. 2) then
      do k=ks,kep1
        do j=js,jep1
          v2b3mv3b2(iep1,j,k) = v2b3mv3b2(ie  ,j,k)*dx1a(ie)*dx1ai(iep1)
          v2b3mv3b2(iep2,j,k) = v2b3mv3b2(ie  ,j,k)*dx1a(ie)*dx1ai(iep2)
        enddo
      enddo
    endif
    if ( noib .eq. 3) then
      do k=ks,kep1
        do j=js,jep1
          v2b3mv3b2(iep1,j,k) = emf1oib1(j,k)
          v2b3mv3b2(iep2,j,k) = emf1oib2(j,k)
        enddo
      enddo
    endif
    if ( noib .eq. 5) then
      do k=ks,kep1
        do j=js,jep1
          v2b3mv3b2(iep1,j,k) =-v2b3mv3b2(ie  ,j,k)
          v2b3mv3b2(iep2,j,k) =-v2b3mv3b2(iem1,j,k)
        enddo
      enddo
    endif
  endif
  !------------------ j - boundary -----------------
  idtop = myid1 + (myid2 + 1)*nproc1 + myid3*(nproc1*nproc2)
  idbot = myid1 + (myid2 - 1)*nproc1 + myid3*(nproc1*nproc2)
  isendtop=2
  isendbot=1
  irecvtop=1
  irecvbot=2
  if(myid2 .eq. 0) idbot = myid1 + (nproc2 - 1)*nproc1 &
                         + myid3*(nproc1*nproc2)
  if(myid2 .eq. nproc2-1) idtop = myid1 + 0*nproc1 &
                                + myid3*(nproc1*nproc2)
  allocate(recvtop(ism2:iep2,jep2:jep3,ks:kep1), &
           recvbot(ism2:iep2,jsm2:jsm1,ks:kep1))
  allocate(sendtop(ism2:iep2,jem1:je,ks:kep1), &
           sendbot(ism2:iep2,jsp1:jsp2,ks:kep1))
  sendtop=v2b3mv3b2(ism2:iep2,jem1:je,ks:kep1)
  sendbot=v2b3mv3b2(ism2:iep2,jsp1:jsp2,ks:kep1)
  call message_pass_nodes(idtop, idbot, &
    iep2-ism2+1, 2, kep1-ks+1, iep2-ism2+1, 2, kep1-ks+1, &
    sendtop, sendbot, recvtop, recvbot, &
    isendtop, isendbot, irecvtop, irecvbot)
  v2b3mv3b2(ism2:iep2,jep2:jep3,ks:kep1) = recvtop
  v2b3mv3b2(ism2:iep2,jsm2:jsm1,ks:kep1) = recvbot
  deallocate(recvtop,recvbot,sendtop,sendbot)
  ! Inner j boundary.
  if(nijb .eq. 3 .or. nojb .eq. 3) then
    call getdiffemf1jb(v2b3mv3b2)
  endif
  if(myid2 .eq. 0) then
    if (nijb .eq. 1) then
      do i=ism2,iep2
        do k=ks,kep1
          v2b3mv3b2(i,js  ,k) = 0.
          v2b3mv3b2(i,jsm1,k) =-v2b3mv3b2(i,jsp1,k)
          v2b3mv3b2(i,jsm2,k) =-v2b3mv3b2(i,jsp2,k)
        enddo
      enddo
    endif
    if (nijb .eq. 2) then
      do i=ism2,iep2
        do k=ks,kep1
          v2b3mv3b2(i,jsm1,k) = v2b3mv3b2(i,js  ,k) &
            + v2b3mv3b2(i,js,k) - v2b3mv3b2(i,jsp1,k)
          v2b3mv3b2(i,jsm2,k) = v2b3mv3b2(i,jsm1  ,k) &
            + v2b3mv3b2(i,jsm1,k) - v2b3mv3b2(i,js,k)
        enddo
      enddo
    endif
    if (nijb .eq. 3) then
      do i=ism2,iep2
        do k=ks,kep1
          v2b3mv3b2(i,js  ,k) = emf1ijb1(k,i)
          v2b3mv3b2(i,jsm1,k) = emf1ijb2(k,i)
          v2b3mv3b2(i,jsm2,k) = emf1ijb3(k,i)
        enddo
      enddo
    endif
    if (nijb .eq. 5) then
      do i=ism2,iep2
        do k=ks,kep1
          v2b3mv3b2(i,js  ,k) = 0.
          v2b3mv3b2(i,jsm1,k) = v2b3mv3b2(i,jsp1,k)
          v2b3mv3b2(i,jsm2,k) = v2b3mv3b2(i,jsp2,k)
        enddo
      enddo
    endif
  endif
  ! Outer j boundary.
  if(myid2 .eq. nproc2-1) then
    if (nojb .eq. 1) then
      do i=ism2,iep2
        do k=ks,kep1
          v2b3mv3b2(i,jep1,k) = 0.
          v2b3mv3b2(i,jep2,k) =-v2b3mv3b2(i,je  ,k)
          v2b3mv3b2(i,jep3,k) =-v2b3mv3b2(i,jem1,k)
        enddo
      enddo
    endif
    if (nojb .eq. 2) then
      do i=ism2,iep2
        do k=ks,kep1
          v2b3mv3b2(i,jep2,k) = v2b3mv3b2(i,jep1,k) &
            + v2b3mv3b2(i,jep1,k) - v2b3mv3b2(i,je,k)
          v2b3mv3b2(i,jep3,k) = v2b3mv3b2(i,jep2,k) &
            + v2b3mv3b2(i,jep2,k) - v2b3mv3b2(i,jep1,k)
        enddo
      enddo
    endif
    if (nojb .eq. 3) then
      do i=ism2,iep2
        do k=ks,kep1
          v2b3mv3b2(i,jep1,k) = emf1ojb1(k,i)
          v2b3mv3b2(i,jep2,k) = emf1ojb2(k,i)
          v2b3mv3b2(i,jep3,k) = emf1ojb3(k,i)
        enddo
      enddo
    endif
    if (nojb .eq. 5) then
      do i=ism2,iep2
        do k=ks,kep1
          v2b3mv3b2(i,jep1,k) = 0.
          v2b3mv3b2(i,jep2,k) = v2b3mv3b2(i,je  ,k)
          v2b3mv3b2(i,jep3,k) = v2b3mv3b2(i,jem1,k)
        enddo
      enddo
    endif
  endif
  !------------------ k - boundary -----------------
  idtop = myid1 + myid2*nproc1 + (myid3+1)*(nproc1*nproc2)
  idbot = myid1 + myid2*nproc1 + (myid3-1)*(nproc1*nproc2)
  isendtop=2
  isendbot=1
  irecvtop=1
  irecvbot=2
  if(myid3==0) idbot = myid1 + myid2*nproc1 + (nproc3-1)*(nproc1*nproc2)
  if(myid3==nproc3-1) idtop = myid1 + myid2*nproc1 + 0*(nproc1*nproc2)
  allocate(recvtop(ism2:iep2,jsm2:jep3,kep2:kep3), &
           recvbot(ism2:iep2,jsm2:jep3,ksm2:ksm1))
  allocate(sendtop(ism2:iep2,jsm2:jep3,kem1:ke), &
           sendbot(ism2:iep2,jsm2:jep3,ksp1:ksp2))
  sendtop=v2b3mv3b2(ism2:iep2,jsm2:jep3,kem1:ke)
  sendbot=v2b3mv3b2(ism2:iep2,jsm2:jep3,ksp1:ksp2)
  call message_pass_nodes(idtop, idbot, &
    iep2-ism2+1, jep3-jsm2+1, 2, iep2-ism2+1, jep3-jsm2+1, 2, &
    sendtop, sendbot, recvtop, recvbot, &
    isendtop, isendbot, irecvtop, irecvbot)
  v2b3mv3b2(ism2:iep2,jsm2:jep3,kep2:kep3) = recvtop
  v2b3mv3b2(ism2:iep2,jsm2:jep3,ksm2:ksm1) = recvbot
  deallocate(recvtop,recvbot,sendtop,sendbot)
  ! Inner k boundary
  if(nikb .eq. 3 .or. nokb .eq. 3) then
    call getdiffemf1kb(v2b3mv3b2)
  endif
  if(myid3 .eq. 0) then
    if ( nikb .eq. 1) then
      do j=jsm2,jep3
        do i=ism2,iep2
          v2b3mv3b2(i,j,ks  ) = 0.
          v2b3mv3b2(i,j,ksm1) =-v2b3mv3b2(i,j,ksp1)
          v2b3mv3b2(i,j,ksm2) =-v2b3mv3b2(i,j,ksp2)
        enddo
      enddo
    endif
    if ( nikb .eq. 2) then
      do j=jsm2,jep3
        do i=ism2,iep2
          v2b3mv3b2(i,j,ksm1) = v2b3mv3b2(i,j,ks  ) &
            + v2b3mv3b2(i,j,ks) - v2b3mv3b2(i,j,ksp1)
          v2b3mv3b2(i,j,ksm2) = v2b3mv3b2(i,j,ksm1  ) &
            + v2b3mv3b2(i,j,ksm1) - v2b3mv3b2(i,j,ks)
        enddo
      enddo
    endif
    if ( nikb .eq. 3) then
      do j=jsm2,jep3
        do i=ism2,iep2
          v2b3mv3b2(i,j,ks  ) = emf1ikb1(i,j)
          v2b3mv3b2(i,j,ksm1) = emf1ikb2(i,j)
          v2b3mv3b2(i,j,ksm2) = emf1ikb3(i,j)
        enddo
      enddo
    endif
    if ( nikb .eq. 5) then
      do j=jsm2,jep3
        do i=ism2,iep2
          v2b3mv3b2(i,j,ks  ) = 0.
          v2b3mv3b2(i,j,ksm1) = v2b3mv3b2(i,j,ksp1)
          v2b3mv3b2(i,j,ksm2) = v2b3mv3b2(i,j,ksp2)
        enddo
      enddo
    endif
  endif
  ! Outer k boundary.
  if(myid3 .eq. nproc3-1) then
    if ( nokb .eq. 1) then
      do j=jsm2,jep3
        do i=ism2,iep2
          v2b3mv3b2(i,j,kep1) = 0.
          v2b3mv3b2(i,j,kep2) =-v2b3mv3b2(i,j,ke  )
          v2b3mv3b2(i,j,kep3) =-v2b3mv3b2(i,j,kem1)
        enddo
      enddo
    endif
    if ( nokb .eq. 2) then
      do j=jsm2,jep3
        do i=ism2,iep2
          v2b3mv3b2(i,j,kep2) = v2b3mv3b2(i,j,kep1) &
            + v2b3mv3b2(i,j,kep1) - v2b3mv3b2(i,j,ke)
          v2b3mv3b2(i,j,kep3) = v2b3mv3b2(i,j,kep2) &
            + v2b3mv3b2(i,j,kep2) - v2b3mv3b2(i,j,kep1)
        enddo
      enddo
    endif
    if ( nokb .eq. 3) then
      do j=jsm2,jep3
        do i=ism2,iep2
          v2b3mv3b2(i,j,kep1) = emf1okb1(i,j)
          v2b3mv3b2(i,j,kep2) = emf1okb2(i,j)
          v2b3mv3b2(i,j,kep3) = emf1okb3(i,j)
        enddo
      enddo
    endif
    if ( nokb .eq. 5) then
      do j=jsm2,jep3
        do i=ism2,iep2
          v2b3mv3b2(i,j,kep1) = 0.
          v2b3mv3b2(i,j,kep2) = v2b3mv3b2(i,j,ke  )
          v2b3mv3b2(i,j,kep3) = v2b3mv3b2(i,j,kem1)
        enddo
      enddo
    endif
  endif

end subroutine bvaldiffemf1

subroutine bvaldiffemf2(v3b1mv1b3)
  use ModPar
  use ModGrid
  use ModBoundary
  use ModUserSetup,   ONLY: getdiffemf2ib, getdiffemf2jb, getdiffemf2kb
  implicit none

  real, intent(inout), dimension(in,jn,kn) :: v3b1mv1b3

  integer idtop, idbot
  integer :: isendtop, isendbot, irecvtop, irecvbot
  real, allocatable :: sendtop(:,:,:), sendbot(:,:,:)
  real, allocatable :: recvtop(:,:,:), recvbot(:,:,:)
  integer i, j, k
  !-----------------------------------------------------------------------------
  !------------------ j - boundary ----------------- 
  idtop = myid1 + (myid2 + 1)*nproc1 + myid3*(nproc1*nproc2)
  idbot = myid1 + (myid2 - 1)*nproc1 + myid3*(nproc1*nproc2)
  isendtop=2
  isendbot=1
  irecvtop=1
  irecvbot=2
  if(myid2 .eq. 0) idbot = myid1 + (nproc2 - 1)*nproc1 &
                         + myid3*(nproc1*nproc2)
  if(myid2 .eq. nproc2-1) idtop = myid1 + 0*nproc1 &
                                + myid3*(nproc1*nproc2)
  allocate(recvtop(is:iep1,jep1:jep2,ks:kep1), &
           recvbot(is:iep1,jsm2:jsm1,ks:kep1))
  allocate(sendtop(is:iep1,jem1:je, ks:kep1), &
           sendbot(is:iep1,js:jsp1, ks:kep1))
  sendtop=v3b1mv1b3(is:iep1,jem1:je, ks:kep1)
  sendbot=v3b1mv1b3(is:iep1,js:jsp1, ks:kep1)
  call message_pass_nodes(idtop, idbot, &
    iep1-is+1, 2, kep1-ks+1, iep1-is+1, 2, kep1-ks+1, &
    sendtop, sendbot, recvtop, recvbot, &
    isendtop, isendbot, irecvtop, irecvbot)
    v3b1mv1b3(is:iep1,jep1:jep2,ks:kep1) = recvtop
    v3b1mv1b3(is:iep1,jsm2:jsm1,ks:kep1) = recvbot
  deallocate(recvtop,recvbot,sendtop,sendbot)
  ! Inner j boundary.
  if(nijb .eq. 3 .or. nojb .eq. 3) then
    call getdiffemf2jb(v3b1mv1b3)
  endif
  if(myid2 .eq. 0) then
    if (nijb .eq. 1) then
      do i=is,iep1
        do k=ks,kep1
          v3b1mv1b3(i,jsm1,k) = v3b1mv1b3(i,js  ,k)
          v3b1mv1b3(i,jsm2,k) = v3b1mv1b3(i,jsp1,k)
        enddo
      enddo
    endif
    if (nijb .eq. 2) then
      do i=is,iep1
        do k=ks,kep1
          v3b1mv1b3(i,jsm1,k) = v3b1mv1b3(i,js  ,k)*dx2a(js)*dx2ai(jsm1)
          v3b1mv1b3(i,jsm2,k) = v3b1mv1b3(i,js  ,k)*dx2a(js)*dx2ai(jsm2)
        enddo
      enddo
    endif
    if (nijb .eq. 3) then
      do i=is,iep1
        do k=ks,kep1
             v3b1mv1b3(i,jsm1,k) = emf2ijb1(k,i)
             v3b1mv1b3(i,jsm2,k) = emf2ijb2(k,i)
        enddo
      enddo
    endif
    if (nijb .eq. 5) then
      do i=is,iep1
        do k=ks,kep1
          v3b1mv1b3(i,jsm1,k) =-v3b1mv1b3(i,js  ,k)
          v3b1mv1b3(i,jsm2,k) =-v3b1mv1b3(i,jsp1,k)
        enddo
      enddo
    endif
  endif
  ! Outer j boundary.
  if(myid2 .eq. nproc2-1) then
    if (nojb .eq. 1) then
      do i=is,iep1
        do k=ks,kep1
          v3b1mv1b3(i,jep1,k) = v3b1mv1b3(i,je  ,k)
          v3b1mv1b3(i,jep2,k) = v3b1mv1b3(i,jem1,k)
        enddo
      enddo
    endif
    if (nojb .eq. 2) then
      do i=is,iep1
        do k=ks,kep1
          v3b1mv1b3(i,jep1,k) = v3b1mv1b3(i,je  ,k)*dx2a(je)*dx2ai(jep1)
          v3b1mv1b3(i,jep2,k) = v3b1mv1b3(i,je  ,k)*dx2a(je)*dx2ai(jep2)
        enddo
      enddo
    endif
    if (nojb .eq. 3) then
      do i=is,iep1
        do k=ks,kep1
          v3b1mv1b3(i,jep1,k) = emf2ojb1(k,i)
          v3b1mv1b3(i,jep2,k) = emf2ojb2(k,i)
        enddo
      enddo
    endif
    if (nojb .eq. 5) then
      do i=is,iep1
        do k=ks,kep1
          v3b1mv1b3(i,jep1,k) =-v3b1mv1b3(i,je  ,k)
          v3b1mv1b3(i,jep2,k) =-v3b1mv1b3(i,jem1,k)
        enddo
      enddo
    endif
  endif
  !------------------ k - boundary -----------------
  idtop = myid1 + myid2*nproc1 + (myid3+1)*(nproc1*nproc2)
  idbot = myid1 + myid2*nproc1 + (myid3-1)*(nproc1*nproc2)
  isendtop=2
  isendbot=1
  irecvtop=1
  irecvbot=2
  if(myid3==0) idbot = myid1 + myid2*nproc1 + (nproc3-1)*(nproc1*nproc2)
  if(myid3==nproc3-1) idtop = myid1 + myid2*nproc1 + 0*(nproc1*nproc2)
  allocate(recvtop(is:iep1,jsm2:jep2,kep2:kep3), &
           recvbot(is:iep1,jsm2:jep2,ksm2:ksm1))
  allocate(sendtop(is:iep1,jsm2:jep2,kem1:ke), &
           sendbot(is:iep1,jsm2:jep2,ksp1:ksp2))
  sendtop=v3b1mv1b3(is:iep1,jsm2:jep2,kem1:ke)
  sendbot=v3b1mv1b3(is:iep1,jsm2:jep2,ksp1:ksp2)
  call message_pass_nodes(idtop, idbot, &
    iep1-is+1, jep2-jsm2+1, 2, iep1-is+1, jep2-jsm2+1, 2, &
    sendtop, sendbot, recvtop, recvbot, &
    isendtop, isendbot, irecvtop, irecvbot)
  v3b1mv1b3(is:iep1,jsm2:jep2,kep2:kep3) = recvtop
  v3b1mv1b3(is:iep1,jsm2:jep2,ksm2:ksm1) = recvbot
  deallocate(recvtop, recvbot, sendtop, sendbot)
  ! Inner k boundary
  if(nikb .eq. 3 .or. nokb .eq. 3) then
    call getdiffemf2kb(v3b1mv1b3)
  endif
  if(myid3 .eq. 0) then
    if (nikb .eq. 1) then
      do j=jsm2,jep2
        do i=is,iep1
          v3b1mv1b3(i,j,ks  ) = 0.
          v3b1mv1b3(i,j,ksm1) =-v3b1mv1b3(i,j,ksp1)
          v3b1mv1b3(i,j,ksm2) =-v3b1mv1b3(i,j,ksp2)
        enddo
      enddo
    endif
    if (nikb .eq. 2) then
      do j=jsm2,jep2
        do i=is,iep1
          v3b1mv1b3(i,j,ksm1) = v3b1mv1b3(i,j,ks  ) &
            + v3b1mv1b3(i,j,ks) - v3b1mv1b3(i,j,ksp1)
          v3b1mv1b3(i,j,ksm2) = v3b1mv1b3(i,j,ksm1  ) &
            + v3b1mv1b3(i,j,ksm1) - v3b1mv1b3(i,j,ks)
        enddo
      enddo
    endif
    if (nikb .eq. 3) then
      do j=jsm2,jep2
        do i=is,iep1
          v3b1mv1b3(i,j,ks  ) = emf2ikb1(i,j)
          v3b1mv1b3(i,j,ksm1) = emf2ikb2(i,j)
          v3b1mv1b3(i,j,ksm2) = emf2ikb3(i,j)
        enddo
      enddo
    endif
    if (nikb .eq. 5) then
      do j=jsm2,jep2
        do i=is,iep1
          v3b1mv1b3(i,j,ks  ) = 0.
          v3b1mv1b3(i,j,ksm1) = v3b1mv1b3(i,j,ksp1)
          v3b1mv1b3(i,j,ksm2) = v3b1mv1b3(i,j,ksp2)
        enddo
      enddo
    endif
  endif
  ! Outer k boundary
  if(myid3 .eq. nproc3-1) then
    if (nokb .eq. 1) then
      do j=jsm2,jep2
        do i=is,iep1
          v3b1mv1b3(i,j,kep1) = 0.
          v3b1mv1b3(i,j,kep2) =-v3b1mv1b3(i,j,ke  )
          v3b1mv1b3(i,j,kep3) =-v3b1mv1b3(i,j,kem1)
        enddo
      enddo
    endif
    if (nokb .eq. 2) then
      do j=jsm2,jep2
        do i=is,iep1
          v3b1mv1b3(i,j,kep2) = v3b1mv1b3(i,j,kep1) &
            + v3b1mv1b3(i,j,kep1) - v3b1mv1b3(i,j,ke)
          v3b1mv1b3(i,j,kep3) = v3b1mv1b3(i,j,kep2) &
            + v3b1mv1b3(i,j,kep2) - v3b1mv1b3(i,j,kep1)
        enddo
      enddo
    endif
    if (nokb .eq. 3) then
      do j=jsm2,jep2
        do i=is,iep1
          v3b1mv1b3(i,j,kep1) = emf2okb1(i,j)
          v3b1mv1b3(i,j,kep2) = emf2okb2(i,j)
          v3b1mv1b3(i,j,kep3) = emf2okb3(i,j)
        enddo
      enddo
    endif
    if (nokb .eq. 5) then
      do j=jsm2,jep2
        do i=is,iep1
          v3b1mv1b3(i,j,kep1) = 0.
          v3b1mv1b3(i,j,kep2) = v3b1mv1b3(i,j,ke  )
          v3b1mv1b3(i,j,kep3) = v3b1mv1b3(i,j,kem1)
        enddo
      enddo
    endif
  endif
  !-------------------------------------------------------
  !------------------ i - boundary ----------------------
  idtop = myid1 + 1 + myid2*nproc1 + myid3*(nproc1*nproc2)
  idbot = myid1 - 1 + myid2*nproc1 + myid3*(nproc1*nproc2)
  isendtop=2
  isendbot=1
  irecvtop=1
  irecvbot=2
  if(myid1 .eq. 0) idbot = nproc1 - 1 + nproc1*myid2 &
                         + (nproc1*nproc2)*myid3
  if(myid1 .eq. nproc1-1) idtop = 0 + nproc1*myid2 &
                                + (nproc1*nproc2)*myid3
  allocate(recvtop(iep2:iep3,jsm2:jep2,ksm2:kep3), &
           recvbot(ism2:ism1,jsm2:jep2,ksm2:kep3))
  allocate(sendtop(iem1:ie,jsm2:jep2,ksm2:kep3), &
           sendbot(isp1:isp2,jsm2:jep2,ksm2:kep3))
  sendtop=v3b1mv1b3(iem1:ie,jsm2:jep2,ksm2:kep3)
  sendbot=v3b1mv1b3(isp1:isp2,jsm2:jep2,ksm2:kep3)
  call message_pass_nodes(idtop, idbot, &
    2, jep2-jsm2+1, kep3-ksm2+1, 2, jep2-jsm2+1, kep3-ksm2+1, &
    sendtop, sendbot, recvtop, recvbot, &
    isendtop, isendbot, irecvtop, irecvbot)
  v3b1mv1b3(iep2:iep3,jsm2:jep2,ksm2:kep3) = recvtop
  v3b1mv1b3(ism2:ism1,jsm2:jep2,ksm2:kep3) = recvbot
  deallocate(recvtop, recvbot, sendtop, sendbot)
  ! Inner i boundary.
  if(niib .eq. 3 .or. noib .eq. 3) then
    call getdiffemf2ib(v3b1mv1b3)
  endif
  if(myid1 .eq. 0) then
    if (niib .eq. 1) then
      do k=ksm2,kep3
        do j=jsm2,jep2
             v3b1mv1b3(is  ,j,k) = 0.
             v3b1mv1b3(ism1,j,k) =-v3b1mv1b3(isp1,j,k)
             v3b1mv1b3(ism2,j,k) =-v3b1mv1b3(isp2,j,k)
        enddo
      enddo
    endif
    if (niib .eq. 2) then
      do k=ksm2,kep3
        do j=jsm2,jep2
          v3b1mv1b3(ism1,j,k) = v3b1mv1b3(is  ,j,k)*g2a(is)*g2ai(ism1) &
            + v3b1mv1b3(is  ,j,k)*g2a(is)*g2ai(ism1) &
            - v3b1mv1b3(isp1,j,k)*g2a(isp1)*g2ai(ism1)
          v3b1mv1b3(ism2,j,k) = v3b1mv1b3(ism1  ,j,k)*g2a(ism1)*g2ai(ism2) &
            + v3b1mv1b3(ism1  ,j,k)*g2a(ism1)*g2ai(ism2) &
            - v3b1mv1b3(is,j,k)*g2a(is)*g2ai(ism2)
        enddo
      enddo
    endif
    if (niib .eq. 3) then
      do k=ksm2,kep3
        do j=jsm2,jep2
          v3b1mv1b3(is  ,j,k) = emf2iib1(j,k)
          v3b1mv1b3(ism1,j,k) = emf2iib2(j,k)
          v3b1mv1b3(ism2,j,k) = emf2iib3(j,k)
        enddo
      enddo
    endif
    if (niib .eq. 5) then
      do k=ksm2,kep3
        do j=jsm2,jep2
          v3b1mv1b3(is  ,j,k) = 0.
          v3b1mv1b3(ism1,j,k) = v3b1mv1b3(isp1,j,k)
          v3b1mv1b3(ism2,j,k) = v3b1mv1b3(isp2,j,k)
        enddo
      enddo
    endif

  endif
  ! Outer i boundary.
  if(myid1 .eq. nproc1-1) then
    if (noib .eq. 1) then
      do k=ksm2,kep3
        do j=jsm2,jep2
          v3b1mv1b3(iep1,j,k) = 0.
          v3b1mv1b3(iep2,j,k) =-v3b1mv1b3(ie  ,j,k)
          v3b1mv1b3(iep3,j,k) =-v3b1mv1b3(iem1,j,k)
        enddo
      enddo
    endif
    if (noib .eq. 2) then
      do k=ksm2,kep3
        do j=jsm2,jep2
          v3b1mv1b3(iep2,j,k) = v3b1mv1b3(iep1,j,k)*g2a(iep1)*g2ai(iep2) &
            + v3b1mv1b3(iep1,j,k)*g2a(iep1)*g2ai(iep2) &
            - v3b1mv1b3(ie,j,k)*g2a(ie)*g2ai(iep2)
          v3b1mv1b3(iep3,j,k) = v3b1mv1b3(iep2,j,k)*g2a(iep2)*g2ai(iep3) &
            + v3b1mv1b3(iep2,j,k)*g2a(iep2)*g2ai(iep3) &
            - v3b1mv1b3(iep1,j,k)*g2a(iep1)*g2ai(iep3)
        enddo
      enddo
    endif
    if (noib .eq. 3) then
      do k=ksm2,kep3
        do j=jsm2,jep2
          v3b1mv1b3(iep1,j,k) = emf2oib1(j,k)
          v3b1mv1b3(iep2,j,k) = emf2oib2(j,k)
          v3b1mv1b3(iep3,j,k) = emf2oib3(j,k)
        enddo
      enddo
    endif
    if (noib .eq. 5) then
      do k=ksm2,kep3
        do j=jsm2,jep2
          v3b1mv1b3(iep1,j,k) = 0.
          v3b1mv1b3(iep2,j,k) = v3b1mv1b3(ie  ,j,k)
          v3b1mv1b3(iep3,j,k) = v3b1mv1b3(iem1,j,k)
        enddo
      enddo
    endif
  endif

end subroutine bvaldiffemf2

subroutine bvaldiffemf3(v1b2mv2b1)
  use ModPar
  use ModGrid
  use ModBoundary
  use ModUserSetup,   ONLY: getdiffemf3ib, getdiffemf3jb, getdiffemf3kb
  implicit none

  real, intent(inout), dimension(in,jn,kn) :: v1b2mv2b1

  integer idtop, idbot
  integer :: isendtop, isendbot, irecvtop, irecvbot
  real, allocatable :: sendtop(:,:,:), sendbot(:,:,:)
  real, allocatable :: recvtop(:,:,:), recvbot(:,:,:)
  integer i, j, k
  !-----------------------------------------------------------------------------
  !------------------ k - boundary -----------------
  idtop = myid1 + myid2*nproc1 + (myid3+1)*(nproc1*nproc2)
  idbot = myid1 + myid2*nproc1 + (myid3-1)*(nproc1*nproc2)
  isendtop=2
  isendbot=1
  irecvtop=1
  irecvbot=2
  if(myid3==0) idbot = myid1 + myid2*nproc1 + (nproc3-1)*(nproc1*nproc2)
  if(myid3==nproc3-1) idtop = myid1 + myid2*nproc1 + 0*(nproc1*nproc2)
  allocate(recvtop(is:iep1,js:jep1,kep1:kep2), &
           recvbot(is:iep1,js:jep1,ksm2:ksm1))
  allocate(sendtop(is:iep1,js:jep1,kem1:ke), &
           sendbot(is:iep1,js:jep1,ks:ksp1))
  sendtop=v1b2mv2b1(is:iep1,js:jep1,kem1:ke)
  sendbot=v1b2mv2b1(is:iep1,js:jep1,ks:ksp1)
  call message_pass_nodes(idtop, idbot, &
    iep1-is+1, jep1-js+1, 2, iep1-is+1, jep1-js+1, 2, &
    sendtop, sendbot, recvtop, recvbot, &
    isendtop, isendbot, irecvtop, irecvbot)
  v1b2mv2b1(is:iep1,js:jep1,kep1:kep2) = recvtop
  v1b2mv2b1(is:iep1,js:jep1,ksm2:ksm1) = recvbot
  deallocate(recvtop, recvbot, sendtop, sendbot)
  ! Inner k boundary
  if(nikb .eq. 3 .or. nokb .eq. 3) then
    call getdiffemf3kb( v1b2mv2b1 )
  endif
  if(myid3 .eq. 0) then
    if (nikb .eq. 1) then
      do j=js,jep1
        do i=is,iep1
          v1b2mv2b1(i,j,ksm1) = v1b2mv2b1(i,j,ks  )
          v1b2mv2b1(i,j,ksm2) = v1b2mv2b1(i,j,ksp1)
        enddo
      enddo
    endif
    if (nikb .eq. 2) then
      do j=js,jep1
        do i=is,iep1
          v1b2mv2b1(i,j,ksm1) = v1b2mv2b1(i,j,ks  )*dx3a(ks)*dx3ai(ksm1)
          v1b2mv2b1(i,j,ksm2) = v1b2mv2b1(i,j,ks  )*dx3a(ks)*dx3ai(ksm2)
        enddo
      enddo
    endif
    if (nikb .eq. 3) then
      do j=js,jep1
        do i=is,iep1
          v1b2mv2b1(i,j,ksm1) = emf3ikb1(i,j)
          v1b2mv2b1(i,j,ksm2) = emf3ikb2(i,j)
        enddo
      enddo
    endif
    if (nikb .eq. 5) then
      do j=js,jep1
        do i=is,iep1
          v1b2mv2b1(i,j,ksm1) =-v1b2mv2b1(i,j,ks  )
          v1b2mv2b1(i,j,ksm2) =-v1b2mv2b1(i,j,ksp1)
        enddo
      enddo
    endif
  endif
  ! Outer k boundary
  if(myid3 .eq. nproc3-1) then
    if (nokb .eq. 1) then
      do j=js,jep1
        do i=is,iep1
          v1b2mv2b1(i,j,kep1) = v1b2mv2b1(i,j,ke  )
          v1b2mv2b1(i,j,kep2) = v1b2mv2b1(i,j,kem1)
        enddo
      enddo
    endif
    if (nokb .eq. 2) then
      do j=js,jep1
        do i=is,iep1
          v1b2mv2b1(i,j,kep1) = v1b2mv2b1(i,j,ke  )*dx3a(ke)*dx3ai(kep1)
          v1b2mv2b1(i,j,kep2) = v1b2mv2b1(i,j,ke  )*dx3a(ke)*dx3ai(kep2)
        enddo
      enddo
    endif
    if (nokb .eq. 3) then
      do j=js,jep1
        do i=is,iep1
          v1b2mv2b1(i,j,kep1) = emf3okb1(i,j)
          v1b2mv2b1(i,j,kep2) = emf3okb2(i,j)
        enddo
      enddo
    endif
    if (nokb .eq. 5) then
      do j=js,jep1
        do i=is,iep1
          v1b2mv2b1(i,j,kep1) =-v1b2mv2b1(i,j,ke  )
          v1b2mv2b1(i,j,kep2) =-v1b2mv2b1(i,j,kem1)
        enddo
      enddo
    endif
  endif
  !------------------ i - boundary -----------------
  idtop = myid1 + 1 + myid2*nproc1 + myid3*(nproc1*nproc2)
  idbot = myid1 - 1 + myid2*nproc1 + myid3*(nproc1*nproc2)
  isendtop=2
  isendbot=1
  irecvtop=1
  irecvbot=2
  if(myid1 .eq. 0) idbot = nproc1 - 1 + nproc1*myid2 &
                         + (nproc1*nproc2)*myid3
  if(myid1 .eq. nproc1-1) idtop = 0 + nproc1*myid2 &
                                + (nproc1*nproc2)*myid3
  allocate(recvtop(iep2:iep3,js:jep1,ksm2:kep2), &
           recvbot(ism2:ism1,js:jep1,ksm2:kep2))
  allocate(sendtop(iem1:ie,js:jep1,ksm2:kep2), &
           sendbot(isp1:isp2,js:jep1,ksm2:kep2))
  sendtop=v1b2mv2b1(iem1:ie,js:jep1,ksm2:kep2)
  sendbot=v1b2mv2b1(isp1:isp2,js:jep1,ksm2:kep2)
  call message_pass_nodes(idtop, idbot, &
    2, jep1-js+1, kep2-ksm2+1, 2, jep1-js+1, kep2-ksm2+1, &
    sendtop, sendbot, recvtop, recvbot, &
    isendtop, isendbot, irecvtop, irecvbot)
  v1b2mv2b1(iep2:iep3,js:jep1,ksm2:kep2) = recvtop
  v1b2mv2b1(ism2:ism1,js:jep1,ksm2:kep2) = recvbot
  deallocate(recvtop, recvbot, sendtop, sendbot)
  ! Inner i boundary
  if(niib .eq. 3 .or. noib .eq. 3) then
    call getdiffemf3ib( v1b2mv2b1 )
  endif
  if(myid1 .eq. 0) then
    if (niib .eq. 1) then
      do k=ksm2,kep2
        do j=js,jep1
          v1b2mv2b1(is  ,j,k) = 0.
          v1b2mv2b1(ism1,j,k) =-v1b2mv2b1(isp1,j,k)
          v1b2mv2b1(ism2,j,k) =-v1b2mv2b1(isp2,j,k)
        enddo
      enddo
    endif
    if (niib .eq. 2) then
      do k=ksm2,kep2
        do j=js,jep1
          v1b2mv2b1(ism1,j,k) = v1b2mv2b1(is  ,j,k)*g31a(is)*g31ai(ism1) &
            + v1b2mv2b1(is  ,j,k)*g31a(is)*g31ai(ism1) &
            - v1b2mv2b1(isp1,j,k)*g31a(isp1)*g31ai(ism1)
          v1b2mv2b1(ism2,j,k) = v1b2mv2b1(ism1  ,j,k)*g31a(ism1)*g31ai(ism2) &
            + v1b2mv2b1(ism1  ,j,k)*g31a(ism1)*g31ai(ism2) &
            - v1b2mv2b1(is,j,k)*g31a(is)*g31ai(ism2)
        enddo
      enddo
    endif
    if (niib .eq. 3) then
      do k=ksm2,kep2
        do j=js,jep1
          v1b2mv2b1(is  ,j,k) = emf3iib1(j,k)
          v1b2mv2b1(ism1,j,k) = emf3iib2(j,k)
          v1b2mv2b1(ism2,j,k) = emf3iib3(j,k)
        enddo
      enddo
    endif
    if (niib .eq. 5) then
      do k=ksm2,kep2
        do j=js,jep1
          v1b2mv2b1(is  ,j,k) = 0.
          v1b2mv2b1(ism1,j,k) = v1b2mv2b1(isp1,j,k)
          v1b2mv2b1(ism2,j,k) = v1b2mv2b1(isp2,j,k)
        enddo
      enddo
    endif
  endif
  ! Outer i boundary
  if(myid1 .eq. nproc1-1) then
    if (noib .eq. 1) then
      do k=ksm2,kep2
        do j=js,jep1
          v1b2mv2b1(iep1,j,k) = 0.
          v1b2mv2b1(iep2,j,k) =-v1b2mv2b1(ie  ,j,k)
          v1b2mv2b1(iep3,j,k) =-v1b2mv2b1(iem1,j,k)
        enddo
      enddo
    endif
    if (noib .eq. 2) then
      do k=ksm2,kep2
        do j=js,jep1
             v1b2mv2b1(iep2,j,k) = v1b2mv2b1(iep1,j,k)*g31a(iep1)*g31ai(iep2) &
               + v1b2mv2b1(iep1,j,k)*g31a(iep1)*g31ai(iep2) &
               - v1b2mv2b1(ie,j,k)*g31a(ie)*g31ai(iep2)
             v1b2mv2b1(iep3,j,k) = v1b2mv2b1(iep2,j,k)*g31a(iep2)*g31ai(iep3) &
               + v1b2mv2b1(iep2,j,k)*g31a(iep2)*g31ai(iep3) &
               - v1b2mv2b1(iep1,j,k)*g31a(iep1)*g31ai(iep3)
        enddo
      enddo
    endif
    if (noib .eq. 3) then
      do k=ksm2,kep2
        do j=js,jep1
          v1b2mv2b1(iep1,j,k) = emf3oib1(j,k)
          v1b2mv2b1(iep2,j,k) = emf3oib2(j,k)
          v1b2mv2b1(iep3,j,k) = emf3oib3(j,k)
        enddo
      enddo
    endif
    if (noib .eq. 5) then
      do k=ksm2,kep2
        do j=js,jep1
          v1b2mv2b1(iep1,j,k) = 0.
          v1b2mv2b1(iep2,j,k) = v1b2mv2b1(ie  ,j,k)
          v1b2mv2b1(iep3,j,k) = v1b2mv2b1(iem1,j,k)
        enddo
      enddo
    endif
  endif
  !------------------ j - boundary -----------------
  idtop = myid1 + (myid2 + 1)*nproc1 + myid3*(nproc1*nproc2)
  idbot = myid1 + (myid2 - 1)*nproc1 + myid3*(nproc1*nproc2)
  isendtop=2
  isendbot=1
  irecvtop=1
  irecvbot=2
  if(myid2 .eq. 0) idbot = myid1 + (nproc2 - 1)*nproc1 &
                         + myid3*(nproc1*nproc2)
  if(myid2 .eq. nproc2-1) idtop = myid1 + 0*nproc1 &
                                + myid3*(nproc1*nproc2)
  allocate(recvtop(ism2:iep3,jep2:jep3,ksm2:kep2), &
           recvbot(ism2:iep3,jsm2:jsm1,ksm2:kep2))
  allocate(sendtop(ism2:iep3,jem1:je,ksm2:kep2), &
           sendbot(ism2:iep3,jsp1:jsp2,ksm2:kep2))
  sendtop=v1b2mv2b1(ism2:iep3,jem1:je,ksm2:kep2)
  sendbot=v1b2mv2b1(ism2:iep3,jsp1:jsp2,ksm2:kep2)
  call message_pass_nodes(idtop, idbot, &
    iep3-ism2+1, 2, kep2-ksm2+1, iep3-ism2+1, 2, kep2-ksm2+1, &
    sendtop, sendbot, recvtop, recvbot, &
    isendtop, isendbot, irecvtop, irecvbot)
  v1b2mv2b1(ism2:iep3,jep2:jep3,ksm2:kep2) = recvtop
  v1b2mv2b1(ism2:iep3,jsm2:jsm1,ksm2:kep2) = recvbot
  deallocate(recvtop,recvbot,sendtop,sendbot)
  ! Inner j boundary.   
  if(nijb .eq. 3 .or. nojb .eq. 3) then
    call getdiffemf3jb( v1b2mv2b1 )
  endif
  if(myid2 .eq. 0) then
    if ( nijb .eq. 1) then
      do i=ism2,iep3
        do k=ksm2,kep2
          v1b2mv2b1(i,js  ,k) = 0.
          v1b2mv2b1(i,jsm1,k) =-v1b2mv2b1(i,jsp1,k)
          v1b2mv2b1(i,jsm2,k) =-v1b2mv2b1(i,jsp2,k)
        enddo
      enddo
    endif
    if ( nijb .eq. 2) then
      do i=ism2,iep3
        do k=ksm2,kep2
             v1b2mv2b1(i,jsm1,k) = v1b2mv2b1(i,js  ,k)*g32a(js)*g32ai(jsm1) &
               + v1b2mv2b1(i,js,k)*g32a(js)*g32ai(jsm1) &
               - v1b2mv2b1(i,jsp1,k)*g32a(jsp1)*g32ai(jsm1)
             v1b2mv2b1(i,jsm2,k) = v1b2mv2b1(i,jsm1  ,k)*g32a(jsm1)*g32ai(jsm2) &
               + v1b2mv2b1(i,jsm1,k)*g32a(jsm1)*g32ai(jsm2) &
               - v1b2mv2b1(i,js,k)*g32a(js)*g32ai(jsm2)
        enddo
      enddo
    endif
    if ( nijb .eq. 3) then
      do i=ism2,iep3
        do k=ksm2,kep2
          v1b2mv2b1(i,js  ,k) = emf3ijb1(k,i)
          v1b2mv2b1(i,jsm1,k) = emf3ijb2(k,i)
          v1b2mv2b1(i,jsm2,k) = emf3ijb3(k,i)
        enddo
      enddo
    endif
    if ( nijb .eq. 5) then
      do i=ism2,iep3
        do k=ksm2,kep2
          v1b2mv2b1(i,js  ,k) = 0.
          v1b2mv2b1(i,jsm1,k) = v1b2mv2b1(i,jsp1,k)
          v1b2mv2b1(i,jsm2,k) = v1b2mv2b1(i,jsp2,k)
        enddo
      enddo
    endif
  endif
  ! Outer j boundary.                      
  if(myid2 .eq. nproc2-1) then
    if (nojb .eq. 1) then
      do i=ism2,iep3
        do k=ksm2,kep2
             v1b2mv2b1(i,jep1,k) = 0.
             v1b2mv2b1(i,jep2,k) =-v1b2mv2b1(i,je  ,k)
             v1b2mv2b1(i,jep3,k) =-v1b2mv2b1(i,jem1,k)
        enddo
      enddo
    endif
    if (nojb .eq. 2) then
      do i=ism2,iep3
        do k=ksm2,kep2
             v1b2mv2b1(i,jep2,k) = v1b2mv2b1(i,jep1,k)*g32a(jep1)*g32ai(jep2) &
               + v1b2mv2b1(i,jep1,k)*g32a(jep1)*g32ai(jep2) &
               - v1b2mv2b1(i,je,k)*g32a(je)*g32ai(jep2)
             v1b2mv2b1(i,jep3,k) = v1b2mv2b1(i,jep2,k)*g32a(jep2)*g32ai(jep3) &
               + v1b2mv2b1(i,jep2,k)*g32a(jep2)*g32ai(jep3) &
               - v1b2mv2b1(i,jep1,k)*g32a(jep1)*g32ai(jep3)
        enddo
      enddo
    endif
    if (nojb .eq. 3) then
      do i=ism2,iep3
        do k=ksm2,kep2
          v1b2mv2b1(i,jep1,k) = emf3ojb1(k,i)
          v1b2mv2b1(i,jep2,k) = emf3ojb2(k,i)
          v1b2mv2b1(i,jep3,k) = emf3ojb3(k,i)
        enddo
      enddo
    endif
    if (nojb .eq. 5) then
      do i=ism2,iep3
        do k=ksm2,kep2
          v1b2mv2b1(i,jep1,k) = 0.
          v1b2mv2b1(i,jep2,k) = v1b2mv2b1(i,je  ,k)
          v1b2mv2b1(i,jep3,k) = v1b2mv2b1(i,jem1,k)
        enddo
      enddo
    endif
  endif

end subroutine bvaldiffemf3

subroutine bvaleta(eta)
  use ModPar
  use ModGrid
  use ModBoundary
  implicit none

  real, intent(inout), dimension(in,jn,kn) :: eta

  integer :: idtop, idbot
  integer :: isendtop, isendbot, irecvtop, irecvbot
  real, allocatable :: sendtop(:,:,:), sendbot(:,:,:)
  real, allocatable :: recvtop(:,:,:), recvbot(:,:,:)
  integer i, j, k
  !---------------------------------------------------------
  !------------------ i - boundary -----------------
  idtop = myid1 + 1 + myid2*nproc1 + myid3*(nproc1*nproc2)
  idbot = myid1 - 1 + myid2*nproc1 + myid3*(nproc1*nproc2)
  isendtop=2
  isendbot=1
  irecvtop=1
  irecvbot=2
  if(myid1 .eq. 0) idbot = nproc1 - 1 + nproc1*myid2 &
                         + (nproc1*nproc2)*myid3
  if(myid1 .eq. nproc1-1) idtop = 0 + nproc1*myid2 &
                                + (nproc1*nproc2)*myid3
  allocate(recvtop(iep1:iep2,js:je,ks:ke), recvbot(ism2:ism1,js:je,ks:ke))
  allocate(sendtop(iem1:ie,js:je,ks:ke), sendbot(is:isp1,js:je,ks:ke))
  sendtop=eta(iem1:ie,js:je,ks:ke)
  sendbot=eta(is:isp1,js:je,ks:ke)
  call message_pass_nodes(idtop, idbot, &
    2, je-js+1, ke-ks+1, 2, je-js+1, ke-ks+1, &
    sendtop, sendbot, recvtop, recvbot, &
    isendtop, isendbot, irecvtop, irecvbot)
  eta(iep1:iep2,js:je,ks:ke) = recvtop
  eta(ism2:ism1,js:je,ks:ke) = recvbot
  deallocate(recvtop, recvbot, sendtop, sendbot)
  ! Inner i boundary
  if(myid1 .eq. 0) then
    if(niib .eq. 1) then
      do k=ks,ke
        do j=js,je
          eta(ism1,j,k) = eta(is  ,j,k)
          eta(ism2,j,k) = eta(isp1,j,k)
        enddo
      enddo
    endif
    if(niib .eq. 2) then
      do k=ks,ke
        do j=js,je
          eta(ism1,j,k) = eta(is  ,j,k)
          eta(ism2,j,k) = eta(ism1,j,k)
        enddo
      enddo
    endif
    if(niib .eq. 3) then
      do k=ks,ke
        do j=js,je
          eta(ism1,j,k) = eta(is  ,j,k)
          eta(ism2,j,k) = eta(ism1,j,k)
        enddo
      enddo
    endif
    if(niib .eq. 5) then
      do k=ks,ke
        do j=js,je
          eta(ism1,j,k) = eta(is  ,j,k)
          eta(ism2,j,k) = eta(isp1,j,k)
        enddo
      enddo
    endif
  endif
  ! Outer i boundary
  if(myid1 .eq. nproc1-1) then
    if(noib .eq. 1) then
      do k=ks,ke
        do j=js,je
          eta(iep1,j,k) = eta(ie,j,k)
          eta(iep2,j,k) = eta(iem1,j,k)
        enddo
      enddo
    endif
    if(noib .eq. 2) then
      do k=ks,ke
        do j=js,je
          eta(iep1,j,k) = eta(ie,j,k)
          eta(iep2,j,k) = eta(iep1,j,k)
        enddo
      enddo
    endif
    if(noib .eq. 3) then
      do k=ks,ke
        do j=js,je
          eta(iep1,j,k) = eta(ie,j,k)
          eta(iep2,j,k) = eta(iep1,j,k)
        enddo
      enddo
    endif
    if(noib .eq. 5) then
      do k=ks,ke
        do j=js,je
          eta(iep1,j,k) = eta(ie  ,j,k)
          eta(iep2,j,k) = eta(iem1,j,k)
        enddo
      enddo
    endif
  endif

  !------------------ j - boundary -----------------
  idtop = myid1 + (myid2 + 1)*nproc1 + myid3*(nproc1*nproc2)
  idbot = myid1 + (myid2 - 1)*nproc1 + myid3*(nproc1*nproc2)
  isendtop=2
  isendbot=1
  irecvtop=1
  irecvbot=2
  if(myid2 .eq. 0) idbot = myid1 + (nproc2 - 1)*nproc1 &
                         + myid3*(nproc1*nproc2)
  if(myid2 .eq. nproc2-1) idtop = myid1 + 0*nproc1 &
                                + myid3*(nproc1*nproc2)
  allocate(recvtop(ism2:iep2,jep1:jep2,ks:ke),recvbot(ism2:iep2,jsm2:jsm1,ks:ke))
  allocate(sendtop(ism2:iep2,jem1:je,ks:ke),sendbot(ism2:iep2,js:jsp1,ks:ke))
  sendtop=eta(ism2:iep2,jem1:je,ks:ke)
  sendbot=eta(ism2:iep2,js:jsp1,ks:ke)
  call message_pass_nodes(idtop, idbot, &
    iep2-ism2+1, 2, ke-ks+1, iep2-ism2+1, 2, ke-ks+1, &
    sendtop, sendbot, recvtop, recvbot, &
    isendtop, isendbot, irecvtop, irecvbot)
  eta(ism2:iep2,jep1:jep2,ks:ke) = recvtop
  eta(ism2:iep2,jsm2:jsm1,ks:ke) = recvbot
  deallocate(recvtop, recvbot, sendtop, sendbot)
  ! Inner j boundary
  if(myid2 .eq. 0) then
    if ( nijb .eq. 1) then
      do i=ism2,iep2
        do k=ks,ke
          eta(i,jsm1,k) = eta(i,js  ,k)
          eta(i,jsm2,k) = eta(i,jsp1,k)
        enddo
      enddo
    endif
    if ( nijb .eq. 2) then
      do i=ism2,iep2
        do k=ks,ke
          eta(i,jsm1,k) = eta(i,js  ,k)
          eta(i,jsm2,k) = eta(i,jsm1,k)
        enddo
      enddo
    endif
    if ( nijb .eq. 3) then
      do i=ism2,iep2
        do k=ks,ke
          eta(i,jsm1,k) = eta(i,js  ,k)
          eta(i,jsm2,k) = eta(i,jsm1,k)
        enddo
      enddo
    endif
    if ( nijb .eq. 5) then
      do i=ism2,iep2
        do k=ks,ke
          eta(i,jsm1,k) = eta(i,js  ,k)
          eta(i,jsm2,k) = eta(i,jsp1,k)
        enddo
      enddo
    endif
  endif
  ! Outer j boundary
  if(myid2 .eq. nproc2-1) then
    if ( nojb .eq. 1) then
      do i=ism2,iep2
        do k=ks,ke
          eta(i,jep1,k) = eta(i,je  ,k)
          eta(i,jep2,k) = eta(i,jem1,k)
        enddo
      enddo
    endif
    if ( nojb .eq. 2) then
      do i=ism2,iep2
        do k=ks,ke
          eta(i,jep1,k) = eta(i,je  ,k)
          eta(i,jep2,k) = eta(i,jep1,k)
        enddo
      enddo
    endif
    if ( nojb .eq. 3) then
      do i=ism2,iep2
        do k=ks,ke
          eta(i,jep1,k) = eta(i,je  ,k)
          eta(i,jep2,k) = eta(i,jep1,k)
        enddo
      enddo
    endif
    if ( nojb .eq. 5) then
      do i=ism2,iep2
        do k=ks,ke
          eta(i,jep1,k) = eta(i,je  ,k)
          eta(i,jep2,k) = eta(i,jem1,k)
        enddo
      enddo
    endif
  endif
  !------------------ k - boundary -----------------
  idtop = myid1 + myid2*nproc1 + (myid3+1)*(nproc1*nproc2)
  idbot = myid1 + myid2*nproc1 + (myid3-1)*(nproc1*nproc2)
  isendtop=2
  isendbot=1
  irecvtop=1
  irecvbot=2
  if(myid3==0) idbot = myid1 + myid2*nproc1 + (nproc3-1)*(nproc1*nproc2)
  if(myid3==nproc3-1) idtop = myid1 + myid2*nproc1 + 0*(nproc1*nproc2)
  allocate(recvtop(ism2:iep2,jsm2:jep2,kep1:kep2),recvbot(ism2:iep2,jsm2:jep2,ksm2:ksm1))
  allocate(sendtop(ism2:iep2,jsm2:jep2,kem1:ke),sendbot(ism2:iep2,jsm2:jep2,ks:ksp1))
  sendtop=eta(ism2:iep2,jsm2:jep2,kem1:ke)
  sendbot=eta(ism2:iep2,jsm2:jep2,ks:ksp1)
  call message_pass_nodes(idtop, idbot, &
    iep2-ism2+1, jep2-jsm2+1, 2, iep2-ism2+1, jep2-jsm2+1, 2, &
    sendtop, sendbot, recvtop, recvbot, &
    isendtop, isendbot, irecvtop, irecvbot)
  eta(ism2:iep2,jsm2:jep2,kep1:kep2) = recvtop
  eta(ism2:iep2,jsm2:jep2,ksm2:ksm1) = recvbot
  deallocate(recvtop, recvbot, sendtop, sendbot)
  ! Inner k boundary
  if(myid3 .eq. 0) then
    if ( nikb .eq. 1) then
      do j=jsm2,jep2
        do i=ism2,iep2
          eta(i,j,ksm1) = eta(i,j,ks  )
          eta(i,j,ksm2) = eta(i,j,ksp1)
        enddo
      enddo
    endif
    if ( nikb .eq. 2) then
      do j=jsm2,jep2
        do i=ism2,iep2
          eta(i,j,ksm1) = eta(i,j,ks  )
          eta(i,j,ksm2) = eta(i,j,ksm1)
        enddo
      enddo
    endif
    if ( nikb .eq. 3) then
      do j=jsm2,jep2
        do i=ism2,iep2
          eta(i,j,ksm1) = eta(i,j,ks  )
          eta(i,j,ksm2) = eta(i,j,ksm1)
        enddo
      enddo
    endif
    if ( nikb .eq. 5) then
      do j=jsm2,jep2
        do i=ism2,iep2
          eta(i,j,ksm1) = eta(i,j,ks  )
          eta(i,j,ksm2) = eta(i,j,ksp1)
        enddo
      enddo
    endif
  endif
  ! Outer k boundary
  if(myid3 .eq. nproc3-1) then
    if ( nokb .eq. 1) then
      do j=jsm2,jep2
        do i=ism2,iep2
          eta(i,j,kep1) = eta(i,j,ke  )
          eta(i,j,kep2) = eta(i,j,kem1)
        enddo
      enddo
    endif
    if ( nokb .eq. 2) then
      do j=jsm2,jep2
        do i=ism2,iep2
          eta(i,j,kep1) = eta(i,j,ke  )
          eta(i,j,kep2) = eta(i,j,kep1)
        enddo
      enddo
    endif
    if ( nokb .eq. 3) then
      do j=jsm2,jep2
        do i=ism2,iep2
          eta(i,j,kep1) = eta(i,j,ke  )
          eta(i,j,kep2) = eta(i,j,kep1)
        enddo
      enddo
    endif
    if ( nokb .eq. 5) then
      do j=jsm2,jep2
        do i=ism2,iep2
          eta(i,j,kep1) = eta(i,j,ke  )
          eta(i,j,kep2) = eta(i,j,kem1)
        enddo
      enddo
    endif
  endif

end subroutine bvaleta

end module ModBval
