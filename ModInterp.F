module ModInterp

  implicit none

  private
  public :: x1tvd
  public :: x1tvd_a
  public :: x2tvd
  public :: x2tvd_a
  public :: x3tvd
  public :: x3tvd_a
  public :: x1zc3d_den
  public :: x2zc3d_den
  public :: x3zc3d_den
  public :: x1del_den
  public :: x2del_den
  public :: x3del_den
  public :: findpos

contains

subroutine x1tvd (qint,ql,qr,dq,qlr)
  use ModPar
  use ModGrid
  implicit none

  real, intent(in) :: qint(ijkn)
  real, intent(out) :: dq(ijkn),ql(ijkn),qr(ijkn),qlr(ijkn)

  integer :: i, j, k
  real :: dqi(in), dqim, dqip
!
!-----------------------------------------------------------
!
         do i=ism1,iep1
           dqim     = ( qint(i  ) - qint(i-1) ) * dx1bi(i  )
           dqip     = ( qint(i+1) - qint(i  ) ) * dx1bi(i+1)
           dqi(i  ) = sign (0.5,dqim) * max ( 0., &
                      min ( abs(dqim), sign(1.,dqim)*dqip ) )
         enddo

         do i=is,iep1
           dqim = qint(i-1) + dx1a(i-1) * dqi(i-1)
           dqip = qint(i  ) - dx1a(i  ) * dqi(i  )
           ql(i) = dqim
           qr(i) = dqip
           dq(i) = dqip - dqim
           qlr(i) = 0.5 * (dqip + dqim)
         enddo

end subroutine x1tvd

subroutine x1tvd_a (qint,ql,qr,dq,qlr)
  use ModPar
  use ModGrid
  implicit none

  real, intent(in) :: qint(ijkn)
  real, intent(out) :: dq(ijkn),ql(ijkn),qr(ijkn),qlr(ijkn)

  integer :: i, j, k
  real :: dqi(in), dqim, dqip
!
!-----------------------------------------------------------
!
         do i=ism1,iep2
           dqim     = ( qint(i  ) - qint(i-1) ) * dx1ai(i-1)
           dqip     = ( qint(i+1) - qint(i  ) ) * dx1ai(i  )
           dqi(i  ) = sign (0.5,dqim) * max ( 0., &
                      min ( abs(dqim), sign(1.,dqim)*dqip ) )
         enddo

         do i=ism1,iep1
           dqim = qint(i  ) + dx1a(i  ) * dqi(i  )
           dqip = qint(i+1) - dx1a(i  ) * dqi(i+1)
           ql(i) = dqim
           qr(i) = dqip
           dq(i) = dqip - dqim
           qlr(i) = 0.5 * (dqip + dqim)
         enddo

end subroutine x1tvd_a

subroutine x2tvd (qint,ql,qr,dq,qlr)
  use ModPar
  use ModGrid
  implicit none

  real, intent(in) :: qint(ijkn)
  real, intent(out) :: dq(ijkn),ql(ijkn),qr(ijkn),qlr(ijkn)

  integer :: i, j, k
  real :: dqj(jn), dqjm, dqjp
!
!-----------------------------------------------------------
!
         do j=jsm1,jep1
           dqjm     = ( qint(j  ) - qint(j-1) ) * dx2bi(j  )
           dqjp     = ( qint(j+1) - qint(j  ) ) * dx2bi(j+1)
           dqj(j  ) = sign (0.5,dqjm) * max ( 0., &
                      min ( abs(dqjm), sign(1.,dqjm)*dqjp ) )
         enddo

         do j=js,jep1
           dqjm = qint(j-1) + dx2a(j-1) * dqj(j-1)
           dqjp = qint(j  ) - dx2a(j  ) * dqj(j  )
           ql(j) = dqjm
           qr(j) = dqjp
           dq(j) = dqjp - dqjm
           qlr(j) = 0.5 * (dqjp + dqjm)
         enddo

end subroutine x2tvd

subroutine x2tvd_a (qint,ql,qr,dq,qlr)
  use ModPar
  use ModGrid
  implicit none

  real, intent(in) :: qint(ijkn)
  real, intent(out) :: dq(ijkn),ql(ijkn),qr(ijkn),qlr(ijkn)

  integer :: i, j, k
  real :: dqj(jn), dqjm, dqjp
!
!-----------------------------------------------------------
!
         do j=jsm1,jep2
           dqjm     = ( qint(j  ) - qint(j-1) ) * dx2ai(j-1)
           dqjp     = ( qint(j+1) - qint(j  ) ) * dx2ai(j  )
           dqj(j  ) = sign (0.5,dqjm) * max ( 0., &
                      min ( abs(dqjm), sign(1.,dqjm)*dqjp ) )
         enddo

         do j=jsm1,jep1
           dqjm = qint(j  ) + dx2a(j  ) * dqj(j  )
           dqjp = qint(j+1) - dx2a(j  ) * dqj(j+1)
           ql(j) = dqjm
           qr(j) = dqjp
           dq(j) = dqjp - dqjm
           qlr(j) = 0.5 * (dqjp + dqjm)
         enddo

end subroutine x2tvd_a

subroutine x3tvd (qint,ql,qr,dq,qlr)
  use ModPar
  use ModGrid
  implicit none

  real, intent(in) :: qint(ijkn)
  real, intent(out) :: dq(ijkn),ql(ijkn),qr(ijkn),qlr(ijkn)

  integer :: i, j, k
  real :: dqk(kn), dqkm, dqkp
!
!-----------------------------------------------------------
!
         do k=ksm1,kep1
           dqkm     = ( qint(k  ) - qint(k-1) ) * dx3bi(k  )
           dqkp     = ( qint(k+1) - qint(k  ) ) * dx3bi(k+1)
           dqk(k  ) = sign (0.5,dqkm) * max ( 0., &
                      min ( abs(dqkm), sign(1.,dqkm)*dqkp ) )
         enddo

         do k=ks,kep1
           dqkm = qint(k-1) + dx3a(k-1) * dqk(k-1)
           dqkp = qint(k  ) - dx3a(k  ) * dqk(k  )
           ql(k) = dqkm
           qr(k) = dqkp
           dq(k) = dqkp - dqkm
           qlr(k) = 0.5 * (dqkp + dqkm)
         enddo

end subroutine x3tvd

subroutine x3tvd_a (qint,ql,qr,dq,qlr)
  use ModPar
  use ModGrid
  implicit none

  real, intent(in) :: qint(ijkn)
  real, intent(out) :: dq(ijkn),ql(ijkn),qr(ijkn),qlr(ijkn)

  integer :: i, j, k
  real :: dqk(kn), dqkm, dqkp
!
!-----------------------------------------------------------
!
         do k=ksm1,kep2
           dqkm     = ( qint(k  ) - qint(k-1) ) * dx3ai(k-1)
           dqkp     = ( qint(k+1) - qint(k  ) ) * dx3ai(k  )
           dqk(k  ) = sign (0.5,dqkm) * max ( 0., &
                      min ( abs(dqkm), sign(1.,dqkm)*dqkp ) )
         enddo

         do k=ksm1,kep1
           dqkm = qint(k  ) + dx3a(k  ) * dqk(k  )
           dqkp = qint(k+1) - dx3a(k  ) * dqk(k+1)
           ql(k) = dqkm
           qr(k) = dqkp
           dq(k) = dqkp - dqkm
           qlr(k) = 0.5 * (dqkp + dqkm)
         enddo

end subroutine x3tvd_a

subroutine x1zc3d_den ( qty, jbeg, kbeg, jend, kend, &
                        ql, qr)
  use ModPar
  use ModGrid
  use ModBoundary
  implicit none

  real, intent(in) :: qty(in,jn,kn)
  real, intent(out) :: ql(in,jn,kn), qr(in,jn,kn)
  integer, intent(in) :: jbeg, kbeg, jend, kend

  integer :: idtop, idbot
  integer :: isendtop, isendbot, irecvtop, irecvbot
  real, allocatable :: recvtop(:,:,:), recvbot(:,:,:)
  real, allocatable :: sendtop(:,:,:), sendbot(:,:,:)

  integer :: i, j, k
  real :: dqim, dqip, dqii(in,jn,kn)
!
!
!-------------- 2nd order interpolated values ---------------
!
!  1.  Evaluate monotonised, slope-limited difference in "q" across the zone.
!
	 do k=kbeg,kend
           do j=jbeg,jend
             do i=ism1,iep1
               dqim        = (qty(i  ,j,k) - qty(i-1,j,k)) * dx1bi(i  )
               dqip        = (qty(i+1,j,k) - qty(i  ,j,k)) * dx1bi(i+1)
               dqii(i,j,k) = sign (0.5,dqim) * max ( 0., &
                           min ( abs(dqim), sign(1.,dqim)*dqip ) )
             enddo
           enddo
         enddo
!
!        Communicate between processors
!
         idtop=myid1+1+nproc1*myid2+(nproc1*nproc2)*myid3
         idbot=myid1-1+nproc1*myid2+(nproc1*nproc2)*myid3
         isendtop=2
         isendbot=1
         irecvtop=1
         irecvbot=2
         if(myid1 .eq. 0) then
           idbot=nproc1-1+nproc1*myid2+(nproc1*nproc2)*myid3
         endif
         if(myid1 .eq. nproc1-1) then
           idtop=0+nproc1*myid2+(nproc1*nproc2)*myid3
         endif
         allocate(recvbot(ism2:ism2,jbeg:jend,kbeg:kend), &
                  recvtop(jep2:jep2,jbeg:jend,kbeg:kend))
         allocate(sendtop(iem1:iem1,jbeg:jend,kbeg:kend), &
                  sendbot(isp1:isp1,jbeg:jend,kbeg:kend))
         sendtop=dqii(iem1:iem1,jbeg:jend,kbeg:kend)
         sendbot=dqii(isp1:isp1,jbeg:jend,kbeg:kend)
         call message_pass_nodes(idtop, idbot,&
           1, jend-jbeg+1, kend-kbeg+1, 1, jend-jbeg+1, kend-kbeg+1, &
           sendtop,sendbot,recvtop,recvbot, &
           isendtop,isendbot,irecvtop,irecvbot)
         dqii(ism2:ism2,jbeg:jend,kbeg:kend)=recvbot
         dqii(iep2:iep2,jbeg:jend,kbeg:kend)=recvtop
         deallocate(recvtop, recvbot, sendtop, sendbot)
!
!        Inner i boundary
!
         if(myid1 .eq. 0) then

           if (niib .ne. 4) then
             do k=kbeg,kend
               do j=jbeg,jend
                 dqii(ism2,j,k)=0.
               enddo
             enddo
           endif
           if (niib .eq. 1) then
             do k=kbeg,kend
               do j=jbeg,jend
                 dqii(ism2,j,k)=-dqii(isp1,j,k)
               enddo
             enddo
           endif
           if (niib .eq. 5) then
             do k=kbeg,kend
               do j=jbeg,jend
                 dqii(ism2,j,k)=-dqii(isp1,j,k)
               enddo
             enddo
           endif

         endif
!
!        Outer i boundary
!
         if(myid1 .eq. nproc1-1) then

           if (noib .ne. 4) then
             do k=kbeg,kend
               do j=jbeg,jend
                 dqii(iep2,j,k)=0.
               enddo
             enddo
           endif
           if (noib .eq. 1) then
             do k=kbeg,kend
               do j=jbeg,jend
                 dqii(iep2,j,k)=-dqii(iem1,j,k)
               enddo
             enddo
           endif
           if (noib .eq. 5) then
             do k=kbeg,kend
               do j=jbeg,jend
                 dqii(iep2,j,k)=-dqii(iem1,j,k)
               enddo
             enddo
           endif

         endif
!
!  2.  compute the l and r interface values.
!
         do k=kbeg,kend
           do j=jbeg,jend
             do i=ism1,iep2
               dqim = qty(i-1,j,k) + dx1a(i-1)*dqii(i-1,j,k)
               dqip = qty(i  ,j,k) - dx1a(i  )*dqii(i  ,j,k)
               ql(i,j,k) = dqim
               qr(i,j,k) = dqip
             enddo
           enddo
         enddo

end subroutine x1zc3d_den

subroutine x2zc3d_den ( qty, kbeg, ibeg, kend, iend, &
                        ql, qr)
  use ModPar
  use ModGrid
  use ModBoundary
  implicit none

  real, intent(in) :: qty(in,jn,kn)
  real, intent(out) :: ql(in,jn,kn), qr(in,jn,kn)
  integer, intent(in) :: kbeg, ibeg, kend, iend

  integer :: idtop, idbot
  integer :: isendtop, isendbot, irecvtop, irecvbot
  real, allocatable :: recvtop(:,:,:), recvbot(:,:,:)
  real, allocatable :: sendtop(:,:,:), sendbot(:,:,:)

  integer :: i, j, k
  real :: dqjm, dqjp, dqji(in,jn,kn)
!
!
!-------------- 2nd order interpolated values ---------------
!
!  1.  Evaluate monotonised, slope-limited difference in "q" across the zone.
!
         do k=kbeg,kend
           do j=jsm1,jep1
             do i=ibeg,iend
               dqjm        = (qty(i,j  ,k) - qty(i,j-1,k)) * dx2bi(j  )
               dqjp        = (qty(i,j+1,k) - qty(i,j  ,k)) * dx2bi(j+1)
               dqji(i,j,k) = sign (0.5,dqjm) * max ( 0., &
                           min ( abs(dqjm), sign(1.,dqjm)*dqjp ) )
             enddo
           enddo
         enddo
!
!        Communicate between processors
!
         idtop=myid1+(myid2+1)*nproc1+myid3*(nproc1*nproc2)
         idbot=myid1+(myid2-1)*nproc1+myid3*(nproc1*nproc2)
         isendtop=2
         isendbot=1
         irecvtop=1
         irecvbot=2
         if(myid2 .eq. 0) then
           idbot=myid1+(nproc2-1)*nproc1+myid3*(nproc1*nproc2)
         endif
         if(myid2 .eq. nproc2-1) then
           idtop=myid1+0*nproc1+myid3*(nproc1*nproc2)
         endif
         allocate(recvbot(ibeg:iend,jsm2:jsm2,kbeg:kend), &
                  recvtop(ibeg:iend,jep2:jep2,kbeg:kend))
         allocate(sendtop(ibeg:iend,jem1:jem1,kbeg:kend), &
                  sendbot(ibeg:iend,jsp1:jsp1,kbeg:kend))
         sendtop=dqji(ibeg:iend,jem1:jem1,kbeg:kend)
         sendbot=dqji(ibeg:iend,jsp1:jsp1,kbeg:kend)
         call message_pass_nodes(idtop, idbot,&
           iend-ibeg+1, 1, kend-kbeg+1, iend-ibeg+1, 1, kend-kbeg+1, &
           sendtop,sendbot,recvtop,recvbot, &
           isendtop,isendbot,irecvtop,irecvbot)           
         dqji(ibeg:iend,jsm2:jsm2,kbeg:kend)=recvbot
         dqji(ibeg:iend,jep2:jep2,kbeg:kend)=recvtop
         deallocate(recvtop, recvbot, sendtop, sendbot)
!
!        Inner j boundary
!
         if(myid2 .eq. 0) then

           if (nijb .ne. 4) then 
             do k=kbeg,kend
               do i=ibeg,iend
                 dqji(i,jsm2,k)=0.
               enddo
             enddo
           endif
           if (nijb .eq. 1) then
             do k=kbeg,kend
               do i=ibeg,iend
                 dqji(i,jsm2,k)=-dqji(i,jsp1,k)
               enddo
             enddo
           endif
           if (nijb .eq. 5) then
             do k=kbeg,kend
               do i=ibeg,iend
                 dqji(i,jsm2,k)=-dqji(i,jsp1,k)
               enddo
             enddo
           endif

         endif
!
!        Outer j boundary
!
         if(myid2 .eq. nproc2-1) then

           if (nojb .ne. 4) then
             do k=kbeg,kend
               do i=ibeg,iend
                 dqji(i,jep2,k)=0.
               enddo
             enddo
           endif
           if (nojb .eq. 1) then
             do k=kbeg,kend
               do i=ibeg,iend
                 dqji(i,jep2,k)=-dqji(i,jem1,k)
               enddo
             enddo
           endif
           if (nojb .eq. 5) then
             do k=kbeg,kend
               do i=ibeg,iend
                 dqji(i,jep2,k)=-dqji(i,jem1,k)
               enddo
             enddo
           endif

         endif
!
!  2.  Choose time averaged, upwinded interface value.
!
         do k=kbeg,kend
           do j=jsm1,jep2
             do i=ibeg,iend
               dqjm = qty(i,j-1,k) + dx2a(j-1)*dqji(i,j-1,k)
               dqjp = qty(i,j  ,k) - dx2a(j  )*dqji(i,j  ,k)
               ql(i,j,k) = dqjm
               qr(i,j,k) = dqjp
             enddo
           enddo
         enddo

end subroutine x2zc3d_den

subroutine x3zc3d_den ( qty, ibeg, jbeg, iend, jend, &
                        ql, qr)
  use ModPar
  use ModGrid
  use ModBoundary
  implicit none

  real, intent(in) :: qty(in,jn,kn)
  real, intent(out) :: ql(in,jn,kn), qr(in,jn,kn)
  integer, intent(in) :: ibeg, jbeg, iend, jend

  integer :: idtop, idbot
  integer :: isendtop, isendbot, irecvtop, irecvbot
  real, allocatable :: recvtop(:,:,:), recvbot(:,:,:)
  real, allocatable :: sendtop(:,:,:), sendbot(:,:,:)

  integer :: i, j, k
  real :: dqkm, dqkp, dqki(in,jn,kn)
!
!
!-------------- 2nd order interpolated values ---------------
!
!  1.  Evaluate monotonised, slope-limited difference in "q" across the zone.
!
	 do k=ksm1,kep1
           do j=jbeg,jend
             do i=ibeg,iend
               dqkm        = (qty(i,j,k  ) - qty(i,j,k-1)) * dx3bi(k  )
               dqkp        = (qty(i,j,k+1) - qty(i,j,k  )) * dx3bi(k+1)
               dqki(i,j,k) = sign (0.5,dqkm) * max ( 0., &
                           min ( abs(dqkm), sign(1.,dqkm)*dqkp ) )
             enddo
           enddo
         enddo
!
!        Communicate between processors
!
         idtop=myid1+myid2*nproc1+(myid3+1)*(nproc1*nproc2)
         idbot=myid1+myid2*nproc1+(myid3-1)*(nproc1*nproc2)
         isendtop=2
         isendbot=1
         irecvtop=1
         irecvbot=2
         if(myid3 .eq. 0) then
           idbot=myid1+myid2*nproc1+(nproc3-1)*(nproc1*nproc2)
         endif
         if(myid3 .eq. nproc3-1) then
           idtop=myid1+myid2*nproc1+0*(nproc1*nproc2)
         endif
         allocate(recvbot(ibeg:iend,jbeg:jend,ksm2:ksm2), &
                  recvtop(ibeg:iend,jbeg:jend,kep2:kep2))
         allocate(sendtop(ibeg:iend,jbeg:jend,kem1:kem1), &
                  sendbot(ibeg:iend,jbeg:jend,ksp1:ksp1))
         sendtop=dqki(ibeg:iend,jbeg:jend,kem1:kem1)
         sendbot=dqki(ibeg:iend,jbeg:jend,ksp1:ksp1)
         call message_pass_nodes(idtop, idbot,&
           iend-ibeg+1, jend-jbeg+1, 1, iend-ibeg+1, jend-jbeg+1, 1, &
           sendtop,sendbot,recvtop,recvbot, &
           isendtop,isendbot,irecvtop,irecvbot)
         dqki(ibeg:iend,jbeg:jend,ksm2:ksm2)=recvbot
         dqki(ibeg:iend,jbeg:jend,kep2:kep2)=recvtop
         deallocate(recvtop, recvbot, sendtop, sendbot)
!
!        Inner k boundary
!
         if(myid3 .eq. 0) then

         if (nikb .ne. 4) then 
           do j=jbeg,jend
             do i=ibeg,iend
               dqki(i,j,ksm2)=0.
             enddo
           enddo
         endif
         if (nikb .eq. 1) then
           do j=jbeg,jend
             do i=ibeg,iend
               dqki(i,j,ksm2)=-dqki(i,j,ksp1)
             enddo
           enddo 
         endif
         if (nikb .eq. 5) then
           do j=jbeg,jend
             do i=ibeg,iend
               dqki(i,j,ksm2)=-dqki(i,j,ksp1)
             enddo
           enddo 
         endif

         endif
!
!        Outer k boundary
!
         if(myid3 .eq. nproc3-1) then

         if (nokb .ne. 4) then
           do j=jbeg,jend
             do i=ibeg,iend
               dqki(i,j,kep2)=0.
             enddo
           enddo
         endif
         if (nokb .eq. 1) then
           do j=jbeg,jend
             do i=ibeg,iend
               dqki(i,j,kep2)=-dqki(i,j,kem1)
             enddo
           enddo
         endif
         if (nokb .eq. 5) then
           do j=jbeg,jend
             do i=ibeg,iend
               dqki(i,j,kep2)=-dqki(i,j,kem1)
             enddo
           enddo
         endif

         endif
!
!  2.  Choose time averaged, upwinded interface value.
!
         do k=ksm1,kep2
           do j=jbeg,jend
             do i=ibeg,iend
               dqkm = qty(i,j,k-1) + dx3a(k-1)*dqki(i,j,k-1)
               dqkp = qty(i,j,k  ) - dx3a(k  )*dqki(i,j,k  )
               ql(i,j,k) = dqkm
               qr(i,j,k) = dqkp
             enddo
           enddo
         enddo

end subroutine x3zc3d_den

subroutine x1del_den(bint, dbint)
  use ModPar
  use ModGrid
  implicit none

  real, intent(in), dimension(in) :: bint
  real, intent(out), dimension(in) :: dbint

  integer :: i
  real :: dqim, dqip, q1, q2
  real, dimension(in) :: dqi
!
!-------------- 2nd order interpolated values ---------------
!
!
!  1.  Evaluate monotonised, slope-limited difference in "q" across the zone.
!
	 do i=ism1,iep1
           dqim     = ( bint(i  ) - bint(i-1) ) * dx1bi(i  )
           dqip     = ( bint(i+1) - bint(i  ) ) * dx1bi(i+1)
           dqi(i  ) = sign (0.5,dqim) * max ( 0., &
                     min ( abs(dqim), sign(1.,dqim)*dqip ) )
         enddo
!
!  2.  Perform an upwinded interpolation of "q" to the time-centred
!      bases of the characteristics.
!
         do i=is,iep1
           q1     = bint(i-1) + dx1a(i-1) * dqi(i-1)
           q2     = bint(i  ) - dx1a(i  ) * dqi(i  )
           dbint(i) = q2 - q1
         enddo

end subroutine x1del_den

subroutine x2del_den(bint, dbint)
  use ModPar
  use ModGrid
  implicit none

  real, intent(in), dimension(jn) :: bint
  real, intent(out), dimension(jn) :: dbint

  integer :: j
  real :: dqjm, dqjp, q1, q2
  real, dimension(jn) :: dqi
!
!-------------- 2nd order interpolated values ---------------
!
!
!  1.  Evaluate monotonised, slope-limited difference in "q" across the zone.
!
	 do j=jsm1,jep1
           dqjm     = ( bint(j  ) - bint(j-1) ) * dx2bi(j  )
           dqjp     = ( bint(j+1) - bint(j  ) ) * dx2bi(j+1)
           dqi(j  ) = sign (0.5,dqjm) * max ( 0., &
                     min ( abs(dqjm), sign(1.,dqjm)*dqjp ) )
         enddo
!
!  2.  Perform an upwinded interpolation of "q" to the time-centred
!      bases of the characteristics.
!
         do j=js,jep1
           q1     = bint(j-1) + dx2a(j-1) * dqi(j-1)
           q2     = bint(j  ) - dx2a(j  ) * dqi(j  )
           dbint(j) = q2 - q1
         enddo

end subroutine x2del_den

subroutine x3del_den(bint, dbint)
  use ModPar
  use ModGrid
  implicit none

  real, intent(in), dimension(kn) :: bint
  real, intent(out), dimension(kn) :: dbint

  integer :: k
  real :: dqkm, dqkp, q1, q2
  real, dimension(kn) :: dqi
!
!-------------- 2nd order interpolated values ---------------
!
!
!  1.  Evaluate monotonised, slope-limited difference in "q" across the zone.
!
	 do k=ksm1,kep1
           dqkm     = ( bint(k  ) - bint(k-1) ) * dx3bi(k  )
           dqkp     = ( bint(k+1) - bint(k  ) ) * dx3bi(k+1)
           dqi(k  ) = sign (0.5,dqkm) * max ( 0., &
                     min ( abs(dqkm), sign(1.,dqkm)*dqkp ) )
         enddo
!
!  2.  Perform an upwinded interpolation of "q" to the time-centred
!      bases of the characteristics.
!
         do k=ks,kep1
           q1     = bint(k-1) + dx3a(k-1) * dqi(k-1)
           q2     = bint(k  ) - dx3a(k  ) * dqi(k  )
           dbint(k) = q2 - q1
         enddo

end subroutine x3del_den

subroutine findpos(nttot,tsub_arr,itsub_arr,itcount_arr, &
  ifilesub_arr,ctime,itsub,itcount,ifilesub)
  implicit none

  include "mpif.h"

  integer, intent(in) :: nttot
  real, intent(in), dimension(nttot) :: tsub_arr
  integer, intent(in), dimension(nttot) :: itsub_arr,itcount_arr, &
    ifilesub_arr
  real, intent(in) :: ctime
  integer, intent(out) :: itsub, itcount, ifilesub

  integer k, klo, khi
  integer ierr

      klo=1
      khi=nttot

      if((ctime.lt.tsub_arr(1)).or.(ctime.ge.tsub_arr(nttot))) then
        write(6,*) 'ctime out of range of timetable'
        write(6,*) 'tsub_arr(1),ctime,tsub_arr(nttot)'
        write(6,*) tsub_arr(1),ctime,tsub_arr(nttot)
        call MPI_ABORT(MPI_COMM_WORLD, 1,ierr)
      endif

      do while ((khi-klo).gt.1)
        k = (khi+klo)/2
        if(tsub_arr(k).gt.ctime) then
          khi = k
        else
          klo = k
        endif
      enddo

      itsub=itsub_arr(khi)
      itcount=itcount_arr(khi)
      ifilesub=ifilesub_arr(khi)

end subroutine findpos

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

end module ModInterp
