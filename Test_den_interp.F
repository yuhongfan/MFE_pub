subroutine x1zc3d_den ( qty, jbeg, kbeg, jend, kend, &
                        ql, qr)
  use ModPar
  use ModGrid
  use ModBoundary
  use ModUserSetup, only: PDMB,fac4th
  implicit none

  real, intent(in) :: qty(in,jn,kn)
  real, intent(out) :: ql(in,jn,kn), qr(in,jn,kn)
  integer, intent(in) :: jbeg, kbeg, jend, kend
  real*8:: f_itp,maxf,minf,f,df0,df1,df2,s0,s1,s2,q0,q1,df_left,df_right

  integer :: idtop, idbot
  integer :: isendtop, isendbot, irecvtop, irecvbot
  real, allocatable :: recvtop(:,:,:), recvbot(:,:,:)
  real, allocatable :: sendtop(:,:,:), sendbot(:,:,:)
  

  integer :: i, j, k
  real :: dqim, dqip, dqii(in,jn,kn)
!
!
!-------------- 4th order interpolated values with PDM limiters---------------
!

         do k=kbeg,kend
           do j=jbeg,jend
             do i=is,iep1
               f_itp = sum(qty(i-2:i+1,j,k)*vol(i-2:i+1,j,k)*fac4th)
                f_itp = f_itp/vol(i,j,k,dir1)
                maxf = max(qty(i-1,j,k),qty(i,j,k))
                minf = min(qty(i-1,j,k),qty(i,j,k))

                f = max(minf,min(f_itp,maxf))

                df0 = PDMB*(qty(i-1,j,k)-qty(i-2,j,k))
                df1 = PDMB*(qty(i,j,k)-qty(i-1,j,k))
                df2 = PDMB*(qty(i+1,j,k)-qty(i,j,k))

                s0 = sign(1.0_8,df0)
                s1 = sign(1.0_8,df1)
                s2 = sign(1.0_8,df2)

                df0 = abs(df0)
                df1 = abs(df1)
                df2 = abs(df2)

                q0 = abs(s0+s1)
                q1 = abs(s1+s2)

                df_left = f - qty(i-1,j,k)
                df_right = qty(i,j,k) - f

                ql(i,j,k) = f - s1*max(0.0_8,abs(df_left) - q0*df0)
                qr(i,j,k) = f + s1*max(0.0_8,abs(df_right) - q1*df2)
             enddo
           enddo
         enddo

        
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


end subroutine x1zc3d_den

subroutine x2zc3d_den ( qty, kbeg, ibeg, kend, iend, &
                        ql, qr)
  use ModPar
  use ModGrid
  use ModBoundary
  use ModUserSetup, only: PDMB,fac4th
  implicit none

  real, intent(in) :: qty(in,jn,kn)
  real, intent(out) :: ql(in,jn,kn), qr(in,jn,kn)
  integer, intent(in) :: kbeg, ibeg, kend, iend

  integer :: idtop, idbot
  integer :: isendtop, isendbot, irecvtop, irecvbot
  real, allocatable :: recvtop(:,:,:), recvbot(:,:,:)
  real, allocatable :: sendtop(:,:,:), sendbot(:,:,:)
  real*8:: f_itp,maxf,minf,f,df0,df1,df2,s0,s1,s2,q0,q1,df_left,df_right

  integer :: i, j, k
  real :: dqjm, dqjp, dqji(in,jn,kn)
!
!
!-------------- 2nd order interpolated values ---------------
!
!  1.  Evaluate monotonised, slope-limited difference in "q" across the zone.
!
         do k=kbeg,kend
           do j=js,jep1
             do i=ibeg,iend
               f_itp = sum(qty(i,j-2:j+1,k)*vol(i,j-2:j+1,k)*fac4th)
                f_itp = f_itp/vol(i,j,k,dir2)
                maxf = max(qty(i,j-1,k),qty(i,j,k))
                minf = min(qty(i,j-1,k),qty(i,j,k))

                f = max(minf,min(f_itp,maxf))

                df0 = PDMB*(qty(i,j-1,k)-qty(i,j-2,k))
                df1 = PDMB*(qty(i,j,k)-qty(i,j-1,k))
                df2 = PDMB*(qty(i,j+1,k)-qty(i,j,k))

                s0 = sign(1.0_8,df0)
                s1 = sign(1.0_8,df1)
                s2 = sign(1.0_8,df2)

                df0 = abs(df0)
                df1 = abs(df1)
                df2 = abs(df2)

                q0 = abs(s0+s1)
                q1 = abs(s1+s2)

                df_left = f - qty(i,j-1,k)
                df_right = qty(i,j,k) - f

                ql(i,j,k) = f - s1*max(0.0_8,abs(df_left) - q0*df0)
                qr(i,j,k) = f + s1*max(0.0_8,abs(df_right) - q1*df2)
             enddo
           enddo
         enddo



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


end subroutine x2zc3d_den

subroutine x3zc3d_den ( qty, ibeg, jbeg, iend, jend, &
                        ql, qr)
  use ModPar
  use ModGrid
  use ModBoundary
  use ModUserSetup, only: PDMB,fac4th
  implicit none

  real, intent(in) :: qty(in,jn,kn)
  real, intent(out) :: ql(in,jn,kn), qr(in,jn,kn)
  integer, intent(in) :: ibeg, jbeg, iend, jend
  real*8:: f_itp,maxf,minf,f,df0,df1,df2,s0,s1,s2,q0,q1,df_left,df_right

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
	 do k=ks,kep1
           do j=jbeg,jend
             do i=ibeg,iend
               f_itp = sum(qty(i,j,k-2:k+1)*vol(i,j,k-2:k+1)*fac4th)
                f_itp = f_itp/vol(i,j,k,dir3)
                maxf = max(qty(i,j,k),qty(i,j,k-1))
                minf = min(qty(i,j,k),qty(i,j,k-1))

                f = max(minf,min(f_itp,maxf))

                df0 = PDMB*(qty(i,j,k-1)-qty(i,j,k-2))
                df1 = PDMB*(qty(i,j,k)-qty(i,j,k-1))
                df2 = PDMB*(qty(i,j,k+1)-qty(i,j,k))

                s0 = sign(1.0_8,df0)
                s1 = sign(1.0_8,df1)
                s2 = sign(1.0_8,df2)

                df0 = abs(df0)
                df1 = abs(df1)
                df2 = abs(df2)

                q0 = abs(s0+s1)
                q1 = abs(s1+s2)

                df_left = f - qty(i,j,k-1)
                df_right = qty(i,j,k) - f

                ql(i,j,k) = f - s1*max(0.0_8,abs(df_left) - q0*df0)
                qr(i,j,k) = f + s1*max(0.0_8,abs(df_right) - q1*df2)
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

end subroutine x3zc3d_den