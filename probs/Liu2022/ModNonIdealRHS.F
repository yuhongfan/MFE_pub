module ModNonIdealRHS

#ifdef THCONDUCT
  implicit none

  private

  public :: derivs_cond
  public :: derivs_radloss
  public :: derivs_heating

contains

subroutine derivs_cond
  use ModPar
  use ModPhysics
  use ModGrid
  use ModField
  use ModDel,             ONLY: delei, delq
  use ModWork
  use ModSundry
  use ModConductivity
  use ModInterp

  implicit none

  real, parameter :: crit_aniso=1.D3, tempmod=5.D5
  integer :: i, j, k
  real :: dr1, dr2, dr3, drmin, dfactor, &
          aval, fsatup, fsatdn, fsatq, tempcond, ovtau, &
          flx1, flx2, flx3, qspitzer, cmax, rr, qq
  real, dimension(in,jn,kn) :: bdelt, bb, qlim, qflx1, qflx2, qflx3
  real, dimension(ijkn) :: flux, fsat, qint, ql, qr, qlr, dq
!
!
!-----------------------------------------------------------------------
!
       aval=1.D0/8.D0/sqrt(2.D0*gamma)*(gamma-1.D0)*sqrt(mpovme) &
           *aval0
!
!      compute grad T components
!
       do k=ksm1,kep2
         do j=jsm1,jep2
           do i=ism1,iep2
             scratch1(i,j,k) = &
               (tfield(i,j,k)-tfield(i-1,j,k))*dx1bi(i)
             scratch2(i,j,k) = (tfield(i,j,k) &
               -tfield(i,j-1,k))*g2bi(i)*dx2bi(j)
             scratch3(i,j,k) = (tfield(i,j,k) &
               -tfield(i,j,k-1))*g31bi(i)*g32bi(j)*dx3bi(k)
           enddo
         enddo
       enddo
!
       do k=ksm1,kep1
         do j=jsm1,jep1
           do i=ism1,iep1
!
!      compute b dot del T
!
             bdelt(i,j,k) = 0.5D0*( b1(i,j,k) &
               *scratch1(i,j,k) &
               +b1(i+1,j,k)*scratch1(i+1,j,k) ) &
               + 0.5D0*( b2(i,j,k)*scratch2(i,j,k) &
               +b2(i,j+1,k)*scratch2(i,j+1,k) ) &
               + 0.5D0*( b3(i,j,k)*scratch3(i,j,k) &
               +b3(i,j,k+1)*scratch3(i,j,k+1) )
!
!      compute b*b strength
!
             bb(i,j,k) = (b1(i,j,k) &
               +b1(i+1,j,k))**2 + (b2(i,j,k) &
               +b2(i,j+1,k))**2 + (b3(i,j,k) &
               +b3(i,j,k+1))**2
             bb(i,j,k)=0.25D0*bb(i,j,k)
!
!      compute diffcond
!
             tempcond=max(tfield(i,j,k),tempmod/unit_temp)
             diffcond(i,j,k)=kappa0 &
               *tempcond**(5.D0/2.D0)
           enddo
         enddo
       enddo
!
!      compute flux limiting ratio
!
       do k=ks,ke
         do j=js,je
           do i=is,iep1
             flx1=b1(i,j,k) &
               *0.5D0*(bdelt(i,j,k)+bdelt(i-1,j,k)) &
               /(0.5D0*(bb(i,j,k)+bb(i-1,j,k))+tiny)
             flux(i)=-0.5D0*(diffcond(i,j,k)+diffcond(i-1,j,k)) &
               *flx1
             qflx1(i,j,k)=b1(i,j,k) &
               /sqrt(0.5D0*(bb(i,j,k)+bb(i-1,j,k))+tiny) &
               *0.5D0*(q(i-1,j,k)+q(i,j,k))
             fsatup=aval &
               *sqrt(gamma*(gamma-1.D0)*eint(i,j,k)/d(i,j,k)) &
               *eint(i,j,k)
             fsatdn=aval &
               *sqrt(gamma*(gamma-1.D0)*eint(i-1,j,k)/d(i-1,j,k)) &
               *eint(i-1,j,k)
             fsatq=0.5D0*(fsatup+fsatdn) &
               -sign(0.5D0,qflx1(i,j,k))*(fsatup-fsatdn)
             qflx1(i,j,k)=sign(1.D0,qflx1(i,j,k))*min(abs(fsatq),abs(qflx1(i,j,k)))
             fsat(i)=0.5D0*(fsatup+fsatdn) &
               -sign(0.5D0,flux(i))*(fsatup-fsatdn)
           enddo
           do i=is,ie
             qlim(i,j,k)=max(1.D0,abs(flux(i))/fsat(i), &
               abs(flux(i+1))/fsat(i+1))
           enddo
         enddo
       enddo
       do k=ks,ke
         do i=is,ie
           do j=js,jep1
             flx2=b2(i,j,k) &
               *0.5D0*(bdelt(i,j,k)+bdelt(i,j-1,k)) &
               /(0.5D0*(bb(i,j,k)+bb(i,j-1,k))+tiny)
             flux(j)=-0.5D0*(diffcond(i,j,k)+diffcond(i,j-1,k)) &
               *flx2
             qflx2(i,j,k)=b2(i,j,k) &
               /sqrt(0.5D0*(bb(i,j,k)+bb(i,j-1,k))+tiny) &
               *0.5D0*(q(i,j-1,k)+q(i,j,k))
             fsatup=aval &
               *sqrt(gamma*(gamma-1.D0)*eint(i,j,k)/d(i,j,k)) &
               *eint(i,j,k)
             fsatdn=aval &
               *sqrt(gamma*(gamma-1.D0)*eint(i,j-1,k)/d(i,j-1,k)) &
               *eint(i,j-1,k)
             fsatq=0.5D0*(fsatup+fsatdn) &
               -sign(0.5D0,qflx2(i,j,k))*(fsatup-fsatdn)
             qflx2(i,j,k)=sign(1.D0,qflx2(i,j,k))*min(abs(fsatq),abs(qflx2(i,j,k)))
             fsat(j)=0.5D0*(fsatup+fsatdn) &
               -sign(0.5D0,flux(j))*(fsatup-fsatdn)
           enddo
           do j=js,je
             qlim(i,j,k)=max(qlim(i,j,k),abs(flux(j))/fsat(j), &
               abs(flux(j+1))/fsat(j+1))
           enddo
         enddo
       enddo
       do j=js,je
         do i=is,ie
           do k=ks,kep1
             flx3=b3(i,j,k) &
               *0.5D0*(bdelt(i,j,k)+bdelt(i,j,k-1)) &
               /(0.5D0*(bb(i,j,k)+bb(i,j,k-1))+tiny)
             flux(k)=-0.5D0*(diffcond(i,j,k)+diffcond(i,j,k-1)) &
               *flx3
             qflx3(i,j,k)=b3(i,j,k) &
               /sqrt(0.5D0*(bb(i,j,k)+bb(i,j,k-1))+tiny) &
               *0.5D0*(q(i,j,k-1)+q(i,j,k))
             fsatup=aval &
               *sqrt(gamma*(gamma-1.D0)*eint(i,j,k)/d(i,j,k)) &
               *eint(i,j,k)
             fsatdn=aval &
               *sqrt(gamma*(gamma-1.D0)*eint(i,j,k-1)/d(i,j,k-1)) &
               *eint(i,j,k-1)
             fsatq=0.5D0*(fsatup+fsatdn) &
               -sign(0.5D0,qflx3(i,j,k))*(fsatup-fsatdn)
             qflx3(i,j,k)=sign(1.D0,qflx3(i,j,k))*min(abs(fsatq),abs(qflx3(i,j,k)))
             fsat(k)=0.5D0*(fsatup+fsatdn) &
               -sign(0.5D0,flux(k))*(fsatup-fsatdn)
           enddo
           do k=ks,ke
             qlim(i,j,k)=max(qlim(i,j,k),abs(flux(k))/fsat(k), &
               abs(flux(k+1))/fsat(k+1))
           enddo
         enddo
       enddo
       call bvalqlim(qlim)
!
       do k=ksm1,kep1
         do j=jsm1,jep1
           do i=ism1,iep1
             diffcond(i,j,k)=diffcond(i,j,k)/qlim(i,j,k)
             dr1 = dx1a(i)
             dr2 = dx2a(j)*g2b(i)
             dr3 = dx3a(k)*g2b(i)*g32b(j)
             drmin = min(dr1,dr2,dr3)
             dfactor=diffcond(i,j,k)/d(i,j,k) &
               *(gamma-1.D0)*unit_v**2/(rstar*unit_temp)
             tau_hc(i,j,k) = max(dt*2.D0,dfactor*dt**2/(drmin*courno)**2)
           enddo
         enddo
       enddo
!
!      compute delq and delei
!
       do k=ks,ke
         do j=js,je
           do i=is,ie
             ovtau=1.D0/tau_hc(i,j,k)
             qspitzer=-diffcond(i,j,k)*bdelt(i,j,k) &
               /sqrt(bb(i,j,k)+tiny)
             delq(i,j,k)=(qspitzer-q(i,j,k))*ovtau
             delei(i,j,k)=delei(i,j,k) &
               - (qflx1(i+1,j,k)*g2a(i+1)*g31a(i+1) &
                 -qflx1(i,j,k)*g2a(i)*g31a(i)) &
               * dx1ai(i) * g2bi(i) * g31bi(i) &
               - (qflx2(i,j+1,k)*g32a(j+1)-qflx2(i,j,k)*g32a(j)) &
               * dx2ai(j) * g2bi(i) * g32bi(j) &
               - (qflx3(i,j,k+1)-qflx3(i,j,k)) &
               * dx3ai(k) * g31bi(i) * g32bi(j)
           enddo
         enddo
       enddo
!
!      add numerical diffusion to delq
!
       do k=ks,ke; do j=js,je
         qint(ism2:iep2) = q(ism2:iep2,j,k)
         call x1tvd(qint,ql,qr,dq,qlr)
         do i = is,iep1
#ifndef ISOMHD
           cs2=gamma*(gamma-1.D0)*0.5D0*(eint(i-1,j,k)+eint(i,j,k)) &
              /(0.5D0*(d(i-1,j,k)+d(i,j,k)))
#endif
           cmax= &
             sqrt( min( clight**2, ( (b1(i,j,k))**2 &
             +(0.25D0*(b2(i,j,k)+b2(i,j+1,k) &
             +b2(i-1,j,k)+b2(i-1,j+1,k)))**2 &
             +(0.25D0*(b3(i,j,k)+b3(i,j,k+1) &
             +b3(i-1,j,k)+b3(i-1,j,k+1)))**2 ) &
             /(0.5D0*(d(i,j,k)+d(i-1,j,k))) ) &
             + cs2 )
           rr = abs(dq(i))/(abs(qint(i)-qint(i-1)) + tiny)
           qq=(min(1.D0,alpha_e*rr))**nlf_e
           flux(i) = (- 0.5D0*cmax*qq &
             *dq(i))*g2a(i) * g31a(i)
         enddo
         do i=is,ie
           delq(i,j,k) = delq(i,j,k) &
             - (flux(i+1) - flux(i))*dvl1ai(i)
         enddo
       enddo; enddo
       do k=ks,ke; do i=is,ie
         qint(jsm2:jep2) = q(i,jsm2:jep2,k)
         call x2tvd (qint,ql,qr,dq,qlr)
         do j=js,jep1
#ifndef ISOMHD
           cs2=gamma*(gamma-1.D0)*0.5D0*(eint(i,j-1,k)+eint(i,j,k)) &
              /(0.5D0*(d(i,j-1,k)+d(i,j,k)))
#endif
           cmax = &
             sqrt( min( clight**2, ( (0.25D0*(b1(i,j,k)+b1(i+1,j,k) &
             +b1(i,j-1,k)+b1(i+1,j-1,k)))**2 &
             +(b2(i,j,k))**2 &
             +(0.25D0*(b3(i,j,k)+b3(i,j,k+1) &
             +b3(i,j-1,k)+b3(i,j-1,k+1)))**2 ) &
             /(0.5D0*(d(i,j,k)+d(i,j-1,k))) ) &
             + cs2 )
           rr = abs(dq(j))/(abs(qint(j)-qint(j-1)) + tiny)
           qq = (min(1.D0, alpha_e*rr))**nlf_e
           flux(j) = (- 0.5D0*cmax*qq*dq(j))*g32a(j)
         enddo
         do j=js,je
           delq(i,j,k) = delq(i,j,k) &
             - (flux(j+1) - flux(j))*g2bi(i)*dvl2ai(j)
         enddo
       enddo; enddo
       do j=js,je; do i=is,ie
         qint(ksm2:kep2) = q(i,j,ksm2:kep2)
         call x3tvd (qint,ql,qr,dq,qlr)
         do k=ks,kep1
#ifndef ISOMHD
           cs2=gamma*(gamma-1.D0)*0.5D0*(eint(i,j,k-1)+eint(i,j,k)) &
              /(0.5D0*(d(i,j,k-1)+d(i,j,k)))
#endif
           cmax = &
             sqrt( min( clight**2, ( (0.25D0*(b1(i,j,k)+b1(i+1,j,k) &
             +b1(i,j,k-1)+b1(i+1,j,k-1)))**2 &
             +(0.25D0*(b2(i,j,k)+b2(i,j+1,k) &
             +b2(i,j,k-1)+b2(i,j+1,k-1)))**2 &
             +(b3(i,j,k))**2 ) &
             /(0.5D0*(d(i,j,k)+d(i,j,k-1))) ) &
             + cs2 )
           rr = abs(dq(k))/(abs(qint(k)-qint(k-1)) + tiny)
           qq = (min(1.D0, alpha_e*rr))**nlf_e
           flux(k) = (- 0.5D0*cmax*qq*dq(k))
         enddo
         do k=ks,ke
            delq(i,j,k) = delq(i,j,k) &
              - (flux(k+1) - flux(k)) &
              *g31bi(i)*g32bi(j)*dvl3ai(k)
         enddo
       enddo; enddo
!

end subroutine derivs_cond

subroutine derivs_radloss

end subroutine derivs_radloss

subroutine derivs_heating

end subroutine derivs_heating

subroutine bvalqlim(qlim)
  use ModPar
  use ModGrid
  use ModBoundary,   ONLY: niib,noib,nijb,nojb,nikb,nokb

  implicit none

  real, intent(inout) :: qlim(in,jn,kn)

  integer :: idtop, idbot
  integer :: isendtop, isendbot, irecvtop, irecvbot
  real, allocatable :: sendtop(:,:,:), sendbot(:,:,:)
  real, allocatable :: recvtop(:,:,:), recvbot(:,:,:)
  integer :: i, j, k
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
  sendtop=qlim(iem1:ie,js:je,ks:ke)
  sendbot=qlim(is:isp1,js:je,ks:ke)
  call message_pass_nodes(idtop, idbot, &
    2, je-js+1, ke-ks+1, 2, je-js+1, ke-ks+1, &
    sendtop, sendbot, recvtop, recvbot, &
    isendtop, isendbot, irecvtop, irecvbot)
  qlim(iep1:iep2,js:je,ks:ke) = recvtop
  qlim(ism2:ism1,js:je,ks:ke) = recvbot
  deallocate(recvtop, recvbot, sendtop, sendbot)
!
! Inner i boundary.
!
  if(myid1 .eq. 0) then
    if(niib .ne. 4) then
    do k=ks,ke
      do j=js,je
        qlim(ism1,j,k) = qlim(is  ,j,k)
        qlim(ism2,j,k) = qlim(is  ,j,k)
      enddo
    enddo
    endif
  endif
!
! Outer i boundary.
!
  if(myid1 .eq. nproc1-1) then
    if(noib .ne. 4) then
    do k=ks,ke
      do j=js,je
        qlim(iep1,j,k) = qlim(ie,j,k)
        qlim(iep2,j,k) = qlim(ie,j,k)
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
  sendtop=qlim(ism2:iep2,jem1:je,ks:ke)
  sendbot=qlim(ism2:iep2,js:jsp1,ks:ke)
  call message_pass_nodes(idtop, idbot, &
    iep2-ism2+1, 2, ke-ks+1, iep2-ism2+1, 2, ke-ks+1, &
    sendtop, sendbot, recvtop, recvbot, &
    isendtop, isendbot, irecvtop, irecvbot)
  qlim(ism2:iep2,jep1:jep2,ks:ke) = recvtop
  qlim(ism2:iep2,jsm2:jsm1,ks:ke) = recvbot
  deallocate(recvtop, recvbot, sendtop, sendbot)
!
! Inner j boundary.
!
  if(myid2 .eq. 0) then
    if(nijb .ne. 4) then
    do i=ism2,iep2
      do k=ks,ke
        qlim(i,jsm1,k) = qlim(i,js  ,k)
        qlim(i,jsm2,k) = qlim(i,js  ,k)
      enddo
    enddo
    endif
  endif
!
! Outer j boundary.
!
  if(myid2 .eq. nproc2-1) then
    if(nojb .ne. 4) then
    do i=ism2,iep2
      do k=ks,ke
        qlim(i,jep1,k) = qlim(i,je  ,k)
        qlim(i,jep2,k) = qlim(i,je  ,k)
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
  sendtop=qlim(ism2:iep2,jsm2:jep2,kem1:ke)
  sendbot=qlim(ism2:iep2,jsm2:jep2,ks:ksp1)
  call message_pass_nodes(idtop, idbot, &
    iep2-ism2+1, jep2-jsm2+1, 2, iep2-ism2+1, jep2-jsm2+1, 2, &
    sendtop, sendbot, recvtop, recvbot, &
    isendtop, isendbot, irecvtop, irecvbot)
  qlim(ism2:iep2,jsm2:jep2,kep1:kep2) = recvtop
  qlim(ism2:iep2,jsm2:jep2,ksm2:ksm1) = recvbot
  deallocate(recvtop, recvbot, sendtop, sendbot)
! Inner k boundary.
  if(myid3 .eq. 0) then
    if(nikb .ne. 4) then
    do j=jsm2,jep2
      do i=ism2,iep2
        qlim(i,j,ksm1) = qlim(i,j,ks  )
        qlim(i,j,ksm2) = qlim(i,j,ks  )
      enddo
    enddo
    endif
  endif
! Outer k boundary.
  if(myid3 .eq. nproc3-1) then
    if(nokb .ne. 4) then
    do j=jsm2,jep2
      do i=ism2,iep2
        qlim(i,j,kep1) = qlim(i,j,ke  )
        qlim(i,j,kep2) = qlim(i,j,ke  )
      enddo
    enddo
    endif
  endif

end subroutine bvalqlim

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

#endif

end module ModNonIdealRHS
