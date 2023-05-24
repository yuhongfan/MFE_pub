module ModIteration

  implicit none

  private

  public :: nudt
  public :: rk3step

contains

subroutine nudt
  use ModPar,        ONLY: tiny, huge
  use ModPhysics,    ONLY: gamma, cs2
  use ModGrid
  use ModField
  use ModSundry,     ONLY: dt, deltat, time, courno, clight, dtmdi2

  implicit none

  include "mpif.h"

  integer ::  ierr, i, j, k, nnow, nnex
  real :: dtnewlc, dtnew
  real :: q1,va,vv1,vv2,vv3,vv,cmax,dr1,dr2,dr3,drmin
  real, parameter :: small=1.D-6
!-----------------------------------------------------------------------
#ifdef ISOMHD
       q1=sqrt(cs2)
#endif
       dtnewlc = huge
       do k=ks,ke
         do j=js,je
           do i=is,ie
#ifndef ISOMHD
             cs2=gamma*(gamma-1.D0)*eint(i,j,k)/d(i,j,k)
             q1=sqrt(cs2)
#endif
             va = sqrt((max(abs(b1(i,j,k)),abs(b1(i+1,j,k))))**2 &
                  + (max(abs(b2(i,j,k)),abs(b2(i,j+1,k))))**2 &
                  + (max(abs(b3(i,j,k)),abs(b3(i,j,k+1))))**2) &
                  / sqrt(d(i,j,k))
             va = min(va, clight)
             vv1 = max( abs(v1(i,j,k)) , abs(v1(i+1,j,k)) )
             vv2 = max( abs(v2(i,j,k)) , abs(v2(i,j+1,k)) )
             vv3 = max( abs(v3(i,j,k)) , abs(v3(i,j,k+1)) )
             vv = sqrt( vv1**2 + vv2**2 + vv3**2)
             dr1 = dx1a(i)
             dr2 = dx2a(j)*g2b(i)
             dr3 = dx3a(k)*g2b(i)*g32b(j)
             drmin = min(dr1,dr2,dr3)
             cmax = vv+va+q1
             dtnewlc = min(dtnewlc, drmin/(2.D0*cmax+tiny), &
               1.D0/sqrt(dtmdi2+tiny))
           enddo
         enddo
       enddo

       call MPI_ALLREDUCE(dtnewlc,dtnew,1, &
         MPI_DOUBLE_PRECISION,MPI_MIN, &
         MPI_COMM_WORLD,ierr)

       dt     =   courno*dtnew

       nnow=int(time/deltat)
       nnex=int((time+dt)/deltat)
       if(nnex > nnow) then
         dt=min(((nnow+1)+small)*deltat-time, (courno+small)*dtnew)
       endif

end subroutine nudt

subroutine rk3step
  use ModGrid
  use ModField
  use ModDel
  use ModRHS
  use ModNonIdealRHS

  implicit none

  integer       i       , j       , k
  real ::       cstp,ccstp,gstp

!-----------------------------------------------------------------------
!
!      step 1
!
       call derivs
#ifdef HEAT_RADLOSS
       call derivs_radloss
       call derivs_heating
#endif
#ifdef HOLLWEG
       call derivs_hollweg
#endif
#ifdef THCONDUCT
       call derivs_cond
#endif
       cstp = 1.D0/3.D0
       ccstp = 1.D0/3.D0
       gstp = -5.D0/9.D0
       call update(cstp,ccstp)
       emf1s_pr(ism2:iep2,jsm2:jep3,ksm2:kep3) = &
         gstp*emf1s(ism2:iep2,jsm2:jep3,ksm2:kep3)
       emf2s_pr(ism2:iep3,jsm2:jep2,ksm2:kep3) = &
         gstp*emf2s(ism2:iep3,jsm2:jep2,ksm2:kep3)
       emf3s_pr(ism2:iep3,jsm2:jep3,ksm2:kep2) = &
         gstp*emf3s(ism2:iep3,jsm2:jep3,ksm2:kep2)
       deld_pr(ism2:iep2,jsm2:jep2,ksm2:kep2) = &
         gstp*deld(ism2:iep2,jsm2:jep2,ksm2:kep2)
       dels1_pr(ism2:iep2,jsm2:jep2,ksm2:kep2) = &
         gstp*dels1(ism2:iep2,jsm2:jep2,ksm2:kep2)
       dels2_pr(ism2:iep2,jsm2:jep2,ksm2:kep2) = &
         gstp*dels2(ism2:iep2,jsm2:jep2,ksm2:kep2)
       dels3_pr(ism2:iep2,jsm2:jep2,ksm2:kep2) = &
         gstp*dels3(ism2:iep2,jsm2:jep2,ksm2:kep2)
#ifndef ISOMHD
       delei_pr(ism2:iep2,jsm2:jep2,ksm2:kep2) = &
         gstp*delei(ism2:iep2,jsm2:jep2,ksm2:kep2)
#endif
#ifdef THCONDUCT
       delq_pr(is:ie,js:je,ks:ke) = &
         gstp*delq(is:ie,js:je,ks:ke)
#endif
!
!      step 2
!
       call derivs
#ifdef HEAT_RADLOSS
       call derivs_radloss
       call derivs_heating
#endif
#ifdef HOLLWEG
       call derivs_hollweg
#endif
#ifdef THCONDUCT
       call derivs_cond
#endif
       emf1s(ism2:iep2,jsm2:jep3,ksm2:kep3) = &
         emf1s_pr(ism2:iep2,jsm2:jep3,ksm2:kep3) + &
         emf1s(ism2:iep2,jsm2:jep3,ksm2:kep3)
       emf2s(ism2:iep3,jsm2:jep2,ksm2:kep3) = &
         emf2s_pr(ism2:iep3,jsm2:jep2,ksm2:kep3) + &
         emf2s(ism2:iep3,jsm2:jep2,ksm2:kep3)
       emf3s(ism2:iep3,jsm2:jep3,ksm2:kep2) = &
         emf3s_pr(ism2:iep3,jsm2:jep3,ksm2:kep2) + &
         emf3s(ism2:iep3,jsm2:jep3,ksm2:kep2)
       deld(ism2:iep2,jsm2:jep2,ksm2:kep2) = &
         deld_pr(ism2:iep2,jsm2:jep2,ksm2:kep2) + &
         deld(ism2:iep2,jsm2:jep2,ksm2:kep2)
       dels1(ism2:iep2,jsm2:jep2,ksm2:kep2) = &
         dels1_pr(ism2:iep2,jsm2:jep2,ksm2:kep2) + &
         dels1(ism2:iep2,jsm2:jep2,ksm2:kep2)
       dels2(ism2:iep2,jsm2:jep2,ksm2:kep2) = &
         dels2_pr(ism2:iep2,jsm2:jep2,ksm2:kep2) + &
         dels2(ism2:iep2,jsm2:jep2,ksm2:kep2)
       dels3(ism2:iep2,jsm2:jep2,ksm2:kep2) = &
         dels3_pr(ism2:iep2,jsm2:jep2,ksm2:kep2) + &
         dels3(ism2:iep2,jsm2:jep2,ksm2:kep2)
#ifndef ISOMHD
       delei(ism2:iep2,jsm2:jep2,ksm2:kep2) = &
         delei_pr(ism2:iep2,jsm2:jep2,ksm2:kep2) + &
         delei(ism2:iep2,jsm2:jep2,ksm2:kep2)
#endif
#ifdef THCONDUCT
       delq(is:ie,js:je,ks:ke) = &
         delq_pr(is:ie,js:je,ks:ke) + &
         delq(is:ie,js:je,ks:ke)
#endif
       cstp=15.D0/16.D0
       ccstp=3.D0/4.D0
       gstp=-153.D0/128.D0
       call update(cstp,ccstp)
       emf1s_pr(ism2:iep2,jsm2:jep3,ksm2:kep3) = &
         gstp*emf1s(ism2:iep2,jsm2:jep3,ksm2:kep3)
       emf2s_pr(ism2:iep3,jsm2:jep2,ksm2:kep3) = &
         gstp*emf2s(ism2:iep3,jsm2:jep2,ksm2:kep3)
       emf3s_pr(ism2:iep3,jsm2:jep3,ksm2:kep2) = &
         gstp*emf3s(ism2:iep3,jsm2:jep3,ksm2:kep2)
       deld_pr(ism2:iep2,jsm2:jep2,ksm2:kep2) = &
         gstp*deld(ism2:iep2,jsm2:jep2,ksm2:kep2)
       dels1_pr(ism2:iep2,jsm2:jep2,ksm2:kep2) = &
         gstp*dels1(ism2:iep2,jsm2:jep2,ksm2:kep2)
       dels2_pr(ism2:iep2,jsm2:jep2,ksm2:kep2) = &
         gstp*dels2(ism2:iep2,jsm2:jep2,ksm2:kep2)
       dels3_pr(ism2:iep2,jsm2:jep2,ksm2:kep2) = &
         gstp*dels3(ism2:iep2,jsm2:jep2,ksm2:kep2)
#ifndef ISOMHD
       delei_pr(ism2:iep2,jsm2:jep2,ksm2:kep2) = &
         gstp*delei(ism2:iep2,jsm2:jep2,ksm2:kep2)
#endif
#ifdef THCONDUCT
       delq_pr(is:ie,js:je,ks:ke) = &
         gstp*delq(is:ie,js:je,ks:ke)
#endif
!
!      step 3
!
       call derivs
#ifdef HEAT_RADLOSS
       call derivs_radloss
       call derivs_heating
#endif
#ifdef HOLLWEG
       call derivs_hollweg
#endif
#ifdef THCONDUCT
       call derivs_cond
#endif
       emf1s(ism2:iep2,jsm2:jep3,ksm2:kep3) = &
         emf1s_pr(ism2:iep2,jsm2:jep3,ksm2:kep3) + &
         emf1s(ism2:iep2,jsm2:jep3,ksm2:kep3)
       emf2s(ism2:iep3,jsm2:jep2,ksm2:kep3) = &
         emf2s_pr(ism2:iep3,jsm2:jep2,ksm2:kep3) + &
         emf2s(ism2:iep3,jsm2:jep2,ksm2:kep3)
       emf3s(ism2:iep3,jsm2:jep3,ksm2:kep2) = &
         emf3s_pr(ism2:iep3,jsm2:jep3,ksm2:kep2) + &
         emf3s(ism2:iep3,jsm2:jep3,ksm2:kep2)
       deld(ism2:iep2,jsm2:jep2,ksm2:kep2) = &
         deld_pr(ism2:iep2,jsm2:jep2,ksm2:kep2) + &
         deld(ism2:iep2,jsm2:jep2,ksm2:kep2)
       dels1(ism2:iep2,jsm2:jep2,ksm2:kep2) = &
         dels1_pr(ism2:iep2,jsm2:jep2,ksm2:kep2) + &
         dels1(ism2:iep2,jsm2:jep2,ksm2:kep2)
       dels2(ism2:iep2,jsm2:jep2,ksm2:kep2) = &
         dels2_pr(ism2:iep2,jsm2:jep2,ksm2:kep2) + &
         dels2(ism2:iep2,jsm2:jep2,ksm2:kep2)
       dels3(ism2:iep2,jsm2:jep2,ksm2:kep2) = &
         dels3_pr(ism2:iep2,jsm2:jep2,ksm2:kep2) + &
         dels3(ism2:iep2,jsm2:jep2,ksm2:kep2)
#ifndef ISOMHD
       delei(ism2:iep2,jsm2:jep2,ksm2:kep2) = &
         delei_pr(ism2:iep2,jsm2:jep2,ksm2:kep2) + &
         delei(ism2:iep2,jsm2:jep2,ksm2:kep2)
#endif
#ifdef THCONDUCT
       delq(is:ie,js:je,ks:ke) = &
         delq_pr(is:ie,js:je,ks:ke) + &
         delq(is:ie,js:je,ks:ke)
#endif
       cstp=8.D0/15.D0
       ccstp=1.D0
       gstp=0.D0
       call update(cstp,ccstp)

end subroutine rk3step

subroutine update(cstp,ccstp)
  use ModPar
  use ModPhysics
  use ModGrid
  use ModField
  use ModDel
  use ModSundry
  use ModBval
  use ModBoundary,     ONLY: niib, noib, nijb, nojb, nikb, nokb
  use ModUserSetup,    ONLY: getbib, getbjb, getbkb
  use ModWork,         ONLY: scratch1, scratch2, scratch3
  implicit none

  include "mpif.h"

  real, intent(in) :: cstp,ccstp
  integer :: ierr
  integer :: i       , ip1     , j       , jp1     , k &
             , kp1
!
!-----------------------------------------------------------------------
!      Update ctime so that boudary conditions are for the updated fields
!-----------------------------------------------------------------------
       ctime=time+ccstp*dt
!
!-----------------------------------------------------------------------
!-----------------------------> Update b1 <-----------------------------
!-----------------------------------------------------------------------
!
       do k=ksm2,kep2
         kp1 = k + kone
         do j=jsm2,jep2
           jp1 = j + jone
           do i=ism2,iep2
             b1(i,j,k) = b1(i,j,k) &
               + cstp * dt * &
               ( emf3s(i,jp1,k)*dx3a(k)*g31a(i)*g32a(jp1) &
                 - emf3s(i,j,k)*dx3a(k)*g31a(i)*g32a(j) &
                 - emf2s(i,j,kp1)*dx2a(j)*g2a(i) &
                 + emf2s(i,j,k)*dx2a(j)*g2a(i) ) &
               * g2ai (i) * g31ai(i) * g32bi(j) &
               * dx2ai(j) * dx3ai(k)
           enddo
         enddo
       enddo
!
!-----------------------------------------------------------------------
!-----------------------------> Update b2 <-----------------------------
!-----------------------------------------------------------------------
!
       do k=ksm2,kep2
         kp1 = k + kone
         do j=jsm2,jep2
           do i=ism2,iep2
             ip1 = i + ione
             b2(i,j,k) = b2(i,j,k) &
               + cstp * dt * &
               ( emf1s(i,j,kp1)*dx1a(i) &
                 - emf1s(i,j,k)*dx1a(i) &
                 - emf3s(ip1,j,k)*dx3a(k)*g31a(ip1)*g32a(j) &
                 + emf3s(i,j,k)*dx3a(k)*g31a(i)*g32a(j) ) &
               * g31bi(i) * g32ai(j) &
               * dx1ai(i) * dx3ai(k)
           enddo
         enddo
       enddo
!
!-----------------------------------------------------------------------
!-----------------------------> Update b3 <-----------------------------
!-----------------------------------------------------------------------
!
       do k=ksm2,kep2
         do j=jsm2,jep2
           jp1 = j + jone
           do i=ism2,iep2
             ip1 = i + ione
             b3(i,j,k) = b3(i,j,k) &
               + cstp * dt * &
               ( emf2s(ip1,j,k)*dx2a(j)*g2a(ip1) &
                - emf2s(i,j,k)*dx2a(j)*g2a(i) &
                - emf1s(i,jp1,k)*dx1a(i) &
                + emf1s(i,j,k)*dx1a(i) ) &
               * g2bi (i) * dx1ai(i) * dx2ai(j)
           enddo
         enddo
       enddo
!
!-----------------------------------------------
!------ user special boundary for b-------------
!-----------------------------------------------
!
       if(niib .eq. 3 .or. noib .eq. 3) then
         call getbib
       endif
       if(nijb .eq. 3 .or. nojb .eq. 3) then
         call getbjb
       endif
       if(nikb .eq. 3 .or. nokb .eq. 3) then
         call getbkb
       endif
!
!-----------------------------------------------
!------ update momenta -------------------------
!-----------------------------------------------
!
       do k=ks,ke
         do j=js,je
           do i=is,iep1
             scratch1(i,j,k) = dels1(i,j,k) * cstp * dt &
                + v1(i,j,k) &
                * 0.5D0 * (d(i-1,j  ,k  )+d(i,j,k))
           enddo
         enddo
       enddo
       do k=ks,ke
         do j=js,jep1
           do i=is,ie
             scratch2(i,j,k) = dels2(i,j,k) * cstp * dt &
                + v2(i,j,k) &
                * 0.5D0 * (d(i  ,j-1,k  )+d(i,j,k))
           enddo
         enddo
       enddo
       do k=ks,kep1
         do j=js,je
           do i=is,ie
             scratch3(i,j,k) = dels3(i,j,k) * cstp * dt &
                + v3(i,j,k) &
                * 0.5D0 * (d(i  ,j  ,k-1)+d(i,j,k))
           enddo
         enddo
       enddo
!
!-----------------------------------------------
!------- update d-------------------------------
!-----------------------------------------------
!
       do k=ks,ke
         do j=js,je
           do i=is,ie
             d(i,j,k) = d(i,j,k) + cstp * dt * deld(i,j,k)
             if(d(i,j,k) .lt. dfloor) then
!               write(6,*) 'd=',d(i,j,k),' at ',x1b(i),x2b(j),x3b(k)
                d(i,j,k)=dfloor
             endif
           enddo
         enddo
       enddo
       call bvald
!
!-----------------------------------------------------
!----- update velocities -----------------------------
!-----------------------------------------------------
!
       do k=ks,ke
         do j=js,je
           do i=is,iep1
             v1(i,j,k) = scratch1(i,j,k) &
                        * 2.D0 / (d(i-1,j  ,k  )+d(i,j,k))
           enddo
         enddo
       enddo
       do k=ks,ke
         do j=js,jep1
           do i=is,ie
             v2(i,j,k) = scratch2(i,j,k) &
                        * 2.D0 / (d(i  ,j-1,k  )+d(i,j,k))
           enddo
         enddo
       enddo
       do k=ks,kep1
         do j=js,je
           do i=is,ie
             v3(i,j,k) = scratch3(i,j,k) &
                        * 2.D0 / (d(i  ,j  ,k-1)+d(i,j,k))
           enddo
         enddo
       enddo
       call bvalv
#ifndef ISOMHD
!
!--------------------------------------------------
!------- update energy ----------------------------
!--------------------------------------------------
!
       do k=ks,ke
         do j=js,je
           do i=is,ie
             eint(i,j,k) = eint(i,j,k) + cstp * dt * delei(i,j,k)
             if(eint(i,j,k) .lt. efloor) then
!               write(6,*) 'eint=',eint(i,j,k),' at ',x1b(i),x2b(j),x3b(k)
                eint(i,j,k)=efloor
             endif
           enddo
         enddo
       enddo
       call bvalei
       do k=ksm2,kep2
         do j=jsm2,jep2
           do i=ism2,iep2
             tfield(i,j,k)=eint(i,j,k)*(gamma-1.D0) &
                /d(i,j,k)*unit_v**2/(rstar*unit_temp)
           enddo
         enddo
       enddo
#endif

#ifdef THCONDUCT
       q(is:ie,js:je,ks:ke)=q(is:ie,js:je,ks:ke) &
          + cstp*dt*delq(is:ie,js:je,ks:ke)
       call bvaleta(q)
#endif
       
end subroutine update

end module ModIteration
