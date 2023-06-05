module ModRHS

  implicit none

  private

  public :: derivs

contains

subroutine derivs
  use ModPar
  use ModGrid
  use ModDel
  use ModBoundary
  use ModUserSetup,     ONLY: getdelsib, getdelsjb, getdelskb
  use ModIO
  implicit none

  integer :: i, j, k
!
! checking
  character(len=80) :: filecheck
!
!-----------------------------------------------------------------------
!
!      initialize the del arrays
       dels1 = 0.
       dels2 = 0.
       dels3 = 0.
       deld = 0.
#ifndef ISOMHD
       delei = 0.
#endif
       emf1s = 0.
       emf2s = 0.
       emf3s = 0.

!      add hydrodynamic source terms to del arrays
       call stv1
       call stv2
       call stv3
!
!      special boundary condition for line-tying
!
       if(myid1 .eq. 0 .or. myid1 .eq. nproc1-1) then
         if(niib .eq. 3 .or. noib .eq. 3) then
           call getdelsib
         endif
       endif
       if(myid2 .eq. 0 .or. myid2 .eq. nproc2-1) then
         if(nijb .eq. 3 .or. nojb .eq. 3) then
           call getdelsjb
         endif
       endif
       if(myid3 .eq. 0 .or. myid3 .eq. nproc3-1) then
         if(nikb .eq. 3 .or. nokb .eq. 3) then
           call getdelskb
         endif
       endif

       call pdv
       call advection
       call lorentz
       call hsmoc
       call diffemf

!      checking
!      filecheck='dels1.dat'
!      call writeqty_mpi(dels1,filecheck)
!      filecheck='dels2.dat'
!      call writeqty_mpi(dels2,filecheck)
!      filecheck='dels3.dat'
!      call writeqty_mpi(dels3,filecheck)
!      filecheck='deld.dat'
!      call writeqty_mpi(deld,filecheck)
!      filecheck='delei.dat'
!      call writeqty_mpi(delei,filecheck)
!      filecheck='emf1s.dat'
!      call writeemf1_mpi(emf1s,filecheck)
!      filecheck='emf2s.dat'
!      call writeemf2_mpi(emf2s,filecheck)
!      filecheck='emf3s.dat'
!      call writeemf3_mpi(emf3s,filecheck)

end subroutine derivs

#ifdef JCROSSB
subroutine lorentz
  use ModPar
  use ModGrid
  use Modfield
  use ModDel
  use ModWork
  use ModBoundary
  use ModSundry,        ONLY: ovc2
  use ModBval,          ONLY: bvaleta
  use ModUserSetup,     ONLY: getdelsib, getdelsjb, getdelskb
  implicit none

  integer ::      i       , j       , k
  real :: q1
  real, dimension(ijkn) :: bave, bstar
  real, dimension(in,jn,kn) :: bdotdels, bb, bdotv, fh1, fh2, fh3
!
!-----------------------------------------------------------------------
!
       do k=ksm1,kep1
         do j=jsm1,jep1
           do i=ism1,iep1
             srd1(i,j,k) = ( (b1(i,j,k))**2 &
               +(0.25D0*(b2(i,j,k)+b2(i,j+1,k) &
               +b2(i-1,j,k)+b2(i-1,j+1,k)))**2 &
               +(0.25D0*(b3(i,j,k)+b3(i,j,k+1) &
               +b3(i-1,j,k)+b3(i-1,j,k+1)))**2 ) &
               /(0.5D0*(d(i,j,k)+d(i-1,j,k)))*ovc2
             srd2(i,j,k) = ( (0.25D0*(b1(i,j,k)+b1(i+1,j,k) &
               +b1(i,j-1,k)+b1(i+1,j-1,k)))**2 &
               +(b2(i,j,k))**2 &
               +(0.25D0*(b3(i,j,k)+b3(i,j,k+1) &
               +b3(i,j-1,k)+b3(i,j-1,k+1)))**2 ) &
               /(0.5D0*(d(i,j,k)+d(i,j-1,k)))*ovc2
             srd3(i,j,k) = ( (0.25D0*(b1(i,j,k)+b1(i+1,j,k) &
               +b1(i,j,k-1)+b1(i+1,j,k-1)))**2 &
               +(0.25D0*(b2(i,j,k)+b2(i,j+1,k) &
               +b2(i,j,k-1)+b2(i,j+1,k-1)))**2 &
               +(b3(i,j,k))**2 ) &
               /(0.5D0*(d(i,j,k)+d(i,j,k-1)))*ovc2
           enddo
         enddo
       enddo

       do k=ks,ke
         do j=js,je
           do i=is,ie
             bdotdels(i,j,k) = 0.5D0*( b1(i,j,k) &
               *dels1(i,j,k) &
               +b1(i+1,j,k)*dels1(i+1,j,k) ) &
               + 0.5D0*( b2(i,j,k)*dels2(i,j,k) &
               +b2(i,j+1,k)*dels2(i,j+1,k) ) &
               + 0.5D0*( b3(i,j,k)*dels3(i,j,k) &
               +b3(i,j,k+1)*dels3(i,j,k+1) )
             bdotv(i,j,k) = 0.5D0*( b1(i,j,k) &
               *v1(i,j,k) &
               +b1(i+1,j,k)*v1(i+1,j,k) ) &
               + 0.5D0*( b2(i,j,k)*v2(i,j,k) &
               +b2(i,j+1,k)*v2(i,j+1,k) ) &
               + 0.5D0*( b3(i,j,k)*v3(i,j,k) &
               +b3(i,j,k+1)*v3(i,j,k+1) )
             bb(i,j,k) = 0.25D0 * ( (b1(i,j,k) &
               +b1(i+1,j,k))**2 + (b2(i,j,k) &
               +b2(i,j+1,k))**2 + (b3(i,j,k) &
               +b3(i,j,k+1))**2 )
           enddo
         enddo
       enddo
       call bvaleta(bdotdels)
       call bvaleta(bdotv)
       call bvaleta(bb)
       call bvaleta(deld)

       do k=ks,ke
         do j=js,je
           do i=is,iep1
             fh1(i,j,k) = srd1(i,j,k)/(1.D0+srd1(i,j,k)) &
               *( (dels1(i,j,k)-0.5D0*(bdotdels(i,j,k)+bdotdels(i-1,j,k)) &
               /(0.5D0*(bb(i,j,k)+bb(i-1,j,k))+tiny)*b1(i,j,k)) &
               - 0.5D0*(deld(i,j,k)+deld(i-1,j,k)) &
               *(v1(i,j,k)-0.5D0*(bdotv(i,j,k)+bdotv(i-1,j,k)) &
               /(0.5D0*(bb(i,j,k)+bb(i-1,j,k))+tiny)*b1(i,j,k)) )
           enddo
         enddo
       enddo
       do k=ks,ke
         do j=js,jep1
           do i=is,ie
             fh2(i,j,k) = srd2(i,j,k)/(1.D0+srd2(i,j,k)) &
               *( (dels2(i,j,k)-0.5D0*(bdotdels(i,j,k)+bdotdels(i,j-1,k)) &
               /(0.5D0*(bb(i,j,k)+bb(i,j-1,k))+tiny)*b2(i,j,k)) &
               - 0.5D0*(deld(i,j,k)+deld(i,j-1,k)) &
               *(v2(i,j,k)-0.5D0*(bdotv(i,j,k)+bdotv(i,j-1,k)) &
               /(0.5D0*(bb(i,j,k)+bb(i,j-1,k))+tiny)*b2(i,j,k)) )
           enddo
         enddo
       enddo
       do k=ks,kep1
         do j=js,je
           do i=is,ie
             fh3(i,j,k) = srd3(i,j,k)/(1.D0+srd3(i,j,k)) &
               *( (dels3(i,j,k)-0.5D0*(bdotdels(i,j,k)+bdotdels(i,j,k-1)) &
               /(0.5D0*(bb(i,j,k)+bb(i,j,k-1))+tiny)*b3(i,j,k)) &
               - 0.5D0*(deld(i,j,k)+deld(i,j,k-1)) &
               *(v3(i,j,k)-0.5D0*(bdotv(i,j,k)+bdotv(i,j,k-1)) &
               /(0.5D0*(bb(i,j,k)+bb(i,j,k-1))+tiny)*b3(i,j,k)) )
           enddo
         enddo
       enddo
       do k=ks,ke
         do j=js,je
           do i=is,iep1
             fh1(i,j,k)=dels1(i,j,k)-fh1(i,j,k)
           enddo
         enddo
       enddo
       do k=ks,ke
         do j=js,jep1
           do i=is,ie
             fh2(i,j,k)=dels2(i,j,k)-fh2(i,j,k)
           enddo
         enddo
       enddo
       do k=ks,kep1
         do j=js,je
           do i=is,ie
             fh3(i,j,k)=dels3(i,j,k)-fh3(i,j,k)
           enddo
         enddo
       enddo
!
!-----------------------------------------------------------------------
!---- 1-force ----------------------------------------------------------
!-----------------------------------------------------------------------
!
       do k=ks,ke
         do i=is,iep1
           do j=js,jep1
             bave (j) = 0.5D0 * ( b2(i,j,k) + b2(i-1,j,k) )
             bstar(j) = ( - b1(i,j-1,k) + b1(i,j,k) ) &
               * g2ai(i) * dx2bi(j)
           enddo
           do j=js,je
             dels1(i,j,k) = ( bave (j+1)*bstar(j+1) &
                        + bave (j)*bstar(j) ) * 0.5D0
           enddo
         enddo
       enddo
!
       do j=js,je
         do i=is,iep1
           do k=ks,kep1
             bave (k) = 0.5D0 * ( b3(i,j,k) + b3(i-1,j,k) )
             bstar(k) = ( - b1(i,j,k-1) + b1(i,j,k) ) &
               * g31ai(i) * g32bi(j) * dx3bi(k)
           enddo
           do k=ks,ke
             dels1(i,j,k) = dels1(i,j,k) &
                        + ( bave (k+1)*bstar(k+1) &
                        + bave (k)*bstar(k) ) * 0.5D0
           enddo
         enddo
       enddo
!
       do k=ks,kep1
         do j=js,jep1
           do i=is,iep1
             scratch1(i,j,k) = (b2(i,j,k)+b2(i-1,j,k)) &
                           * (g2b(i)*b2(i,j,k)-g2b(i-1)*b2(i-1,j,k)) &
                           * g2ai (i)
             scratch2(i,j,k) = (b3(i,j,k)+b3(i-1,j,k)) &
                           * (g31b(i)*b3(i,j,k)-g31b(i-1)*b3(i-1,j,k)) &
                           * g31ai(i)
           enddo
         enddo
       enddo
       do k=ks,ke
         do j=js,je
           do i=is,iep1
             dels1(i,j,k)  = dels1(i,j,k) &
                     - ( scratch1(i,j,k) + scratch1(i,j+1,k  ) &
                       + scratch2(i,j,k) + scratch2(i,j  ,k+1) ) &
                     * dx1bi(i) * 0.25D0
           enddo
         enddo
       enddo
!
!-----------------------------------------------------------------------
!---- 2-force ----------------------------------------------------------
!-----------------------------------------------------------------------
!
!
       do j=js,jep1
         do i=is,ie
           do k=ks,kep1
             bave (k) = 0.5D0 * ( b3(i,j,k) + b3(i,j-1,k) )
             bstar(k) = ( - b2(i,j,k-1) + b2(i,j,k) ) &
               * g31bi(i) * g32ai(j) * dx3bi(k)
           enddo
           do k=ks,ke
             dels2(i,j,k) = ( bave (k+1)*bstar(k+1) &
                        + bave (k)*bstar(k) ) * 0.5D0
           enddo
         enddo
       enddo
!
       do k=ks,ke
         do j=js,jep1
           do i=is,iep1
             bave (i) = 0.5D0 * ( b1(i,j,k) + b1(i,j-1,k) )
             bstar(i) = ( - g2b(i-1)*b2(i-1,j,k) + g2b(i)*b2(i,j,k) ) &
               * g2ai(i) * dx1bi(i)
           enddo
           do i=is,ie
             dels2(i,j,k) = dels2(i,j,k) &
                        + ( bave (i+1)*bstar(i+1) &
                        + bave (i)*bstar(i) ) * 0.5D0
           enddo
         enddo
       enddo
!
       do k=ks,kep1
         do j=js,jep1
           do i=is,iep1
             scratch1(i,j,k) = (b3(i,j,k)+b3(i,j-1,k)) &
                           * (g32b(j)*b3(i,j,k)-g32b(j-1)*b3(i,j-1,k)) &
                           * g32ai(j) * g2bi(i)
             scratch2(i,j,k) = (b1(i,j,k)+b1(i,j-1,k)) &
                           * (b1(i,j,k)-b1(i,j-1,k)) * g2ai(i)
           enddo
         enddo
       enddo
       do k=ks,ke
         do j=js,jep1
           do i=is,ie
             dels2(i,j,k)  = dels2(i,j,k) &
                     - ( scratch1(i,j,k) + scratch1(i  ,j,k+1) &
                       + scratch2(i,j,k) + scratch2(i+1,j,k  ) ) &
                     * dx2bi(j) * 0.25D0

           enddo
         enddo
       enddo
!
!-----------------------------------------------------------------------
!---- 3-force ----------------------------------------------------------
!-----------------------------------------------------------------------
!
!
       do k=ks,kep1
         do j=js,je
           do i=is,iep1
             bave (i) = 0.5D0 * ( b1(i,j,k) + b1(i,j,k-1) )
             bstar(i) = ( -g31b(i-1)*b3(i-1,j,k) + g31b(i)*b3(i,j,k) ) &
               * g31ai(i) * dx1bi(i)
           enddo
           do i=is,ie
             dels3(i,j,k) = ( bave (i+1)*bstar(i+1) &
                        + bave (i)*bstar(i) ) * 0.5D0
           enddo
         enddo
         do i=is,ie
           do j=js,jep1
             bave (j) = 0.5D0 * ( b2(i,j,k) + b2(i,j,k-1) )
             bstar(j) = ( - g32b(j-1)*b3(i,j-1,k) + g32b(j)*b3(i,j,k) ) &
               * g2bi(i) * g32ai(j) * dx2bi(j)
           enddo
           do j=js,je
             dels3(i,j,k) = dels3(i,j,k) &
                        + ( bave (j+1)*bstar(j+1) &
                        + bave (j)*bstar(j) ) * 0.5D0
           enddo
         enddo
       enddo

       do k=ks,kep1
         do j=js,jep1
           do i=is,iep1
             scratch1(i,j,k) = (b1(i,j,k  ) + b1(i,j,k-1)) &
                           * (b1(i,j,k) - b1(i,j,k-1)) &
                           * g31ai(i) * g32bi(j)
             scratch2(i,j,k) = (b2(i,j,k  ) + b2(i,j,k-1)) &
                           * (b2(i,j,k) - b2(i,j,k-1)) &
                           * g31bi(i) * g32ai(j)
           enddo
         enddo
       enddo
       do k=ks,kep1
         do j=js,je
           do i=is,ie
             dels3(i,j,k)  = dels3(i,j,k) &
                     - ( scratch1(i,j,k) + scratch1(i+1,j  ,k) &
                       + scratch2(i,j,k) + scratch2(i  ,j+1,k) ) &
                     * dx3bi(k) &
                     * 0.25D0
           enddo
         enddo
       enddo
!
!      special boundary condition for line-tying
!
       if(myid1 .eq. 0 .or. myid1 .eq. nproc1-1) then
         if(niib .eq. 3 .or. noib .eq. 3) then
           call getdelsib
         endif
       endif
       if(myid2 .eq. 0 .or. myid2 .eq. nproc2-1) then
         if(nijb .eq. 3 .or. nojb .eq. 3) then
           call getdelsjb
         endif
       endif
       if(myid3 .eq. 0 .or. myid3 .eq. nproc3-1) then
         if(nikb .eq. 3 .or. nokb .eq. 3) then
           call getdelskb
         endif
       endif

       do k=ks,ke
         do j=js,je
           do i=is,iep1
             dels1(i,j,k) = dels1(i,j,k)/(1.D0+srd1(i,j,k)) &
               +fh1(i,j,k)
           enddo
         enddo
       enddo
       do k=ks,ke
         do j=js,jep1
           do i=is,ie
             dels2(i,j,k) = dels2(i,j,k)/(1.D0+srd2(i,j,k)) &
               +fh2(i,j,k)
           enddo
         enddo
       enddo
       do k=ks,kep1
         do j=js,je
           do i=is,ie
             dels3(i,j,k) = dels3(i,j,k)/(1.D0+srd3(i,j,k)) &
               +fh3(i,j,k)
           enddo
         enddo
       enddo
!
!      checking dels
!
!      write(idcpu,'(i3.3)') myid
!      open(unit=13, file='/ptmp/yfan/hinode/dels_lorentz.cpu'//idcpu,
!    &  form='unformatted')
!      write(13) in,jn,kn
!      write(13) (((dels1(i,j,k),i=is,ie),j=js,je),k=ks,ke)
!      write(13) (((dels2(i,j,k),i=is,ie),j=js,je),k=ks,ke)
!      write(13) (((dels3(i,j,k),i=is,ie),j=js,je),k=ks,ke)
!      close(13)
!

end subroutine lorentz
#else
subroutine lorentz
  use ModPar
  use ModGrid
  use Modfield
  use ModDel
  use ModWork
  use ModBoundary
  use ModSundry,        ONLY: ovc2
  use ModBval,          ONLY: bvaleta
  use ModUserSetup,     ONLY: getdelsib, getdelsjb, getdelskb
  implicit none

  integer ::       i       , j       , k
  real :: q1
  real, dimension(1:ijkn) :: scra1, scra2, scra3
  real, dimension(in,jn,kn) :: bdotdels, bb, bdotv, fh1, fh2, fh3
!
!-----------------------------------------------------------------------
!
       do k=ksm1,kep1
         do j=jsm1,jep1
           do i=ism1,iep1
             srd1(i,j,k) = ( (b1(i,j,k))**2 &
               +(0.25D0*(b2(i,j,k)+b2(i,j+1,k) &
               +b2(i-1,j,k)+b2(i-1,j+1,k)))**2 &
               +(0.25D0*(b3(i,j,k)+b3(i,j,k+1) &
               +b3(i-1,j,k)+b3(i-1,j,k+1)))**2 ) &
               /(0.5D0*(d(i,j,k)+d(i-1,j,k)))*ovc2
             srd2(i,j,k) = ( (0.25D0*(b1(i,j,k)+b1(i+1,j,k) &
               +b1(i,j-1,k)+b1(i+1,j-1,k)))**2 &
               +(b2(i,j,k))**2 &
               +(0.25D0*(b3(i,j,k)+b3(i,j,k+1) &
               +b3(i,j-1,k)+b3(i,j-1,k+1)))**2 ) &
               /(0.5D0*(d(i,j,k)+d(i,j-1,k)))*ovc2
             srd3(i,j,k) = ( (0.25D0*(b1(i,j,k)+b1(i+1,j,k) &
               +b1(i,j,k-1)+b1(i+1,j,k-1)))**2 &
               +(0.25D0*(b2(i,j,k)+b2(i,j+1,k) &
               +b2(i,j,k-1)+b2(i,j+1,k-1)))**2 &
               +(b3(i,j,k))**2 ) &
               /(0.5D0*(d(i,j,k)+d(i,j,k-1)))*ovc2
           enddo
         enddo
       enddo

       do k=ks,ke
         do j=js,je
           do i=is,ie
             bdotdels(i,j,k) = 0.5D0*( b1(i,j,k) &
               *dels1(i,j,k) &
               +b1(i+1,j,k)*dels1(i+1,j,k) ) &
               + 0.5D0*( b2(i,j,k)*dels2(i,j,k) &
               +b2(i,j+1,k)*dels2(i,j+1,k) ) &
               + 0.5D0*( b3(i,j,k)*dels3(i,j,k) &
               +b3(i,j,k+1)*dels3(i,j,k+1) )
             bdotv(i,j,k) = 0.5D0*( b1(i,j,k) &
               *v1(i,j,k) &
               +b1(i+1,j,k)*v1(i+1,j,k) ) &
               + 0.5D0*( b2(i,j,k)*v2(i,j,k) &
               +b2(i,j+1,k)*v2(i,j+1,k) ) &
               + 0.5D0*( b3(i,j,k)*v3(i,j,k) &
               +b3(i,j,k+1)*v3(i,j,k+1) )
             bb(i,j,k) = 0.25D0 * ( (b1(i,j,k) &
               +b1(i+1,j,k))**2 + (b2(i,j,k) &
               +b2(i,j+1,k))**2 + (b3(i,j,k) &
               +b3(i,j,k+1))**2 )
           enddo
         enddo
       enddo
       call bvaleta(bdotdels)
       call bvaleta(bdotv)
       call bvaleta(bb)
       call bvaleta(deld)

       do k=ks,ke
         do j=js,je
           do i=is,iep1
             fh1(i,j,k) = srd1(i,j,k)/(1.D0+srd1(i,j,k)) &
               *( (dels1(i,j,k)-0.5D0*(bdotdels(i,j,k)+bdotdels(i-1,j,k)) &
               /(0.5D0*(bb(i,j,k)+bb(i-1,j,k))+tiny)*b1(i,j,k)) &
               - 0.5D0*(deld(i,j,k)+deld(i-1,j,k)) &
               *(v1(i,j,k)-0.5D0*(bdotv(i,j,k)+bdotv(i-1,j,k)) &
               /(0.5D0*(bb(i,j,k)+bb(i-1,j,k))+tiny)*b1(i,j,k)) )
           enddo
         enddo
       enddo
       do k=ks,ke
         do j=js,jep1
           do i=is,ie
             fh2(i,j,k) = srd2(i,j,k)/(1.D0+srd2(i,j,k)) &
               *( (dels2(i,j,k)-0.5D0*(bdotdels(i,j,k)+bdotdels(i,j-1,k)) &
               /(0.5D0*(bb(i,j,k)+bb(i,j-1,k))+tiny)*b2(i,j,k)) &
               - 0.5D0*(deld(i,j,k)+deld(i,j-1,k)) &
               *(v2(i,j,k)-0.5D0*(bdotv(i,j,k)+bdotv(i,j-1,k)) &
               /(0.5D0*(bb(i,j,k)+bb(i,j-1,k))+tiny)*b2(i,j,k)) )
           enddo
         enddo
       enddo
       do k=ks,kep1
         do j=js,je
           do i=is,ie
             fh3(i,j,k) = srd3(i,j,k)/(1.D0+srd3(i,j,k)) &
               *( (dels3(i,j,k)-0.5D0*(bdotdels(i,j,k)+bdotdels(i,j,k-1)) &
               /(0.5D0*(bb(i,j,k)+bb(i,j,k-1))+tiny)*b3(i,j,k)) &
               - 0.5D0*(deld(i,j,k)+deld(i,j,k-1)) &
               *(v3(i,j,k)-0.5D0*(bdotv(i,j,k)+bdotv(i,j,k-1)) &
               /(0.5D0*(bb(i,j,k)+bb(i,j,k-1))+tiny)*b3(i,j,k)) )
           enddo
         enddo
       enddo
       do k=ks,ke
         do j=js,je
           do i=is,iep1
             fh1(i,j,k)=dels1(i,j,k)-fh1(i,j,k)
           enddo
         enddo
       enddo
       do k=ks,ke
         do j=js,jep1
           do i=is,ie
             fh2(i,j,k)=dels2(i,j,k)-fh2(i,j,k)
           enddo
         enddo
       enddo
       do k=ks,kep1
         do j=js,je
           do i=is,ie
             fh3(i,j,k)=dels3(i,j,k)-fh3(i,j,k)
           enddo
         enddo
       enddo

    !----------------------------------------------------------------------------------

    do k=ks,kep1; do i=is,iep1; do j=js,jep1
       scratch1(i,j,k) = 0.25d0*(b1(i,j-1,k) + b1(i,j,k))*(b2(i,j,k) + b2(i-1,j,k))
       scratch2(i,j,k) = 0.25d0*(b3(i,j,k) + b3(i-1,j,k))*(b1(i,j,k-1) + b1(i,j,k))
       scratch3(i,j,k) = 0.25d0*(b3(i,j,k) + b3(i,j-1,k))*(b2(i,j,k-1) + b2(i,j,k))
    end do; enddo; enddo

    ! 1-force
    do k=ks,ke; do j=js,je; do i=is,iep1
       q1 = g31ai(i)*g32bi(j)
       dels1(i,j,k) = (g32a(j+1)*scratch1(i,j+1,k) - &
            g32a(j)*scratch1(i,j,k))*g2ai(i)*g32bi(j)*dx2ai(j) + &
            (scratch2(i,j,k+1) - scratch2(i,j,k))*q1*dx3ai(k)
    end do; end do; end do
    do k=ks,ke; do j=js,je
       do i=ism1,iep1
          scra1(i) = g2b(i)**2*(b2(i,j,k)**2 + b2(i,j+1,k)**2)
          scra2(i) = g31b(i)**2*(b3(i,j,k)**2 + b3(i,j,k+1)**2)
          scra3(i) = (g2b(i)*g31b(i)*(b1(i+1,j,k) + b1(i,j,k)))**2*0.5d0
       enddo
       do i=is,iep1
          dels1(i,j,k) = dels1(i,j,k) - 0.25d0*dx1bi(i)*( &
               g2ai(i)**2*(scra1(i) - scra1(i-1)) + &
               g31ai(i)**2*(scra2(i) -scra2(i-1)) - &
               (g2ai(i)*g31ai(i))**2*(scra3(i) - scra3(i-1)))
       enddo
    enddo; enddo

    ! 2-force
    do k=ks,ke; do j=js,jep1; do i=is,ie
       q1 = g31bi(i)*g32ai(j)
       dels2(i,j,k) = (scratch3(i,j,k+1) - scratch3(i,j,k))*q1*dx3ai(k) + &
            (g2a(i+1)*g31a(i+1)*scratch1(i+1,j,k)*g2a(i+1) - &
            g2a(i)*g2a(i)*g31a(i)*scratch1(i,j,k))*g2bi(i)*g2bi(i)*g31bi(i)*dx1ai(i)
    end do; end do; end do
    do k=ks,ke; do i=is,ie
       do j=jsm1,jep1
          scra1(j) = g32b(j)**2*(b3(i,j,k)**2 + b3(i,j,k+1)**2)
          scra2(j) = b1(i,j,k)**2 + b1(i+1,j,k)**2
          scra3(j) = (g32b(j)*(b2(i,j,k) + b2(i,j+1,k)))**2*0.5d0
       enddo
       do j=js,jep1
          dels2(i,j,k) = dels2(i,j,k) - 0.25d0*dx2bi(j)*g2bi(i)*( &
               g32ai(j)**2*(scra1(j) - scra1(j-1)) + &
               (scra2(j) - scra2(j-1)) - &
               g32ai(j)**2*(scra3(j) - scra3(j-1)))
       enddo
    enddo; enddo

    ! 3-force
    do k=ks,kep1; do j=js,je; do i=is,ie
       dels3(i,j,k) = (g2a(i+1)*g31a(i+1)*scratch2(i+1,j,k)*g31a(i+1) - &
            g31a(i)*g2a(i)*g31a(i)*scratch2(i,j,k))*g2bi(i)*g31bi(i)*g31bi(i)*dx1ai(i) + &
            (g32a(j+1)*scratch3(i,j+1,k)*g32a(j+1) - g32a(j)*g32a(j)*scratch3(i,j,k)) &
            *g32bi(j)*g32bi(j)*g2bi(i)*dx2ai(j)
    end do; end do; end do
    do j=js,je; do i=is,ie
       do k=ksm1,kep1
          scra1(k) = b1(i,j,k)**2 + b1(i+1,j,k)**2
          scra2(k) = b2(i,j,k)**2 + b2(i,j+1,k)**2
          scra3(k) = (b3(i,j,k) + b3(i,j,k+1))**2*0.5d0
       enddo
       do k=ks,kep1
          dels3(i,j,k) = dels3(i,j,k) - 0.25d0*dx3bi(k)*g31bi(i)*g32bi(j)*( &
               scra1(k) - scra1(k-1) + scra2(k) - scra2(k-1) - scra3(k)+scra3(k-1))
       enddo
    enddo; enddo
!
!      special boundary condition for line-tying
!
       if(myid1 .eq. 0 .or. myid1 .eq. nproc1-1) then
         if(niib .eq. 3 .or. noib .eq. 3) then
           call getdelsib
         endif
       endif
       if(myid2 .eq. 0 .or. myid2 .eq. nproc2-1) then
         if(nijb .eq. 3 .or. nojb .eq. 3) then
           call getdelsjb
         endif
       endif
       if(myid3 .eq. 0 .or. myid3 .eq. nproc3-1) then
         if(nikb .eq. 3 .or. nokb .eq. 3) then
           call getdelskb
         endif
       endif

       do k=ks,ke
         do j=js,je
           do i=is,iep1
             dels1(i,j,k) = dels1(i,j,k)/(1.D0+srd1(i,j,k)) &
               +fh1(i,j,k)
           enddo
         enddo
       enddo
       do k=ks,ke
         do j=js,jep1
           do i=is,ie
             dels2(i,j,k) = dels2(i,j,k)/(1.D0+srd2(i,j,k)) &
               +fh2(i,j,k)
           enddo
         enddo
       enddo
       do k=ks,kep1
         do j=js,je
           do i=is,ie
             dels3(i,j,k) = dels3(i,j,k)/(1.D0+srd3(i,j,k)) &
               +fh3(i,j,k)
           enddo
         enddo
       enddo

end subroutine lorentz
#endif

subroutine stv1
  use ModPar
  use ModPhysics
  use ModGrid
  use ModField
  use ModDel
  use ModUserSetup,       ONLY: etasqrt
  use ModSundry,          ONLY: ctime
  implicit none

  integer ::      i       , j       , k
  real :: rho, rs
!
!-----------------------------------------------------------------------
!
       rs=(etasqrt*ctime)**(7.D0/6.D0)
!
!      Add source terms to "dels1".
!
       do k=ks,ke
         do j=js,je
           do i=is,iep1
             rho = 0.5 * ( d(i-1,j,k) + d(i,j,k) )
!
!  1.  pressure gradient
!
             dels1(i,j,k)  = dels1(i,j,k) &
#ifndef ISOMHD
                     - (gamma-1.)*(eint(i,j,k) - eint(i-1,j,k)) * dx1bi(i)
#else
                     - cs2 * ( d (i,j,k) - d (i-1,j,k) ) * dx1bi(i)
#endif
#ifdef RTPGR
!
!  2.  gravity
!
             if(x1a(i) .le. rs) then
             dels1(i,j,k) = dels1(i,j,k) &
                     - rho*gacc/x1a(i)**2
             endif
#endif
#ifdef XYZGX
!
!  2.  gravity
!
             dels1(i,j,k) = dels1(i,j,k) - rho*gacc
#endif
!
!  3.  rotational pseudo-forces
!
             dels1(i,j,k)  = dels1(i,j,k) &
                  + ( v2(i  ,j,k) + v2(i  ,j+1,k  ) &
                    + v2(i-1,j,k) + v2(i-1,j+1,k  ) )**2 &
                  /16. * dg2bd1(i) * g2ai(i) &
                  * rho
             dels1(i,j,k)  = dels1(i,j,k) &
                  + ( v3(i  ,j,k) + v3(i  ,j  ,k+1) &
                    + v3(i-1,j,k) + v3(i-1,j  ,k+1) )**2 &
                  /16. * dg31bd1(i) * g31ai(i) &
                  * rho
!
           enddo
         enddo
       enddo

end subroutine stv1

subroutine stv2
  use ModPar
  use ModPhysics
  use ModGrid
  use ModField
  use ModDel
  implicit none

  integer ::      i       , j       , k
  real :: rho
!
!-----------------------------------------------------------------------
!
!      Add source terms to "dels2".
!
       do k=ks,ke
         do j=js,jep1
           do i=is,ie
             rho = 0.5 * ( d(i,j-1,k) + d(i,j,k) )
!
!  1.  pressure gradient
!
             dels2(i,j,k)  = dels2(i,j,k) &
#ifndef ISOMHD
                     - (gamma-1.)*(eint(i,j,k) - eint(i,j-1,k)) &
#else
                     - cs2 * ( d (i,j,k) - d (i,j-1,k) ) &
#endif
                     * dx2bi(j) * g2bi(i)
#ifdef XYZGY
!
!  2.  gravity
!
             dels2(i,j,k) = dels2(i,j,k) - rho*gacc
#endif
!
!  3.  rotational pseudo-force
!
             dels2(i,j,k)  = dels2(i,j,k) &
                  + ( v3(i,j  ,k) + v3(i,j  ,k+1) &
                    + v3(i,j-1,k) + v3(i,j-1,k+1) )**2 &
                  /16. * g2bi(i) * dg32bd2(j) * g32ai(j) &
                  * rho
!
           enddo
         enddo
       enddo

end subroutine stv2

subroutine stv3
  use ModPar
  use ModPhysics
  use ModGrid
  use ModField
  use ModDel
  implicit none

  integer ::       i       , j       , k
  real :: rho
!
!-----------------------------------------------------------------------
!
!      Add source terms to dels3.
!
       do k=ks,kep1
         do j=js,je
           do i=is,ie
             rho = 0.5 * ( d(i,j,k-1) + d(i,j,k) )
!
!  1.  pressure gradient
!
             dels3(i,j,k)  = dels3(i,j,k) &
#ifndef ISOMHD
                     - (gamma-1.)*(eint(i,j,k) - eint(i,j,k-1)) &
#else
                     - cs2 * ( d (i,j,k) - d (i,j,k-1) ) &
#endif
                     * dx3bi(k) * g31bi(i) * g32bi(j)
#ifdef XYZGZ
!
!  2.  gravity
!
             dels3(i,j,k)  = dels3(i,j,k) - rho*gacc
#endif
!
           enddo
         enddo
       enddo

end subroutine stv3

subroutine pdv
#ifndef ISOMHD
  use ModPar
  use ModPhysics
  use ModGrid
  use ModField
  use ModDel,           ONLY:delei
  use ModWork,          ONLY:scratch1
  implicit none

  integer ::       i       , j       , k
  real :: q1
!
!-----------------------------------------------------------------------
!
!      Compute divergence of velocity field.
!
       do k=ks,ke
         do j=js,je
           do i=is,ie
             scratch1(i,j,k) = ( g2a(i+1) * g31a(i+1) * v1(i+1,j,k) &
                        - g2a(i  ) * g31a(i  ) * v1(i  ,j,k) ) &
                        * dx1ai(i)*g2bi(i)*g31bi(i) &
                        + ( g32a(j+1) * v2(i,j+1,k) &
                        - g32a(j  ) * v2(i,j  ,k) ) &
                        * dx2ai(j) * g2bi(i) * g32bi(j) &
                        + ( v3(i,j,k+1) - v3(i,j,k) ) &
                        * dx3ai(k) * g31bi(i) * g32bi(j)
           enddo
         enddo
       enddo
!
!      Add pdv source term to energy density
!
       q1 = (gamma-1.)
       do k=ks,ke
         do j=js,je
           do i=is,ie
             delei(i,j,k) = delei(i,j,k) &
               - q1 * eint(i,j,k) * scratch1(i,j,k)
           enddo
         enddo
       enddo

#endif
end subroutine pdv

subroutine advection
  use ModPar
  use ModPhysics
  use ModGrid
  use ModField
  use ModSundry
  use ModDel
  use ModWork
  use ModInterp
  implicit none

  integer :: i, j, k, ierr
  real :: cmax, rr, qq, q1
  real, dimension(1:ijkn) :: flux, qint, ql, qr, qlr, dq
  real, dimension(in,jn,kn) :: cmaxc, cmax1, cmax2, cmax3

!
!-----------------------------------------------------------------------
!
!      For the purposes of consistent advection, construct the mass
!  flux across each interface.  The mass flux will be used to create
!  the fluxes of all variables, including the momenta


       call x1zc3d_den ( d  , jsm1, ksm1, jep1, kep1, &
                         scratch2, scratch3)

       do k=ksm1,kep1
         do j=jsm1,jep1
           do i=ism1,iep2
#ifndef ISOMHD
             cs2=gamma*(gamma-1.)*0.5*(eint(i-1,j,k)+eint(i,j,k)) &
                /(0.5*(d(i-1,j,k)+d(i,j,k)))
#endif
             cmax= &
               sqrt( min( clight**2, ( (b1(i,j,k))**2 &
               +(0.25*(b2(i,j,k)+b2(i,j+1,k) &
               +b2(i-1,j,k)+b2(i-1,j+1,k)))**2 &
               +(0.25*(b3(i,j,k)+b3(i,j,k+1) &
               +b3(i-1,j,k)+b3(i-1,j,k+1)))**2 ) &
               /(0.5*(d(i,j,k)+d(i-1,j,k))) ) &
               + cs2 )
             rr=abs(scratch3(i,j,k)-scratch2(i,j,k)) &
                /(abs(d(i,j,k)-d(i-1,j,k))+tiny)
             q1=(min(1.,alpha_d*rr))**nlf_d
             srd1(i,j,k) = 0.5*(scratch3(i,j,k)+scratch2(i,j,k)) &
               *v1(i,j,k) - 0.5*(abs(v1(i,j,k))+cmax*q1) &
               *(scratch3(i,j,k)-scratch2(i,j,k))
             cmax1(i,j,k) = cmax
           enddo
         enddo
       enddo

       call x2zc3d_den ( d  , ksm1, ism1, kep1, iep1, &
                         scratch2, scratch3)

       do k=ksm1,kep1
         do j=jsm1,jep2
           do i=ism1,iep1
#ifndef ISOMHD
             cs2=gamma*(gamma-1.)*0.5*(eint(i,j-1,k)+eint(i,j,k)) &
                /(0.5*(d(i,j-1,k)+d(i,j,k)))
#endif
             cmax = &
               sqrt( min( clight**2, ( (0.25*(b1(i,j,k)+b1(i+1,j,k) &
               +b1(i,j-1,k)+b1(i+1,j-1,k)))**2 &
               +(b2(i,j,k))**2 &
               +(0.25*(b3(i,j,k)+b3(i,j,k+1) &
               +b3(i,j-1,k)+b3(i,j-1,k+1)))**2 ) &
               /(0.5*(d(i,j,k)+d(i,j-1,k))) ) &
               + cs2 )
             rr=abs(scratch3(i,j,k)-scratch2(i,j,k)) &
                /(abs(d(i,j,k)-d(i,j-1,k))+tiny)
             q1=(min(1.,alpha_d*rr))**nlf_d
             srd2(i,j,k) = 0.5*(scratch3(i,j,k)+scratch2(i,j,k)) &
               *v2(i,j,k) - 0.5*(abs(v2(i,j,k))+cmax*q1) &
               *(scratch3(i,j,k)-scratch2(i,j,k))
             cmax2(i,j,k) = cmax
           enddo
         enddo
       enddo

       call x3zc3d_den ( d  , ism1, jsm1, iep1, jep1, &
                         scratch2, scratch3)

       do k=ksm1,kep2
         do j=jsm1,jep1
           do i=ism1,iep1
#ifndef ISOMHD
             cs2=gamma*(gamma-1.)*0.5*(eint(i,j,k-1)+eint(i,j,k)) &
                /(0.5*(d(i,j,k-1)+d(i,j,k)))
#endif
             cmax = &
               sqrt( min( clight**2, ( (0.25*(b1(i,j,k)+b1(i+1,j,k) &
               +b1(i,j,k-1)+b1(i+1,j,k-1)))**2 &
               +(0.25*(b2(i,j,k)+b2(i,j+1,k) &
               +b2(i,j,k-1)+b2(i,j+1,k-1)))**2 &
               +(b3(i,j,k))**2 ) &
               /(0.5*(d(i,j,k)+d(i,j,k-1))) ) &
               + cs2 )
             rr=abs(scratch3(i,j,k)-scratch2(i,j,k)) &
                /(abs(d(i,j,k)-d(i,j,k-1))+tiny)
             q1=(min(1.,alpha_d*rr))**nlf_d
             srd3(i,j,k) = 0.5*(scratch3(i,j,k)+scratch2(i,j,k)) &
               *v3(i,j,k) - 0.5*(abs(v3(i,j,k))+cmax*q1) &
               *(scratch3(i,j,k)-scratch2(i,j,k))
             cmax3(i,j,k) = cmax
           enddo
         enddo
       enddo
!
!
!      Perform mass density advection.
!
       do k=ks,ke
         do j=js,je
           do i=is,ie
             deld(i,j,k) = - (g2a(i+1)*g31a(i+1)*srd1(i+1,j,k) &
               - g2a(i)*g31a(i)*srd1(i,j,k)) &
               *dx1ai(i)*g2bi(i)*g31bi(i) &
               - (g32a(j+1)*srd2(i,j+1,k) &
               - g32a(j)*srd2(i,j,k)) &
               *dx2ai(j)*g2bi(i)*g32bi(j) &
               - (srd3(i,j,k+1) - srd3(i,j,k)) &
               *dx3ai(k)*g31bi(i)*g32bi(j)
           enddo
         enddo
       enddo

#ifndef ISOMHD
!
!      Compute specific internal energy for adiabatic evolution
!
       do k=ks,ke; do j=js,je
         qint(ism2:iep2) = eint(ism2:iep2,j,k) / d(ism2:iep2,j,k)
         call x1tvd(qint,ql,qr,dq,qlr)
         do i = is,iep1
           cmax = cmax1(i,j,k)
           rr = abs(dq(i))/(abs(qint(i)-qint(i-1)) + tiny)
           qq=(min(1.,alpha_e*rr))**nlf_e
           flux(i) = (srd1(i,j,k)*qlr(i) - 0.5*(abs(srd1(i,j,k)) &
             + cmax*qq*0.5*(d(i-1,j,k)+d(i,j,k))) &
             *dq(i))*g2a(i)*g31a(i)
         enddo
         do i=is,ie
           delei(i,j,k) = delei(i,j,k) &
             - (flux(i+1) - flux(i))*dx1ai(i)*g2bi(i)*g31bi(i)
         enddo
       enddo; enddo
       do k=ks,ke; do i=is,ie
         qint(jsm2:jep2) = eint(i,jsm2:jep2,k) / d(i,jsm2:jep2,k)
         call x2tvd (qint,ql,qr,dq,qlr)
         do j=js,jep1
           cmax = cmax2(i,j,k)
           rr = abs(dq(j))/(abs(qint(j)-qint(j-1)) + tiny)
           qq = (min(1., alpha_e*rr))**nlf_e
           flux(j) = (srd2(i,j,k)*qlr(j) - 0.5*(abs(srd2(i,j,k)) &
             + cmax*qq*0.5*(d(i,j-1,k)+d(i,j,k)))*dq(j))*g32a(j)
         enddo
         do j=js,je
           delei(i,j,k) = delei(i,j,k) &
             - (flux(j+1) - flux(j))*dx2ai(j)*g2bi(i)*g32bi(j)
         enddo
       enddo; enddo
       do j=js,je; do i=is,ie
         qint(ksm2:kep2) = eint(i,j,ksm2:kep2) / d(i,j,ksm2:kep2)
         call x3tvd (qint,ql,qr,dq,qlr)
         do k=ks,kep1
           cmax = cmax3(i,j,k)
           rr = abs(dq(k))/(abs(qint(k)-qint(k-1)) + tiny)
           qq = (min(1., alpha_e*rr))**nlf_e
           flux(k) = (srd3(i,j,k)*qlr(k) - 0.5*(abs(srd3(i,j,k)) &
             + cmax*qq*0.5*(d(i,j,k-1)+d(i,j,k)))*dq(k))
         enddo
         do k=ks,ke
            delei(i,j,k) = delei(i,j,k) &
              - (flux(k+1) - flux(k)) &
              *dx3ai(k)*g31bi(i)*g32bi(j)
         enddo
       enddo; enddo
#endif
!
!      advection of momenta
!
       do k=ksm1,kep1; do j=jsm1,jep1; do i=ism1,iep1
#ifndef ISOMHD
             cs2=gamma*(gamma-1.)*eint(i,j,k)/d(i,j,k)
#endif
!            cmax at cell center
             cmaxc(i,j,k) =  sqrt( min( clight**2, ( &
                 (0.5*(b1(i,j,k)+b1(i+1,j,k)))**2 &
               + (0.5*(b2(i,j,k)+b2(i,j+1,k)))**2 &
               + (0.5*(b3(i,j,k)+b3(i,j,k+1)))**2) &
               /(d(i,j,k)) ) &
               + cs2 )
#ifndef ISOMHD
             cs2=gamma*(gamma-1.) &
                *0.25*(eint(i,j,k)+eint(i-1,j,k)+eint(i-1,j-1,k)+eint(i,j-1,k)) &
                /(0.25*(d(i,j,k)+d(i-1,j,k)+d(i-1,j-1,k)+d(i,j-1,k)))
#endif
!            cmax at 3-edge
             cmax3(i,j,k) = &
               sqrt( min( clight**2, ( (0.5*(b1(i,j-1,k)+b1(i,j,k)))**2 &
               + (0.5*(b2(i,j,k)+b2(i-1,j,k)))**2 &
               + (0.125*(b3(i,j,k)+b3(i,j,k+1) &
                        +b3(i,j-1,k)+b3(i,j-1,k+1) &
                        +b3(i-1,j-1,k)+b3(i-1,j-1,k+1) &
                        +b3(i-1,j,k)+b3(i-1,j,k+1)))**2 ) &
               /(0.25*(d(i,j,k)+d(i-1,j,k) &
                       +d(i-1,j-1,k)+d(i,j-1,k))) ) &
               + cs2 )
#ifndef ISOMHD
             cs2=gamma*(gamma-1.) &
                *0.25*(eint(i,j,k)+eint(i-1,j,k)+eint(i-1,j,k-1)+eint(i,j,k-1)) &
                /(0.25*(d(i,j,k)+d(i-1,j,k)+d(i-1,j,k-1)+d(i,j,k-1)))
#endif
!            cmax at 2-edge
             cmax2(i,j,k) = &
               sqrt( min( clight**2, ( (0.5*(b1(i,j,k-1)+b1(i,j,k)))**2 &
               + (0.5*(b3(i,j,k)+b3(i-1,j,k)))**2 &
               + (0.125*(b2(i,j,k)+b2(i,j+1,k) &
                        +b2(i,j,k-1)+b2(i,j+1,k-1) &
                        +b2(i-1,j,k-1)+b2(i-1,j+1,k-1) &
                        +b2(i-1,j,k)+b2(i-1,j+1,k)))**2 ) &
               /(0.25*(d(i,j,k)+d(i-1,j,k) &
                       +d(i-1,j,k-1)+d(i,j,k-1))) ) &
               + cs2 )
#ifndef ISOMHD
             cs2=gamma*(gamma-1.) &
                *0.25*(eint(i,j,k)+eint(i,j-1,k)+eint(i,j-1,k-1)+eint(i,j,k-1)) &
                /(0.25*(d(i,j,k)+d(i,j-1,k)+d(i,j-1,k-1)+d(i,j,k-1)))
#endif
!            cmax at 1-edge
             cmax1(i,j,k) = &
               sqrt( min( clight**2, ( (0.5*(b3(i,j,k)+b3(i,j-1,k)))**2 &
               + (0.5*(b2(i,j,k-1)+b2(i,j,k)))**2 &
               + (0.125*(b1(i,j,k)+b1(i+1,j,k) &
                        +b1(i,j,k-1)+b1(i+1,j,k-1) &
                        +b1(i,j-1,k-1)+b1(i+1,j-1,k-1) &
                        +b1(i,j-1,k)+b1(i+1,j-1,k)))**2) &
               /(0.25*(d(i,j,k)+d(i,j-1,k) &
                       +d(i,j-1,k-1)+d(i,j,k-1))) ) &
               + cs2 )
       enddo; enddo; enddo
!      advection for dels1
       do k=ks,ke; do j=js,je
         qint(ism2:iep3) = v1(ism2:iep3,j,k)
         call x1tvd_a(qint,ql,qr,dq,qlr)
         do i=ism1,iep1
           cmax=cmaxc(i,j,k)
           rr=abs(dq(i))/(abs(qint(i+1)-qint(i)) + tiny)
           qq=(min(1.,alpha_v*rr))**nlf_v
           q1=0.5*abs(0.5*(srd1(i,j,k)+srd1(i+1,j,k)))*dq(i) &
             +0.5*cmax*qq*d(i,j,k)*dq(i)
           flux(i)=0.5*(srd1(i,j,k)+srd1(i+1,j,k))*qlr(i) - q1
           scratch1(i,j,k)=q1*(qint(i+1)-qint(i))*dx1ai(i)
         enddo
         do i=is,iep1
           dels1(i,j,k) = dels1(i,j,k) &
             - ( g2b(i)*g31b(i)*flux(i) &
                - g2b(i-1)*g31b(i-1)*flux(i-1) ) &
             * dx1bi(i)*g2ai(i)*g31ai(i)
         enddo
       enddo; enddo
       do k=ks,ke; do i=is,iep1
         qint(jsm2:jep2) = v1(i,jsm2:jep2,k)
         call x2tvd(qint,ql,qr,dq,qlr)
         do j=js,jep1
           cmax = cmax3(i,j,k)
           rr = abs(dq(j))/(abs(qint(j)-qint(j-1)) + tiny)
           qq=(min(1.,alpha_v*rr))**nlf_v
           q1=0.5*abs(0.5*(srd2(i-1,j,k)+srd2(i,j,k)))*dq(j) &
             +0.5*cmax*qq*0.25*(d(i,j,k)+d(i-1,j,k) &
             +d(i-1,j-1,k)+d(i,j-1,k))*dq(j)
           flux(j) = 0.5*(srd2(i-1,j,k)+srd2(i,j,k))*qlr(j) - q1
           scratch2(i,j,k)=q1*(qint(j)-qint(j-1)) &
             *g2ai(i)*dx2bi(j)
         enddo
         do j=js,je
           dels1(i,j,k) = dels1(i,j,k) &
             - (g32a(j+1)*flux(j+1) - g32a(j)*flux(j)) &
               * dx2ai(j) * g2ai(i) * g32bi(j)
         enddo
       enddo; enddo
       do j=js,je; do i=is,iep1
         qint(ksm2:kep2) = v1(i,j,ksm2:kep2)
         call x3tvd(qint,ql,qr,dq,qlr)
         do k=ks,kep1
           cmax = cmax2(i,j,k)
           rr = abs(dq(k))/(abs(qint(k)-qint(k-1)) + tiny)
           qq=(min(1.,alpha_v*rr))**nlf_v
           q1=0.5*abs(0.5*(srd3(i-1,j,k)+srd3(i,j,k)))*dq(k) &
             +0.5*cmax*qq*0.25*(d(i,j,k)+d(i-1,j,k) &
             +d(i-1,j,k-1)+d(i,j,k-1))*dq(k)
           flux(k) = 0.5*(srd3(i-1,j,k)+srd3(i,j,k))*qlr(k) - q1
           scratch3(i,j,k)=q1*(qint(k)-qint(k-1)) &
             *g31ai(i)*g32bi(j)*dx3bi(k)
         enddo
         do k=ks,ke
           dels1(i,j,k) = dels1(i,j,k) &
             - ( flux(k+1) - flux(k) ) &
             * dx3ai(k) * g31ai(i) * g32bi(j)
         enddo
       enddo; enddo
#ifndef ISOMHD
!      add viscous heating due to the upwinding
       do k=ks,ke
         do j=js,je
           do i=is,ie
             delei(i,j,k)=delei(i,j,k) &
                + scratch1(i,j,k) &
                  +(scratch2(i,j,k)+scratch2(i+1,j,k) &
                  +scratch2(i,j+1,k)+scratch2(i+1,j+1,k))*0.25 &
                  +(scratch3(i,j,k)+scratch3(i+1,j,k) &
                  +scratch3(i,j,k+1)+scratch3(i+1,j,k+1))*0.25
           enddo
         enddo
       enddo
#endif
!      advection for dels2
       do k=ks,ke; do j=js,jep1
         qint(ism2:iep2) = v2(ism2:iep2,j,k)*g2bi(ism2:iep2)
         call x1tvd(qint,ql,qr,dq,qlr)
         do i=is,iep1
           cmax = cmax3(i,j,k)
           rr = abs(dq(i))/(abs(qint(i)-qint(i-1)) + tiny)
           qq=(min(1.,alpha_v*rr))**nlf_v
           q1=0.5*abs(0.5*(srd1(i,j-1,k)+srd1(i,j,k)))*dq(i)*(g2a(i))**2 &
             +(0.5*cmax*qq*0.25*(d(i,j,k)+d(i-1,j,k) &
             +d(i-1,j-1,k)+d(i,j-1,k))*dq(i))*(g2a(i))**2
           flux(i) = 0.5*(srd1(i,j-1,k)+srd1(i,j,k))*qlr(i)*(g2a(i))**2 - q1
           scratch1(i,j,k)=q1*(qint(i)-qint(i-1))*dx1bi(i)
         enddo
         do i=is,ie
           dels2(i,j,k) = dels2(i,j,k) &
             - (g2a(i+1)*g31a(i+1)*flux(i+1) &
             - g2a(i)*g31a(i)*flux(i)) &
             * dx1ai(i) * g2bi(i) * g31bi(i) * g2bi(i)
         enddo
       enddo; enddo
       do k=ks,ke; do i=is,ie
         qint(jsm2:jep3) = v2(i,jsm2:jep3,k)*g2bi(i)
         call x2tvd_a(qint,ql,qr,dq,qlr)
         do j=jsm1,jep1
           cmax = cmaxc(i,j,k)
           rr = abs(dq(j))/(abs(qint(j+1)-qint(j)) + tiny)
           qq=(min(1.,alpha_v*rr))**nlf_v
           q1=0.5*abs(0.5*(srd2(i,j,k)+srd2(i,j+1,k)))*dq(j) * (g2b(i))**2 &
             +0.5*cmax*qq*d(i,j,k)*dq(j) * (g2b(i))**2
           flux(j) = 0.5*(srd2(i,j,k)+srd2(i,j+1,k))*qlr(j)*(g2b(i))**2 - q1
           scratch2(i,j,k)=q1*(qint(j+1)-qint(j)) &
             *g2bi(i)*dx2ai(j)
         enddo
         do j=js,jep1
           dels2(i,j,k) = dels2(i,j,k) &
             - (g32b(j)*flux(j) - g32b(j-1)*flux(j-1)) &
               * dx2bi(j) * g2bi(i) * g32ai(j) * g2bi(i)
         enddo
       enddo; enddo
       do j=js,jep1; do i=is,ie
         qint(ksm2:kep2) = v2(i,j,ksm2:kep2)*g2bi(i)
         call x3tvd(qint,ql,qr,dq,qlr)
         do k=ks,kep1
           cmax = cmax1(i,j,k)
           rr = abs(dq(k))/(abs(qint(k)-qint(k-1)) + tiny)
           qq=(min(1.,alpha_v*rr))**nlf_v
           q1=0.5*abs(0.5*(srd3(i,j-1,k)+srd3(i,j,k)))*dq(k) * (g2b(i))**2 &
             +0.5*cmax*qq*0.25*(d(i,j,k)+d(i,j-1,k) &
                     +d(i,j-1,k-1)+d(i,j,k-1))*dq(k) * (g2b(i))**2
           flux(k) = 0.5*(srd3(i,j-1,k)+srd3(i,j,k))*qlr(k)*(g2b(i))**2 - q1
           scratch3(i,j,k)=q1*(qint(k)-qint(k-1)) &
             *g31bi(i)*g32ai(j)*dx3bi(k)
         enddo
         do k=ks,ke
           dels2(i,j,k) = dels2(i,j,k) &
             - (flux(k+1) - flux(k)) &
             * dx3ai(k) * g31bi(i) * g32ai(j) * g2bi(i)
         enddo
       enddo; enddo
#ifndef ISOMHD
!      add viscous heating due to the upwinding
       do k=ks,ke
         do j=js,je
           do i=is,ie
             delei(i,j,k)=delei(i,j,k) &
               +(scratch1(i,j,k)+scratch1(i,j+1,k) &
               +scratch1(i+1,j,k)+scratch1(i+1,j+1,k))*0.25 &
               +scratch2(i,j,k) &
               +(scratch3(i,j,k)+scratch3(i,j+1,k) &
               +scratch3(i,j,k+1)+scratch3(i,j+1,k+1))*0.25
           enddo
         enddo
       enddo
#endif
!      advection for dels3   
       do k=ks,kep1; do j=js,je
         qint(ism2:iep2) = v3(ism2:iep2,j,k)*g31bi(ism2:iep2)*g32bi(j)
         call x1tvd(qint,ql,qr,dq,qlr)
         do i = is,iep1
           cmax = cmax2(i,j,k)
           rr = abs(dq(i))/(abs(qint(i)-qint(i-1)) + tiny)
           qq=(min(1.,alpha_v*rr))**nlf_v
           q1=0.5*abs(0.5*(srd1(i,j,k-1)+srd1(i,j,k)))*dq(i) &
             *(g31a(i)*g32b(j))**2 &
             +(0.5*cmax*qq*0.25*(d(i,j,k)+d(i-1,j,k) &
                     +d(i-1,j,k-1)+d(i,j,k-1))*dq(i)) &
                     *(g31a(i)*g32b(j))**2
           flux(i) = 0.5*(srd1(i,j,k-1)+srd1(i,j,k))*qlr(i) &
             *(g31a(i)*g32b(j))**2 - q1
           scratch1(i,j,k)=q1*(qint(i)-qint(i-1)) &
             * dx1bi(i)
         enddo
         do i=is,ie
           dels3(i,j,k) = dels3(i,j,k) &
             - (g2a(i+1)*g31a(i+1)*flux(i+1) &
                - g2a(i)*g31a(i)*flux(i)) &
             * dx1ai(i) * g2bi(i) * g31bi(i) * g31bi(i) * g32bi(j)
         enddo
       enddo; enddo
       do k=ks,kep1; do i=is,ie
         qint(jsm2:jep2) = v3(i,jsm2:jep2,k)*g31bi(i)*g32bi(jsm2:jep2)
         call x2tvd(qint,ql,qr,dq,qlr)
         do j=js,jep1
           cmax = cmax1(i,j,k)
           rr = abs(dq(j))/(abs(qint(j)-qint(j-1)) + tiny)
           qq=(min(1.,alpha_v*rr))**nlf_v
           q1=0.5*abs(0.5*(srd2(i,j,k)+srd2(i,j,k-1)))*dq(j) &
             *(g31b(i)* g32a(j))**2 &
             +0.5*cmax*qq*0.25*(d(i,j,k)+d(i,j-1,k) &
                     +d(i,j-1,k-1)+d(i,j,k-1))*dq(j) &
                     *(g31b(i)*g32a(j))**2
           flux(j)=0.5*(srd2(i,j,k)+srd2(i,j,k-1))*qlr(j) &
             *(g31b(i)*g32a(j))**2 - q1
           scratch2(i,j,k)=q1*(qint(j)-qint(j-1)) &
             * g2bi(i) * dx2bi(j)
         enddo
         do j=js,je
           dels3(i,j,k) = dels3(i,j,k) &
             - (g32a(j+1)*flux(j+1) - g32a(j)*flux(j)) &
             * dx2ai(j) * g2bi(i) * g32bi(j) * g31bi(i) * g32bi(j)
         enddo
       enddo; enddo
       do j=js,je; do i=is,ie
         qint(ksm2:kep3) = v3(i,j,ksm2:kep3)*g31bi(i)*g32bi(j)
         call x3tvd_a(qint,ql,qr,dq,qlr)
         do k=ksm1,kep1
           cmax = cmaxc(i,j,k)
           rr = abs(dq(k))/(abs(qint(k+1)-qint(k)) + tiny)
           qq=(min(1.,alpha_v*rr))**nlf_v
           q1=0.5*abs(0.5*(srd3(i,j,k)+srd3(i,j,k+1)))*dq(k) &
             *(g31b(i)*g32b(j))**2 &
             +0.5*cmax*qq*d(i,j,k)*dq(k) &
             *(g31b(i)*g32b(j))**2
           flux(k) = 0.5*(srd3(i,j,k)+srd3(i,j,k+1))*qlr(k) &
             *(g31b(i)*g32b(j))**2 - q1
           scratch3(i,j,k)=q1*(qint(k+1)-qint(k)) &
             * g31bi(i) * g32bi(j) * dx3ai(k)
         enddo
         do k=ks,kep1
           dels3(i,j,k) = dels3(i,j,k) &
             - (flux(k) - flux(k-1)) &
             * dx3bi(k) * g31bi(i) * g32bi(j) * g31bi(i) * g32bi(j)
         enddo
       enddo; enddo
#ifndef ISOMHD
!      add viscous heating due to the upwinding
       do k=ks,ke
         do j=js,je
           do i=is,ie
             delei(i,j,k)=delei(i,j,k) &
               +(scratch1(i,j,k)+scratch1(i,j,k+1) &
               +scratch1(i+1,j,k)+scratch1(i+1,j,k+1))*0.25 &
               +(scratch2(i,j,k)+scratch2(i,j,k+1) &
               +scratch2(i,j+1,k)+scratch2(i,j+1,k+1))*0.25 &
               +scratch3(i,j,k)
           enddo
         enddo
       enddo
#endif

end subroutine advection

subroutine hsmoc
  use ModPar
  use ModGrid
  use ModField
  use ModWork
  use ModDel
  use ModBval
  implicit none

  integer :: i, j, k
!
!-----------------------------------------------------------------------
!
  do k=ks,kep1
     do j=js,jep1
        do i=is,ie
           scratch1(i,j,k)=0.25D0*(v2(i,j,k)+v2(i,j,k-1))*(b3(i,j,k)+b3(i,j-1,k))
           scratch2(i,j,k)=0.25D0*(b2(i,j,k)+b2(i,j,k-1))*(v3(i,j,k)+v3(i,j-1,k))
           srd1(i,j,k)=scratch1(i,j,k) - scratch2(i,j,k)
        enddo
     enddo
  enddo
!
!-----------------------------------------------------------------------
!---- 2.  emf2 ---------------------------------------------------------
!-----------------------------------------------------------------------
!
  do k=ks,kep1
     do j=js,je
        do i=is,iep1
           scratch1(i,j,k)=0.25D0*(v3(i,j,k)+v3(i-1,j,k))*(b1(i,j,k)+b1(i,j,k-1))
           scratch2(i,j,k)=0.25D0*(b3(i,j,k)+b3(i-1,j,k))*(v1(i,j,k)+v1(i,j,k-1))
           srd2(i,j,k)=scratch1(i,j,k) - scratch2(i,j,k)
        enddo
     enddo
  enddo
!
!-----------------------------------------------------------------------
!---- 3.  emf3 ---------------------------------------------------------
!-----------------------------------------------------------------------
!
  do k=ks,ke
     do j=js,jep1
        do i=is,iep1
           scratch1(i,j,k)=0.25D0*(v1(i,j,k)+v1(i,j-1,k))*(b2(i,j,k)+b2(i-1,j,k))
           scratch2(i,j,k)=0.25D0*(b1(i,j,k)+b1(i,j-1,k))*(v2(i,j,k)+v2(i-1,j,k))
           srd3(i,j,k)=scratch1(i,j,k) - scratch2(i,j,k)
        enddo
     enddo
  enddo
!
!-----------------------------------------------------------------------
!
!      Set boundary values for "term1" and "term2".
!
       call bvalemf1 (srd1)
!
!      Compute "emf1" for all 1-edges, including the ghost zones.
!
       do k=ksm2,kep3
         do j=jsm2,jep3
           do i=ism2,iep2
             emf1s(i,j,k) = emf1s(i,j,k) &
               + srd1(i,j,k)
           enddo
         enddo
       enddo
!
!-----------------------------------------------------------------------
!
!      Set boundary values for "term1" and "term2".
!
       call bvalemf2 (srd2)
!
!      Compute "emf2" for all 2-edges, including the ghost zones.
!
       do k=ksm2,kep3
         do j=jsm2,jep2
           do i=ism2,iep3
             emf2s(i,j,k) = emf2s(i,j,k) &
               + srd2(i,j,k)
           enddo
         enddo
       enddo
!
!-----------------------------------------------------------------------
!
!      Set boundary values for "term1" and "term2".
!
       call bvalemf3 (srd3)
!
!      Compute "emf3" for all 3-edges, including the ghost zones.
!
       do k=ksm2,kep2
         do j=jsm2,jep3
           do i=ism2,iep3
             emf3s(i,j,k) = emf3s(i,j,k) &
               + srd3(i,j,k)
           enddo
         enddo
       enddo

end subroutine hsmoc

subroutine diffemf
  use ModPar
  use ModGrid
  use ModField
  use ModWork
  use ModDel
  use ModBoundary,      ONLY: isplt
  use ModSundry,        ONLY: dtmdi2, nlf_d, alpha_d, dt, clight
  use ModBval
  use ModInterp,        ONLY: x1del_den, x2del_den, x3del_den

  implicit none

  include "mpif.h"

  integer :: i, j, k, ierr, nlf
  real, dimension(ijkn) :: bint, dbint
  real :: dr1, dr2, dr3, drmin, current, dtmdi2tmp, dtmdi2lc, &
          alpha, vp, fq, q1, vv1, vv2, vv3, etatot
  real, parameter :: ww=1.D-5

  real, dimension(in,jn,kn) :: emf1, emf2, emf3, eta, &
                               curl1, curl2, curl3
!
!-----------------------------------------------------------------------
!
       nlf=nlf_d
       alpha=alpha_d
!
!-----------------------------------------------------------------------
!---- 1.  emf1 ---------------------------------------------------------
!-----------------------------------------------------------------------
!
!      calculate del2(b3 h3 dx3)
!
       do k=ks,kep1
         do i=is,ie
           do j=jsm2,jep2
             bint(j)=b3(i,j,k)*g31b(i)*g32b(j)*dx3b(k)
           enddo
!          if js and jep1 is at x2a=0 and pi/2
!          if (isplt .eq. -1) then
!           if (myid2 .eq. 0) then
!            bint(jsm1)=-bint(js)
!            bint(jsm2)=-bint(jsp1)
!           endif
!           if (myid2 .eq. nproc2-1) then
!            bint(jep1)=-bint(je)
!            bint(jep2)=-bint(jem1)
!           endif
!          endif
           call x2del_den(bint,dbint)
           do j=js,jep1
             scratch1(i,j,k)=-dbint(j)
           enddo
         enddo
       enddo
!
!      calculate del3(b2 h2 dx2)
!
       do j=js,jep1
         do i=is,ie
           do k=ksm2,kep2
             bint(k)=b2(i,j,k)*g2b(i)*dx2b(j)
           enddo
           call x3del_den(bint,dbint)
           do k=ks,kep1
             scratch2(i,j,k)=-dbint(k)
           enddo
         enddo
       enddo
!
       do k=ks,kep1
         do j=js,jep1
           do i=is,ie
             emf1(i,j,k) = scratch1(i,j,k) - scratch2(i,j,k)
           enddo
         enddo
       enddo
!
!-----------------------------------------------------------------------
!---- 2.  emf2 ---------------------------------------------------------
!-----------------------------------------------------------------------
!
!      calculate del3(b1 h1 dx1)
!
       do j=js,je
         do i=is,iep1
           do k=ksm2,kep2
             bint(k)=b1(i,j,k)*dx1b(i)
           enddo
           call x3del_den(bint,dbint)
           do k=ks,kep1
             scratch1(i,j,k)=-dbint(k)
           enddo
         enddo
       enddo
!
!      calculate del1(b3 h3 dx3)
!
       do k=ks,kep1
         do j=js,je
           do i=ism2,iep2
             bint(i)=b3(i,j,k)*g31b(i)*g32b(j)*dx3b(k)
           enddo
           call x1del_den(bint,dbint)
           do i=is,iep1
             scratch2(i,j,k)=-dbint(i)
           enddo
         enddo
       enddo
!
       do k=ks,kep1
         do j=js,je
           do i=is,iep1
             emf2(i,j,k) = scratch1(i,j,k) - scratch2(i,j,k)
           enddo
         enddo
       enddo
!
!-----------------------------------------------------------------------
!---- 3.  emf3 ---------------------------------------------------------
!-----------------------------------------------------------------------
!
!      calculate del1(b2 h2 dx2)
!
       do k=ks,ke
         do j=js,jep1
           do i=ism2,iep2
             bint(i)=b2(i,j,k)*g2b(i)*dx2b(j)
           enddo
           call x1del_den(bint,dbint)
           do i=is,iep1
             scratch1(i,j,k)=-dbint(i)
           enddo
         enddo
       enddo
!
!      calculate del2(b1 h1 dx1)
!
       do k=ks,ke
         do i=is,iep1
           do j=jsm2,jep2
             bint(j)=b1(i,j,k)*dx1b(i)
           enddo
           call x2del_den(bint,dbint)
           do j=js,jep1
             scratch2(i,j,k)=-dbint(j)
           enddo
         enddo
       enddo
!
       do k=ks,ke
         do j=js,jep1
           do i=is,iep1
             emf3(i,j,k) = scratch1(i,j,k) - scratch2(i,j,k)
           enddo
         enddo
       enddo
!
!      calculate eta
!
       dtmdi2lc=0.D0
       do k=ks,ke
         do j=js,je
           do i=is,ie
             current=sqrt((0.25D0*(emf1(i,j,k)+emf1(i,j+1,k) &
               +emf1(i,j,k+1)+emf1(i,j+1,k+1)))**2 &
               +(0.25D0*(emf2(i,j,k)+emf2(i+1,j,k) &
               +emf2(i,j,k+1)+emf2(i+1,j,k+1)))**2 &
               +(0.25D0*(emf3(i,j,k)+emf3(i+1,j,k) &
               +emf3(i,j+1,k)+emf3(i+1,j+1,k)))**2)
             vv1 = max( abs(v1(i,j,k)) , abs(v1(i+1,j,k)) )
             vv2 = max( abs(v2(i,j,k)) , abs(v2(i,j+1,k)) )
             vv3 = max( abs(v3(i,j,k)) , abs(v3(i,j,k+1)) )
             dr1=dx1a(i)
             dr2=dx2a(j)*g2b(i)
             dr3=dx3a(k)*g31b(i)*g32b(j)
             drmin=min(dr1,dr2,dr3)
             eta(i,j,k)=min(0.5D0*current/sqrt(d(i,j,k)), &
               0.25D0*drmin*clight)
             etatot=eta(i,j,k)+0.5D0*(vv1*dr1+vv2*dr2+vv3*dr3)
             dtmdi2tmp=(4.D0*etatot/drmin**2)**2
             dtmdi2lc=max(dtmdi2lc,dtmdi2tmp)
           enddo
         enddo
       enddo
       call bvaleta(eta)
       call MPI_ALLREDUCE(dtmdi2lc,dtmdi2,1, &
         MPI_DOUBLE_PRECISION,MPI_MAX, &
         MPI_COMM_WORLD,ierr)
!
!-----------------------------------------------------------------------
!---- 1.  emf1 ---------------------------------------------------------
!-----------------------------------------------------------------------
!
!      calculate eta on emf1 edge
!
       do k=ks,kep1
         do j=js,jep1
           do i=is,ie
             scratch3(i,j,k)=0.25D0*(eta(i,j,k) &
               +eta(i,j-1,k)+eta(i,j-1,k-1) &
               +eta(i,j,k-1))
             vp=0.5D0*(v2(i,j,k)+v2(i,j,k-1))
             fq=tanh(vp*dt*g2bi(i)*dx2bi(j)/ww)
             scratch3(i,j,k)=scratch3(i,j,k) &
               +0.5D0*vp*fq*g2b(i)*dx2b(j)
             vp=0.5D0*(v3(i,j,k)+v3(i,j-1,k))
             fq=tanh(vp*dt*g31bi(i)*g32ai(j)*dx3bi(k)/ww)
             scratch3(i,j,k)=scratch3(i,j,k) &
               +0.5D0*vp*fq*g31b(i)*g32a(j)*dx3b(k)
           enddo
         enddo
       enddo
!
!      calculate del2(b3 h3 dx3)
!
       do k=ks,kep1
         do i=is,ie
           do j=jsm2,jep2
             bint(j)=b3(i,j,k)*g31b(i)*g32b(j)*dx3b(k)
           enddo
!          if js and jep1 is at x2a=0 and pi/2
!          if (isplt .eq. -1) then
!           if (myid2 .eq. 0) then
!            bint(jsm1)=-bint(js)
!            bint(jsm2)=-bint(jsp1)
!           endif
!           if (myid2 .eq. nproc2-1) then
!            bint(jep1)=-bint(je)
!            bint(jep2)=-bint(jem1)
!           endif
!          endif
           call x2del_den(bint,dbint)
           do j=js,jep1
             scratch1(i,j,k)= 0.
             emf1(i,j,k)=dbint(j)*g2bi(i)*dx2bi(j) &
               *g31bi(i)*g32ai(j)*dx3bi(k)
             curl1(i,j,k) = (bint(j) - bint(j-1))*g2bi(i)*dx2bi(j) &
                *g31bi(i)*g32ai(j)*dx3bi(k)
           enddo
         enddo
       enddo
!
!      calculate del3(b2 h2 dx2)
!
       do j=js,jep1
         do i=is,ie
           do k=ksm2,kep2
             bint(k)=b2(i,j,k)*g2b(i)*dx2b(j)
           enddo
           call x3del_den(bint,dbint)
           do k=ks,kep1
             scratch2(i,j,k)= 0.
             emf1(i,j,k) = emf1(i,j,k) &
              - dbint(k)*g2bi(i)*dx2bi(j) &
               *g31bi(i)*g32ai(j)*dx3bi(k)
             curl1(i,j,k) = curl1(i,j,k) &
                - (bint(k) - bint(k-1))*g2bi(i)*dx2bi(j) &
                *g31bi(i)*g32ai(j)*dx3bi(k)
             scratch1(i,j,k)=scratch1(i,j,k) &
               + (-scratch3(i,j,k))*emf1(i,j,k)
           enddo
         enddo
       enddo
!
       do k=ks,kep1
         do j=js,jep1
           do i=is,ie
             srd1(i,j,k) = scratch1(i,j,k) - scratch2(i,j,k)
           enddo
         enddo
       enddo
!
!-----------------------------------------------------------------------
!---- 2.  emf2 ---------------------------------------------------------
!-----------------------------------------------------------------------
!
!      calculate eta and on emf2 edge
!
       do k=ks,kep1
         do j=js,je
           do i=is,iep1
             scratch3(i,j,k)=0.25D0*(eta(i,j,k) &
               +eta(i-1,j,k)+eta(i-1,j,k-1) &
               +eta(i,j,k-1))
             vp=0.5D0*(v3(i,j,k)+v3(i-1,j,k))
             fq=tanh(vp*dt*g31ai(i)*g32bi(j)*dx3bi(k)/ww)
             scratch3(i,j,k)=scratch3(i,j,k) &
               +0.5D0*vp*fq*g31a(i)*g32b(j)*dx3b(k)
             vp=0.5D0*(v1(i,j,k)+v1(i,j,k-1))
             fq=tanh(vp*dt*dx1bi(i)/ww)
             scratch3(i,j,k)=scratch3(i,j,k) &
               +0.5D0*vp*fq*dx1b(i)
           enddo
         enddo
       enddo
!
!      calculate del3(b1 h1 dx1)
!
       do j=js,je
         do i=is,iep1
           do k=ksm2,kep2
             bint(k)=b1(i,j,k)*dx1b(i)
           enddo
           call x3del_den(bint,dbint)
           do k=ks,kep1
             scratch1(i,j,k)= 0.
             emf2(i,j,k)=dbint(k) &
               *g31ai(i)*g32bi(j)*dx3bi(k)*dx1bi(i)
             curl2(i,j,k) = (bint(k) - bint(k-1)) &
               *g31ai(i)*g32bi(j)*dx3bi(k)*dx1bi(i)
           enddo
         enddo
       enddo
!
!      calculate del1(b3 h3 dx3)
!
       do k=ks,kep1
         do j=js,je
           do i=ism2,iep2
             bint(i)=b3(i,j,k)*g31b(i)*g32b(j)*dx3b(k)
           enddo
           call x1del_den(bint,dbint)
           do i=is,iep1
             scratch2(i,j,k)= 0.
             curl2(i,j,k) = curl2(i,j,k) &
               - (bint(i) - bint(i-1))*g31ai(i)*g32bi(j)*dx3bi(k) &
               *dx1bi(i)
             emf2(i,j,k)=emf2(i,j,k) &
               -dbint(i)*g31ai(i)*g32bi(j)*dx3bi(k)*dx1bi(i)
             scratch1(i,j,k)=scratch1(i,j,k) &
               + (-scratch3(i,j,k))*emf2(i,j,k)
           enddo
         enddo
       enddo
!
       do k=ks,kep1
         do j=js,je
           do i=is,iep1
             srd2(i,j,k) = scratch1(i,j,k) - scratch2(i,j,k)
           enddo
         enddo
       enddo
!
!-----------------------------------------------------------------------
!---- 3.  emf3 ---------------------------------------------------------
!-----------------------------------------------------------------------
!
!      calculate eta on emf3 edge
       do k=ks,ke
         do j=js,jep1
           do i=is,iep1
             scratch3(i,j,k)=0.25D0*(eta(i,j,k) &
               +eta(i-1,j,k)+eta(i-1,j-1,k) &
               +eta(i,j-1,k))
             vp=0.5D0*(v1(i,j,k)+v1(i,j-1,k))
             fq=tanh(vp*dt*dx1bi(i)/ww)
             scratch3(i,j,k)=scratch3(i,j,k) &
               +0.5D0*vp*fq*dx1b(i)
             vp=0.5D0*(v2(i,j,k)+v2(i-1,j,k))
             fq=tanh(vp*dt*g2ai(i)*dx2bi(j)/ww)
             scratch3(i,j,k)=scratch3(i,j,k) &
               +0.5D0*vp*fq*g2a(i)*dx2b(j)
           enddo
         enddo
       enddo
!
!      calculate del1(b2 h2 dx2)
!
       do k=ks,ke
         do j=js,jep1
           do i=ism2,iep2
             bint(i)=b2(i,j,k)*g2b(i)*dx2b(j)
           enddo
           call x1del_den(bint,dbint)
           do i=is,iep1
             scratch1(i,j,k)= 0.
             emf3(i,j,k)=dbint(i) &
               *dx1bi(i)*g2ai(i)*dx2bi(j)
             curl3(i,j,k) = (bint(i) - bint(i-1)) &
               *dx1bi(i)*g2ai(i)*dx2bi(j)
           enddo
         enddo
       enddo
!
!      calculate del2(b1 h1 dx1)
!
       do k=ks,ke
         do i=is,iep1
           do j=jsm2,jep2
             bint(j)=b1(i,j,k)*dx1b(i)
           enddo
           call x2del_den(bint,dbint)
           do j=js,jep1
             scratch2(i,j,k)= 0.
             curl3(i,j,k) = curl3(i,j,k) &
               - (bint(j)-bint(j-1))*dx1bi(i)*g2ai(i)*dx2bi(j)
             emf3(i,j,k)=emf3(i,j,k) &
               -dbint(j)*dx1bi(i)*g2ai(i)*dx2bi(j)
             scratch1(i,j,k)=scratch1(i,j,k) &
               + (-scratch3(i,j,k))*emf3(i,j,k)
           enddo
         enddo
       enddo
!
       do k=ks,ke
         do j=js,jep1
           do i=is,iep1
             srd3(i,j,k) = scratch1(i,j,k) - scratch2(i,j,k)
           enddo
         enddo
       enddo
!
!-----------------------------------------------------------------------
!
!      Set boundary values for "term1" and "term2".
!
       call bvaldiffemf1 (srd1)
!
!      Compute "emf1" for all 1-edges, including the ghost zones.
!
       do k=ksm2,kep3
         do j=jsm2,jep3
           do i=ism2,iep2
             emf1(i,j,k) =  srd1(i,j,k)
             emf1s(i,j,k) = emf1s(i,j,k) + emf1(i,j,k)
           enddo
         enddo
       enddo
!
!-----------------------------------------------------------------------
!
!      Set boundary values for "term1" and "term2".
!
       call bvaldiffemf2 (srd2)
!
!      Compute "emf2" for all 2-edges, including the ghost zones.
!
       do k=ksm2,kep3
         do j=jsm2,jep2
           do i=ism2,iep3
             emf2(i,j,k) = srd2(i,j,k)
             emf2s(i,j,k) = emf2s(i,j,k) + emf2(i,j,k)
           enddo
         enddo
       enddo
!
!-----------------------------------------------------------------------
!
!      Set boundary values for "term1" and "term2".
!
       call bvaldiffemf3 (srd3)
!
!      Compute "emf3" for all 3-edges, including the ghost zones.
!
       do k=ksm2,kep2
         do j=jsm2,jep3
           do i=ism2,iep3
             emf3(i,j,k) = srd3(i,j,k)
             emf3s(i,j,k) = emf3s(i,j,k) + emf3(i,j,k)
           enddo
         enddo
       enddo
!
#ifndef ISOMHD
  do k=ks,ke
    do j=js,je
      do i=is,ie
        q1 = -(emf1(i,j,k)*curl1(i,j,k) &
          + emf1(i,j+1,k)*curl1(i,j+1,k) &
          + emf1(i,j,k+1)*curl1(i,j,k+1) &
          + emf1(i,j+1,k+1)*curl1(i,j+1,k+1)) * 0.25D0 &
          - (emf2(i,j,k)*curl2(i,j,k) &
          + emf2(i+1,j,k)*curl2(i+1,j,k) &
          + emf2(i,j,k+1)*curl2(i,j,k+1) &
          + emf2(i+1,j,k+1)*curl2(i+1,j,k+1)) * 0.25D0 &
          - (emf3(i,j,k)*curl3(i,j,k) &
          + emf3(i+1,j,k)*curl3(i+1,j,k) &
          + emf3(i,j+1,k)*curl3(i,j+1,k) &
          + emf3(i+1,j+1,k)*curl3(i+1,j+1,k)) * 0.25D0
!       q1 = - ((curl1(i,j,k) + curl1(i,j+1,k) + curl1(i,j,k+1) + curl1(i,j+1,k+1)) &
!         *(emf1(i,j,k) + emf1(i,j+1,k) + emf1(i,j,k+1) + emf1(i,j+1,k+1)) + &
!          (curl2(i,j,k) + curl2(i+1,j,k) + curl2(i,j,k+1) + curl2(i+1,j,k+1)) &
!         *(emf2(i,j,k) + emf2(i+1,j,k) + emf2(i,j,k+1) + emf2(i+1,j,k+1)) + &
!          (curl3(i,j,k) + curl3(i+1,j,k) + curl3(i,j+1,k) + curl3(i+1,j+1,k)) &
!         *(emf3(i,j,k) + emf3(i+1,j,k) + emf3(i,j+1,k) + emf3(i+1,j+1,k))) &
!         /16.D0
        delei(i,j,k) = delei(i,j,k) + max(0.D0,q1)
      enddo
    enddo
  enddo
#endif

end subroutine diffemf

end module ModRHS
