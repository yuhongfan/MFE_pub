module ModGetQtys
  private
  
  public :: get_global

contains
  
#ifndef ISOMHD
  subroutine get_global(em, ek, eth)
    use ModPar
    use ModGrid
    use ModField,    ONLY: v1, v2, v3, b1, b2, b3, d, eint
    implicit none
    
    include "mpif.h"
    
    real, intent(out) :: em, ek, eth
    
    integer :: i, j, k, ierr
    real :: dVol, energylc(1:3), energy(1:3)
    !----------------------------------------------------------------------------------
    
    energy   = 0.d0
    energylc = 0.d0
    
    do k=ks,ke; do j=js,je; do i=is,ie
       dVol = dvl1a(i)*dvl2a(j)*dvl3a(k)
       energylc(1) = energylc(1) + 0.125d0*dVol*((b1(i,j,k) + b1(i+1,j,k))**2 + &
            (b2(i,j,k) + b2(i,j+1,k))**2 + (b3(i,j,k) + b3(i,j,k+1))**2)
       energylc(2) = energylc(2) + 0.125d0*dVol*d(i,j,k) &
         *((v1(i,j,k)+v1(i+1,j,k))**2 + (v2(i,j,k)+v2(i,j+1,k))**2 &
         + (v3(i,j,k)+v3(i,j,k+1))**2)
       energylc(3) = energylc(3) + dVol*eint(i,j,k)
    enddo; enddo; enddo
    call MPI_ALLREDUCE(energylc, energy, 3, MPI_DOUBLE_PRECISION, MPI_SUM, &
         MPI_COMM_WORLD,ierr)
    
    em  = energy(1)
    ek  = energy(2)
    eth  = energy(3)
    
  end subroutine get_global
#else
  subroutine get_global(em, ek)
    use ModPar
    use ModGrid
    use ModField,    ONLY: v1, v2, v3, b1, b2, b3, d
    implicit none
    
    include "mpif.h"
    
    real, intent(out) :: em, ek
    
    integer :: i, j, k, ierr
    real :: dVol, energylc(1:2), energy(1:2)
    !----------------------------------------------------------------------------------
    
    energy   = 0.d0
    energylc = 0.d0
    
    do k=ks,ke; do j=js,je; do i=is,ie
       dVol = dvl1a(i)*dvl2a(j)*dvl3a(k)
       energylc(1) = energylc(1) + 0.125d0*dVol*((b1(i,j,k) + b1(i+1,j,k))**2 + &
            (b2(i,j,k) + b2(i,j+1,k))**2 + (b3(i,j,k) + b3(i,j,k+1))**2)
       energylc(2) = energylc(2) + 0.125d0*dVol*d(i,j,k) &
         *((v1(i,j,k)+v1(i+1,j,k))**2 + (v2(i,j,k)+v2(i,j+1,k))**2 &
         + (v3(i,j,k)+v3(i,j,k+1))**2)
    enddo; enddo; enddo
    call MPI_ALLREDUCE(energylc, energy, 2, MPI_DOUBLE_PRECISION, MPI_SUM, &
         MPI_COMM_WORLD,ierr)
    
    em  = energy(1)
    ek  = energy(2)
    
  end subroutine get_global
#endif

end module ModGetQtys
