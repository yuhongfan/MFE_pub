module ModPar
!      Primary parameters:

  integer, parameter :: inmax = 425, jnmax =125, knmax = 13, &
                        ijknmax = max(inmax,jnmax,knmax)
  integer, parameter :: nproc1 = 30, nproc2 = 6, nproc3=1, &
                        nproc=nproc1*nproc2*nproc3
  integer, parameter :: in = (inmax-5)/nproc1+5, &
                        jn = (jnmax-5)/nproc2+5, &
                        kn = (knmax-5)/nproc3+5, &
                        ijkn = max(in,jn,kn)
  integer :: myid, myid1, myid2, myid3

  real(8), parameter :: tiny = 1.0D-99, huge = 1.0D+99

end module ModPar
