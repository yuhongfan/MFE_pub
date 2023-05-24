module ModPhysics
!
!      Physical constants
!
  real, parameter :: g_const = 6.672D-8, &
                     rgas = 8.3D7, &
                     kboltz = 1.381D-16, &
                     mproton = 1.673D-24, &
                     pi   = 3.1415926535898D0, &
                     mpovme = 1836.15267245D0
!
!      Physical parameters and units specified in input_param
!
  real :: gamma = 1.6666666666667D0, muconst = 1.D0
! fundamental units
  real :: unit_rho, unit_len, unit_b, unit_temp
! derived units: NOTE it is assumed that unit_v=unit_b/sqrt(4*pi*unit_rho)
  real :: unit_v, unit_time
!
!      Physical parameters specified in field_init
!
  real :: cs2, gacc, rstar
#ifdef THCONDUCT
  real :: kappa0_dim = 0.D0, kappa0 = 0.D0, &
          diffcoef0=0.D0, aval0 = 1.D0
#endif

end module ModPhysics
