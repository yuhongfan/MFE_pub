module ModSundry
! initialized in field_init

  integer :: itnow, ntcond
  real :: dt, time, ctime
  real :: courno = 0.5D0

  real :: dtmdi2 = 0.D0
  integer :: nlf_d = 4, nlf_v = 4, nlf_e = 4
  real :: alpha_d = 1.D0, alpha_v = 1.D0, alpha_e = 1.D0

  real :: ovc2 = 1.D-99, clight=1.D99

  real :: deltat = 1.D99

#ifndef ISOMHD
  real :: dfloor = 1.D-99, efloor = 1.D-99
#else
  real :: dfloor = 1.D-99
#endif

end module ModSundry
