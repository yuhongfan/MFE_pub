module ModConductivity
#ifdef THCONDUCT
  use ModPar
  implicit none

  real, dimension(in,jn,kn) :: diffcond, tau_hc
#endif
end module ModConductivity
