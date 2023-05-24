module ModField
  use ModPar,   ONLY: in, jn, kn
  implicit none

  real, dimension(1:in,1:jn,1:kn) :: v1, v2, v3, b1, b2, b3, d
#ifndef ISOMHD
  real, dimension(1:in,1:jn,1:kn) :: eint, tfield
#endif
#ifdef THCONDUCT
  real, dimension(1:in,1:jn,1:kn) :: q
#endif

end module ModField
