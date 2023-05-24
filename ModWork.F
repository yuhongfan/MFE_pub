module ModWork
  use ModPar,   ONLY: in, jn, kn
  implicit none

  real, dimension(in,jn,kn) :: srd1, srd2, srd3, &
                               scratch1, scratch2, scratch3

end module ModWork
