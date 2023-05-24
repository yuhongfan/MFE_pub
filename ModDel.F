module ModDel
  use ModPar,   ONLY: in, jn, kn
  implicit none

  real, dimension(1:in,1:jn,1:kn) :: dels1, dels2, dels3, deld
  real, dimension(1:in,1:jn,1:kn) :: emf1s, emf2s, emf3s
  real, dimension(1:in,1:jn,1:kn) :: dels1_pr, dels2_pr, dels3_pr, deld_pr
  real, dimension(1:in,1:jn,1:kn) :: emf1s_pr, emf2s_pr, emf3s_pr
#ifndef ISOMHD
  real, dimension(1:in,1:jn,1:kn) :: delei
  real, dimension(1:in,1:jn,1:kn) :: delei_pr
#endif
#ifdef THCONDUCT
  real, dimension(1:in,1:jn,1:kn) :: delq
  real, dimension(1:in,1:jn,1:kn) :: delq_pr
#endif

end module ModDel
