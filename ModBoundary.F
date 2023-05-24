module ModBoundary
  use ModPar,   ONLY: in, jn, kn
  implicit none

  integer :: niib, noib, nijb, nojb, nikb, nokb
  integer :: isplt = 1

  real :: diib1(jn,kn),diib2(jn,kn),doib1(jn,kn),doib2(jn,kn)
  real :: dijb1(kn,in),dijb2(kn,in),dojb1(kn,in),dojb2(kn,in)
  real :: dikb1(in,jn),dikb2(in,jn),dokb1(in,jn),dokb2(in,jn)

#ifndef ISOMHD
  real :: eint_iib1(jn,kn),eint_iib2(jn,kn), &
          eint_oib1(jn,kn),eint_oib2(jn,kn)
  real :: eint_ijb1(kn,in),eint_ijb2(kn,in), &
          eint_ojb1(kn,in),eint_ojb2(kn,in)
  real :: eint_ikb1(in,jn),eint_ikb2(in,jn), &
          eint_okb1(in,jn),eint_okb2(in,jn)
#endif

  real :: v1iib1(jn,kn),v1iib2(jn,kn),v1iib3(jn,kn), &
          v1oib1(jn,kn),v1oib2(jn,kn),v1oib3(jn,kn)
  real :: v1ijb1(kn,in),v1ijb2(kn,in),v1ijb0(kn,in), &
          v1ojb1(kn,in),v1ojb2(kn,in)
  real :: v1ikb1(in,jn),v1ikb2(in,jn),v1ikb0(in,jn), &
          v1okb1(in,jn),v1okb2(in,jn)

  real :: v2iib1(jn,kn),v2iib2(jn,kn),v2iib0(jn,kn), &
          v2oib1(jn,kn),v2oib2(jn,kn)
  real :: v2ijb1(kn,in),v2ijb2(kn,in),v2ijb3(kn,in), &
          v2ojb1(kn,in),v2ojb2(kn,in),v2ojb3(kn,in)
  real :: v2ikb1(in,jn),v2ikb2(in,jn),v2ikb0(in,jn), &
          v2okb1(in,jn),v2okb2(in,jn)

  real :: v3iib1(jn,kn),v3iib2(jn,kn),v3iib0(jn,kn), &
          v3oib1(jn,kn),v3oib2(jn,kn)
  real :: v3ijb1(kn,in),v3ijb2(kn,in),v3ijb0(kn,in), &
          v3ojb1(kn,in),v3ojb2(kn,in)
  real :: v3ikb1(in,jn),v3ikb2(in,jn),v3ikb3(in,jn), &
          v3okb1(in,jn),v3okb2(in,jn),v3okb3(in,jn)

  real :: emf1iib1(jn,kn),emf1iib2(jn,kn), &
          emf1oib1(jn,kn),emf1oib2(jn,kn)
  real :: emf1ijb1(kn,in),emf1ijb2(kn,in),emf1ijb3(kn,in), &
          emf1ojb1(kn,in),emf1ojb2(kn,in),emf1ojb3(kn,in)
  real :: emf1ikb1(in,jn),emf1ikb2(in,jn),emf1ikb3(in,jn), &
          emf1okb1(in,jn),emf1okb2(in,jn),emf1okb3(in,jn)

  real :: emf2iib1(jn,kn),emf2iib2(jn,kn),emf2iib3(jn,kn), &
          emf2oib1(jn,kn),emf2oib2(jn,kn),emf2oib3(jn,kn)
  real :: emf2ijb1(kn,in),emf2ijb2(kn,in), &
          emf2ojb1(kn,in),emf2ojb2(kn,in)
  real :: emf2ikb1(in,jn),emf2ikb2(in,jn),emf2ikb3(in,jn), &
          emf2okb1(in,jn),emf2okb2(in,jn),emf2okb3(in,jn)

  real :: emf3iib1(jn,kn),emf3iib2(jn,kn),emf3iib3(jn,kn), &
          emf3oib1(jn,kn),emf3oib2(jn,kn),emf3oib3(jn,kn)
  real :: emf3ijb1(kn,in),emf3ijb2(kn,in),emf3ijb3(kn,in), &
          emf3ojb1(kn,in),emf3ojb2(kn,in),emf3ojb3(kn,in)
  real :: emf3ikb1(in,jn),emf3ikb2(in,jn), &
          emf3okb1(in,jn),emf3okb2(in,jn)

  real :: e2ibs(1:jn)=0., e2ibe(1:jn)=0., e3ibs(1:kn)=0., e3ibe(1:kn)=0.
  real :: e3jbs(1:kn)=0., e3jbe(1:kn)=0., e1jbs(1:in)=0., e1jbe(1:in)=0.
  real :: e1kbs(1:in)=0., e1kbe(1:in)=0., e2kbs(1:jn)=0., e2kbe(1:jn)=0.

end module ModBoundary
