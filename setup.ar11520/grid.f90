Module grid
      implicit none
      save
      integer, public ::       is  , js  , ks  , ie  , je  , ke
      integer, public ::       ism1, jsm1, ksm1, ism2, jsm2, ksm2
      integer, public ::       isp1, jsp1, ksp1, isp2, jsp2, ksp2
      integer, public ::       iem1, jem1, kem1, iem2, jem2, kem2
      integer, public ::       iep1, jep1, kep1, iep2, jep2, kep2
      integer, public ::       iep3, jep3, kep3
      integer, public ::       ione, jone, kone
!
      real(8), allocatable, public :: dx1a(:), dx2a(:), dx3a(:)
      real(8), allocatable, public :: dx1b(:), dx2b(:), dx3b(:)
      real(8), allocatable, public :: dx1ai(:), dx2ai(:), dx3ai(:)
      real(8), allocatable, public :: dx1bi(:), dx2bi(:), dx3bi(:)
      real(8), allocatable, public :: x1a(:), x2a(:), x3a(:)
      real(8), allocatable, public :: x1b(:), x2b(:), x3b(:)
!
      real(8), allocatable, public :: g2a(:), dg2ad1(:), g2ai(:)
      real(8), allocatable, public :: g2b(:), dg2bd1(:), g2bi(:)
      real(8), allocatable, public :: g31a(:), dg31ad1(:), g31ai(:)
      real(8), allocatable, public :: g31b(:), dg31bd1(:), g31bi(:)
      real(8), allocatable, public :: g32a(:), dg32ad2(:), g32ai(:)
      real(8), allocatable, public :: g32b(:), dg32bd2(:), g32bi(:)
!
      real(8), allocatable, public :: dvl1a(:), dvl1ai(:), dvl1b(:), dvl1bi(:)
      real(8), allocatable, public :: dvl2a(:), dvl2ai(:), dvl2b(:), dvl2bi(:)
      real(8), allocatable, public :: dvl3a(:), dvl3ai(:), dvl3b(:), dvl3bi(:)
End Module grid
