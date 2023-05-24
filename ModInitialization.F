module ModInitialization

  implicit none

  private

  public :: initialize

contains

subroutine initialize
  use ModPhysics
  use ModGrid
  use ModField
  use ModUserSetup
  use ModBval

  implicit none
  integer :: i,j,k

!--------------------------------------------------------
       call input_param

       call grid

       call field_init
       call bvald
       call bvalei
       call bvalv
#ifndef ISOMHD
       do k=ksm2,kep2
         do j=jsm2,jep2
           do i=ism2,iep2
             tfield(i,j,k)=eint(i,j,k)*(gamma-1.D0) &
                /d(i,j,k)*unit_v**2/(rstar*unit_temp)
           enddo
         enddo
       enddo
#endif

end subroutine initialize

subroutine input_param
  use ModPhysics
  use ModControl

  implicit none

  integer :: itmax_in,itintv_in,itout_in,ifile_in,irst_in
  real :: tout_in,tend_in,wtlimit_in
  character(len=80) infile_in,rstfile_in, &
    b1file_in,b2file_in,b3file_in, &
    v1file_in,v2file_in,v3file_in, dfile_in, efile_in, &
    engfile_in,gridfile_in,qfile_in
  namelist / control_namelist / infile_in,rstfile_in, &
    b1file_in,b2file_in,b3file_in, &
    v1file_in,v2file_in,v3file_in,dfile_in,efile_in, &
    engfile_in,gridfile_in, &
    qfile_in, &
    itmax_in, itintv_in, itout_in, ifile_in, irst_in, &
    tout_in, tend_in, wtlimit_in

  real :: gamma_in, muconst_in
  real ::  unit_rho_in, unit_len_in, unit_b_in, unit_temp_in
  namelist / physicalparams_namelist / gamma_in, muconst_in, &
    unit_rho_in, unit_len_in, unit_b_in, unit_temp_in

!--------------------------------------------------------

       open(unit=20, file="inparam", status="old", &
         position="rewind")
       read(unit=20, nml=control_namelist)
       read(unit=20, nml=physicalparams_namelist)
       close(unit=20)

       infile=infile_in
       rstfile=rstfile_in
       b1file=b1file_in
       b2file=b2file_in
       b3file=b3file_in
       v1file=v1file_in
       v2file=v2file_in
       v3file=v3file_in
       dfile=dfile_in
       efile=efile_in
       engfile=engfile_in
       gridfile=gridfile_in
       qfile=qfile_in
       itmax=itmax_in
       itintv=itintv_in
       itout=itout_in
       ifile=ifile_in
       irst=irst_in
       tout=tout_in
       tend=tend_in
       wtlimit=wtlimit_in

       gamma=gamma_in
       muconst=muconst_in
       unit_rho=unit_rho_in
       unit_len=unit_len_in
       unit_b=unit_b_in
       unit_temp=unit_temp_in
       unit_v=unit_b/sqrt(4.D0*pi*unit_rho)
       unit_time=unit_len/unit_v

end subroutine input_param

end module ModInitialization
