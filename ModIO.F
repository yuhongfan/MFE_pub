module ModIO

  implicit none

  private

  public :: readdata_mpi
  public :: readrst_mpi
  public :: writedata_mpi
  public :: writerst_mpi
  public :: writegrid
  public :: writedata_cpu
  public :: writegrid_cpu
  public :: writeqty_mpi
  public :: readqty_mpi
  public :: writeemf1_mpi
  public :: writeemf2_mpi
  public :: writeemf3_mpi
  public :: read12qty_mpi
  public :: read12v1_mpi
  public :: read12v2_mpi
  public :: read12emf1_mpi
  public :: read12emf2_mpi
  public :: read23qty_mpi
  public :: read23emf2_mpi
  public :: read23emf3_mpi
  public :: read23emf1_mpi

contains
!-----------------------------------------------------------------------
subroutine readdata_mpi
  use ModPar
  use ModControl
  use ModSundry
  use ModField
  implicit none

  include "mpif.h"

  integer i,j,k,counter
  integer sizes(3), subsizes(3), starts(3)
  integer filetype,fileinfo,fhandl,mysize,ierr
  integer(kind=MPI_OFFSET_KIND) disp
  integer mpi_status(MPI_STATUS_SIZE)
  character(len=4) idout
  integer inmax_in,jnmax_in,knmax_in
  real time_in
  real, allocatable :: buf(:,:,:)
!--------------------------------------------------------------
       allocate(buf(in-1,jn-1,kn-1))
!
       write(idout,'(i4.4)') ifile
!
       sizes(1)=inmax-1
       sizes(2)=jnmax-1
       sizes(3)=knmax-1
!
       subsizes(1)=in-1
       starts(1)=myid1*(in-5)
!
       subsizes(2)=jn-1
       starts(2)=myid2*(jn-5)
!
       subsizes(3)=kn-1
       starts(3)=myid3*(kn-5)
!
       call MPI_Type_create_subarray(3, sizes, subsizes, starts, &
         MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, filetype, ierr)
       call MPI_Type_commit(filetype, ierr)
       call MPI_Info_create(fileinfo, ierr)
!
       fhandl=13
       call MPI_File_open(MPI_COMM_WORLD, &
         trim(v1file)//idout//'.dat', &
         MPI_MODE_RDONLY, &
         fileinfo, fhandl, ierr)
!
       if(myid .eq. 0) then
         call MPI_File_read(fhandl, time_in, 1, &
           MPI_DOUBLE_PRECISION, mpi_status, ierr)
         call MPI_File_read(fhandl, inmax_in, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_read(fhandl, jnmax_in, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_read(fhandl, knmax_in, 1, MPI_INTEGER, &
           mpi_status, ierr)
       endif
       call MPI_BCAST(time_in,1,MPI_DOUBLE_PRECISION,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(inmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(jnmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(knmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       if(inmax_in .ne. inmax .or.  &
         jnmax_in .ne. jnmax .or.  &
         knmax_in .ne. knmax) then
          write(6,*) 'inmax,jnmax,knmax in v1file', &
            inmax_in,jnmax_in,knmax_in, &
            ' do not match inmax, jnmax, knmax in par.h'
          call MPI_ABORT(MPI_COMM_WORLD, 1,ierr)
       endif
       time=time_in
!
       disp=sizeof(time)+3*sizeof(inmax)
       call MPI_FILE_SET_VIEW(fhandl, disp, &
         MPI_DOUBLE_PRECISION, filetype, 'native', &
         fileinfo, ierr)
!
       mysize=subsizes(1)*subsizes(2)*subsizes(3)
       call MPI_File_read_all(fhandl, buf, &
         mysize, MPI_DOUBLE_PRECISION, mpi_status, ierr)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call MPI_File_close(fhandl, ierr)
       do k=1,kn-1
         do j=1,jn-1
           do i=1,in-1
             v1(i,j,k)=buf(i,j,k)
           enddo
         enddo
       enddo
!
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       fhandl=13
       call MPI_File_open(MPI_COMM_WORLD, &
         trim(v2file)//idout//'.dat', &
         MPI_MODE_RDONLY, &
         fileinfo, fhandl, ierr)
!
       if(myid .eq. 0) then
         call MPI_File_read(fhandl, time_in, 1, &
           MPI_DOUBLE_PRECISION, mpi_status, ierr)
         call MPI_File_read(fhandl, inmax_in, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_read(fhandl, jnmax_in, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_read(fhandl, knmax_in, 1, MPI_INTEGER, &
           mpi_status, ierr)
       endif
       call MPI_BCAST(time_in,1,MPI_DOUBLE_PRECISION,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(inmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(jnmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(knmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       if(inmax_in .ne. inmax .or.  &
         jnmax_in .ne. jnmax .or.  &
         knmax_in .ne. knmax) then
          write(6,*) 'inmax,jnmax,knmax in v2file', &
            inmax_in,jnmax_in,knmax_in, &
            ' do not match inmax, jnmax, knmax in par.h'
          call MPI_ABORT(MPI_COMM_WORLD, 1,ierr)
       endif
!
       disp=sizeof(time)+3*sizeof(inmax)
       call MPI_FILE_SET_VIEW(fhandl, disp, &
         MPI_DOUBLE_PRECISION, filetype, 'native', &
         fileinfo, ierr)
!
       mysize=subsizes(1)*subsizes(2)*subsizes(3)
       call MPI_File_read_all(fhandl, buf, &
         mysize, MPI_DOUBLE_PRECISION, mpi_status, ierr)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call MPI_File_close(fhandl, ierr)
       do k=1,kn-1
         do j=1,jn-1
           do i=1,in-1
             v2(i,j,k)=buf(i,j,k)
           enddo
         enddo
       enddo
!
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       fhandl=13
       call MPI_File_open(MPI_COMM_WORLD, &
         trim(v3file)//idout//'.dat', &
         MPI_MODE_RDONLY, &
         fileinfo, fhandl, ierr)
!
       if(myid .eq. 0) then
         call MPI_File_read(fhandl, time_in, 1, &
           MPI_DOUBLE_PRECISION, mpi_status, ierr)
         call MPI_File_read(fhandl, inmax_in, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_read(fhandl, jnmax_in, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_read(fhandl, knmax_in, 1, MPI_INTEGER, &
           mpi_status, ierr)
       endif
       call MPI_BCAST(time_in,1,MPI_DOUBLE_PRECISION,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(inmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(jnmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(knmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       if(inmax_in .ne. inmax .or.  &
         jnmax_in .ne. jnmax .or.  &
         knmax_in .ne. knmax) then
          write(6,*) 'inmax,jnmax,knmax in v3file', &
            inmax_in,jnmax_in,knmax_in, &
            ' do not match inmax, jnmax, knmax in par.h'
          call MPI_ABORT(MPI_COMM_WORLD, 1,ierr)
       endif
!
       disp=sizeof(time)+3*sizeof(inmax)
       call MPI_FILE_SET_VIEW(fhandl, disp, &
         MPI_DOUBLE_PRECISION, filetype, 'native', &
         fileinfo, ierr)
!
       mysize=subsizes(1)*subsizes(2)*subsizes(3)
       call MPI_File_read_all(fhandl, buf, &
         mysize, MPI_DOUBLE_PRECISION, mpi_status, ierr)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call MPI_File_close(fhandl, ierr)
       do k=1,kn-1
         do j=1,jn-1
           do i=1,in-1
             v3(i,j,k)=buf(i,j,k)
           enddo
         enddo
       enddo
!
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       fhandl=13
       call MPI_File_open(MPI_COMM_WORLD, &
         trim(b1file)//idout//'.dat', &
         MPI_MODE_RDONLY, &
         fileinfo, fhandl, ierr)
!
       if(myid .eq. 0) then
         call MPI_File_read(fhandl, time_in, 1, &
           MPI_DOUBLE_PRECISION, mpi_status, ierr)
         call MPI_File_read(fhandl, inmax_in, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_read(fhandl, jnmax_in, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_read(fhandl, knmax_in, 1, MPI_INTEGER, &
           mpi_status, ierr)
       endif
       call MPI_BCAST(time_in,1,MPI_DOUBLE_PRECISION,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(inmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(jnmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(knmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       if(inmax_in .ne. inmax .or.  &
         jnmax_in .ne. jnmax .or.  &
         knmax_in .ne. knmax) then
          write(6,*) 'inmax,jnmax,knmax in b1file', &
            inmax_in,jnmax_in,knmax_in, &
            ' do not match inmax, jnmax, knmax in par.h'
          call MPI_ABORT(MPI_COMM_WORLD, 1,ierr)
       endif
!
       disp=sizeof(time)+3*sizeof(inmax)
       call MPI_FILE_SET_VIEW(fhandl, disp, &
         MPI_DOUBLE_PRECISION, filetype, 'native', &
         fileinfo, ierr)
!
       mysize=subsizes(1)*subsizes(2)*subsizes(3)
       call MPI_File_read_all(fhandl, buf, &
         mysize, MPI_DOUBLE_PRECISION, mpi_status, ierr)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call MPI_File_close(fhandl, ierr)
       do k=1,kn-1
         do j=1,jn-1
           do i=1,in-1
             b1(i,j,k)=buf(i,j,k)
           enddo
         enddo
       enddo
!
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       fhandl=13
       call MPI_File_open(MPI_COMM_WORLD, &
         trim(b2file)//idout//'.dat', &
         MPI_MODE_RDONLY, &
         fileinfo, fhandl, ierr)
!
       if(myid .eq. 0) then
         call MPI_File_read(fhandl, time_in, 1, &
           MPI_DOUBLE_PRECISION, mpi_status, ierr)
         call MPI_File_read(fhandl, inmax_in, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_read(fhandl, jnmax_in, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_read(fhandl, knmax_in, 1, MPI_INTEGER, &
           mpi_status, ierr)
       endif
       call MPI_BCAST(time_in,1,MPI_DOUBLE_PRECISION,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(inmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(jnmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(knmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       if(inmax_in .ne. inmax .or.  &
         jnmax_in .ne. jnmax .or.  &
         knmax_in .ne. knmax) then
          write(6,*) 'inmax,jnmax,knmax in b2file', &
            inmax_in,jnmax_in,knmax_in, &
            ' do not match inmax, jnmax, knmax in par.h'
          call MPI_ABORT(MPI_COMM_WORLD, 1,ierr)
       endif
!
       disp=sizeof(time)+3*sizeof(inmax)
       call MPI_FILE_SET_VIEW(fhandl, disp, &
         MPI_DOUBLE_PRECISION, filetype, 'native', &
         fileinfo, ierr)
!
       mysize=subsizes(1)*subsizes(2)*subsizes(3)
       call MPI_File_read_all(fhandl, buf, &
         mysize, MPI_DOUBLE_PRECISION, mpi_status, ierr)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call MPI_File_close(fhandl, ierr)
       do k=1,kn-1
         do j=1,jn-1
           do i=1,in-1
             b2(i,j,k)=buf(i,j,k)
           enddo
         enddo
       enddo
!
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       fhandl=13
       call MPI_File_open(MPI_COMM_WORLD, &
         trim(b3file)//idout//'.dat', &
         MPI_MODE_RDONLY, &
         fileinfo, fhandl, ierr)
!
       if(myid .eq. 0) then
         call MPI_File_read(fhandl, time_in, 1, &
           MPI_DOUBLE_PRECISION, mpi_status, ierr)
         call MPI_File_read(fhandl, inmax_in, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_read(fhandl, jnmax_in, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_read(fhandl, knmax_in, 1, MPI_INTEGER, &
           mpi_status, ierr)
       endif
       call MPI_BCAST(time_in,1,MPI_DOUBLE_PRECISION,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(inmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(jnmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(knmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       if(inmax_in .ne. inmax .or.  &
         jnmax_in .ne. jnmax .or.  &
         knmax_in .ne. knmax) then
          write(6,*) 'inmax,jnmax,knmax in b3file', &
            inmax_in,jnmax_in,knmax_in, &
            ' do not match inmax, jnmax, knmax in par.h'
          call MPI_ABORT(MPI_COMM_WORLD, 1,ierr)
       endif
!
       disp=sizeof(time)+3*sizeof(inmax)
       call MPI_FILE_SET_VIEW(fhandl, disp, &
         MPI_DOUBLE_PRECISION, filetype, 'native', &
         fileinfo, ierr)
!
       mysize=subsizes(1)*subsizes(2)*subsizes(3)
       call MPI_File_read_all(fhandl, buf, &
         mysize, MPI_DOUBLE_PRECISION, mpi_status, ierr)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call MPI_File_close(fhandl, ierr)
       do k=1,kn-1
         do j=1,jn-1
           do i=1,in-1
             b3(i,j,k)=buf(i,j,k)
           enddo
         enddo
       enddo
!
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       fhandl=13
       call MPI_File_open(MPI_COMM_WORLD, &
         trim(dfile)//idout//'.dat', &
         MPI_MODE_RDONLY, &
         fileinfo, fhandl, ierr)
!
       if(myid .eq. 0) then
         call MPI_File_read(fhandl, time_in, 1, &
           MPI_DOUBLE_PRECISION, mpi_status, ierr)
         call MPI_File_read(fhandl, inmax_in, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_read(fhandl, jnmax_in, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_read(fhandl, knmax_in, 1, MPI_INTEGER, &
           mpi_status, ierr)
       endif
       call MPI_BCAST(time_in,1,MPI_DOUBLE_PRECISION,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(inmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(jnmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(knmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       if(inmax_in .ne. inmax .or.  &
         jnmax_in .ne. jnmax .or.  &
         knmax_in .ne. knmax) then
          write(6,*) 'inmax,jnmax,knmax in dfile', &
            inmax_in,jnmax_in,knmax_in, &
            ' do not match inmax, jnmax, knmax in par.h'
          call MPI_ABORT(MPI_COMM_WORLD, 1,ierr)
       endif
!
       disp=sizeof(time)+3*sizeof(inmax)
       call MPI_FILE_SET_VIEW(fhandl, disp, &
         MPI_DOUBLE_PRECISION, filetype, 'native', &
         fileinfo, ierr)
!
       mysize=subsizes(1)*subsizes(2)*subsizes(3)
       call MPI_File_read_all(fhandl, buf, &
         mysize, MPI_DOUBLE_PRECISION, mpi_status, ierr)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call MPI_File_close(fhandl, ierr)
       do k=1,kn-1
         do j=1,jn-1
           do i=1,in-1
             d(i,j,k)=buf(i,j,k)
           enddo
         enddo
       enddo
!
#ifndef ISOMHD
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       fhandl=13
       call MPI_File_open(MPI_COMM_WORLD, &
         trim(efile)//idout//'.dat', &
         MPI_MODE_RDONLY, &
         fileinfo, fhandl, ierr)
!
       if(myid .eq. 0) then
         call MPI_File_read(fhandl, time_in, 1, &
           MPI_DOUBLE_PRECISION, mpi_status, ierr)
         call MPI_File_read(fhandl, inmax_in, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_read(fhandl, jnmax_in, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_read(fhandl, knmax_in, 1, MPI_INTEGER, &
           mpi_status, ierr)
       endif
       call MPI_BCAST(time_in,1,MPI_DOUBLE_PRECISION,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(inmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(jnmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(knmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       if(inmax_in .ne. inmax .or.  &
         jnmax_in .ne. jnmax .or.  &
         knmax_in .ne. knmax) then
          write(6,*) 'inmax,jnmax,knmax in efile', &
            inmax_in,jnmax_in,knmax_in, &
            ' do not match inmax, jnmax, knmax in par.h'
          call MPI_ABORT(MPI_COMM_WORLD, 1,ierr)
       endif
!
       disp=sizeof(time)+3*sizeof(inmax)
       call MPI_FILE_SET_VIEW(fhandl, disp, &
         MPI_DOUBLE_PRECISION, filetype, 'native', &
         fileinfo, ierr)
!
       mysize=subsizes(1)*subsizes(2)*subsizes(3)
       call MPI_File_read_all(fhandl, buf, &
         mysize, MPI_DOUBLE_PRECISION, mpi_status, ierr)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call MPI_File_close(fhandl, ierr)
       do k=1,kn-1
         do j=1,jn-1
           do i=1,in-1
             eint(i,j,k)=buf(i,j,k)
           enddo
         enddo
       enddo
#endif
!
#ifdef THCONDUCT
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       fhandl=13
       call MPI_File_open(MPI_COMM_WORLD, &
         trim(qfile)//idout//'.dat', &
         MPI_MODE_RDONLY, &
         fileinfo, fhandl, ierr)
!
       if(myid .eq. 0) then
         call MPI_File_read(fhandl, time_in, 1, &
           MPI_DOUBLE_PRECISION, mpi_status, ierr)
         call MPI_File_read(fhandl, inmax_in, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_read(fhandl, jnmax_in, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_read(fhandl, knmax_in, 1, MPI_INTEGER, &
           mpi_status, ierr)
       endif
       call MPI_BCAST(time_in,1,MPI_DOUBLE_PRECISION,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(inmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(jnmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(knmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       if(inmax_in .ne. inmax .or.  &
         jnmax_in .ne. jnmax .or.  &
         knmax_in .ne. knmax) then
          write(6,*) 'inmax,jnmax,knmax in qfile', &
            inmax_in,jnmax_in,knmax_in, &
            ' do not match inmax, jnmax, knmax in par.h'
          call MPI_ABORT(MPI_COMM_WORLD, 1,ierr)
       endif
       time=time_in
!
       disp=sizeof(time)+3*sizeof(inmax)
       call MPI_FILE_SET_VIEW(fhandl, disp, &
         MPI_DOUBLE_PRECISION, filetype, 'native', &
         fileinfo, ierr)
!
       mysize=subsizes(1)*subsizes(2)*subsizes(3)
       call MPI_File_read_all(fhandl, buf, &
         mysize, MPI_DOUBLE_PRECISION, mpi_status, ierr)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call MPI_File_close(fhandl, ierr)
       do k=1,kn-1
         do j=1,jn-1
           do i=1,in-1
             q(i,j,k)=buf(i,j,k)
           enddo
         enddo
       enddo
!
#endif
       call MPI_Barrier(MPI_COMM_WORLD,ierr)

       deallocate(buf)

end subroutine readdata_mpi


subroutine readrst_mpi
  use ModPar
  use ModControl
  use ModSundry
  use ModField
  implicit none

  include "mpif.h"

  integer i,j,k,counter
  integer sizes(3), subsizes(3), starts(3)
  integer filetype,fileinfo,fhandl,mysize,ierr
  integer(kind=MPI_OFFSET_KIND) disp
  integer mpi_status(MPI_STATUS_SIZE)
  character(len=3) idout
  integer inmax_in,jnmax_in,knmax_in
  real time_in
  real, allocatable :: buf(:,:,:)
!--------------------------------------------------------------
       allocate(buf(in-1,jn-1,kn-1))
!
       idout='rst'
!
       sizes(1)=inmax-1
       sizes(2)=jnmax-1
       sizes(3)=knmax-1
!
       subsizes(1)=in-1
       starts(1)=myid1*(in-5)
!
       subsizes(2)=jn-1
       starts(2)=myid2*(jn-5)
!
       subsizes(3)=kn-1
       starts(3)=myid3*(kn-5)
!
       call MPI_Type_create_subarray(3, sizes, subsizes, starts, &
         MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, filetype, ierr)
       call MPI_Type_commit(filetype, ierr)
       call MPI_Info_create(fileinfo, ierr)
!
       fhandl=13
       call MPI_File_open(MPI_COMM_WORLD, &
         trim(v1file)//idout//'.dat', &
         MPI_MODE_RDONLY, &
         fileinfo, fhandl, ierr)
!
       if(myid .eq. 0) then
         call MPI_File_read(fhandl, time_in, 1, &
           MPI_DOUBLE_PRECISION, mpi_status, ierr)
         call MPI_File_read(fhandl, inmax_in, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_read(fhandl, jnmax_in, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_read(fhandl, knmax_in, 1, MPI_INTEGER, &
           mpi_status, ierr)
       endif
       call MPI_BCAST(time_in,1,MPI_DOUBLE_PRECISION,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(inmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(jnmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(knmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       if(inmax_in .ne. inmax .or.  &
         jnmax_in .ne. jnmax .or.  &
         knmax_in .ne. knmax) then
          write(6,*) 'inmax,jnmax,knmax in v1file', &
            inmax_in,jnmax_in,knmax_in, &
            ' do not match inmax, jnmax, knmax in par.h'
          call MPI_ABORT(MPI_COMM_WORLD, 1,ierr)
       endif
       time=time_in
!
       disp=sizeof(time)+3*sizeof(inmax)
       call MPI_FILE_SET_VIEW(fhandl, disp, &
         MPI_DOUBLE_PRECISION, filetype, 'native', &
         fileinfo, ierr)
!
       mysize=subsizes(1)*subsizes(2)*subsizes(3)
       call MPI_File_read_all(fhandl, buf, &
         mysize, MPI_DOUBLE_PRECISION, mpi_status, ierr)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call MPI_File_close(fhandl, ierr)
       do k=1,kn-1
         do j=1,jn-1
           do i=1,in-1
             v1(i,j,k)=buf(i,j,k)
           enddo
         enddo
       enddo
!
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       fhandl=13
       call MPI_File_open(MPI_COMM_WORLD, &
         trim(v2file)//idout//'.dat', &
         MPI_MODE_RDONLY, &
         fileinfo, fhandl, ierr)
!
       if(myid .eq. 0) then
         call MPI_File_read(fhandl, time_in, 1, &
           MPI_DOUBLE_PRECISION, mpi_status, ierr)
         call MPI_File_read(fhandl, inmax_in, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_read(fhandl, jnmax_in, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_read(fhandl, knmax_in, 1, MPI_INTEGER, &
           mpi_status, ierr)
       endif
       call MPI_BCAST(time_in,1,MPI_DOUBLE_PRECISION,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(inmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(jnmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(knmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       if(inmax_in .ne. inmax .or.  &
         jnmax_in .ne. jnmax .or.  &
         knmax_in .ne. knmax) then
          write(6,*) 'inmax,jnmax,knmax in v2file', &
            inmax_in,jnmax_in,knmax_in, &
            ' do not match inmax, jnmax, knmax in par.h'
          call MPI_ABORT(MPI_COMM_WORLD, 1,ierr)
       endif
!
       disp=sizeof(time)+3*sizeof(inmax)
       call MPI_FILE_SET_VIEW(fhandl, disp, &
         MPI_DOUBLE_PRECISION, filetype, 'native', &
         fileinfo, ierr)
!
       mysize=subsizes(1)*subsizes(2)*subsizes(3)
       call MPI_File_read_all(fhandl, buf, &
         mysize, MPI_DOUBLE_PRECISION, mpi_status, ierr)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call MPI_File_close(fhandl, ierr)
       do k=1,kn-1
         do j=1,jn-1
           do i=1,in-1
             v2(i,j,k)=buf(i,j,k)
           enddo
         enddo
       enddo
!
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       fhandl=13
       call MPI_File_open(MPI_COMM_WORLD, &
         trim(v3file)//idout//'.dat', &
         MPI_MODE_RDONLY, &
         fileinfo, fhandl, ierr)
!
       if(myid .eq. 0) then
         call MPI_File_read(fhandl, time_in, 1, &
           MPI_DOUBLE_PRECISION, mpi_status, ierr)
         call MPI_File_read(fhandl, inmax_in, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_read(fhandl, jnmax_in, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_read(fhandl, knmax_in, 1, MPI_INTEGER, &
           mpi_status, ierr)
       endif
       call MPI_BCAST(time_in,1,MPI_DOUBLE_PRECISION,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(inmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(jnmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(knmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       if(inmax_in .ne. inmax .or.  &
         jnmax_in .ne. jnmax .or.  &
         knmax_in .ne. knmax) then
          write(6,*) 'inmax,jnmax,knmax in v3file', &
            inmax_in,jnmax_in,knmax_in, &
            ' do not match inmax, jnmax, knmax in par.h'
          call MPI_ABORT(MPI_COMM_WORLD, 1,ierr)
       endif
!
       disp=sizeof(time)+3*sizeof(inmax)
       call MPI_FILE_SET_VIEW(fhandl, disp, &
         MPI_DOUBLE_PRECISION, filetype, 'native', &
         fileinfo, ierr)
!
       mysize=subsizes(1)*subsizes(2)*subsizes(3)
       call MPI_File_read_all(fhandl, buf, &
         mysize, MPI_DOUBLE_PRECISION, mpi_status, ierr)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call MPI_File_close(fhandl, ierr)
       do k=1,kn-1
         do j=1,jn-1
           do i=1,in-1
             v3(i,j,k)=buf(i,j,k)
           enddo
         enddo
       enddo
!
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       fhandl=13
       call MPI_File_open(MPI_COMM_WORLD, &
         trim(b1file)//idout//'.dat', &
         MPI_MODE_RDONLY, &
         fileinfo, fhandl, ierr)
!
       if(myid .eq. 0) then
         call MPI_File_read(fhandl, time_in, 1, &
           MPI_DOUBLE_PRECISION, mpi_status, ierr)
         call MPI_File_read(fhandl, inmax_in, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_read(fhandl, jnmax_in, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_read(fhandl, knmax_in, 1, MPI_INTEGER, &
           mpi_status, ierr)
       endif
       call MPI_BCAST(time_in,1,MPI_DOUBLE_PRECISION,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(inmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(jnmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(knmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       if(inmax_in .ne. inmax .or.  &
         jnmax_in .ne. jnmax .or.  &
         knmax_in .ne. knmax) then
          write(6,*) 'inmax,jnmax,knmax in b1file', &
            inmax_in,jnmax_in,knmax_in, &
            ' do not match inmax, jnmax, knmax in par.h'
          call MPI_ABORT(MPI_COMM_WORLD, 1,ierr)
       endif
!
       disp=sizeof(time)+3*sizeof(inmax)
       call MPI_FILE_SET_VIEW(fhandl, disp, &
         MPI_DOUBLE_PRECISION, filetype, 'native', &
         fileinfo, ierr)
!
       mysize=subsizes(1)*subsizes(2)*subsizes(3)
       call MPI_File_read_all(fhandl, buf, &
         mysize, MPI_DOUBLE_PRECISION, mpi_status, ierr)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call MPI_File_close(fhandl, ierr)
       do k=1,kn-1
         do j=1,jn-1
           do i=1,in-1
             b1(i,j,k)=buf(i,j,k)
           enddo
         enddo
       enddo
!
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       fhandl=13
       call MPI_File_open(MPI_COMM_WORLD, &
         trim(b2file)//idout//'.dat', &
         MPI_MODE_RDONLY, &
         fileinfo, fhandl, ierr)
!
       if(myid .eq. 0) then
         call MPI_File_read(fhandl, time_in, 1, &
           MPI_DOUBLE_PRECISION, mpi_status, ierr)
         call MPI_File_read(fhandl, inmax_in, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_read(fhandl, jnmax_in, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_read(fhandl, knmax_in, 1, MPI_INTEGER, &
           mpi_status, ierr)
       endif
       call MPI_BCAST(time_in,1,MPI_DOUBLE_PRECISION,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(inmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(jnmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(knmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       if(inmax_in .ne. inmax .or.  &
         jnmax_in .ne. jnmax .or.  &
         knmax_in .ne. knmax) then
          write(6,*) 'inmax,jnmax,knmax in b2file', &
            inmax_in,jnmax_in,knmax_in, &
            ' do not match inmax, jnmax, knmax in par.h'
          call MPI_ABORT(MPI_COMM_WORLD, 1,ierr)
       endif
!
       disp=sizeof(time)+3*sizeof(inmax)
       call MPI_FILE_SET_VIEW(fhandl, disp, &
         MPI_DOUBLE_PRECISION, filetype, 'native', &
         fileinfo, ierr)
!
       mysize=subsizes(1)*subsizes(2)*subsizes(3)
       call MPI_File_read_all(fhandl, buf, &
         mysize, MPI_DOUBLE_PRECISION, mpi_status, ierr)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call MPI_File_close(fhandl, ierr)
       do k=1,kn-1
         do j=1,jn-1
           do i=1,in-1
             b2(i,j,k)=buf(i,j,k)
           enddo
         enddo
       enddo
!
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       fhandl=13
       call MPI_File_open(MPI_COMM_WORLD, &
         trim(b3file)//idout//'.dat', &
         MPI_MODE_RDONLY, &
         fileinfo, fhandl, ierr)
!
       if(myid .eq. 0) then
         call MPI_File_read(fhandl, time_in, 1, &
           MPI_DOUBLE_PRECISION, mpi_status, ierr)
         call MPI_File_read(fhandl, inmax_in, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_read(fhandl, jnmax_in, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_read(fhandl, knmax_in, 1, MPI_INTEGER, &
           mpi_status, ierr)
       endif
       call MPI_BCAST(time_in,1,MPI_DOUBLE_PRECISION,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(inmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(jnmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(knmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       if(inmax_in .ne. inmax .or.  &
         jnmax_in .ne. jnmax .or.  &
         knmax_in .ne. knmax) then
          write(6,*) 'inmax,jnmax,knmax in b3file', &
            inmax_in,jnmax_in,knmax_in, &
            ' do not match inmax, jnmax, knmax in par.h'
          call MPI_ABORT(MPI_COMM_WORLD, 1,ierr)
       endif
!
       disp=sizeof(time)+3*sizeof(inmax)
       call MPI_FILE_SET_VIEW(fhandl, disp, &
         MPI_DOUBLE_PRECISION, filetype, 'native', &
         fileinfo, ierr)
!
       mysize=subsizes(1)*subsizes(2)*subsizes(3)
       call MPI_File_read_all(fhandl, buf, &
         mysize, MPI_DOUBLE_PRECISION, mpi_status, ierr)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call MPI_File_close(fhandl, ierr)
       do k=1,kn-1
         do j=1,jn-1
           do i=1,in-1
             b3(i,j,k)=buf(i,j,k)
           enddo
         enddo
       enddo
!
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       fhandl=13
       call MPI_File_open(MPI_COMM_WORLD, &
         trim(dfile)//idout//'.dat', &
         MPI_MODE_RDONLY, &
         fileinfo, fhandl, ierr)
!
       if(myid .eq. 0) then
         call MPI_File_read(fhandl, time_in, 1, &
           MPI_DOUBLE_PRECISION, mpi_status, ierr)
         call MPI_File_read(fhandl, inmax_in, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_read(fhandl, jnmax_in, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_read(fhandl, knmax_in, 1, MPI_INTEGER, &
           mpi_status, ierr)
       endif
       call MPI_BCAST(time_in,1,MPI_DOUBLE_PRECISION,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(inmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(jnmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(knmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       if(inmax_in .ne. inmax .or.  &
         jnmax_in .ne. jnmax .or.  &
         knmax_in .ne. knmax) then
          write(6,*) 'inmax,jnmax,knmax in dfile', &
            inmax_in,jnmax_in,knmax_in, &
            ' do not match inmax, jnmax, knmax in par.h'
          call MPI_ABORT(MPI_COMM_WORLD, 1,ierr)
       endif
!
       disp=sizeof(time)+3*sizeof(inmax)
       call MPI_FILE_SET_VIEW(fhandl, disp, &
         MPI_DOUBLE_PRECISION, filetype, 'native', &
         fileinfo, ierr)
!
       mysize=subsizes(1)*subsizes(2)*subsizes(3)
       call MPI_File_read_all(fhandl, buf, &
         mysize, MPI_DOUBLE_PRECISION, mpi_status, ierr)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call MPI_File_close(fhandl, ierr)
       do k=1,kn-1
         do j=1,jn-1
           do i=1,in-1
             d(i,j,k)=buf(i,j,k)
           enddo
         enddo
       enddo
!
#ifndef ISOMHD
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       fhandl=13
       call MPI_File_open(MPI_COMM_WORLD, &
         trim(efile)//idout//'.dat', &
         MPI_MODE_RDONLY, &
         fileinfo, fhandl, ierr)
!
       if(myid .eq. 0) then
         call MPI_File_read(fhandl, time_in, 1, &
           MPI_DOUBLE_PRECISION, mpi_status, ierr)
         call MPI_File_read(fhandl, inmax_in, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_read(fhandl, jnmax_in, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_read(fhandl, knmax_in, 1, MPI_INTEGER, &
           mpi_status, ierr)
       endif
       call MPI_BCAST(time_in,1,MPI_DOUBLE_PRECISION,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(inmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(jnmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(knmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       if(inmax_in .ne. inmax .or.  &
         jnmax_in .ne. jnmax .or.  &
         knmax_in .ne. knmax) then
          write(6,*) 'inmax,jnmax,knmax in efile', &
            inmax_in,jnmax_in,knmax_in, &
            ' do not match inmax, jnmax, knmax in par.h'
          call MPI_ABORT(MPI_COMM_WORLD, 1,ierr)
       endif
!
       disp=sizeof(time)+3*sizeof(inmax)
       call MPI_FILE_SET_VIEW(fhandl, disp, &
         MPI_DOUBLE_PRECISION, filetype, 'native', &
         fileinfo, ierr)
!
       mysize=subsizes(1)*subsizes(2)*subsizes(3)
       call MPI_File_read_all(fhandl, buf, &
         mysize, MPI_DOUBLE_PRECISION, mpi_status, ierr)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call MPI_File_close(fhandl, ierr)
       do k=1,kn-1
         do j=1,jn-1
           do i=1,in-1
             eint(i,j,k)=buf(i,j,k)
           enddo
         enddo
       enddo
#endif
!
#ifdef THCONDUCT
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       fhandl=13
       call MPI_File_open(MPI_COMM_WORLD, &
         trim(qfile)//idout//'.dat', &
         MPI_MODE_RDONLY, &
         fileinfo, fhandl, ierr)
!
       if(myid .eq. 0) then
         call MPI_File_read(fhandl, time_in, 1, &
           MPI_DOUBLE_PRECISION, mpi_status, ierr)
         call MPI_File_read(fhandl, inmax_in, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_read(fhandl, jnmax_in, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_read(fhandl, knmax_in, 1, MPI_INTEGER, &
           mpi_status, ierr)
       endif
       call MPI_BCAST(time_in,1,MPI_DOUBLE_PRECISION,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(inmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(jnmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(knmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       if(inmax_in .ne. inmax .or.  &
         jnmax_in .ne. jnmax .or.  &
         knmax_in .ne. knmax) then
          write(6,*) 'inmax,jnmax,knmax in qfile', &
            inmax_in,jnmax_in,knmax_in, &
            ' do not match inmax, jnmax, knmax in par.h'
          call MPI_ABORT(MPI_COMM_WORLD, 1,ierr)
       endif
       time=time_in
!
       disp=sizeof(time)+3*sizeof(inmax)
       call MPI_FILE_SET_VIEW(fhandl, disp, &
         MPI_DOUBLE_PRECISION, filetype, 'native', &
         fileinfo, ierr)
!
       mysize=subsizes(1)*subsizes(2)*subsizes(3)
       call MPI_File_read_all(fhandl, buf, &
         mysize, MPI_DOUBLE_PRECISION, mpi_status, ierr)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call MPI_File_close(fhandl, ierr)
       do k=1,kn-1
         do j=1,jn-1
           do i=1,in-1
             q(i,j,k)=buf(i,j,k)
           enddo
         enddo
       enddo
!
#endif

       call MPI_Barrier(MPI_COMM_WORLD,ierr)

       deallocate(buf)

end subroutine readrst_mpi


subroutine writedata_mpi
  use ModPar
  use ModGrid
  use ModControl
  use ModSundry
  use ModField
  implicit none

  include "mpif.h"

  integer i,j,k,counter
  integer sizes(3), subsizes(3), starts(3)
  integer filetype,fileinfo,fhandl,mysize,ierr
  integer(kind=MPI_OFFSET_KIND) disp
  integer mpi_status(MPI_STATUS_SIZE)
  character(len=4) idout
  real, allocatable :: buf(:)
!-----------------------------------------------------
       allocate(buf(in*jn*kn))
!
       write(idout,'(i4.4)') ifile
!
       sizes(1)=inmax-1
       sizes(2)=jnmax-1
       sizes(3)=knmax-1
!
       subsizes(1)=in-5
       starts(1)=is-1+myid1*(in-5)
       if(myid1 .eq. 0) then
         subsizes(1)=subsizes(1)+2
         starts(1)=0
       endif
       if(myid1 .eq. nproc1-1) then
         subsizes(1)=subsizes(1)+2
       endif
!
       subsizes(2)=jn-5
       starts(2)=js-1+myid2*(jn-5)
       if(myid2 .eq. 0) then
         subsizes(2)=subsizes(2)+2
         starts(2)=0
       endif
       if(myid2 .eq. nproc2-1) then
         subsizes(2)=subsizes(2)+2
       endif
!
       subsizes(3)=kn-5
       starts(3)=ks-1+myid3*(kn-5)
       if(myid3 .eq. 0) then
         subsizes(3)=subsizes(3)+2
         starts(3)=0
       endif
       if(myid3 .eq. nproc3-1) then
         subsizes(3)=subsizes(3)+2
       endif
!
       call MPI_Type_create_subarray(3, sizes, subsizes, starts, &
         MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, filetype, ierr)
       call MPI_Type_commit(filetype, ierr)
       call MPI_Info_create(fileinfo, ierr)
!
       fhandl=13
       call MPI_File_open(MPI_COMM_WORLD, &
         trim(v1file)//idout//'.dat', &
         IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), &
         fileinfo, fhandl, ierr)
!
       if(myid .eq. 0) then
         call MPI_File_write(fhandl, time, 1, &
           MPI_DOUBLE_PRECISION, mpi_status, ierr)
         call MPI_File_write(fhandl, inmax, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_write(fhandl, jnmax, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_write(fhandl, knmax, 1, MPI_INTEGER, &
           mpi_status, ierr)
       endif
!
       disp=sizeof(time)+3*sizeof(inmax)
       call MPI_FILE_SET_VIEW(fhandl, disp, &
         MPI_DOUBLE_PRECISION, filetype, 'native', &
         fileinfo, ierr)
!
       mysize=subsizes(1)*subsizes(2)*subsizes(3)
       do k=starts(3)+1,starts(3)+subsizes(3)
         do j=starts(2)+1,starts(2)+subsizes(2)
           do i=starts(1)+1,starts(1)+subsizes(1)
             counter=1+i-starts(1)-1 &
               +(j-starts(2)-1)*subsizes(1) &
               +(k-starts(3)-1)*subsizes(1)*subsizes(2)
             buf(counter)=v1(i-myid1*(in-5),j-myid2*(jn-5),k-myid3*(kn-5))
           enddo
         enddo
       enddo
       call MPI_File_write_all(fhandl, buf, &
         mysize, MPI_DOUBLE_PRECISION, mpi_status, ierr)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call MPI_File_close(fhandl, ierr)
!
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       fhandl=13
       call MPI_File_open(MPI_COMM_WORLD, &
         trim(v2file)//idout//'.dat', &
         IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), &
         fileinfo, fhandl, ierr)
!
       if(myid .eq. 0) then
         call MPI_File_write(fhandl, time, 1, &
           MPI_DOUBLE_PRECISION, mpi_status, ierr)
         call MPI_File_write(fhandl, inmax, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_write(fhandl, jnmax, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_write(fhandl, knmax, 1, MPI_INTEGER, &
           mpi_status, ierr)
       endif
!
       disp=sizeof(time)+3*sizeof(inmax)
       call MPI_FILE_SET_VIEW(fhandl, disp, &
         MPI_DOUBLE_PRECISION, filetype, 'native', &
         fileinfo, ierr)
!
       mysize=subsizes(1)*subsizes(2)*subsizes(3)
       do k=starts(3)+1,starts(3)+subsizes(3)
         do j=starts(2)+1,starts(2)+subsizes(2)
           do i=starts(1)+1,starts(1)+subsizes(1)
             counter=1+i-starts(1)-1 &
               +(j-starts(2)-1)*subsizes(1) &
               +(k-starts(3)-1)*subsizes(1)*subsizes(2)
             buf(counter)=v2(i-myid1*(in-5),j-myid2*(jn-5),k-myid3*(kn-5))
           enddo
         enddo
       enddo
       call MPI_File_write_all(fhandl, buf, &
         mysize, MPI_DOUBLE_PRECISION, mpi_status, ierr)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call MPI_File_close(fhandl, ierr)
!
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       fhandl=13
       call MPI_File_open(MPI_COMM_WORLD, &
         trim(v3file)//idout//'.dat', &
         IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), &
         fileinfo, fhandl, ierr)
!
       if(myid .eq. 0) then
         call MPI_File_write(fhandl, time, 1, &
           MPI_DOUBLE_PRECISION, mpi_status, ierr)
         call MPI_File_write(fhandl, inmax, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_write(fhandl, jnmax, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_write(fhandl, knmax, 1, MPI_INTEGER, &
           mpi_status, ierr)
       endif
!
       disp=sizeof(time)+3*sizeof(inmax)
       call MPI_FILE_SET_VIEW(fhandl, disp, &
         MPI_DOUBLE_PRECISION, filetype, 'native', &
         fileinfo, ierr)
!
       mysize=subsizes(1)*subsizes(2)*subsizes(3)
       do k=starts(3)+1,starts(3)+subsizes(3)
         do j=starts(2)+1,starts(2)+subsizes(2)
           do i=starts(1)+1,starts(1)+subsizes(1)
             counter=1+i-starts(1)-1 &
               +(j-starts(2)-1)*subsizes(1) &
               +(k-starts(3)-1)*subsizes(1)*subsizes(2)
             buf(counter)=v3(i-myid1*(in-5),j-myid2*(jn-5),k-myid3*(kn-5))
           enddo
         enddo
       enddo
       call MPI_File_write_all(fhandl, buf, &
         mysize, MPI_DOUBLE_PRECISION, mpi_status, ierr)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call MPI_File_close(fhandl, ierr)
!
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       fhandl=13
       call MPI_File_open(MPI_COMM_WORLD, &
         trim(b1file)//idout//'.dat', &
         IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), &
         fileinfo, fhandl, ierr)
!
       if(myid .eq. 0) then
         call MPI_File_write(fhandl, time, 1, &
           MPI_DOUBLE_PRECISION, mpi_status, ierr)
         call MPI_File_write(fhandl, inmax, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_write(fhandl, jnmax, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_write(fhandl, knmax, 1, MPI_INTEGER, &
           mpi_status, ierr)
       endif
!
       disp=sizeof(time)+3*sizeof(inmax)
       call MPI_FILE_SET_VIEW(fhandl, disp, &
         MPI_DOUBLE_PRECISION, filetype, 'native', &
         fileinfo, ierr)
!
       mysize=subsizes(1)*subsizes(2)*subsizes(3)
       do k=starts(3)+1,starts(3)+subsizes(3)
         do j=starts(2)+1,starts(2)+subsizes(2)
           do i=starts(1)+1,starts(1)+subsizes(1)
             counter=1+i-starts(1)-1 &
               +(j-starts(2)-1)*subsizes(1) &
               +(k-starts(3)-1)*subsizes(1)*subsizes(2)
             buf(counter)=b1(i-myid1*(in-5),j-myid2*(jn-5),k-myid3*(kn-5))
           enddo
         enddo
       enddo
       call MPI_File_write_all(fhandl, buf, &
         mysize, MPI_DOUBLE_PRECISION, mpi_status, ierr)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call MPI_File_close(fhandl, ierr)
!
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       fhandl=13
       call MPI_File_open(MPI_COMM_WORLD, &
         trim(b2file)//idout//'.dat', &
         IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), &
         fileinfo, fhandl, ierr)
!
       if(myid .eq. 0) then
         call MPI_File_write(fhandl, time, 1, &
           MPI_DOUBLE_PRECISION, mpi_status, ierr)
         call MPI_File_write(fhandl, inmax, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_write(fhandl, jnmax, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_write(fhandl, knmax, 1, MPI_INTEGER, &
           mpi_status, ierr)
       endif
!
       disp=sizeof(time)+3*sizeof(inmax)
       call MPI_FILE_SET_VIEW(fhandl, disp, &
         MPI_DOUBLE_PRECISION, filetype, 'native', &
         fileinfo, ierr)
!
       mysize=subsizes(1)*subsizes(2)*subsizes(3)
       do k=starts(3)+1,starts(3)+subsizes(3)
         do j=starts(2)+1,starts(2)+subsizes(2)
           do i=starts(1)+1,starts(1)+subsizes(1)
             counter=1+i-starts(1)-1 &
               +(j-starts(2)-1)*subsizes(1) &
               +(k-starts(3)-1)*subsizes(1)*subsizes(2)
             buf(counter)=b2(i-myid1*(in-5),j-myid2*(jn-5),k-myid3*(kn-5))
           enddo
         enddo
       enddo
       call MPI_File_write_all(fhandl, buf, &
         mysize, MPI_DOUBLE_PRECISION, mpi_status, ierr)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call MPI_File_close(fhandl, ierr)
!
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       fhandl=13
       call MPI_File_open(MPI_COMM_WORLD, &
         trim(b3file)//idout//'.dat', &
         IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), &
         fileinfo, fhandl, ierr)
!
       if(myid .eq. 0) then
         call MPI_File_write(fhandl, time, 1, &
           MPI_DOUBLE_PRECISION, mpi_status, ierr)
         call MPI_File_write(fhandl, inmax, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_write(fhandl, jnmax, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_write(fhandl, knmax, 1, MPI_INTEGER, &
           mpi_status, ierr)
       endif
!
       disp=sizeof(time)+3*sizeof(inmax)
       call MPI_FILE_SET_VIEW(fhandl, disp, &
         MPI_DOUBLE_PRECISION, filetype, 'native', &
         fileinfo, ierr)
!
       mysize=subsizes(1)*subsizes(2)*subsizes(3)
       do k=starts(3)+1,starts(3)+subsizes(3)
         do j=starts(2)+1,starts(2)+subsizes(2)
           do i=starts(1)+1,starts(1)+subsizes(1)
             counter=1+i-starts(1)-1 &
               +(j-starts(2)-1)*subsizes(1) &
               +(k-starts(3)-1)*subsizes(1)*subsizes(2)
             buf(counter)=b3(i-myid1*(in-5),j-myid2*(jn-5),k-myid3*(kn-5))
           enddo
         enddo
       enddo
       call MPI_File_write_all(fhandl, buf, &
         mysize, MPI_DOUBLE_PRECISION, mpi_status, ierr)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call MPI_File_close(fhandl, ierr)
!
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       fhandl=13
       call MPI_File_open(MPI_COMM_WORLD, &
         trim(dfile)//idout//'.dat', &
         IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), &
         fileinfo, fhandl, ierr)
!
       if(myid .eq. 0) then
         call MPI_File_write(fhandl, time, 1, &
           MPI_DOUBLE_PRECISION, mpi_status, ierr)
         call MPI_File_write(fhandl, inmax, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_write(fhandl, jnmax, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_write(fhandl, knmax, 1, MPI_INTEGER, &
           mpi_status, ierr)
       endif
!
       disp=sizeof(time)+3*sizeof(inmax)
       call MPI_FILE_SET_VIEW(fhandl, disp, &
         MPI_DOUBLE_PRECISION, filetype, 'native', &
         fileinfo, ierr)
!
       mysize=subsizes(1)*subsizes(2)*subsizes(3)
       do k=starts(3)+1,starts(3)+subsizes(3)
         do j=starts(2)+1,starts(2)+subsizes(2)
           do i=starts(1)+1,starts(1)+subsizes(1)
             counter=1+i-starts(1)-1 &
               +(j-starts(2)-1)*subsizes(1) &
               +(k-starts(3)-1)*subsizes(1)*subsizes(2)
             buf(counter)=d(i-myid1*(in-5),j-myid2*(jn-5),k-myid3*(kn-5))
           enddo
         enddo
       enddo
       call MPI_File_write_all(fhandl, buf, &
         mysize, MPI_DOUBLE_PRECISION, mpi_status, ierr)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call MPI_File_close(fhandl, ierr)
!
#ifndef ISOMHD
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       fhandl=13
       call MPI_File_open(MPI_COMM_WORLD, &
         trim(efile)//idout//'.dat', &
         IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), &
         fileinfo, fhandl, ierr)
!
       if(myid .eq. 0) then
         call MPI_File_write(fhandl, time, 1, &
           MPI_DOUBLE_PRECISION, mpi_status, ierr)
         call MPI_File_write(fhandl, inmax, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_write(fhandl, jnmax, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_write(fhandl, knmax, 1, MPI_INTEGER, &
           mpi_status, ierr)
       endif
!
       disp=sizeof(time)+3*sizeof(inmax)
       call MPI_FILE_SET_VIEW(fhandl, disp, &
         MPI_DOUBLE_PRECISION, filetype, 'native', &
         fileinfo, ierr)
!
       mysize=subsizes(1)*subsizes(2)*subsizes(3)
       do k=starts(3)+1,starts(3)+subsizes(3)
         do j=starts(2)+1,starts(2)+subsizes(2)
           do i=starts(1)+1,starts(1)+subsizes(1)
             counter=1+i-starts(1)-1 &
               +(j-starts(2)-1)*subsizes(1) &
               +(k-starts(3)-1)*subsizes(1)*subsizes(2)
             buf(counter)=eint(i-myid1*(in-5),j-myid2*(jn-5),k-myid3*(kn-5))
           enddo
         enddo
       enddo
       call MPI_File_write_all(fhandl, buf, &
         mysize, MPI_DOUBLE_PRECISION, mpi_status, ierr)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call MPI_File_close(fhandl, ierr)
#endif
!
#ifdef THCONDUCT
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       fhandl=13
       call MPI_File_open(MPI_COMM_WORLD, &
         trim(qfile)//idout//'.dat', &
         IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), &
         fileinfo, fhandl, ierr)
!
       if(myid .eq. 0) then
         call MPI_File_write(fhandl, time, 1, &
           MPI_DOUBLE_PRECISION, mpi_status, ierr)
         call MPI_File_write(fhandl, inmax, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_write(fhandl, jnmax, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_write(fhandl, knmax, 1, MPI_INTEGER, &
           mpi_status, ierr)
       endif
!
       disp=sizeof(time)+3*sizeof(inmax)
       call MPI_FILE_SET_VIEW(fhandl, disp, &
         MPI_DOUBLE_PRECISION, filetype, 'native', &
         fileinfo, ierr)
!
       mysize=subsizes(1)*subsizes(2)*subsizes(3)
       do k=starts(3)+1,starts(3)+subsizes(3)
         do j=starts(2)+1,starts(2)+subsizes(2)
           do i=starts(1)+1,starts(1)+subsizes(1)
             counter=1+i-starts(1)-1 &
               +(j-starts(2)-1)*subsizes(1) &
               +(k-starts(3)-1)*subsizes(1)*subsizes(2)
             buf(counter)=q(i-myid1*(in-5),j-myid2*(jn-5),k-myid3*(kn-5))
           enddo
         enddo
       enddo
       call MPI_File_write_all(fhandl, buf, &
         mysize, MPI_DOUBLE_PRECISION, mpi_status, ierr)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call MPI_File_close(fhandl, ierr)
#endif
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
       deallocate(buf)

end subroutine writedata_mpi


subroutine writerst_mpi
  use ModPar
  use ModGrid
  use ModControl
  use ModSundry
  use ModField
  implicit none

  include "mpif.h"

  integer i,j,k,counter
  integer sizes(3), subsizes(3), starts(3)
  integer filetype,fileinfo,fhandl,mysize,ierr
  integer(kind=MPI_OFFSET_KIND) disp
  integer mpi_status(MPI_STATUS_SIZE)
  character(len=3) idout
  real, allocatable :: buf(:)
!-----------------------------------------------------
       allocate(buf(in*jn*kn))
!
       idout='rst'
!
       sizes(1)=inmax-1
       sizes(2)=jnmax-1
       sizes(3)=knmax-1
!
       subsizes(1)=in-5
       starts(1)=is-1+myid1*(in-5)
       if(myid1 .eq. 0) then
         subsizes(1)=subsizes(1)+2
         starts(1)=0
       endif
       if(myid1 .eq. nproc1-1) then
         subsizes(1)=subsizes(1)+2
       endif
!
       subsizes(2)=jn-5
       starts(2)=js-1+myid2*(jn-5)
       if(myid2 .eq. 0) then
         subsizes(2)=subsizes(2)+2
         starts(2)=0
       endif
       if(myid2 .eq. nproc2-1) then
         subsizes(2)=subsizes(2)+2
       endif
!
       subsizes(3)=kn-5
       starts(3)=ks-1+myid3*(kn-5)
       if(myid3 .eq. 0) then
         subsizes(3)=subsizes(3)+2
         starts(3)=0
       endif
       if(myid3 .eq. nproc3-1) then
         subsizes(3)=subsizes(3)+2
       endif
!
       call MPI_Type_create_subarray(3, sizes, subsizes, starts, &
         MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, filetype, ierr)
       call MPI_Type_commit(filetype, ierr)
       call MPI_Info_create(fileinfo, ierr)
!
       fhandl=13
       call MPI_File_open(MPI_COMM_WORLD, &
         trim(v1file)//idout//'.dat', &
         IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), &
         fileinfo, fhandl, ierr)
!
       if(myid .eq. 0) then
         call MPI_File_write(fhandl, time, 1, &
           MPI_DOUBLE_PRECISION, mpi_status, ierr)
         call MPI_File_write(fhandl, inmax, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_write(fhandl, jnmax, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_write(fhandl, knmax, 1, MPI_INTEGER, &
           mpi_status, ierr)
       endif
!
       disp=sizeof(time)+3*sizeof(inmax)
       call MPI_FILE_SET_VIEW(fhandl, disp, &
         MPI_DOUBLE_PRECISION, filetype, 'native', &
         fileinfo, ierr)
!
       mysize=subsizes(1)*subsizes(2)*subsizes(3)
       do k=starts(3)+1,starts(3)+subsizes(3)
         do j=starts(2)+1,starts(2)+subsizes(2)
           do i=starts(1)+1,starts(1)+subsizes(1)
             counter=1+i-starts(1)-1 &
               +(j-starts(2)-1)*subsizes(1) &
               +(k-starts(3)-1)*subsizes(1)*subsizes(2)
             buf(counter)=v1(i-myid1*(in-5),j-myid2*(jn-5),k-myid3*(kn-5))
           enddo
         enddo
       enddo
       call MPI_File_write_all(fhandl, buf, &
         mysize, MPI_DOUBLE_PRECISION, mpi_status, ierr)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call MPI_File_close(fhandl, ierr)
!
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       fhandl=13
       call MPI_File_open(MPI_COMM_WORLD, &
         trim(v2file)//idout//'.dat', &
         IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), &
         fileinfo, fhandl, ierr)
!
       if(myid .eq. 0) then
         call MPI_File_write(fhandl, time, 1, &
           MPI_DOUBLE_PRECISION, mpi_status, ierr)
         call MPI_File_write(fhandl, inmax, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_write(fhandl, jnmax, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_write(fhandl, knmax, 1, MPI_INTEGER, &
           mpi_status, ierr)
       endif
!
       disp=sizeof(time)+3*sizeof(inmax)
       call MPI_FILE_SET_VIEW(fhandl, disp, &
         MPI_DOUBLE_PRECISION, filetype, 'native', &
         fileinfo, ierr)
!
       mysize=subsizes(1)*subsizes(2)*subsizes(3)
       do k=starts(3)+1,starts(3)+subsizes(3)
         do j=starts(2)+1,starts(2)+subsizes(2)
           do i=starts(1)+1,starts(1)+subsizes(1)
             counter=1+i-starts(1)-1 &
               +(j-starts(2)-1)*subsizes(1) &
               +(k-starts(3)-1)*subsizes(1)*subsizes(2)
             buf(counter)=v2(i-myid1*(in-5),j-myid2*(jn-5),k-myid3*(kn-5))
           enddo
         enddo
       enddo
       call MPI_File_write_all(fhandl, buf, &
         mysize, MPI_DOUBLE_PRECISION, mpi_status, ierr)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call MPI_File_close(fhandl, ierr)
!
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       fhandl=13
       call MPI_File_open(MPI_COMM_WORLD, &
         trim(v3file)//idout//'.dat', &
         IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), &
         fileinfo, fhandl, ierr)
!
       if(myid .eq. 0) then
         call MPI_File_write(fhandl, time, 1, &
           MPI_DOUBLE_PRECISION, mpi_status, ierr)
         call MPI_File_write(fhandl, inmax, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_write(fhandl, jnmax, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_write(fhandl, knmax, 1, MPI_INTEGER, &
           mpi_status, ierr)
       endif
!
       disp=sizeof(time)+3*sizeof(inmax)
       call MPI_FILE_SET_VIEW(fhandl, disp, &
         MPI_DOUBLE_PRECISION, filetype, 'native', &
         fileinfo, ierr)
!
       mysize=subsizes(1)*subsizes(2)*subsizes(3)
       do k=starts(3)+1,starts(3)+subsizes(3)
         do j=starts(2)+1,starts(2)+subsizes(2)
           do i=starts(1)+1,starts(1)+subsizes(1)
             counter=1+i-starts(1)-1 &
               +(j-starts(2)-1)*subsizes(1) &
               +(k-starts(3)-1)*subsizes(1)*subsizes(2)
             buf(counter)=v3(i-myid1*(in-5),j-myid2*(jn-5),k-myid3*(kn-5))
           enddo
         enddo
       enddo
       call MPI_File_write_all(fhandl, buf, &
         mysize, MPI_DOUBLE_PRECISION, mpi_status, ierr)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call MPI_File_close(fhandl, ierr)
!
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       fhandl=13
       call MPI_File_open(MPI_COMM_WORLD, &
         trim(b1file)//idout//'.dat', &
         IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), &
         fileinfo, fhandl, ierr)
!
       if(myid .eq. 0) then
         call MPI_File_write(fhandl, time, 1, &
           MPI_DOUBLE_PRECISION, mpi_status, ierr)
         call MPI_File_write(fhandl, inmax, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_write(fhandl, jnmax, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_write(fhandl, knmax, 1, MPI_INTEGER, &
           mpi_status, ierr)
       endif
!
       disp=sizeof(time)+3*sizeof(inmax)
       call MPI_FILE_SET_VIEW(fhandl, disp, &
         MPI_DOUBLE_PRECISION, filetype, 'native', &
         fileinfo, ierr)
!
       mysize=subsizes(1)*subsizes(2)*subsizes(3)
       do k=starts(3)+1,starts(3)+subsizes(3)
         do j=starts(2)+1,starts(2)+subsizes(2)
           do i=starts(1)+1,starts(1)+subsizes(1)
             counter=1+i-starts(1)-1 &
               +(j-starts(2)-1)*subsizes(1) &
               +(k-starts(3)-1)*subsizes(1)*subsizes(2)
             buf(counter)=b1(i-myid1*(in-5),j-myid2*(jn-5),k-myid3*(kn-5))
           enddo
         enddo
       enddo
       call MPI_File_write_all(fhandl, buf, &
         mysize, MPI_DOUBLE_PRECISION, mpi_status, ierr)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call MPI_File_close(fhandl, ierr)
!
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       fhandl=13
       call MPI_File_open(MPI_COMM_WORLD, &
         trim(b2file)//idout//'.dat', &
         IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), &
         fileinfo, fhandl, ierr)
!
       if(myid .eq. 0) then
         call MPI_File_write(fhandl, time, 1, &
           MPI_DOUBLE_PRECISION, mpi_status, ierr)
         call MPI_File_write(fhandl, inmax, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_write(fhandl, jnmax, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_write(fhandl, knmax, 1, MPI_INTEGER, &
           mpi_status, ierr)
       endif
!
       disp=sizeof(time)+3*sizeof(inmax)
       call MPI_FILE_SET_VIEW(fhandl, disp, &
         MPI_DOUBLE_PRECISION, filetype, 'native', &
         fileinfo, ierr)
!
       mysize=subsizes(1)*subsizes(2)*subsizes(3)
       do k=starts(3)+1,starts(3)+subsizes(3)
         do j=starts(2)+1,starts(2)+subsizes(2)
           do i=starts(1)+1,starts(1)+subsizes(1)
             counter=1+i-starts(1)-1 &
               +(j-starts(2)-1)*subsizes(1) &
               +(k-starts(3)-1)*subsizes(1)*subsizes(2)
             buf(counter)=b2(i-myid1*(in-5),j-myid2*(jn-5),k-myid3*(kn-5))
           enddo
         enddo
       enddo
       call MPI_File_write_all(fhandl, buf, &
         mysize, MPI_DOUBLE_PRECISION, mpi_status, ierr)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call MPI_File_close(fhandl, ierr)
!
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       fhandl=13
       call MPI_File_open(MPI_COMM_WORLD, &
         trim(b3file)//idout//'.dat', &
         IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), &
         fileinfo, fhandl, ierr)
!
       if(myid .eq. 0) then
         call MPI_File_write(fhandl, time, 1, &
           MPI_DOUBLE_PRECISION, mpi_status, ierr)
         call MPI_File_write(fhandl, inmax, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_write(fhandl, jnmax, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_write(fhandl, knmax, 1, MPI_INTEGER, &
           mpi_status, ierr)
       endif
!
       disp=sizeof(time)+3*sizeof(inmax)
       call MPI_FILE_SET_VIEW(fhandl, disp, &
         MPI_DOUBLE_PRECISION, filetype, 'native', &
         fileinfo, ierr)
!
       mysize=subsizes(1)*subsizes(2)*subsizes(3)
       do k=starts(3)+1,starts(3)+subsizes(3)
         do j=starts(2)+1,starts(2)+subsizes(2)
           do i=starts(1)+1,starts(1)+subsizes(1)
             counter=1+i-starts(1)-1 &
               +(j-starts(2)-1)*subsizes(1) &
               +(k-starts(3)-1)*subsizes(1)*subsizes(2)
             buf(counter)=b3(i-myid1*(in-5),j-myid2*(jn-5),k-myid3*(kn-5))
           enddo
         enddo
       enddo
       call MPI_File_write_all(fhandl, buf, &
         mysize, MPI_DOUBLE_PRECISION, mpi_status, ierr)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call MPI_File_close(fhandl, ierr)
!
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       fhandl=13
       call MPI_File_open(MPI_COMM_WORLD, &
         trim(dfile)//idout//'.dat', &
         IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), &
         fileinfo, fhandl, ierr)
!
       if(myid .eq. 0) then
         call MPI_File_write(fhandl, time, 1, &
           MPI_DOUBLE_PRECISION, mpi_status, ierr)
         call MPI_File_write(fhandl, inmax, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_write(fhandl, jnmax, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_write(fhandl, knmax, 1, MPI_INTEGER, &
           mpi_status, ierr)
       endif
!
       disp=sizeof(time)+3*sizeof(inmax)
       call MPI_FILE_SET_VIEW(fhandl, disp, &
         MPI_DOUBLE_PRECISION, filetype, 'native', &
         fileinfo, ierr)
!
       mysize=subsizes(1)*subsizes(2)*subsizes(3)
       do k=starts(3)+1,starts(3)+subsizes(3)
         do j=starts(2)+1,starts(2)+subsizes(2)
           do i=starts(1)+1,starts(1)+subsizes(1)
             counter=1+i-starts(1)-1 &
               +(j-starts(2)-1)*subsizes(1) &
               +(k-starts(3)-1)*subsizes(1)*subsizes(2)
             buf(counter)=d(i-myid1*(in-5),j-myid2*(jn-5),k-myid3*(kn-5))
           enddo
         enddo
       enddo
       call MPI_File_write_all(fhandl, buf, &
         mysize, MPI_DOUBLE_PRECISION, mpi_status, ierr)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call MPI_File_close(fhandl, ierr)
!
#ifndef ISOMHD
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       fhandl=13
       call MPI_File_open(MPI_COMM_WORLD, &
         trim(efile)//idout//'.dat', &
         IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), &
         fileinfo, fhandl, ierr)
!
       if(myid .eq. 0) then
         call MPI_File_write(fhandl, time, 1, &
           MPI_DOUBLE_PRECISION, mpi_status, ierr)
         call MPI_File_write(fhandl, inmax, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_write(fhandl, jnmax, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_write(fhandl, knmax, 1, MPI_INTEGER, &
           mpi_status, ierr)
       endif
!
       disp=sizeof(time)+3*sizeof(inmax)
       call MPI_FILE_SET_VIEW(fhandl, disp, &
         MPI_DOUBLE_PRECISION, filetype, 'native', &
         fileinfo, ierr)
!
       mysize=subsizes(1)*subsizes(2)*subsizes(3)
       do k=starts(3)+1,starts(3)+subsizes(3)
         do j=starts(2)+1,starts(2)+subsizes(2)
           do i=starts(1)+1,starts(1)+subsizes(1)
             counter=1+i-starts(1)-1 &
               +(j-starts(2)-1)*subsizes(1) &
               +(k-starts(3)-1)*subsizes(1)*subsizes(2)
             buf(counter)=eint(i-myid1*(in-5),j-myid2*(jn-5),k-myid3*(kn-5))
           enddo
         enddo
       enddo
       call MPI_File_write_all(fhandl, buf, &
         mysize, MPI_DOUBLE_PRECISION, mpi_status, ierr)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call MPI_File_close(fhandl, ierr)
#endif
!
#ifdef THCONDUCT
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       fhandl=13
       call MPI_File_open(MPI_COMM_WORLD, &
         trim(qfile)//idout//'.dat', &
         IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), &
         fileinfo, fhandl, ierr)
!
       if(myid .eq. 0) then
         call MPI_File_write(fhandl, time, 1, &
           MPI_DOUBLE_PRECISION, mpi_status, ierr)
         call MPI_File_write(fhandl, inmax, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_write(fhandl, jnmax, 1, MPI_INTEGER, &
           mpi_status, ierr)
         call MPI_File_write(fhandl, knmax, 1, MPI_INTEGER, &
           mpi_status, ierr)
       endif
!
       disp=sizeof(time)+3*sizeof(inmax)
       call MPI_FILE_SET_VIEW(fhandl, disp, &
         MPI_DOUBLE_PRECISION, filetype, 'native', &
         fileinfo, ierr)
!
       mysize=subsizes(1)*subsizes(2)*subsizes(3)
       do k=starts(3)+1,starts(3)+subsizes(3)
         do j=starts(2)+1,starts(2)+subsizes(2)
           do i=starts(1)+1,starts(1)+subsizes(1)
             counter=1+i-starts(1)-1 &
               +(j-starts(2)-1)*subsizes(1) &
               +(k-starts(3)-1)*subsizes(1)*subsizes(2)
             buf(counter)=q(i-myid1*(in-5),j-myid2*(jn-5),k-myid3*(kn-5))
           enddo
         enddo
       enddo
       call MPI_File_write_all(fhandl, buf, &
         mysize, MPI_DOUBLE_PRECISION, mpi_status, ierr)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call MPI_File_close(fhandl, ierr)
#endif

       call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
       deallocate(buf)

end subroutine writerst_mpi

subroutine writegrid
  use ModPar
  use ModGrid
  use ModControl
  implicit none

  integer i,j,k
  real subvect1(in),vect1(inmax)
  real subvect2(jn),vect2(jnmax)
  real subvect3(kn),vect3(knmax)
!------------------------------------------------------------------
!
       if(myid .eq. 0) then
         open(unit=13, file=trim(gridfile)//'.dat', &
           form='unformatted',access='stream')
         write(13) inmax,jnmax,knmax
       endif
!
!      write x1a
!
       do i=1,in
         subvect1(i)=x1a(i)
       enddo
       call combine_nproc1(subvect1,vect1)
       if(myid .eq. 0) then       
       write(13) (vect1(i),i=1,inmax)
       endif
!
!      write x1b
!
       do i=1,in
         subvect1(i)=x1b(i)
       enddo
       call combine_nproc1(subvect1,vect1)
       if(myid .eq. 0) then
       write(13) (vect1(i),i=1,inmax-1)
       endif
!
!      write x2a
!
       do j=1,jn
         subvect2(j)=x2a(j)
       enddo
       call combine_nproc2(subvect2,vect2)
       if(myid .eq. 0) then
       write(13) (vect2(j),j=1,jnmax)
       endif
!
!      write x2b
!
       do j=1,jn
         subvect2(j)=x2b(j)
       enddo
       call combine_nproc2(subvect2,vect2)
       if(myid .eq. 0) then
       write(13) (vect2(j),j=1,jnmax-1)
       endif
!
!      write x3a
!
       do k=1,kn
         subvect3(k)=x3a(k)
       enddo
       call combine_nproc3(subvect3,vect3)
       if(myid .eq. 0) then
       write(13) (vect3(k),k=1,knmax)
       endif
!
!      write x3b
!
       do k=1,kn
         subvect3(k)=x3b(k)
       enddo
       call combine_nproc3(subvect3,vect3)
       if(myid .eq. 0) then
       write(13) (vect3(k),k=1,knmax-1)
       endif
!
!      write dx1a
!
       do i=1,in
         subvect1(i)=dx1a(i)
       enddo
       call combine_nproc1(subvect1,vect1)
       if(myid .eq. 0) then
       write(13) (vect1(i),i=1,inmax-1)
       endif
!
!      write dx1b
!
       do i=1,in
         subvect1(i)=dx1b(i)
       enddo
       call combine_nproc1(subvect1,vect1)
       if(myid .eq. 0) then
       write(13) (vect1(i),i=1,inmax)
       endif
!
!      write dx2a
!
       do j=1,jn
         subvect2(j)=dx2a(j)
       enddo
       call combine_nproc2(subvect2,vect2)
       if(myid .eq. 0) then       
       write(13) (vect2(j),j=1,jnmax-1)
       endif
!
!      write dx2b
!
       do j=1,jn
         subvect2(j)=dx2b(j)
       enddo
       call combine_nproc2(subvect2,vect2)
       if(myid .eq. 0) then
       write(13) (vect2(j),j=1,jnmax)
       endif
!
!      write dx3a
!
       do k=1,kn
         subvect3(k)=dx3a(k)
       enddo
       call combine_nproc3(subvect3,vect3)
       if(myid .eq. 0) then
       write(13) (vect3(k),k=1,knmax-1)
       endif
!
!      write dx3b
!
       do k=1,kn
         subvect3(k)=dx3b(k)
       enddo
       call combine_nproc3(subvect3,vect3)
       if(myid .eq. 0) then
       write(13) (vect3(k),k=1,knmax)
       endif
!
!     write g2a
!
       do i=1,in
         subvect1(i)=g2a(i)
       enddo
       call combine_nproc1(subvect1,vect1)
       if(myid .eq. 0) then       
       write(13) (vect1(i),i=1,inmax)
       endif
!
!      write g2b
!
       do i=1,in
         subvect1(i)=g2b(i)
       enddo
       call combine_nproc1(subvect1,vect1)
       if(myid .eq. 0) then
       write(13) (vect1(i),i=1,inmax-1)
       endif
!
!      write g31a
!
       do i=1,in
         subvect1(i)=g31a(i)
       enddo
       call combine_nproc1(subvect1,vect1)
       if(myid .eq. 0) then       
       write(13) (vect1(i),i=1,inmax)
       endif
!
!      write g31b
!
       do i=1,in
         subvect1(i)=g31b(i)
       enddo
       call combine_nproc1(subvect1,vect1)
       if(myid .eq. 0) then
       write(13) (vect1(i),i=1,inmax-1)
       endif
!
!      write g32a
!
       do j=1,jn
         subvect2(j)=g32a(j)
       enddo
       call combine_nproc2(subvect2,vect2)
       if(myid .eq. 0) then
       write(13) (vect2(j),j=1,jnmax)
       endif
!
!      write g32b
!
       do j=1,jn
         subvect2(j)=g32b(j)
       enddo
       call combine_nproc2(subvect2,vect2)
       if(myid .eq. 0) then
       write(13) (vect2(j),j=1,jnmax-1)
       endif
!
!      write dg2bd1
!
       do i=1,in
         subvect1(i)=dg2bd1(i)
       enddo
       call combine_nproc1(subvect1,vect1)
       if(myid .eq. 0) then       
       write(13) (vect1(i),i=1,inmax)
       endif
!
!      write dg2ad1
!
       do i=1,in
         subvect1(i)=dg2ad1(i)
       enddo
       call combine_nproc1(subvect1,vect1)
       if(myid .eq. 0) then
       write(13) (vect1(i),i=1,inmax-1)
       endif
!
!      write dg31bd1
!
       do i=1,in
         subvect1(i)=dg31bd1(i)
       enddo
       call combine_nproc1(subvect1,vect1)
       if(myid .eq. 0) then       
       write(13) (vect1(i),i=1,inmax)
       endif
!
!      write dg31ad1
!
       do i=1,in
         subvect1(i)=dg31ad1(i)
       enddo
       call combine_nproc1(subvect1,vect1)
       if(myid .eq. 0) then
       write(13) (vect1(i),i=1,inmax-1)
       endif
!
!      write dg32bd2
!
       do j=1,jn
         subvect2(j)=dg32bd2(j)
       enddo
       call combine_nproc2(subvect2,vect2)
       if(myid .eq. 0) then
       write(13) (vect2(j),j=1,jnmax)
       endif
!
!      write dg32ad2
!
       do j=1,jn
         subvect2(j)=dg32ad2(j)
       enddo
       call combine_nproc2(subvect2,vect2)
       if(myid .eq. 0) then
       write(13) (vect2(j),j=1,jnmax-1)
       endif
!
!      write dvl1a
!
       do i=1,in
         subvect1(i)=dvl1a(i)
       enddo
       call combine_nproc1(subvect1,vect1)
       if(myid .eq. 0) then
       write(13) (vect1(i),i=1,inmax-1)
       endif
!
!      write dvl1b
!
       do i=1,in
         subvect1(i)=dvl1b(i)
       enddo
       call combine_nproc1(subvect1,vect1)
       if(myid .eq. 0) then       
       write(13) (vect1(i),i=1,inmax)
       endif
!
!      write dvl2a
!
       do j=1,jn
         subvect2(j)=dvl2a(j)
       enddo
       call combine_nproc2(subvect2,vect2)
       if(myid .eq. 0) then
       write(13) (vect2(j),j=1,jnmax-1)
       endif
!
!      write dvl2b
!
       do j=1,jn
         subvect2(j)=dvl2b(j)
       enddo
       call combine_nproc2(subvect2,vect2)
       if(myid .eq. 0) then
       write(13) (vect2(j),j=1,jnmax)
       endif
!
!      write dvl3a
!
       do k=1,kn
         subvect3(k)=dvl3a(k)
       enddo
       call combine_nproc3(subvect3,vect3)
       if(myid .eq. 0) then
       write(13) (vect3(k),k=1,knmax-1)
       endif
!
!      write dvl3b
!
       do k=1,kn
         subvect3(k)=dvl3b(k)
       enddo
       call combine_nproc3(subvect3,vect3)
       if(myid .eq. 0) then
       write(13) (vect3(k),k=1,knmax)
       endif
!
!      write dx1ai
!
       do i=1,in
         subvect1(i)=dx1ai(i)
       enddo
       call combine_nproc1(subvect1,vect1)
       if(myid .eq. 0) then
       write(13) (vect1(i),i=1,inmax-1)
       endif
!
!      write dx1bi
!
       do i=1,in
         subvect1(i)=dx1bi(i)
       enddo
       call combine_nproc1(subvect1,vect1)
       if(myid .eq. 0) then       
       write(13) (vect1(i),i=1,inmax)
       endif
!
!      write dx2ai
!
       do j=1,jn
         subvect2(j)=dx2ai(j)
       enddo
       call combine_nproc2(subvect2,vect2)
       if(myid .eq. 0) then
       write(13) (vect2(j),j=1,jnmax-1)
       endif
!
!      write dx2bi
!
       do j=1,jn
         subvect2(j)=dx2bi(j)
       enddo
       call combine_nproc2(subvect2,vect2)
       if(myid .eq. 0) then
       write(13) (vect2(j),j=1,jnmax)
       endif
!
!      write dx3ai
!
       do k=1,kn
         subvect3(k)=dx3ai(k)
       enddo
       call combine_nproc3(subvect3,vect3)
       if(myid .eq. 0) then
       write(13) (vect3(k),k=1,knmax-1)
       endif
!
!      write dx3bi
!
       do k=1,kn
         subvect3(k)=dx3bi(k)
       enddo
       call combine_nproc3(subvect3,vect3)
       if(myid .eq. 0) then
       write(13) (vect3(k),k=1,knmax)
       endif
!
!      write g2ai
!
       do i=1,in
         subvect1(i)=g2ai(i)
       enddo
       call combine_nproc1(subvect1,vect1)
       if(myid .eq. 0) then       
       write(13) (vect1(i),i=1,inmax)
       endif
!
!      write g2bi
!
       do i=1,in
         subvect1(i)=g2bi(i)
       enddo
       call combine_nproc1(subvect1,vect1)
       if(myid .eq. 0) then
       write(13) (vect1(i),i=1,inmax-1)
       endif
!
!      write g31ai
!
       do i=1,in
         subvect1(i)=g31ai(i)
       enddo
       call combine_nproc1(subvect1,vect1)
       if(myid .eq. 0) then       
       write(13) (vect1(i),i=1,inmax)
       endif
!
!      write g31bi
!
       do i=1,in
         subvect1(i)=g31bi(i)
       enddo
       call combine_nproc1(subvect1,vect1)
       if(myid .eq. 0) then
       write(13) (vect1(i),i=1,inmax-1)
       endif
!
!      write g32ai
!
       do j=1,jn
         subvect2(j)=g32ai(j)
       enddo
       call combine_nproc2(subvect2,vect2)
       if(myid .eq. 0) then
       write(13) (vect2(j),j=1,jnmax)
       endif
!
!      write g32bi
!
       do j=1,jn
         subvect2(j)=g32bi(j)
       enddo
       call combine_nproc2(subvect2,vect2)
       if(myid .eq. 0) then
       write(13) (vect2(j),j=1,jnmax-1)
       endif
!
!      write dvl1ai
!
       do i=1,in
         subvect1(i)=dvl1ai(i)
       enddo
       call combine_nproc1(subvect1,vect1)
       if(myid .eq. 0) then
       write(13) (vect1(i),i=1,inmax-1)
       endif
!
!      write dvl1bi
!
       do i=1,in
         subvect1(i)=dvl1bi(i)
       enddo
       call combine_nproc1(subvect1,vect1)
       if(myid .eq. 0) then       
       write(13) (vect1(i),i=1,inmax)
       endif
!
!      write dvl2ai
!
       do j=1,jn
         subvect2(j)=dvl2ai(j)
       enddo
       call combine_nproc2(subvect2,vect2)
       if(myid .eq. 0) then
       write(13) (vect2(j),j=1,jnmax-1)
       endif
!
!      write dvl2bi
!
       do j=1,jn
         subvect2(j)=dvl2bi(j)
       enddo
       call combine_nproc2(subvect2,vect2)
       if(myid .eq. 0) then
       write(13) (vect2(j),j=1,jnmax)
       endif
!
!      write dvl3ai
!
       do k=1,kn
         subvect3(k)=dvl3ai(k)
       enddo
       call combine_nproc3(subvect3,vect3)
       if(myid .eq. 0) then
       write(13) (vect3(k),k=1,knmax-1)
       endif
!
!      write dvl3bi
!
       do k=1,kn
         subvect3(k)=dvl3bi(k)
       enddo
       call combine_nproc3(subvect3,vect3)
       if(myid .eq. 0) then
       write(13) (vect3(k),k=1,knmax)
       endif
!
       if(myid .eq. 0) then
       close(13)
       endif
!
end subroutine writegrid

subroutine combine_nproc1(subvect,vect)
  use ModPar
  use ModGrid
  implicit none

  include "mpif.h"

  integer i
  real subvect(in),vect(inmax)
  real sendbuf1(in),recvbuf1(in)
  integer cputo,cpufrom,icpu
  integer istatus(MPI_STATUS_SIZE), ierr
!
       if(myid .eq. 0) then       
         do i=ism2,iep3
           vect(i)=subvect(i)
         enddo
!        if necessary receive data
         if(nproc1 .gt. 1) then
           do icpu=1,nproc1-1
             cpufrom=icpu
             call MPI_Recv(recvbuf1,in,MPI_DOUBLE_PRECISION, &
               cpufrom, icpu, MPI_COMM_WORLD, istatus, ierr)
             do i=1,in
               vect(i+icpu*(in-5))=recvbuf1(i)
             enddo
           enddo
         endif
       else
         if(myid2 .eq. 0 .and. myid3 .eq. 0) then
           cputo=0
           do i=1,in
             sendbuf1(i)=subvect(i)
           enddo
           call MPI_Send(sendbuf1,in,MPI_DOUBLE_PRECISION, &
             cputo, myid1, MPI_COMM_WORLD, ierr)
         endif
       endif
!
       call MPI_BARRIER(MPI_COMM_WORLD, ierr)
       call MPI_BCAST(vect,inmax,MPI_DOUBLE_PRECISION,0, &
         MPI_COMM_WORLD, ierr)
!
end subroutine combine_nproc1

subroutine combine_nproc2(subvect,vect)
  use ModPar
  use ModGrid
  implicit none

  include "mpif.h"

  integer j
  real subvect(jn),vect(jnmax)
  real sendbuf2(jn),recvbuf2(jn)
  integer cputo,cpufrom,jcpu
  integer istatus(MPI_STATUS_SIZE), ierr
!
       if(myid .eq. 0) then       
         do j=jsm2,jep3
           vect(j)=subvect(j)
         enddo
!        if necessary receive data
         if(nproc2 .gt. 1) then
           do jcpu=1,nproc2-1
             cpufrom=jcpu*nproc1
             call MPI_Recv(recvbuf2,jn,MPI_DOUBLE_PRECISION, &
               cpufrom, jcpu, MPI_COMM_WORLD, istatus, ierr)
             do j=1,jn
               vect(j+jcpu*(jn-5))=recvbuf2(j)
             enddo
           enddo
         endif
       else
         if(myid1 .eq. 0 .and. myid3 .eq. 0) then
           cputo=0
           do j=1,jn
             sendbuf2(j)=subvect(j)
           enddo
           call MPI_Send(sendbuf2,jn,MPI_DOUBLE_PRECISION, &
             cputo, myid2, MPI_COMM_WORLD, ierr)
         endif
       endif
!
       call MPI_BARRIER(MPI_COMM_WORLD, ierr)
       call MPI_BCAST(vect,jnmax,MPI_DOUBLE_PRECISION,0, &
         MPI_COMM_WORLD, ierr)
!
end subroutine combine_nproc2

subroutine combine_nproc3(subvect,vect)
  use ModPar
  use ModGrid
  implicit none

  include "mpif.h"

  integer k
  real subvect(kn),vect(knmax)
  real sendbuf3(kn),recvbuf3(kn)
  integer cputo,cpufrom,kcpu
  integer istatus(MPI_STATUS_SIZE), ierr
!
       if(myid .eq. 0) then       
         do k=ksm2,kep3
           vect(k)=subvect(k)
         enddo
!        if necessary receive data
         if(nproc3 .gt. 1) then
           do kcpu=1,nproc3-1
             cpufrom=kcpu*(nproc1*nproc2)
             call MPI_Recv(recvbuf3,kn,MPI_DOUBLE_PRECISION, &
               cpufrom, kcpu, MPI_COMM_WORLD, istatus, ierr)
             do k=1,kn
               vect(k+kcpu*(kn-5))=recvbuf3(k)
             enddo
           enddo
         endif
       else
         if(myid1 .eq. 0 .and. myid2 .eq. 0) then
           cputo=0
           do k=1,kn
             sendbuf3(k)=subvect(k)
           enddo
           call MPI_Send(sendbuf3,kn,MPI_DOUBLE_PRECISION, &
             cputo, myid3, MPI_COMM_WORLD, ierr)
         endif
       endif
!
       call MPI_BARRIER(MPI_COMM_WORLD, ierr)
       call MPI_BCAST(vect,knmax,MPI_DOUBLE_PRECISION,0, &
         MPI_COMM_WORLD, ierr)
!
end subroutine combine_nproc3

subroutine writedata_cpu
  use ModPar
  use ModGrid
  use ModSundry
  use ModField
  use ModControl
  implicit none

  character(len=4) :: idout, idcpu
  integer i,j,k
!------------------------------------------------------
       write(idcpu,'(i4.4)') myid
       write(idout,'(i4.4)') ifile

       open(unit=13, file=trim(rstfile)//idout//'.cpu'//idcpu, &
         form='unformatted', access='stream')
       write(13) in,jn,kn
       write(13) (((b1(i,j,k),i=ism2,iep2),j=jsm2,jep2),k=ksm2,kep2)
       write(13) (((b2(i,j,k),i=ism2,iep2),j=jsm2,jep2),k=ksm2,kep2)
       write(13) (((b3(i,j,k),i=ism2,iep2),j=jsm2,jep2),k=ksm2,kep2)
       write(13) (((v1(i,j,k),i=ism2,iep2),j=jsm2,jep2),k=ksm2,kep2)
       write(13) (((v2(i,j,k),i=ism2,iep2),j=jsm2,jep2),k=ksm2,kep2)
       write(13) (((v3(i,j,k),i=ism2,iep2),j=jsm2,jep2),k=ksm2,kep2)
       write(13) (((d(i,j,k),i=ism2,iep2),j=jsm2,jep2),k=ksm2,kep2)
#ifndef ISOMHD
       write(13) (((eint(i,j,k),i=ism2,iep2),j=jsm2,jep2),k=ksm2,kep2)
#endif
       write(13) time
       close(13)

end subroutine writedata_cpu

subroutine writegrid_cpu
  use ModPar
  use ModGrid
  use ModControl
  implicit none

  integer i,j,k
  character(len=4) idcpu
!---------------------------------------------------------
       write(idcpu,'(i4.4)') myid
       open(unit=13, file=trim(gridfile)//'.cpu'//idcpu, &
         form='unformatted',access='stream')
       write(13) in,jn,kn
       write(13) (x1a(i),i=ism2,iep3)
       write(13) (x1b(i),i=ism2,iep2)
       write(13) (x2a(j),j=jsm2,jep3)
       write(13) (x2b(j),j=jsm2,jep2)
       write(13) (x3a(k),k=ksm2,kep3)
       write(13) (x3b(k),k=ksm2,kep2)

       write(13) (dx1a(i),i=ism2,iep2)
       write(13) (dx1b(i),i=ism2,iep3)
       write(13) (dx2a(j),j=jsm2,jep2)
       write(13) (dx2b(j),j=jsm2,jep3)
       write(13) (dx3a(k),k=ksm2,kep2)
       write(13) (dx3b(k),k=ksm2,kep3)

       write(13) (g2a(i),i=ism2,iep3)
       write(13) (g2b(i),i=ism2,iep2)
       write(13) (g31a(i),i=ism2,iep3)
       write(13) (g31b(i),i=ism2,iep2)
       write(13) (g32a(j),j=jsm2,jep3)
       write(13) (g32b(j),j=jsm2,jep2)

       write(13) (dg2bd1(i),i=ism2,iep3)
       write(13) (dg2ad1(i),i=ism2,iep2)
       write(13) (dg31bd1(i),i=ism2,iep3)
       write(13) (dg31ad1(i),i=ism2,iep2)
       write(13) (dg32bd2(j),j=jsm2,jep3)
       write(13) (dg32ad2(j),j=jsm2,jep2)

       write(13) (dvl1a(i),i=ism2,iep2)
       write(13) (dvl1b(i),i=ism2,iep3)
       write(13) (dvl2a(j),j=jsm2,jep2)
       write(13) (dvl2b(j),j=jsm2,jep3)
       write(13) (dvl3a(k),k=ksm2,kep2)
       write(13) (dvl3b(k),k=ksm2,kep3)

       write(13) (dx1ai(i),i=ism2,iep2)
       write(13) (dx1bi(i),i=ism2,iep3)
       write(13) (dx2ai(j),j=jsm2,jep2)
       write(13) (dx2bi(j),j=jsm2,jep3)
       write(13) (dx3ai(k),k=ksm2,kep2)
       write(13) (dx3bi(k),k=ksm2,kep3)

       write(13) (g2ai(i),i=ism2,iep3)
       write(13) (g2bi(i),i=ism2,iep2)
       write(13) (g31ai(i),i=ism2,iep3)
       write(13) (g31bi(i),i=ism2,iep2)
       write(13) (g32ai(j),j=jsm2,jep3)
       write(13) (g32bi(j),j=jsm2,jep2)

       write(13) (dvl1ai(i),i=ism2,iep2)
       write(13) (dvl1bi(i),i=ism2,iep3)
       write(13) (dvl2ai(j),j=jsm2,jep2)
       write(13) (dvl2bi(j),j=jsm2,jep3)
       write(13) (dvl3ai(k),k=ksm2,kep2)
       write(13) (dvl3bi(k),k=ksm2,kep3)
       close(13)

end subroutine writegrid_cpu

subroutine writeqty_mpi(qty,filename)
  use ModPar
  use ModGrid
  use ModSundry
  implicit none
  include "mpif.h"

  real, intent(in) :: qty(in,jn,kn)
  character(len=80), intent(in) :: filename

  integer i,j,k,counter
  integer sizes(3), subsizes(3), starts(3)
  integer filetype,fileinfo,fhandl,mysize,ierr
  integer(kind=MPI_OFFSET_KIND) disp
  integer mpi_status(MPI_STATUS_SIZE)
  real, allocatable :: buf(:)
!-----------------------------------------------------------
       allocate(buf(in*jn*kn))

       sizes(1)=inmax-1
       sizes(2)=jnmax-1
       sizes(3)=knmax-1
!
       subsizes(1)=in-5
       starts(1)=is-1+myid1*(in-5)
       if(myid1 .eq. 0) then
         subsizes(1)=subsizes(1)+2
         starts(1)=0
       endif
       if(myid1 .eq. nproc1-1) then
         subsizes(1)=subsizes(1)+2
       endif
!
       subsizes(2)=jn-5
       starts(2)=js-1+myid2*(jn-5)
       if(myid2 .eq. 0) then
         subsizes(2)=subsizes(2)+2
         starts(2)=0
       endif
       if(myid2 .eq. nproc2-1) then
         subsizes(2)=subsizes(2)+2
       endif
!
       subsizes(3)=kn-5
       starts(3)=ks-1+myid3*(kn-5)
       if(myid3 .eq. 0) then
         subsizes(3)=subsizes(3)+2
         starts(3)=0
       endif
       if(myid3 .eq. nproc3-1) then
         subsizes(3)=subsizes(3)+2
       endif
!
       call MPI_Type_create_subarray(3, sizes, subsizes, starts, &
            MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, filetype, ierr)
       call MPI_Type_commit(filetype, ierr)
       call MPI_Info_create(fileinfo, ierr)
!
       fhandl=13
       call MPI_File_open(MPI_COMM_WORLD, &
            trim(filename), &
            IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), &
            fileinfo, fhandl, ierr)
!
       if(myid .eq. 0) then
         call MPI_File_write(fhandl, time, 1, &
              MPI_DOUBLE_PRECISION, mpi_status, ierr)
         call MPI_File_write(fhandl, inmax, 1, MPI_INTEGER, &
              mpi_status, ierr)
         call MPI_File_write(fhandl, jnmax, 1, MPI_INTEGER, &
              mpi_status, ierr)
         call MPI_File_write(fhandl, knmax, 1, MPI_INTEGER, &
              mpi_status, ierr)
       endif
!
       disp=sizeof(time)+3*sizeof(inmax)
       call MPI_FILE_SET_VIEW(fhandl, disp, &
            MPI_DOUBLE_PRECISION, filetype, 'native', &
            fileinfo, ierr)
!
       mysize=subsizes(1)*subsizes(2)*subsizes(3)
       do k=starts(3)+1,starts(3)+subsizes(3)
         do j=starts(2)+1,starts(2)+subsizes(2)
           do i=starts(1)+1,starts(1)+subsizes(1)
             counter=1+i-starts(1)-1 &
               +(j-starts(2)-1)*subsizes(1) &
               +(k-starts(3)-1)*subsizes(1)*subsizes(2)
             buf(counter)=qty(i-myid1*(in-5),j-myid2*(jn-5),k-myid3*(kn-5))
           enddo
         enddo
       enddo
       call MPI_File_write_all(fhandl, buf, &
            mysize, MPI_DOUBLE_PRECISION, mpi_status, ierr)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call MPI_File_close(fhandl, ierr)
!
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
       deallocate(buf)

end subroutine writeqty_mpi

subroutine readqty_mpi(qty,filename)
  use ModPar
  use ModGrid
  implicit none
  include "mpif.h"

  character(len=80), intent(in) :: filename
  real, intent(out) :: qty(in,jn,kn)

  integer i,j,k
  integer sizes(3), subsizes(3), starts(3)
  integer filetype,fileinfo,fhandl,mysize,ierr
  integer(kind=MPI_OFFSET_KIND) disp
  integer mpi_status(MPI_STATUS_SIZE)
  integer inmax_in, jnmax_in, knmax_in
  real time_in
  real, allocatable :: buf(:,:,:)
!----------------------------------------------
       allocate(buf(in-1,jn-1,kn-1))
!
       do k=1,kn
         do j=1,jn
           do i=1,in
             qty(i,j,k)=0.D0
           enddo
         enddo
       enddo
!
       sizes(1)=inmax-1
       sizes(2)=jnmax-1
       sizes(3)=knmax-1
!
       subsizes(1)=in-1
       starts(1)=myid1*(in-5)
!
       subsizes(2)=jn-1
       starts(2)=myid2*(jn-5)
!
       subsizes(3)=kn-1
       starts(3)=myid3*(kn-5)
!
       call MPI_Type_create_subarray(3, sizes, subsizes, starts, &
            MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, filetype, ierr)
       call MPI_Type_commit(filetype, ierr)
       call MPI_Info_create(fileinfo, ierr)
!
       fhandl=13
       call MPI_File_open(MPI_COMM_WORLD, &
            trim(filename), &
            MPI_MODE_RDONLY, &
            fileinfo, fhandl, ierr)
!
       if(myid .eq. 0) then
         call MPI_File_read(fhandl, time_in, 1, &
              MPI_DOUBLE_PRECISION, mpi_status, ierr)
         call MPI_File_read(fhandl, inmax_in, 1, MPI_INTEGER, &
              mpi_status, ierr)
         call MPI_File_read(fhandl, jnmax_in, 1, MPI_INTEGER, &
              mpi_status, ierr)
         call MPI_File_read(fhandl, knmax_in, 1, MPI_INTEGER, &
              mpi_status, ierr)
       endif
       call MPI_BCAST(time_in,1,MPI_DOUBLE_PRECISION,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(inmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(jnmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(knmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       if(inmax_in .ne. inmax .or.  &
          jnmax_in .ne. jnmax .or.  &
          knmax_in .ne. knmax) then
          write(6,*) 'inmax,jnmax,knmax in file', &
            inmax_in,jnmax_in,knmax_in, &
            ' do not match inmax, jnmax, knmax in par.h and timetable'
          call MPI_ABORT(MPI_COMM_WORLD, 1,ierr)
       endif
!
       disp=sizeof(time_in)+3*sizeof(inmax)
       call MPI_FILE_SET_VIEW(fhandl, disp, &
            MPI_DOUBLE_PRECISION, filetype, 'native', &
            fileinfo, ierr)
!
       mysize=subsizes(1)*subsizes(2)*subsizes(3)
       call MPI_File_read_all(fhandl, buf, &
            mysize, MPI_DOUBLE_PRECISION, mpi_status, ierr)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call MPI_File_close(fhandl, ierr)
       do k=1,kn-1
         do j=1,jn-1
           do i=1,in-1
             qty(i,j,k)=buf(i,j,k)
           enddo
         enddo
       enddo
!
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       deallocate(buf)

end subroutine readqty_mpi

subroutine writeemf1_mpi(qty,filename)
  use ModPar
  use ModGrid
  use ModSundry
  implicit none
  include "mpif.h"

  real, intent(in) :: qty(in,jn,kn)
  character(len=80), intent(in) :: filename

  integer i,j,k,counter
  integer sizes(3), subsizes(3), starts(3)
  integer filetype,fileinfo,fhandl,mysize,ierr
  integer(kind=MPI_OFFSET_KIND) disp
  integer mpi_status(MPI_STATUS_SIZE)
  real, allocatable :: buf(:)
!-----------------------------------------------------------
       allocate(buf(in*jn*kn))

       sizes(1)=inmax-1
       sizes(2)=jnmax
       sizes(3)=knmax
!
       subsizes(1)=in-5
       starts(1)=is-1+myid1*(in-5)
       if(myid1 .eq. 0) then
         subsizes(1)=subsizes(1)+2
         starts(1)=0
       endif
       if(myid1 .eq. nproc1-1) then
         subsizes(1)=subsizes(1)+2
       endif
!
       subsizes(2)=jn-5
       starts(2)=js-1+myid2*(jn-5)
       if(myid2 .eq. 0) then
         subsizes(2)=subsizes(2)+2
         starts(2)=0
       endif
       if(myid2 .eq. nproc2-1) then
         subsizes(2)=subsizes(2)+3
       endif
!
       subsizes(3)=kn-5
       starts(3)=ks-1+myid3*(kn-5)
       if(myid3 .eq. 0) then
         subsizes(3)=subsizes(3)+2
         starts(3)=0
       endif
       if(myid3 .eq. nproc3-1) then
         subsizes(3)=subsizes(3)+3
       endif
!
       call MPI_Type_create_subarray(3, sizes, subsizes, starts, &
            MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, filetype, ierr)
       call MPI_Type_commit(filetype, ierr)
       call MPI_Info_create(fileinfo, ierr)
!
       fhandl=13
       call MPI_File_open(MPI_COMM_WORLD, &
            trim(filename), &
            IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), &
            fileinfo, fhandl, ierr)
!
       if(myid .eq. 0) then
         call MPI_File_write(fhandl, time, 1, &
              MPI_DOUBLE_PRECISION, mpi_status, ierr)
         call MPI_File_write(fhandl, inmax, 1, MPI_INTEGER, &
              mpi_status, ierr)
         call MPI_File_write(fhandl, jnmax, 1, MPI_INTEGER, &
              mpi_status, ierr)
         call MPI_File_write(fhandl, knmax, 1, MPI_INTEGER, &
              mpi_status, ierr)
       endif
!
       disp=sizeof(time)+3*sizeof(inmax)
       call MPI_FILE_SET_VIEW(fhandl, disp, &
            MPI_DOUBLE_PRECISION, filetype, 'native', &
            fileinfo, ierr)
!
       mysize=subsizes(1)*subsizes(2)*subsizes(3)
       do k=starts(3)+1,starts(3)+subsizes(3)
         do j=starts(2)+1,starts(2)+subsizes(2)
           do i=starts(1)+1,starts(1)+subsizes(1)
             counter=1+i-starts(1)-1 &
               +(j-starts(2)-1)*subsizes(1) &
               +(k-starts(3)-1)*subsizes(1)*subsizes(2)
             buf(counter)=qty(i-myid1*(in-5),j-myid2*(jn-5),k-myid3*(kn-5))
           enddo
         enddo
       enddo
       call MPI_File_write_all(fhandl, buf, &
            mysize, MPI_DOUBLE_PRECISION, mpi_status, ierr)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call MPI_File_close(fhandl, ierr)
!
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
       deallocate(buf)

end subroutine writeemf1_mpi

subroutine writeemf2_mpi(qty,filename)
  use ModPar
  use ModGrid
  use ModSundry
  implicit none
  include "mpif.h"

  real, intent(in) :: qty(in,jn,kn)
  character(len=80), intent(in) :: filename

  integer i,j,k,counter
  integer sizes(3), subsizes(3), starts(3)
  integer filetype,fileinfo,fhandl,mysize,ierr
  integer(kind=MPI_OFFSET_KIND) disp
  integer mpi_status(MPI_STATUS_SIZE)
  real, allocatable :: buf(:)
!-----------------------------------------------------------
       allocate(buf(in*jn*kn))

       sizes(1)=inmax
       sizes(2)=jnmax-1
       sizes(3)=knmax
!
       subsizes(1)=in-5
       starts(1)=is-1+myid1*(in-5)
       if(myid1 .eq. 0) then
         subsizes(1)=subsizes(1)+2
         starts(1)=0
       endif
       if(myid1 .eq. nproc1-1) then
         subsizes(1)=subsizes(1)+3
       endif
!
       subsizes(2)=jn-5
       starts(2)=js-1+myid2*(jn-5)
       if(myid2 .eq. 0) then
         subsizes(2)=subsizes(2)+2
         starts(2)=0
       endif
       if(myid2 .eq. nproc2-1) then
         subsizes(2)=subsizes(2)+2
       endif
!
       subsizes(3)=kn-5
       starts(3)=ks-1+myid3*(kn-5)
       if(myid3 .eq. 0) then
         subsizes(3)=subsizes(3)+2
         starts(3)=0
       endif
       if(myid3 .eq. nproc3-1) then
         subsizes(3)=subsizes(3)+3
       endif
!
       call MPI_Type_create_subarray(3, sizes, subsizes, starts, &
            MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, filetype, ierr)
       call MPI_Type_commit(filetype, ierr)
       call MPI_Info_create(fileinfo, ierr)
!
       fhandl=13
       call MPI_File_open(MPI_COMM_WORLD, &
            trim(filename), &
            IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), &
            fileinfo, fhandl, ierr)
!
       if(myid .eq. 0) then
         call MPI_File_write(fhandl, time, 1, &
              MPI_DOUBLE_PRECISION, mpi_status, ierr)
         call MPI_File_write(fhandl, inmax, 1, MPI_INTEGER, &
              mpi_status, ierr)
         call MPI_File_write(fhandl, jnmax, 1, MPI_INTEGER, &
              mpi_status, ierr)
         call MPI_File_write(fhandl, knmax, 1, MPI_INTEGER, &
              mpi_status, ierr)
       endif
!
       disp=sizeof(time)+3*sizeof(inmax)
       call MPI_FILE_SET_VIEW(fhandl, disp, &
            MPI_DOUBLE_PRECISION, filetype, 'native', &
            fileinfo, ierr)
!
       mysize=subsizes(1)*subsizes(2)*subsizes(3)
       do k=starts(3)+1,starts(3)+subsizes(3)
         do j=starts(2)+1,starts(2)+subsizes(2)
           do i=starts(1)+1,starts(1)+subsizes(1)
             counter=1+i-starts(1)-1 &
               +(j-starts(2)-1)*subsizes(1) &
               +(k-starts(3)-1)*subsizes(1)*subsizes(2)
             buf(counter)=qty(i-myid1*(in-5),j-myid2*(jn-5),k-myid3*(kn-5))
           enddo
         enddo
       enddo
       call MPI_File_write_all(fhandl, buf, &
            mysize, MPI_DOUBLE_PRECISION, mpi_status, ierr)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call MPI_File_close(fhandl, ierr)
!
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
       deallocate(buf)

end subroutine writeemf2_mpi

subroutine writeemf3_mpi(qty,filename)
  use ModPar
  use ModGrid
  use ModSundry
  implicit none
  include "mpif.h"

  real, intent(in) :: qty(in,jn,kn)
  character(len=80), intent(in) :: filename

  integer i,j,k,counter
  integer sizes(3), subsizes(3), starts(3)
  integer filetype,fileinfo,fhandl,mysize,ierr
  integer(kind=MPI_OFFSET_KIND) disp
  integer mpi_status(MPI_STATUS_SIZE)
  real, allocatable :: buf(:)
!-----------------------------------------------------------
       allocate(buf(in*jn*kn))

       sizes(1)=inmax
       sizes(2)=jnmax
       sizes(3)=knmax-1
!
       subsizes(1)=in-5
       starts(1)=is-1+myid1*(in-5)
       if(myid1 .eq. 0) then
         subsizes(1)=subsizes(1)+2
         starts(1)=0
       endif
       if(myid1 .eq. nproc1-1) then
         subsizes(1)=subsizes(1)+3
       endif
!
       subsizes(2)=jn-5
       starts(2)=js-1+myid2*(jn-5)
       if(myid2 .eq. 0) then
         subsizes(2)=subsizes(2)+2
         starts(2)=0
       endif
       if(myid2 .eq. nproc2-1) then
         subsizes(2)=subsizes(2)+3
       endif
!
       subsizes(3)=kn-5
       starts(3)=ks-1+myid3*(kn-5)
       if(myid3 .eq. 0) then
         subsizes(3)=subsizes(3)+2
         starts(3)=0
       endif
       if(myid3 .eq. nproc3-1) then
         subsizes(3)=subsizes(3)+2
       endif
!
       call MPI_Type_create_subarray(3, sizes, subsizes, starts, &
            MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, filetype, ierr)
       call MPI_Type_commit(filetype, ierr)
       call MPI_Info_create(fileinfo, ierr)
!
       fhandl=13
       call MPI_File_open(MPI_COMM_WORLD, &
            trim(filename), &
            IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), &
            fileinfo, fhandl, ierr)
!
       if(myid .eq. 0) then
         call MPI_File_write(fhandl, time, 1, &
              MPI_DOUBLE_PRECISION, mpi_status, ierr)
         call MPI_File_write(fhandl, inmax, 1, MPI_INTEGER, &
              mpi_status, ierr)
         call MPI_File_write(fhandl, jnmax, 1, MPI_INTEGER, &
              mpi_status, ierr)
         call MPI_File_write(fhandl, knmax, 1, MPI_INTEGER, &
              mpi_status, ierr)
       endif
!
       disp=sizeof(time)+3*sizeof(inmax)
       call MPI_FILE_SET_VIEW(fhandl, disp, &
            MPI_DOUBLE_PRECISION, filetype, 'native', &
            fileinfo, ierr)
!
       mysize=subsizes(1)*subsizes(2)*subsizes(3)
       do k=starts(3)+1,starts(3)+subsizes(3)
         do j=starts(2)+1,starts(2)+subsizes(2)
           do i=starts(1)+1,starts(1)+subsizes(1)
             counter=1+i-starts(1)-1 &
               +(j-starts(2)-1)*subsizes(1) &
               +(k-starts(3)-1)*subsizes(1)*subsizes(2)
             buf(counter)=qty(i-myid1*(in-5),j-myid2*(jn-5),k-myid3*(kn-5))
           enddo
         enddo
       enddo
       call MPI_File_write_all(fhandl, buf, &
            mysize, MPI_DOUBLE_PRECISION, mpi_status, ierr)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call MPI_File_close(fhandl, ierr)
!
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
!
       deallocate(buf)

end subroutine writeemf3_mpi

subroutine read12qty_mpi(ntsub,qty,filename)
  use ModPar
  use ModGrid
  implicit none
  include "mpif.h"

  integer, intent(in) :: ntsub
  character(len=80), intent(in) :: filename
  real, intent(out) :: qty(in,jn,ntsub)

  integer i,j,k
  integer sizes(3), subsizes(3), starts(3)
  integer filetype,fileinfo,fhandl,mysize,ierr
  integer(kind=MPI_OFFSET_KIND) disp
  integer mpi_status(MPI_STATUS_SIZE)
  integer inmax_in, jnmax_in, ntsub_in
  real, allocatable :: buf(:,:,:)
!----------------------------------------------
       allocate(buf(in-1,jn-1,ntsub))
!
       do k=1,ntsub
         do j=1,jn
           do i=1,in
             qty(i,j,k)=0.D0
           enddo
         enddo
       enddo
!
       sizes(1)=inmax-1
       sizes(2)=jnmax-1
       sizes(3)=ntsub
!
       subsizes(1)=in-1
       starts(1)=myid1*(in-5)
!
       subsizes(2)=jn-1
       starts(2)=myid2*(jn-5)
!
       subsizes(3)=ntsub
       starts(3)=0
!
       call MPI_Type_create_subarray(3, sizes, subsizes, starts, &
            MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, filetype, ierr)
       call MPI_Type_commit(filetype, ierr)
       call MPI_Info_create(fileinfo, ierr)
!
       fhandl=13
       call MPI_File_open(MPI_COMM_WORLD, &
            trim(filename), &
            MPI_MODE_RDONLY, &
            fileinfo, fhandl, ierr)
!
       if(myid .eq. 0) then
         call MPI_File_read(fhandl, inmax_in, 1, MPI_INTEGER, &
              mpi_status, ierr)
         call MPI_File_read(fhandl, jnmax_in, 1, MPI_INTEGER, &
              mpi_status, ierr)
         call MPI_File_read(fhandl, ntsub_in, 1, MPI_INTEGER, &
              mpi_status, ierr)
       endif
       call MPI_BCAST(inmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(jnmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(ntsub_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       if(inmax_in .ne. inmax .or.  &
          jnmax_in .ne. jnmax .or.  &
          ntsub_in .ne. ntsub) then
          write(6,*) 'inmax,jnmax,ntsub in file', &
            inmax_in,jnmax_in,ntsub_in, &
            ' do not match inmax, jnmax, ntsub in par.h and timetable'
          call MPI_ABORT(MPI_COMM_WORLD, 1,ierr)
       endif
!
       disp=2*sizeof(inmax)+sizeof(ntsub)
       call MPI_FILE_SET_VIEW(fhandl, disp, &
            MPI_DOUBLE_PRECISION, filetype, 'native', &
            fileinfo, ierr)
!
       mysize=subsizes(1)*subsizes(2)*subsizes(3)
       call MPI_File_read_all(fhandl, buf, &
            mysize, MPI_DOUBLE_PRECISION, mpi_status, ierr)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call MPI_File_close(fhandl, ierr)
       do k=1,ntsub
         do j=1,jn-1
           do i=1,in-1
             qty(i,j,k)=buf(i,j,k)
           enddo
         enddo
       enddo
!
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       deallocate(buf)

end subroutine read12qty_mpi

subroutine read12v1_mpi(ntsub,qty,filename)
  use ModPar
  use ModGrid
  implicit none
  include "mpif.h"

  integer, intent(in) :: ntsub
  character(len=80), intent(in) :: filename
  real, intent(out) :: qty(in,jn,ntsub)

  integer i,j,k
  integer sizes(3), subsizes(3), starts(3)
  integer filetype,fileinfo,fhandl,mysize,ierr
  integer(kind=MPI_OFFSET_KIND) disp
  integer mpi_status(MPI_STATUS_SIZE)
  integer inmax_in, jnmax_in, ntsub_in
  real, allocatable :: buf(:,:,:)
!----------------------------------------------
       allocate(buf(in,jn-1,ntsub))
!
       do k=1,ntsub
         do j=1,jn
           do i=1,in
             qty(i,j,k)=0.D0
           enddo
         enddo
       enddo
!
       sizes(1)=inmax
       sizes(2)=jnmax-1
       sizes(3)=ntsub
!
       subsizes(1)=in
       starts(1)=myid1*(in-5)
!
       subsizes(2)=jn-1
       starts(2)=myid2*(jn-5)
!
       subsizes(3)=ntsub
       starts(3)=0
!
       call MPI_Type_create_subarray(3, sizes, subsizes, starts, &
            MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, filetype, ierr)
       call MPI_Type_commit(filetype, ierr)
       call MPI_Info_create(fileinfo, ierr)
!
       fhandl=13
       call MPI_File_open(MPI_COMM_WORLD, &
            trim(filename), &
            MPI_MODE_RDONLY, &
            fileinfo, fhandl, ierr)
!
       if(myid .eq. 0) then
         call MPI_File_read(fhandl, inmax_in, 1, MPI_INTEGER, &
              mpi_status, ierr)
         call MPI_File_read(fhandl, jnmax_in, 1, MPI_INTEGER, &
              mpi_status, ierr)
         call MPI_File_read(fhandl, ntsub_in, 1, MPI_INTEGER, &
              mpi_status, ierr)
       endif
       call MPI_BCAST(inmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(jnmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(ntsub_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       if(inmax_in .ne. inmax .or.  &
          jnmax_in .ne. jnmax .or.  &
          ntsub_in .ne. ntsub) then
          write(6,*) 'inmax,jnmax,ntsub in file', &
            inmax_in,jnmax_in,ntsub_in, &
            ' do not match inmax, jnmax, ntsub in par.h and timetable'
          call MPI_ABORT(MPI_COMM_WORLD, 1,ierr)
       endif
!
       disp=2*sizeof(inmax)+sizeof(ntsub)
       call MPI_FILE_SET_VIEW(fhandl, disp, &
            MPI_DOUBLE_PRECISION, filetype, 'native', &
            fileinfo, ierr)
!
       mysize=subsizes(1)*subsizes(2)*subsizes(3)
       call MPI_File_read_all(fhandl, buf, &
            mysize, MPI_DOUBLE_PRECISION, mpi_status, ierr)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call MPI_File_close(fhandl, ierr)
       do k=1,ntsub
         do j=1,jn-1
           do i=1,in
             qty(i,j,k)=buf(i,j,k)
           enddo
         enddo
       enddo
!
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       deallocate(buf)

end subroutine read12v1_mpi

subroutine read12v2_mpi(ntsub,qty,filename)
  use ModPar
  use ModGrid
  implicit none
  include "mpif.h"

  integer, intent(in) :: ntsub
  character(len=80), intent(in) :: filename
  real, intent(out) :: qty(in,jn,ntsub)

  integer i,j,k
  integer sizes(3), subsizes(3), starts(3)
  integer filetype,fileinfo,fhandl,mysize,ierr
  integer(kind=MPI_OFFSET_KIND) disp
  integer mpi_status(MPI_STATUS_SIZE)
  integer inmax_in, jnmax_in, ntsub_in
  real, allocatable :: buf(:,:,:)
!----------------------------------------------
       allocate(buf(in-1,jn,ntsub))
!
       do k=1,ntsub
         do j=1,jn
           do i=1,in
             qty(i,j,k)=0.D0
           enddo
         enddo
       enddo
!
       sizes(1)=inmax-1
       sizes(2)=jnmax
       sizes(3)=ntsub
!
       subsizes(1)=in-1
       starts(1)=myid1*(in-5)
!
       subsizes(2)=jn
       starts(2)=myid2*(jn-5)
!
       subsizes(3)=ntsub
       starts(3)=0
!
       call MPI_Type_create_subarray(3, sizes, subsizes, starts, &
            MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, filetype, ierr)
       call MPI_Type_commit(filetype, ierr)
       call MPI_Info_create(fileinfo, ierr)
!
       fhandl=13
       call MPI_File_open(MPI_COMM_WORLD, &
            trim(filename), &
            MPI_MODE_RDONLY, &
            fileinfo, fhandl, ierr)
!
       if(myid .eq. 0) then
         call MPI_File_read(fhandl, inmax_in, 1, MPI_INTEGER, &
              mpi_status, ierr)
         call MPI_File_read(fhandl, jnmax_in, 1, MPI_INTEGER, &
              mpi_status, ierr)
         call MPI_File_read(fhandl, ntsub_in, 1, MPI_INTEGER, &
              mpi_status, ierr)
       endif
       call MPI_BCAST(inmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(jnmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(ntsub_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       if(inmax_in .ne. inmax .or.  &
          jnmax_in .ne. jnmax .or.  &
          ntsub_in .ne. ntsub) then
          write(6,*) 'inmax,jnmax,ntsub in file', &
            inmax_in,jnmax_in,ntsub_in, &
            ' do not match inmax, jnmax, ntsub in par.h and timetable'
          call MPI_ABORT(MPI_COMM_WORLD, 1,ierr)
       endif
!
       disp=2*sizeof(inmax)+sizeof(ntsub)
       call MPI_FILE_SET_VIEW(fhandl, disp, &
            MPI_DOUBLE_PRECISION, filetype, 'native', &
            fileinfo, ierr)
!
       mysize=subsizes(1)*subsizes(2)*subsizes(3)
       call MPI_File_read_all(fhandl, buf, &
            mysize, MPI_DOUBLE_PRECISION, mpi_status, ierr)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call MPI_File_close(fhandl, ierr)
       do k=1,ntsub
         do j=1,jn
           do i=1,in-1
             qty(i,j,k)=buf(i,j,k)
           enddo
         enddo
       enddo
!
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       deallocate(buf)

end subroutine read12v2_mpi

subroutine read12emf1_mpi(ntsub,qty,filename)
  use ModPar
  use ModGrid
  implicit none
  include "mpif.h"

  integer, intent(in) :: ntsub
  character(len=80), intent(in) :: filename
  real, intent(out) :: qty(in,jn,ntsub)

  integer i,j,k
  integer sizes(3), subsizes(3), starts(3)
  integer filetype,fileinfo,fhandl,mysize,ierr
  integer(kind=MPI_OFFSET_KIND) disp
  integer mpi_status(MPI_STATUS_SIZE)
  integer inmax_in, jnmax_in, ntsub_in
  real, allocatable :: buf(:,:,:)
!----------------------------------------------
       allocate(buf(in-1,jn,ntsub))
!
       do k=1,ntsub
         do j=1,jn
           do i=1,in
             qty(i,j,k)=0.D0
           enddo
         enddo
       enddo
!
       sizes(1)=inmax-1
       sizes(2)=jnmax
       sizes(3)=ntsub
!
       subsizes(1)=in-1
       starts(1)=myid1*(in-5)
!
       subsizes(2)=jn
       starts(2)=myid2*(jn-5)
!
       subsizes(3)=ntsub
       starts(3)=0
!
       call MPI_Type_create_subarray(3, sizes, subsizes, starts, &
            MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, filetype, ierr)
       call MPI_Type_commit(filetype, ierr)
       call MPI_Info_create(fileinfo, ierr)
!
       fhandl=13
       call MPI_File_open(MPI_COMM_WORLD, &
            trim(filename), &
            MPI_MODE_RDONLY, &
            fileinfo, fhandl, ierr)
!
       if(myid .eq. 0) then
         call MPI_File_read(fhandl, inmax_in, 1, MPI_INTEGER, &
              mpi_status, ierr)
         call MPI_File_read(fhandl, jnmax_in, 1, MPI_INTEGER, &
              mpi_status, ierr)
         call MPI_File_read(fhandl, ntsub_in, 1, MPI_INTEGER, &
              mpi_status, ierr)
       endif
       call MPI_BCAST(inmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(jnmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(ntsub_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       if(inmax_in .ne. inmax .or.  &
          jnmax_in .ne. jnmax .or.  &
          ntsub_in .ne. ntsub) then
          write(6,*) 'inmax,jnmax,ntsub in file', &
            inmax_in,jnmax_in,ntsub_in, &
            ' do not match inmax, jnmax, ntsub in par.h and timetable'
          call MPI_ABORT(MPI_COMM_WORLD, 1,ierr)
       endif
!
       disp=2*sizeof(inmax)+sizeof(ntsub)
       call MPI_FILE_SET_VIEW(fhandl, disp, &
            MPI_DOUBLE_PRECISION, filetype, 'native', &
            fileinfo, ierr)
!
       mysize=subsizes(1)*subsizes(2)*subsizes(3)
       call MPI_File_read_all(fhandl, buf, &
            mysize, MPI_DOUBLE_PRECISION, mpi_status, ierr)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call MPI_File_close(fhandl, ierr)
       do k=1,ntsub
         do j=1,jn
           do i=1,in-1
             qty(i,j,k)=buf(i,j,k)
           enddo
         enddo
       enddo
!
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       deallocate(buf)

end subroutine read12emf1_mpi

subroutine read12emf2_mpi(ntsub,qty,filename)
  use ModPar
  use ModGrid
  implicit none
  include "mpif.h"

  integer, intent(in) :: ntsub
  character(len=80), intent(in) :: filename
  real, intent(out) :: qty(in,jn,ntsub)

  integer i,j,k
  integer sizes(3), subsizes(3), starts(3)
  integer filetype,fileinfo,fhandl,mysize,ierr
  integer(kind=MPI_OFFSET_KIND) disp
  integer mpi_status(MPI_STATUS_SIZE)
  integer inmax_in, jnmax_in, ntsub_in
  real, allocatable :: buf(:,:,:)
!----------------------------------------------
       allocate(buf(in,jn-1,ntsub))
!
       do k=1,ntsub
         do j=1,jn
           do i=1,in
             qty(i,j,k)=0.D0
           enddo
         enddo
       enddo
!
       sizes(1)=inmax
       sizes(2)=jnmax-1
       sizes(3)=ntsub
!
       subsizes(1)=in
       starts(1)=myid1*(in-5)
!
       subsizes(2)=jn-1
       starts(2)=myid2*(jn-5)
!
       subsizes(3)=ntsub
       starts(3)=0
!
       call MPI_Type_create_subarray(3, sizes, subsizes, starts, &
            MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, filetype, ierr)
       call MPI_Type_commit(filetype, ierr)
       call MPI_Info_create(fileinfo, ierr)
!
       fhandl=13
       call MPI_File_open(MPI_COMM_WORLD, &
            trim(filename), &
            MPI_MODE_RDONLY, &
            fileinfo, fhandl, ierr)
!
       if(myid .eq. 0) then
         call MPI_File_read(fhandl, inmax_in, 1, MPI_INTEGER, &
              mpi_status, ierr)
         call MPI_File_read(fhandl, jnmax_in, 1, MPI_INTEGER, &
              mpi_status, ierr)
         call MPI_File_read(fhandl, ntsub_in, 1, MPI_INTEGER, &
              mpi_status, ierr)
       endif
       call MPI_BCAST(inmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(jnmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(ntsub_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       if(inmax_in .ne. inmax .or.  &
          jnmax_in .ne. jnmax .or.  &
          ntsub_in .ne. ntsub) then
          write(6,*) 'inmax,jnmax,ntsub in file', &
            inmax_in,jnmax_in,ntsub_in, &
            ' do not match inmax, jnmax, ntsub in par.h and timetable'
          call MPI_ABORT(MPI_COMM_WORLD, 1,ierr)
       endif
!
       disp=2*sizeof(inmax)+sizeof(ntsub)
       call MPI_FILE_SET_VIEW(fhandl, disp, &
            MPI_DOUBLE_PRECISION, filetype, 'native', &
            fileinfo, ierr)
!
       mysize=subsizes(1)*subsizes(2)*subsizes(3)
       call MPI_File_read_all(fhandl, buf, &
            mysize, MPI_DOUBLE_PRECISION, mpi_status, ierr)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call MPI_File_close(fhandl, ierr)
       do k=1,ntsub
         do j=1,jn-1
           do i=1,in
             qty(i,j,k)=buf(i,j,k)
           enddo
         enddo
       enddo
!
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       deallocate(buf)

end subroutine read12emf2_mpi

subroutine read23qty_mpi(ntsub,qty,filename)
  use ModPar
  use ModGrid
  implicit none
  include "mpif.h"

  integer, intent(in) :: ntsub
  character(len=80), intent(in) :: filename
  real, intent(out) :: qty(ntsub,jn,kn)

  integer i,j,k
  integer sizes(3), subsizes(3), starts(3)
  integer filetype,fileinfo,fhandl,mysize,ierr
  integer(kind=MPI_OFFSET_KIND) disp
  integer mpi_status(MPI_STATUS_SIZE)
  integer ntsub_in, jnmax_in, knmax_in
  real, allocatable :: buf(:,:,:)
!----------------------------------------------
       allocate(buf(ntsub,jn,kn))
!
       do k=1,kn
         do j=1,jn
           do i=1,ntsub
             qty(i,j,k)=0.D0
           enddo
         enddo
       enddo
!
       sizes(1)=ntsub
       sizes(2)=jnmax
       sizes(3)=knmax
!
       subsizes(1)=ntsub
       starts(1)=0
!
       subsizes(2)=jn
       starts(2)=myid2*(jn-5)
!
       subsizes(3)=kn
       starts(3)=myid3*(kn-5)
!
       call MPI_Type_create_subarray(3, sizes, subsizes, starts, &
            MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, filetype, ierr)
       call MPI_Type_commit(filetype, ierr)
       call MPI_Info_create(fileinfo, ierr)
!
       fhandl=13
       call MPI_File_open(MPI_COMM_WORLD, &
            trim(filename), &
            MPI_MODE_RDONLY, &
            fileinfo, fhandl, ierr)
!
       if(myid .eq. 0) then
         call MPI_File_read(fhandl, ntsub_in, 1, MPI_INTEGER, &
              mpi_status, ierr)
         call MPI_File_read(fhandl, jnmax_in, 1, MPI_INTEGER, &
              mpi_status, ierr)
         call MPI_File_read(fhandl, knmax_in, 1, MPI_INTEGER, &
              mpi_status, ierr)
       endif
       call MPI_BCAST(ntsub_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(jnmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(knmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       if(ntsub_in .ne. ntsub .or.  &
          jnmax_in .ne. jnmax .or.  &
          knmax_in .ne. knmax) then
          write(6,*) 'ntsub,jnmax,knmax in file', &
            ntsub_in,jnmax_in,knmax_in, &
            ' do not match ntsub, jnmax, knmax in par.h and timetable'
          call MPI_ABORT(MPI_COMM_WORLD, 1,ierr)
       endif
!
       disp=2*sizeof(jnmax)+sizeof(ntsub)
       call MPI_FILE_SET_VIEW(fhandl, disp, &
            MPI_DOUBLE_PRECISION, filetype, 'native', &
            fileinfo, ierr)
!
       mysize=subsizes(1)*subsizes(2)*subsizes(3)
       call MPI_File_read_all(fhandl, buf, &
            mysize, MPI_DOUBLE_PRECISION, mpi_status, ierr)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call MPI_File_close(fhandl, ierr)
       do k=1,kn
         do j=1,jn
           do i=1,ntsub
             qty(i,j,k)=buf(i,j,k)
           enddo
         enddo
       enddo
!
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       deallocate(buf)

end subroutine read23qty_mpi

subroutine read23emf2_mpi(qty,filename)
  use ModPar
  use ModGrid
  implicit none
  include "mpif.h"

  character(len=80), intent(in) :: filename
  real, intent(out) :: qty(1,jn,kn)

  integer i,j,k
  integer sizes(3), subsizes(3), starts(3)
  integer filetype,fileinfo,fhandl,mysize,ierr
  integer(kind=MPI_OFFSET_KIND) disp
  integer mpi_status(MPI_STATUS_SIZE)
  integer jnmax_in, knmax_in, icount
  real, allocatable :: buf(:,:,:)
!----------------------------------------------
       allocate(buf(1,jn-1,kn))
!
       qty=0.
!
       sizes(1)=1
       sizes(2)=jnmax-1
       sizes(3)=knmax
!
       subsizes(1)=1
       starts(1)=0
!
       subsizes(2)=jn-1
       starts(2)=myid2*(jn-5)
!
       subsizes(3)=kn
       starts(3)=myid3*(kn-5)
!
       call MPI_Type_create_subarray(3, sizes, subsizes, starts, &
            MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, filetype, ierr)
       call MPI_Type_commit(filetype, ierr)
       call MPI_Info_create(fileinfo, ierr)
!
       fhandl=13
       call MPI_File_open(MPI_COMM_WORLD, &
            trim(filename), &
            MPI_MODE_RDONLY, &
            fileinfo, fhandl, ierr)
!
       if(myid .eq. 0) then
         call MPI_File_read(fhandl, icount, 1, MPI_INTEGER, &
              mpi_status, ierr)
         call MPI_File_read(fhandl, jnmax_in, 1, MPI_INTEGER, &
              mpi_status, ierr)
         call MPI_File_read(fhandl, knmax_in, 1, MPI_INTEGER, &
              mpi_status, ierr)
       endif
       call MPI_BCAST(jnmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(knmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       if(jnmax_in .ne. jnmax .or.  &
          knmax_in .ne. knmax) then
          write(6,*) 'jnmax,knmax,in file', &
            jnmax_in,knmax_in, &
            ' do not match jnmax, knmax in par.h and timetable'
          call MPI_ABORT(MPI_COMM_WORLD, 1,ierr)
       endif
!
       disp=2*sizeof(jnmax)+sizeof(icount)
       call MPI_FILE_SET_VIEW(fhandl, disp, &
            MPI_DOUBLE_PRECISION, filetype, 'native', &
            fileinfo, ierr)
!
       mysize=subsizes(1)*subsizes(2)*subsizes(3)
       call MPI_File_read_all(fhandl, buf, &
            mysize, MPI_DOUBLE_PRECISION, mpi_status, ierr)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call MPI_File_close(fhandl, ierr)
       do k=1,kn
         do j=1,jn-1
             qty(1,j,k)=buf(1,j,k)
         enddo
       enddo
!
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       deallocate(buf)

end subroutine read23emf2_mpi

subroutine read23emf3_mpi(qty,filename)
  use ModPar
  use ModGrid
  implicit none
  include "mpif.h"

  character(len=80), intent(in) :: filename
  real, intent(out) :: qty(1,jn,kn)

  integer i,j,k
  integer sizes(3), subsizes(3), starts(3)
  integer filetype,fileinfo,fhandl,mysize,ierr
  integer(kind=MPI_OFFSET_KIND) disp
  integer mpi_status(MPI_STATUS_SIZE)
  integer jnmax_in, knmax_in, icount
  real, allocatable :: buf(:,:,:)
!----------------------------------------------
       allocate(buf(1,jn,kn-1))
!
       qty=0.
!
       sizes(1)=1
       sizes(2)=jnmax
       sizes(3)=knmax-1
!
       subsizes(1)=1
       starts(1)=0
!
       subsizes(2)=jn
       starts(2)=myid2*(jn-5)
!
       subsizes(3)=kn-1
       starts(3)=myid3*(kn-5)
!
       call MPI_Type_create_subarray(3, sizes, subsizes, starts, &
            MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, filetype, ierr)
       call MPI_Type_commit(filetype, ierr)
       call MPI_Info_create(fileinfo, ierr)
!
       fhandl=13
       call MPI_File_open(MPI_COMM_WORLD, &
            trim(filename), &
            MPI_MODE_RDONLY, &
            fileinfo, fhandl, ierr)
!
       if(myid .eq. 0) then
         call MPI_File_read(fhandl, icount, 1, MPI_INTEGER, &
              mpi_status, ierr)
         call MPI_File_read(fhandl, jnmax_in, 1, MPI_INTEGER, &
              mpi_status, ierr)
         call MPI_File_read(fhandl, knmax_in, 1, MPI_INTEGER, &
              mpi_status, ierr)
       endif
       call MPI_BCAST(jnmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(knmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       if(jnmax_in .ne. jnmax .or.  &
          knmax_in .ne. knmax) then
          write(6,*) 'jnmax,knmax,in file', &
            jnmax_in,knmax_in, &
            ' do not match jnmax, knmax in par.h and timetable'
          call MPI_ABORT(MPI_COMM_WORLD, 1,ierr)
       endif
!
       disp=2*sizeof(jnmax)+sizeof(icount)
       call MPI_FILE_SET_VIEW(fhandl, disp, &
            MPI_DOUBLE_PRECISION, filetype, 'native', &
            fileinfo, ierr)
!
       mysize=subsizes(1)*subsizes(2)*subsizes(3)
       call MPI_File_read_all(fhandl, buf, &
            mysize, MPI_DOUBLE_PRECISION, mpi_status, ierr)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call MPI_File_close(fhandl, ierr)
       do k=1,kn-1
         do j=1,jn
             qty(1,j,k)=buf(1,j,k)
         enddo
       enddo
!
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       deallocate(buf)

end subroutine read23emf3_mpi

subroutine read23emf1_mpi(qty,filename)
  use ModPar
  use ModGrid
  implicit none
  include "mpif.h"

  character(len=80), intent(in) :: filename
  real, intent(out) :: qty(1,jn,kn)

  integer i,j,k
  integer sizes(3), subsizes(3), starts(3)
  integer filetype,fileinfo,fhandl,mysize,ierr
  integer(kind=MPI_OFFSET_KIND) disp
  integer mpi_status(MPI_STATUS_SIZE)
  integer jnmax_in, knmax_in, icount
  real, allocatable :: buf(:,:,:)
!----------------------------------------------
       allocate(buf(1,jn,kn))
!
       qty=0.
!
       sizes(1)=1
       sizes(2)=jnmax
       sizes(3)=knmax
!
       subsizes(1)=1
       starts(1)=0
!
       subsizes(2)=jn
       starts(2)=myid2*(jn-5)
!
       subsizes(3)=kn
       starts(3)=myid3*(kn-5)
!
       call MPI_Type_create_subarray(3, sizes, subsizes, starts, &
            MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, filetype, ierr)
       call MPI_Type_commit(filetype, ierr)
       call MPI_Info_create(fileinfo, ierr)
!
       fhandl=13
       call MPI_File_open(MPI_COMM_WORLD, &
            trim(filename), &
            MPI_MODE_RDONLY, &
            fileinfo, fhandl, ierr)
!
       if(myid .eq. 0) then
         call MPI_File_read(fhandl, icount, 1, MPI_INTEGER, &
              mpi_status, ierr)
         call MPI_File_read(fhandl, jnmax_in, 1, MPI_INTEGER, &
              mpi_status, ierr)
         call MPI_File_read(fhandl, knmax_in, 1, MPI_INTEGER, &
              mpi_status, ierr)
       endif
       call MPI_BCAST(jnmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       call MPI_BCAST(knmax_in,1,MPI_INTEGER,0, &
         MPI_COMM_WORLD,ierr)
       if(jnmax_in .ne. jnmax .or.  &
          knmax_in .ne. knmax) then
          write(6,*) 'jnmax,knmax in file', &
            jnmax_in,knmax_in, &
            ' do not match jnmax, knmax in par.h and timetable'
          call MPI_ABORT(MPI_COMM_WORLD, 1,ierr)
       endif
!
       disp=2*sizeof(jnmax)+sizeof(icount)
       call MPI_FILE_SET_VIEW(fhandl, disp, &
            MPI_DOUBLE_PRECISION, filetype, 'native', &
            fileinfo, ierr)
!
       mysize=subsizes(1)*subsizes(2)*subsizes(3)
       call MPI_File_read_all(fhandl, buf, &
            mysize, MPI_DOUBLE_PRECISION, mpi_status, ierr)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       call MPI_File_close(fhandl, ierr)
       do k=1,kn
         do j=1,jn
             qty(1,j,k)=buf(1,j,k)
         enddo
       enddo
!
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
       deallocate(buf)

end subroutine read23emf1_mpi

end module ModIO
