program xmfe
  use ModPar
  use ModPhysics
  use ModControl
  use ModSundry
  use ModGetQtys
  use ModIO
  use ModInitialization
  use ModIteration
  use ModUserSetup,        ONLY: userstep

  implicit none

  include "mpif.h"

  integer :: isize, ierr, update_rate
  real :: wt0, wtnow, wtelpslc, wtelps, wtprev
  real :: em, ek, eth

!-------------------------------------------------------
!     initialize
!-------------------------------------------------------

      call MPI_Init(ierr)
      call MPI_Comm_rank(MPI_COMM_WORLD,myid,ierr)
      call MPI_Comm_size(MPI_COMM_WORLD,isize,ierr)

!     start wall clock
      wt0=MPI_WTIME()

!     check processors
      if(isize .ne. nproc) then
        write(6, *) 'wrong number of processors'
        write(6, *) nproc, ' specified in par.h, ', isize, ' allocated'
        call MPI_ABORT(MPI_COMM_WORLD, 1,ierr)
      endif

      if(myid==0) then
        write(6,'(a)') '-------------------- MFE starts ------------- '
        write(6,'(a,i4,a)') ' Running MFE on ', nproc, ' CPUs'
        write(6,'(3(a,i4))') ' inmax = ', inmax, ' , jnmax = ', &
          jnmax, ', knmax = ', knmax
        write(6,'(3(a,i4))') ' in = ', in, ', jn = ', jn, ', &
          kn = ', kn
        write(6,'(a)') '--------------------------------------------- '
        write(6,*) ' '
        call flush(6)
      endif

      myid3 = myid/(nproc1*nproc2)
      myid2 = (myid-nproc1*nproc2*myid3)/nproc1
      myid1 = myid-(nproc1*nproc2)*myid3-nproc1*myid2

      call initialize

      if(myid==0) then
        write(6,'(a)') '------------ Physical Parameters ------------- '
        write(6,'(a, e23.14)') 'g_const = ', g_const
        write(6,'(a, e23.14)') 'rgas = ', rgas
        write(6,'(a, e23.14)') 'kboltz = ', kboltz
        write(6,'(a, e23.14)') 'mproton = ', mproton
        write(6,'(a, e23.14)') 'pi = ', pi
        write(6,'(a, e23.14)') 'mpovme = ', mpovme
        write(6,'(a, e23.14)') 'gamma = ', gamma
        write(6,'(a, e23.14)') 'muconst = ', muconst
        write(6,'(a, e23.14)') 'unit_rho = ', unit_rho
        write(6,'(a, e23.14)') 'unit_len = ', unit_len
        write(6,'(a, e23.14)') 'unit_b = ', unit_b
        write(6,'(a, e23.14)') 'unit_temp = ', unit_temp
        write(6,'(a, e23.14)') 'unit_v = ', unit_v
        write(6,'(a, e23.14)') 'unit_time = ', unit_time
        write(6,'(a, e23.14)') 'gacc = ', gacc
        write(6,'(a, e23.14)') 'rstar = ', rstar
#ifdef ISOMHD
        write(6,'(a, e23.14)') 'cs2 = ', cs2
#endif
#ifdef THCONDUCT
        write(6,'(a, e23.14)') 'kappa0_dim = ', kappa0_dim
        write(6,'(a, e23.14)') 'kappa0 = ', kappa0
        write(6,'(a, e23.14)') 'diffcoef0 = ', diffcoef0
        write(6,'(a, e23.14)') 'aval0 = ', aval0
#endif
        write(6,'(a)') '--------------------------------------------- '
        write(6,*) ' '
        write(6,'(a)') '------------ sundry parameters ------------- '
        write(6,'(a, i8)') 'itnow = ', itnow
        write(6,'(a, i8)') 'ntcond = ', ntcond
        write(6,'(a, e23.14)') 'dt = ', dt
        write(6,'(a, e23.14)') 'time = ', time
        write(6,'(a, e23.14)') 'ctime = ', ctime
        write(6,'(a, e23.14)') 'courno = ', courno
        write(6,'(a, e23.14)') 'dtmdi2 = ', dtmdi2
        write(6,'(a, i8)') 'nlf_d = ', nlf_d
        write(6,'(a, i8)') 'nlf_v = ', nlf_v
        write(6,'(a, e23.14)') 'alpha_d = ', alpha_d
        write(6,'(a, e23.14)') 'alpha_v = ', alpha_v
        write(6,'(a, e23.14)') 'dfloor = ', dfloor
#ifndef ISOMHD
        write(6,'(a, e23.14)') 'efloor = ', efloor
#endif
        write(6,'(a)') '--------------------------------------------- '
        write(6,*) ' '
      endif

!     open monitoring file
      if(myid .eq. 0) then
        open(unit=16, file=engfile, position="append")
      endif
#ifndef ISOMHD
      call get_global(em,ek,eth)
#else
      call get_global(em,ek)
#endif
      if(myid .eq. 0) then
#ifndef ISOMHD
        write(16,'(2i8, 5e23.14)') itnow, ntcond, time, dt, em, ek, eth
#else
        write(16,'(2i8, 4e23.14)') itnow, ntcond, time, dt, em, ek
#endif
        call flush(16)
      endif

!     output initial field
      if(ifile .lt. 0) then
        ifile=time/tout
        tfile=tout*ifile+tout
      else
        call writedata_mpi
        tfile=time+tout
      endif
      ifile=ifile+1

!     output the grid

      call writegrid
!
!------------------------------------------------------------
!      main iteration loop
!------------------------------------------------------------
!
!      start iteration
!
       irst=0
       wtnow=MPI_WTIME()
       wtprev=wtnow
       wtelpslc=wtnow-wt0
       call MPI_ALLREDUCE(wtelpslc,wtelps,1, &
         MPI_DOUBLE_PRECISION,MPI_MAX, &
         MPI_COMM_WORLD,ierr)

       do while (time .lt. tend .and. itnow .lt. itmax &
                 .and. wtelps .lt. wtlimit)

       call nudt
       call rk3step
       call userstep
       time=time+dt
       ctime=time
       itnow=itnow+1
       wtnow=MPI_WTIME()
       wtelpslc=wtnow-wt0
       call MPI_ALLREDUCE(wtelpslc,wtelps,1, &
         MPI_DOUBLE_PRECISION,MPI_MAX, &
         MPI_COMM_WORLD,ierr)
!
!      monitoring global quantities
!
       if(mod(itnow,itintv) .eq. 0) then
#ifndef ISOMHD
         call get_global(em,ek,eth)
#else
         call get_global(em,ek)
#endif
         if(myid .eq. 0) then
#ifndef ISOMHD
           write(16,'(2i8, 5e23.14)') itnow, ntcond, time, dt, em, ek, eth
#else
           write(16,'(2i8, 4e23.14)') itnow, ntcond, time, dt, em, ek
#endif
           call flush(16)
           update_rate = itintv*(inmax-5)*(jnmax-5)*(knmax-5)/nproc/(wtnow - wtprev)
           write(6,'(a,I7,a,I9)') ' Grid cell updates per CPU sec = ', &
             update_rate, ' at iteration = ', itnow
           wtprev = wtnow
           call flush(6)
         endif
       endif
!
!      check for output time
!
       if(mod(itnow,itout) .eq. itout-1) then
         call writerst_mpi
       endif
!
       if(time .ge. tfile) then
         call writedata_mpi
         ifile=ifile+1
         tfile=tfile+tout
       endif

       enddo
!
!-------------------------------------------------
!      finalize
!-------------------------------------------------
!
#ifndef ISOMHD
       call get_global(em,ek,eth)
#else
       call get_global(em,ek)
#endif
       if(myid .eq. 0) then
#ifndef ISOMHD
         write(16,'(2i8, 5e23.14)') itnow, ntcond, time, dt, em, ek, eth
#else
         write(16,'(2i8, 4e23.14)') itnow, ntcond, time, dt, em, ek
#endif
         close(16)
       endif
 
       call writedata_mpi
       call writerst_mpi
 
       call MPI_Finalize(ierr)

end program xmfe
