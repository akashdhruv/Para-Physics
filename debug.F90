program debug


        use, intrinsic :: ISO_C_BINDING, ONLY: C_PTR, C_F_POINTER

        implicit none

        include "mpif.h"

        integer :: solver_comm, shared_comm, myid, procs
        integer :: shared_id,shared_procs,ierr
        integer :: resultlen
        character*(MPI_MAX_PROCESSOR_NAME) :: nameproc
        integer(kind=MPI_ADDRESS_KIND) :: sze
        integer :: disp_unit

        type(C_PTR) :: baseptr
        integer :: win

        real,pointer,dimension(:,:) :: shared_data

        integer :: Nx,Ny

        real :: a

        Nx = 10
        Ny = 1

        solver_comm = MPI_COMM_WORLD

        call MPI_INIT(ierr)
        call MPI_COMM_RANK(solver_comm, myid, ierr)
        call MPI_COMM_SIZE(solver_comm, procs, ierr)

        call MPI_GET_PROCESSOR_NAME(nameproc,resultlen,ierr)

        call MPI_COMM_SPLIT_TYPE(solver_comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, shared_comm, ierr)

        call MPI_COMM_RANK(shared_comm,shared_id,ierr)
        call MPI_COMM_SIZE(shared_comm,shared_procs,ierr)

        !print *,"Global rank: ",myid,"Shared rank: ",shared_id,"Proc: ",nameproc(1:resultlen)
 
        sze = sizeof(a)
        disp_unit = sizeof(a)

        call MPI_WIN_ALLOCATE_SHARED(sze,disp_unit,MPI_INFO_NULL,shared_comm,baseptr,win,ierr)

        call MPI_WIN_SHARED_QUERY(win, shared_id, sze, disp_unit, baseptr)

        call MPI_BARRIER(shared_comm,ierr)

        !call C_F_POINTER(baseptr, shared_data,[Nx,Ny])

        call MPI_WIN_LOCK_ALL(0,win,ierr)

        !shared_data(:,1) = 0

        if(shared_id == 1) shared_data(1,1) = 5 !shared_data(4,1) = 5

        call MPI_BARRIER(shared_comm,ierr)

        if(shared_id == 5) shared_data(1,1) = 2*shared_data(1,1) !shared_data(4,1) = 2*shared_data(4,1)
      
        call MPI_WIN_SYNC(win,ierr)
        call MPI_BARRIER(shared_comm,ierr)

        call MPI_WIN_UNLOCK_ALL(win,ierr)

        if(shared_id == 0) print *,"rank = ",shared_id,"memory = ",shared_data(:,1),sze,disp_unit
    

        call MPI_WIN_FREE(win,ierr)

        call MPI_COMM_FREE(shared_comm,ierr)
        call MPI_FINALIZE(ierr)

end program debug
