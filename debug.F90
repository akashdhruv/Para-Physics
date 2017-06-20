program debug


        implicit none

        include "mpif.h"

        integer :: solver_comm, shared_comm, myid, procs
        integer :: shared_id,shared_procs,ierr
        integer :: resultlen
        character*(MPI_MAX_PROCESSOR_NAME) :: nameproc
        integer(kind=MPI_ADDRESS_KIND) :: sze, baseptr
        integer :: disp_unit,info,win


        sze  = 9000
        info = MPI_INFO_NULL

        solver_comm = MPI_COMM_WORLD

        call MPI_INIT(ierr)
        call MPI_COMM_RANK(solver_comm, myid, ierr)
        call MPI_COMM_SIZE(solver_comm, procs, ierr)

        call MPI_GET_PROCESSOR_NAME(nameproc,resultlen,ierr)

        call MPI_COMM_SPLIT_TYPE(solver_comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, shared_comm, ierr)

        call MPI_COMM_RANK(shared_comm,shared_id,ierr)
        call MPI_COMM_SIZE(shared_comm,shared_procs,ierr)

        print *,"Global rank: ",myid,"Shared rank: ",shared_id,"Proc: ",nameproc(1:resultlen)

        !disp_unit = sze/shared_procs

        !call MPI_WIN_ALLOCATE_SHARED(sze,disp_unit,info,shared_comm,baseptr,win,ierr)

        !print *,"Global rank: ",myid,"Shared rank: ",shared_id,"Proc: ",nameproc(1:resultlen),disp_unit

        !call MPI_WIN_SHARED_QUERY(win,shared_id,sze,disp_unit,baseptr,ierr)

        !print *,"Global rank: ",myid,"Shared rank: ",shared_id,"Proc: ",nameproc(1:resultlen),disp_unit

        call MPI_COMM_FREE(shared_comm,ierr)
        call MPI_FINALIZE(ierr)



end program debug
