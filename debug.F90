program debug


        implicit none

        include "mpif.h"

        integer :: solver_comm, shared_comm, myid, procs
        integer :: shared_id,shared_procs,ierr
        integer :: resultlen
        character*(MPI_MAX_PROCESSOR_NAME) :: nameproc


        solver_comm = MPI_COMM_WORLD

        call MPI_INIT(ierr)
        call MPI_COMM_RANK(solver_comm, myid, ierr)
        call MPI_COMM_SIZE(solver_comm, procs, ierr)

        call MPI_GET_PROCESSOR_NAME(nameproc,resultlen,ierr)

        call MPI_COMM_SPLIT_TYPE(solver_comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, shared_comm, ierr)
        !call MPI_COMM_SPLIT(solver_comm,myid/3,myid,shared_comm,ierr)


        call MPI_COMM_RANK(shared_comm,shared_id,ierr)
        call MPI_COMM_SIZE(shared_comm,shared_procs,ierr)

        print *,"Global rank: ",myid,"Shared rank: ",shared_id,"Proc: ",nameproc(1:resultlen)

        call MPI_FINALIZE(ierr)



end program debug
