subroutine MPIsolver_init()

#include "Solver.h"

    !$ use omp_lib
    use MPI_data
    use Grid_data, only: blockCount,blockID
    implicit none

    integer :: status = 0

    include "mpif.h"

    solver_comm = MPI_COMM_WORLD
    
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(solver_comm, myid, ierr)
    call MPI_COMM_SIZE(solver_comm, procs, ierr)    

    if (nblockx*nblocky /= procs) then
       
        call MPI_FINALIZE(ierr)
        if (myid == 0) print *,"RUNTIME ERROR: The number of MPI processes must be equal to total number of blocks."
        call exit(status)

    end if
 
    call MPI_COMM_SPLIT(solver_comm,myid/nblockx,myid,x_comm,ierr)
    call MPI_COMM_SPLIT(solver_comm,mod(myid,nblockx),myid,y_comm,ierr)

    call MPI_COMM_RANK(x_comm,x_id,ierr)
    call MPI_COMM_SIZE(x_comm,x_procs,ierr)

    call MPI_COMM_RANK(y_comm,y_id,ierr)
    call MPI_COMM_size(y_comm,y_procs,ierr)

    !call cpu_time(start)
    start = omp_get_wtime()

end subroutine MPIsolver_init 
