subroutine MPIsolver_init()

#include "Solver.h"

    !$ use omp_lib

    use MPI_data
    use physicaldata, only: blockCount,blockID,blockLC
    use morton_interface, only: morton_sort

    implicit none

    integer :: checkSumMPI
    integer :: status = 0
    integer :: i

    include "mpif.h"

    solver_comm = MPI_COMM_WORLD
    
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(solver_comm, myid, ierr)
    call MPI_COMM_SIZE(solver_comm, procs, ierr)    

    blockCount = ((nblockx*nblocky)/procs)
    checkSumMPI = blockCount*procs

    if (checkSumMPI /= nblockx*nblocky) then
       
        call MPI_FINALIZE(ierr)

        if (myid == 0) then 
        print *,"RUNTIME ERROR: The number blocks should be greater than and exaclty divisible by total number of MPI processes."
        end if

        call exit(status)

    end if

    allocate(blockID(blockCount))
    allocate(blockLC(nblockx*nblocky,2))

    blockID = (/(I,I=1,blockCount)/)
 
    call morton_sort(blockCount,myid,procs,blockLC)

    call MPI_COMM_SPLIT(solver_comm,myid/nblockx,myid,x_comm,ierr)
    call MPI_COMM_SPLIT(solver_comm,mod(myid,nblockx),myid,y_comm,ierr)

    call MPI_COMM_RANK(x_comm,x_id,ierr)
    call MPI_COMM_SIZE(x_comm,x_procs,ierr)

    call MPI_COMM_RANK(y_comm,y_id,ierr)
    call MPI_COMM_size(y_comm,y_procs,ierr)

    !call cpu_time(start)
    start = omp_get_wtime()

end subroutine MPIsolver_init 
