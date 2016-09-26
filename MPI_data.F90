module MPI_data

#include "Solver.h"

      implicit none


      integer, save :: ierr, myid, procs, solver_comm, x_id, x_procs, x_comm
      integer, save :: y_id, y_procs, y_comm

      integer, save :: mpi_dir_flag

      double precision, save :: start, finish, exec_time

      integer,save :: blockID(MAX_BLOCKS),blockLC(nblockx*nblocky,2)
      integer,save :: blockCount

end module MPI_data
