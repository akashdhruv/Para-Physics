module MPI_data

#include "Solver.h"

      use, intrinsic :: ISO_C_BINDING,  ONLY: C_PTR, C_F_POINTER

      implicit none

      include "mpif.h"

      integer, save :: ierr, myid, procs, solver_comm, x_id, x_procs,x_comm,shared_comm
      integer, save :: y_id, y_procs,y_comm,shared_procs,shared_id,world_group,shared_group

      integer, allocatable, save, dimension(:) :: world_part, shared_part

      integer, save :: mpi_dir_flag

      double precision, save :: start, finish, exec_time

      integer,save :: blockID(MAX_BLOCKS),blockLC(nblockx*nblocky,2)
      integer,save :: blockCount

      integer(kind=MPI_ADDRESS_KIND), save :: center_size,facex_size,facey_size,north_size,east_size
      integer, save :: disp_unit
      type(C_PTR), save :: center_ptr,facex_ptr, facey_ptr
      integer, save :: center_win,facex_win,facey_win,north_win,south_win,east_win,west_win

      integer,save :: mpi_info_key

end module MPI_data
