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

      integer, allocatable, save, dimension(:) :: blockID,blockLC
      integer, allocatable, save, dimension(:) :: xLC, yLC
      integer,save :: blockCount,blockOffset

      integer(kind=MPI_ADDRESS_KIND), save :: center_size,facex_size,facey_size,RMA_size
      integer, save :: disp_unit
      type(C_PTR), save :: center_ptr,facex_ptr, facey_ptr, RMA_ptr
      integer, save :: center_win,facex_win,facey_win,RMA_win

      integer,save :: mpi_info_key

      integer, allocatable, save, dimension(:) :: reqs, req_stat

end module MPI_data
