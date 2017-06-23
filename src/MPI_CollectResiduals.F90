subroutine MPI_CollectResiduals(res,res1,collect_type)

#include "Solver.h"

      use MPI_data

      implicit none
     
      real, intent(inout) :: res,res1
      integer, intent(in) :: collect_type
      integer :: status(MPI_STATUS_SIZE)
      
      if(collect_type == SUM_DATA) then

         call MPI_ALLREDUCE(res,res1,1,MPI_REAL,MPI_SUM,solver_comm,ierr)

      else if(collect_type == MIN_DATA) then

         call MPI_ALLREDUCE(res,res1,1,MPI_REAL,MPI_MIN,solver_comm,ierr)

      else if(collect_type == MAX_DATA) then

         call MPI_ALLREDUCE(res,res1,1,MPI_REAL,MPI_MAX,solver_comm,ierr)

      end if


end subroutine MPI_CollectResiduals
