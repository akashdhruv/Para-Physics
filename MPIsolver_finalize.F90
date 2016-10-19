subroutine MPIsolver_finalize(sim_Complete)

      !$ use omp_lib
      use MPI_data

      implicit none

      logical,intent(in) :: sim_Complete

      include "mpif.h"

      call MPI_COMM_FREE(x_comm,ierr)
      call MPI_COMM_FREE(y_comm,ierr)
  
      !call cpu_time(finish)
      !finish = omp_get_wtime()
      finish = MPI_Wtime()

      exec_time = (finish - start) 

      if (sim_Complete) print '("Execution time: ",f20.10," seconds")',exec_time

      call MPI_FINALIZE(ierr)

end subroutine MPIsolver_finalize
