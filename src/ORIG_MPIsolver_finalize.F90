subroutine MPIsolver_finalize(sim_Complete)

      !$ use omp_lib
      use MPI_data
      use IncompNS_data, only: ins_timePoisson

      implicit none

      logical,intent(in) :: sim_Complete

      call MPI_WIN_FREE(center_win,ierr)
      call MPI_WIN_FREE(facex_win,ierr)
      call MPI_WIN_FREE(facey_win,ierr)

      call MPI_COMM_FREE(x_comm,ierr)
      call MPI_COMM_FREE(y_comm,ierr)
      call MPI_COMM_FREE(shared_comm,ierr)  

      finish = MPI_Wtime()

      exec_time = (finish - start) 

      if (sim_Complete) print '("Execution time: ",f20.10," seconds")',exec_time

      print '("Poisson time:   ",f20.10," seconds")',ins_timePoisson

      call MPI_FINALIZE(ierr)

end subroutine MPIsolver_finalize
