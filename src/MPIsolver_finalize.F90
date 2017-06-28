subroutine MPIsolver_finalize(sim_Complete)

      !$ use omp_lib
      use MPI_data
      use IncompNS_data, only: ins_timePoisson

      implicit none

      logical,intent(in) :: sim_Complete


      finish = MPI_Wtime()

      exec_time = (finish - start)

#ifdef MPI_DIST
      nullify(localCENTER,localFACEX,localFACEY)
      nullify(world_part)
#endif

#ifdef MPI_SHRD
      call MPI_WIN_FREE(center_win,ierr)
      call MPI_WIN_FREE(facex_win,ierr)
      call MPI_WIN_FREE(facey_win,ierr)
      call MPI_COMM_FREE(shared_comm,ierr)

      deallocate(world_part)
      deallocate(shared_part)
#endif

#ifdef MPI_RMA
      call MPI_WIN_FREE(RMA_win,ierr)
      nullify(localCENTER,localFACEX,localFACEY)
      nullify(eastORIGIN,westORIGIN,northORIGIN,southORIGIN)
      deallocate(world_part)
      deallocate(shared_part)
#endif

      call MPI_COMM_FREE(x_comm,ierr)
      call MPI_COMM_FREE(y_comm,ierr)

      if (sim_Complete) print '("Execution time: ",f20.10," seconds")',exec_time

      print '("Poisson time:   ",f20.10," seconds")',ins_timePoisson

      call MPI_FINALIZE(ierr)

end subroutine MPIsolver_finalize
