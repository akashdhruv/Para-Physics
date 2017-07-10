! Main program for the solver which calls three separate units,
! Solver_init     - Initialization of unit specific and physical data
! Solver_evovle   - Time evolution of the solution
! Solver_finalize - Data cleanup

program Solver

   use Solver_interface, only: Solver_init, Solver_evolve, Solver_finalize

   implicit none

    call Solver_init()
    call Solver_evolve()
    call Solver_finalize()

end program Solver
