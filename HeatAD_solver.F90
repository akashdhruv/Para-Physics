subroutine HeatAD_solver(tstep)

#include "Solver.h"

   use HeatAD_interface, only: HeatAD_SolveTemp
 
   implicit none
   
   integer, intent(in) :: tstep

   call HeatAD_SolveTemp(tstep)

end subroutine
