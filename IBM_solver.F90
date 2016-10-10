subroutine IBM_solver(tstep)

#include "Solver.h"

   use IBM_interface, only: ibm_evolve
   use Driver_data, only: dr_dt
   use physicaldata

   implicit none
   
   integer, intent(in) :: tstep

   real :: x0,y0,x1,y1,r0,r1
   real :: pi = 3.14

   r0 = 0.15
   r1 = 0.10

   y0 = sin(tstep*dr_dt*pi)*r0
   y1 = sin(tstep*dr_dt*pi)*r1

   x0 = -0.40
   x1 = -0.10

   call ibm_evolve(x0,y0,r0,x1,y1,r1)

end subroutine   
