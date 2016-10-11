subroutine IBM_solver(tstep)

#include "Solver.h"

   use IBM_interface, only: ibm_evolve
   use Driver_data, only: dr_dt
   use physicaldata
   use IBM_data

   implicit none
   
   integer, intent(in) :: tstep

   real :: x0,y0,x1,y1,r0,r1
   real :: pi = 3.14

   ibm_y0 = sin(tstep*dr_dt*pi)*ibm_r0
   ibm_y1 = sin(tstep*dr_dt*pi)*ibm_r1

   call ibm_evolve(ibm_x0,ibm_y0,ibm_r0,ibm_x1,ibm_y1,ibm_r1)

end subroutine   
