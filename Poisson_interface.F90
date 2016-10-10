module Poisson_interface

#include "Solver.h"
 
       implicit none

       interface
             subroutine Poisson_solver(ps_RHS,ps,ps_res,ps_counter,ps_quant)
             implicit none
             real, dimension(:,:), intent(in) :: ps_RHS
             real, dimension(:,:), intent(inout) :: ps
             real, intent(out) :: ps_res
             integer, intent(out) :: ps_counter
             integer, intent(in) :: ps_quant
             end subroutine Poisson_solver
       end interface

       interface
             subroutine Poisson_direct
             implicit none
             end subroutine
       end interface

       interface 
             subroutine Poisson_analytical
             implicit none
             end subroutine
       end interface


end module Poisson_interface
