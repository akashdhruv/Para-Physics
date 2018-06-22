module Poisson_interface

#include "Solver.h"
 
       implicit none

       interface
             subroutine Poisson_solver(rvar,ivar,ps_counter)
             implicit none
             integer, intent(in)  :: ivar,rvar
             integer, intent(out) :: ps_counter
             end subroutine Poisson_solver
       end interface

       interface
             subroutine Poisson_solver_VC(rvar,ivar,ps_counter,cvar,dvar)
             implicit none
             integer, intent(in)  :: ivar,rvar,cvar
             integer, optional, intent(in) :: dvar
             integer, intent(out) :: ps_counter
             end subroutine Poisson_solver_VC
       end interface

       interface
             subroutine Poisson_direct
             implicit none
             end subroutine
       end interface

       interface
             subroutine Poisson_test(tstep,p_counter,ext_mean,num_mean,error_min,error_max)
             implicit none
             integer, intent(in) :: tstep
             integer, intent(out) :: p_counter
             real, intent(out) :: ext_mean,num_mean,error_min,error_max
             end subroutine
       end interface

       interface 
             subroutine Poisson_analytical
             implicit none
             end subroutine
       end interface


end module Poisson_interface
