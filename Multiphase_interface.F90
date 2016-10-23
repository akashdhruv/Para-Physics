module Multiphase_interface

#include "Solver.h"

implicit none


  interface
    subroutine Multiphase_init()
    end subroutine Multiphase_init
  end interface

  interface 
    subroutine Multiphase_solver(tstep,solnX)
     implicit none
     integer, intent(in) :: tstep
     real, intent(out) :: solnX
    end subroutine Multiphase_solver
  end interface

  interface
    subroutine mph_FillVars(s,pf,thco,cprs,visc,rhox,rhoy,alpx,alpy,T,T_old,beta)
    implicit none
    real, intent(inout), dimension(Nxb+2,Nyb+2) :: s,pf,thco,cprs,visc,rhox,rhoy,alpx,alpy
    real, intent(in), dimension(Nxb+2,Nyb+2) :: T,T_old
    real, intent(in) :: beta
    end subroutine mph_FillVars
  end interface

end module Multiphase_interface
