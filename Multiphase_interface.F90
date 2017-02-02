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
    subroutine mph_FillVars(s,pf,crv,thco,cprs,visc,rho1x,rho1y,rho2x,rho2y,al1x,al1y,al2x,al2y,nrmx,nrmy,T,T_old,beta)
    implicit none
    real, intent(inout), dimension(Nxb+2,Nyb+2) :: s,pf,thco,cprs,visc,crv,rho1x,rho1y,rho2x,rho2y,al1x,al1y,al2x,al2y,nrmx,nrmy
    real, intent(in), dimension(Nxb+2,Nyb+2) :: T,T_old
    real, intent(in) :: beta
    end subroutine mph_FillVars
  end interface

  interface
    subroutine mph_PressureJumps(s,pf,crv,rho1x,rho1y,rho2x,rho2y,w,sigx,sigy,mdot)
    implicit none
    real, dimension(:,:), intent(inout):: s,crv,rho1x,rho2x,rho1y, &
                                            rho2y,pf,w,sigx,sigy
    real, dimension(:,:), intent(in) :: mdot
    end subroutine
  end interface

end module Multiphase_interface
