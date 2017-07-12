module Multiphase_interface

#include "Solver.h"

implicit none


  interface
    subroutine Multiphase_init()
    end subroutine Multiphase_init
  end interface

  interface 
    subroutine Multiphase_solver(tstep,solnX,jump_flag)
     implicit none
     integer, intent(in) :: tstep
     real, intent(out) :: solnX
     logical, intent(in) :: jump_flag
    end subroutine Multiphase_solver
  end interface

  interface
    subroutine mph_FillVars(s,pf,crv,thco,cprs,visc,rho1x,rho1y,rho2x,rho2y,al1x,al1y,al2x,al2y,nrmx,nrmy,smhv,smrh)
    implicit none
    real, intent(inout), dimension(Nxb+2,Nyb+2) :: pf,thco,cprs,visc,crv,rho1x,rho1y,rho2x,rho2y,&
                                                   al1x,al1y,al2x,al2y,nrmx,nrmy,smhv,smrh
    real, intent(in), dimension(:,:) :: s
    end subroutine mph_FillVars
  end interface

  interface
    subroutine mph_PressureJumps(s,pf,crv,rho1x,rho1y,rho2x,rho2y,w,sigx,sigy,mdot)
    implicit none
    real, dimension(:,:), intent(inout):: crv,rho1x,rho2x,rho1y, &
                                            rho2y,pf,w,sigx,sigy
    real, dimension(:,:), intent(in) :: mdot,s
    end subroutine mph_PressureJumps
  end interface

  interface
    subroutine mph_getInterfaceVelocity(u,v,u_int,v_int,smrh,mdot,nrmx,nrmy)
    implicit none
    real, dimension(:,:), intent(in) :: u, v, smrh, mdot, nrmx, nrmy
    real, dimension(:,:), intent(out) :: u_int,v_int
    end subroutine mph_getInterfaceVelocity
  end interface

 interface
  subroutine mph_FillVars_ibm(s,pf,thco,cprs,visc,rhox,rhoy,alpx,alpy,T,T_old,beta,St)
    implicit none
    real,intent(inout),dimension(Nxb+2,Nyb+2) :: s,pf,thco,cprs,visc,rhox,rhoy,alpx,alpy
    real,intent(in),dimension(Nxb+2,Nyb+2) :: T,T_old
    real,intent(in) :: beta,St
  end subroutine mph_FillVars_ibm
 end interface

end module Multiphase_interface
