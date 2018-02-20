subroutine ins_predictor(tstep,u,v,ut,vt,g1_old,g2_old,temp)

       use Grid_data
       use Driver_data
       use IncompNS_data
       use HeatAD_data, only: ht_alpha

#include "Solver.h"

       implicit none

       !-Arguments
       integer, intent(in) :: tstep
       real, intent(inout), dimension(:,:) :: u, v
       real, intent(inout), dimension(:,:) :: ut,vt,g1_old,g2_old
       real, intent(in), dimension(:,:) :: temp

       !-Local Variables
       real, dimension(Nxb,Nyb) :: C1,G1,D1,C2,G2,D2

       ! Predictor Step

       call Convective_U(u,v,gr_dx,gr_dy,C1)
       call Diffusive_U(u,gr_dx,gr_dy,ins_inRe,D1)

       G1 = C1 + D1 + ins_gravX

       if (tstep == 0) then

              ut(2:Nxb+1,2:Nyb+1)=u(2:Nxb+1,2:Nyb+1)+(dr_dt/1)*(G1)
              g1_old(2:Nxb+1,2:Nyb+1) = G1
       else

              ut(2:Nxb+1,2:Nyb+1)=u(2:Nxb+1,2:Nyb+1)+(dr_dt/2)*(3*g1_old(2:Nxb+1,2:Nyb+1)-G1)
              g1_old(2:Nxb+1,2:Nyb+1) = G1
       endif


       call Convective_V(u,v,gr_dx,gr_dy,C2)
       call Diffusive_V(v,gr_dx,gr_dy,ins_inRe,D2)

       G2 = C2 + D2 - ins_gravY*ht_alpha*(temp(2:Nxb+1,2:Nyb+1)+temp(2:Nxb+1,3:Nyb+2))/2.0d0

       if (tstep == 0) then

              vt(2:Nxb+1,2:Nyb+1)=v(2:Nxb+1,2:Nyb+1)+(dr_dt/1)*(G2)
              g2_old(2:Nxb+1,2:Nyb+1) = G2
       else

              vt(2:Nxb+1,2:Nyb+1)=v(2:Nxb+1,2:Nyb+1)+(dr_dt/2)*(3*g2_old(2:Nxb+1,2:Nyb+1)-G2)
              g2_old(2:Nxb+1,2:Nyb+1) = G2
       endif

end subroutine ins_predictor


!! CONVECTIVE U !!
subroutine Convective_U(ut_cal,vt_cal,dx,dy,C1)

#include "Solver.h"

      implicit none

      real,dimension(Nxb+2,Nyb+2), intent(in) :: ut_cal
      real,dimension(Nxb+2,Nyb+2), intent(in) :: vt_cal

      real,intent(in) :: dx
      real,intent(in) :: dy

      !real, allocatable, dimension(:,:) :: ue
      !real, allocatable, dimension(:,:) :: uw
      !real, allocatable, dimension(:,:) :: us
      !real, allocatable, dimension(:,:) :: un
      !real, allocatable, dimension(:,:) :: vs
      !real, allocatable, dimension(:,:) :: vn

      real, dimension(Nxb,Nyb) :: ue,uw,us,un,vs,vn
      real, dimension(Nxb,Nyb), intent(out) :: C1

      !allocate(ue(Nxb,Nyb))
      !allocate(uw(Nxb,Nyb))
      !allocate(us(Nxb,Nyb))
      !allocate(un(Nxb,Nyb))
      !allocate(vs(Nxb,Nyb))
      !allocate(vn(Nxb,Nyb))

      ue = (ut_cal(2:Nxb+1,2:Nyb+1)+ut_cal(3:Nxb+2,2:Nyb+1))/2
      uw = (ut_cal(2:Nxb+1,2:Nyb+1)+ut_cal(1:Nxb,2:Nyb+1))/2
      us = (ut_cal(2:Nxb+1,2:Nyb+1)+ut_cal(2:Nxb+1,1:Nyb))/2
      un = (ut_cal(2:Nxb+1,2:Nyb+1)+ut_cal(2:Nxb+1,3:Nyb+2))/2
      vs = (vt_cal(2:Nxb+1,1:Nyb)+vt_cal(3:Nxb+2,1:Nyb))/2
      vn = (vt_cal(2:Nxb+1,2:Nyb+1)+vt_cal(3:Nxb+2,2:Nyb+1))/2

      C1 = -((ue**2)-(uw**2))/dx - ((un*vn)-(us*vs))/dy

      !deallocate(ue,uw,us,un,vs,vn)

end subroutine Convective_U

!! CONVECTIVE V !!
subroutine Convective_V(ut_cal,vt_cal,dx,dy,C2)

#include "Solver.h"

      implicit none

      real,dimension(Nxb+2,Nyb+2), intent(in) :: ut_cal
      real,dimension(Nxb+2,Nyb+2), intent(in) :: vt_cal

      real, intent(in) :: dx
      real, intent(in) :: dy

      !real, allocatable,dimension(:,:) :: vn, vs, ve, vw, ue, uw
      real,dimension(Nxb,Nyb) :: vn,vs,ve,vw,ue,uw
      real, dimension(Nxb,Nyb), intent(out) :: C2

      !allocate(vn(Nxb,Nyb))
      !allocate(vs(Nxb,Nyb))
      !allocate(ve(Nxb,Nyb))
      !allocate(vw(Nxb,Nyb))
      !allocate(ue(Nxb,Nyb))
      !allocate(uw(Nxb,Nyb))

      vs = (vt_cal(2:Nxb+1,2:Nyb+1)+vt_cal(2:Nxb+1,1:Nyb))/2
      vn = (vt_cal(2:Nxb+1,2:Nyb+1)+vt_cal(2:Nxb+1,3:Nyb+2))/2
      ve = (vt_cal(2:Nxb+1,2:Nyb+1)+vt_cal(3:Nxb+2,2:Nyb+1))/2
      vw = (vt_cal(2:Nxb+1,2:Nyb+1)+vt_cal(1:Nxb,2:Nyb+1))/2
      ue = (ut_cal(2:Nxb+1,2:Nyb+1)+ut_cal(2:Nxb+1,3:Nyb+2))/2
      uw = (ut_cal(1:Nxb,2:Nyb+1)+ut_cal(1:Nxb,3:Nyb+2))/2

      C2 = -((ue*ve)-(uw*vw))/dx - ((vn**2)-(vs**2))/dy

      !deallocate(vn,vs,ve,vw,ue,uw)

end subroutine Convective_V

!! DIFFUSIVE U !!
subroutine Diffusive_U(ut_cal,dx,dy,inRe,D1)

#include "Solver.h"

      implicit none

      real,dimension(Nxb+2,Nyb+2), intent(in) :: ut_cal

      real, intent(in) :: dx
      real, intent(in) :: dy

      real, intent(in) :: inRe

      !real, allocatable, dimension(:,:) :: uP
      !real, allocatable, dimension(:,:) :: uN
      !real, allocatable, dimension(:,:) :: uS
      !real, allocatable, dimension(:,:) :: uE
      !real, allocatable, dimension(:,:) :: uW

      real,dimension(Nxb,Nyb) :: uP,uN,uS,uE,uW
      real, dimension(Nxb,Nyb), intent(out) :: D1

      !allocate(uP(Nxb,Nyb))
      !allocate(uN(Nxb,Nyb))
      !allocate(uS(Nxb,Nyb))
      !allocate(uE(Nxb,Nyb))
      !allocate(uW(Nxb,Nyb))

      uP = ut_cal(2:Nxb+1,2:Nyb+1)
      uE = ut_cal(3:Nxb+2,2:Nyb+1)
      uW = ut_cal(1:Nxb,2:Nyb+1)
      uN = ut_cal(2:Nxb+1,3:Nyb+2)
      uS = ut_cal(2:Nxb+1,1:Nyb)

      !D1 = (inRe/dx)*(((uE-uP)/dx)-((uP-uW)/dx)) + (inRe/dy)*(((uN-uP)/dy)-((uP-uS)/dy))

      !D1 = (inRe/dx)*((uE-uP)/dx)&
      !    -(inRe/dx)*((uP-uW)/dx)&
      !    +(inRe/dy)*((uN-uP)/dy)&
      !    -(inRe/dy)*((uP-uS)/dy)

      D1 = ((inRe*(uE-uP)/dx)-(inRe*(uP-uW)/dx))/dx + &
           ((inRe*(uN-uP)/dy)-(inRe*(uP-uS)/dy))/dy

      !deallocate(uP,uN,uS,uE,uW)

end subroutine Diffusive_U

!! DIFFUSIVE V !!
subroutine Diffusive_V(vt_cal,dx,dy,inRe,D2)

#include "Solver.h"

      implicit none

      real,dimension(Nxb+2,Nyb+2), intent(in) :: vt_cal

      real, intent(in) :: dx
      real, intent(in) :: dy

      real, intent(in) :: inRe

      !real, allocatable, dimension(:,:) :: vP,vE,vW,vN,vS
      real,dimension(Nxb,Nyb) :: vP,vE,vW,vN,vS

      real, dimension(Nxb,Nyb), intent(out) :: D2

      !allocate(vP(Nxb,Nyb))
      !allocate(vE(Nxb,Nyb))
      !allocate(vW(Nxb,Nyb)) 
      !allocate(vN(Nxb,Nyb))
      !allocate(vS(Nxb,Nyb))

      vP = vt_cal(2:Nxb+1,2:Nyb+1)
      vE = vt_cal(3:Nxb+2,2:Nyb+1)
      vW = vt_cal(1:Nxb,2:Nyb+1)
      vN = vt_cal(2:Nxb+1,3:Nyb+2)
      vS = vt_cal(2:Nxb+1,1:Nyb)

      !D2 = (inRe/dx)*(((vE-vP)/dx)-((vP-vW)/dx)) + (inRe/dy)*(((vN-vP)/dy)-((vP-vS)/dy))

      !D2 = (inRe/dx)*((vE-vP)/dx)&
      !    -(inRe/dx)*((vP-vW)/dx)&
      !    +(inRe/dy)*((vN-vP)/dy)&
      !    -(inRe/dy)*((vP-vS)/dy)

      D2 = ((inRe*(vE-vP)/dx)-(inRe*(vP-vW)/dx))/dx + &
           ((inRe*(vN-vP)/dy)-(inRe*(vP-vS)/dy))/dy

      !deallocate(vP,vE,vW,vN,vS)

end subroutine Diffusive_V
