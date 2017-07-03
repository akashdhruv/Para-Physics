subroutine ins_predictor_VD(tstep,u,v,ut,vt,g1_old,g2_old,visc,rho1x,rho1y,rho2x,rho2y)

       use Grid_data
       use Driver_data
       use IncompNS_data

#include "Solver.h"

       implicit none

       !-Arguments
       integer, intent(in) :: tstep
       real, intent(inout), dimension(:,:) :: u,v,visc,rho1x,rho1y,rho2x,rho2y,ut,vt,g1_old,g2_old

       !-Local Variables
       real, dimension(Nxb,Nyb) :: C1,G1,D1,C2,G2,D2

       ! Predictor Step

       call Convective_U_VD(u,v,gr_dx,gr_dy,C1)
       call Diffusive_U_VD(u,visc,rho1x,rho2x,gr_dx,gr_dy,ins_inRe,D1)

       G1 = C1 + D1 + ins_gravX

       if (tstep == 0) then

              ut(2:Nxb+1,2:Nyb+1)=u(2:Nxb+1,2:Nyb+1)+(dr_dt/1)*(G1)
              g1_old(2:Nxb+1,2:Nyb+1) = G1
       else

              ut(2:Nxb+1,2:Nyb+1)=u(2:Nxb+1,2:Nyb+1)+(dr_dt/2)*(3*g1_old(2:Nxb+1,2:Nyb+1)-G1)
              g1_old(2:Nxb+1,2:Nyb+1) = G1
       endif


       call Convective_V_VD(u,v,gr_dx,gr_dy,C2)
       call Diffusive_V_VD(v,visc,rho1y,rho2y,gr_dx,gr_dy,ins_inRe,D2)

       G2 = C2 + D2 + ins_gravY

       if (tstep == 0) then

              vt(2:Nxb+1,2:Nyb+1)=v(2:Nxb+1,2:Nyb+1)+(dr_dt/1)*(G2)
              g2_old(2:Nxb+1,2:Nyb+1) = G2
       else

              vt(2:Nxb+1,2:Nyb+1)=v(2:Nxb+1,2:Nyb+1)+(dr_dt/2)*(3*g2_old(2:Nxb+1,2:Nyb+1)-G2)
              g2_old(2:Nxb+1,2:Nyb+1) = G2
       endif

end subroutine ins_predictor_VD


!! CONVECTIVE U !!
subroutine Convective_U_VD(ut,vt,dx,dy,C1)

#include "Solver.h"
       
      implicit none

      real,dimension(Nxb+2,Nyb+2), intent(in) :: ut
      real,dimension(Nxb+2,Nyb+2), intent(in) :: vt

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

      ue = (ut(2:Nxb+1,2:Nyb+1)+ut(3:Nxb+2,2:Nyb+1))/2
      uw = (ut(2:Nxb+1,2:Nyb+1)+ut(1:Nxb,2:Nyb+1))/2
      us = (ut(2:Nxb+1,2:Nyb+1)+ut(2:Nxb+1,1:Nyb))/2
      un = (ut(2:Nxb+1,2:Nyb+1)+ut(2:Nxb+1,3:Nyb+2))/2
      vs = (vt(2:Nxb+1,1:Nyb)+vt(3:Nxb+2,1:Nyb))/2
      vn = (vt(2:Nxb+1,2:Nyb+1)+vt(3:Nxb+2,2:Nyb+1))/2

      C1 = -((ue**2)-(uw**2))/dx - ((un*vn)-(us*vs))/dy

      !deallocate(ue,uw,us,un,vs,vn)

end subroutine Convective_U_VD

!! CONVECTIVE V !!
subroutine Convective_V_VD(ut,vt,dx,dy,C2)

#include "Solver.h"

      implicit none

      real,dimension(Nxb+2,Nyb+2), intent(in) :: ut
      real,dimension(Nxb+2,Nyb+2), intent(in) :: vt

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

      vs = (vt(2:Nxb+1,2:Nyb+1)+vt(2:Nxb+1,1:Nyb))/2
      vn = (vt(2:Nxb+1,2:Nyb+1)+vt(2:Nxb+1,3:Nyb+2))/2
      ve = (vt(2:Nxb+1,2:Nyb+1)+vt(3:Nxb+2,2:Nyb+1))/2
      vw = (vt(2:Nxb+1,2:Nyb+1)+vt(1:Nxb,2:Nyb+1))/2
      ue = (ut(2:Nxb+1,2:Nyb+1)+ut(2:Nxb+1,3:Nyb+2))/2
      uw = (ut(1:Nxb,2:Nyb+1)+ut(1:Nxb,3:Nyb+2))/2

      C2 = -((ue*ve)-(uw*vw))/dx - ((vn**2)-(vs**2))/dy

      !deallocate(vn,vs,ve,vw,ue,uw)

end subroutine Convective_V_VD

!! DIFFUSIVE U !!
subroutine Diffusive_U_VD(ut,visc,rho1x,rho2x,dx,dy,inRe,D1)

#include "Solver.h"

      implicit none

      real,dimension(Nxb+2,Nyb+2), intent(in) :: ut
      real,dimension(Nxb+2,Nyb+2), intent(in) :: visc,rho1x,rho2x
      real, intent(in) :: dx
      real, intent(in) :: dy

      real, intent(in) :: inRe

      !real, allocatable, dimension(:,:) :: uP
      !real, allocatable, dimension(:,:) :: uN
      !real, allocatable, dimension(:,:) :: uS
      !real, allocatable, dimension(:,:) :: uE
      !real, allocatable, dimension(:,:) :: uW

      real,dimension(Nxb,Nyb) :: uP,uN,uS,uE,uW,txxp,txxm,tyyp,tyym
      real, dimension(Nxb,Nyb), intent(out) :: D1

      !allocate(uP(Nxb,Nyb))
      !allocate(uN(Nxb,Nyb))
      !allocate(uS(Nxb,Nyb))
      !allocate(uE(Nxb,Nyb))
      !allocate(uW(Nxb,Nyb))

      uP = ut(2:Nxb+1,2:Nyb+1)
      uE = ut(3:Nxb+2,2:Nyb+1)
      uW = ut(1:Nxb,2:Nyb+1)
      uN = ut(2:Nxb+1,3:Nyb+2)
      uS = ut(2:Nxb+1,1:Nyb)

      txxp = (inRe*(uE-uP)/dx)*visc(2:Nxb+1,2:Nyb+1)
      txxm = (inRe*(uP-uW)/dx)*visc(1:Nxb,2:Nyb+1)

      tyyp = inRe*(uN-uP)/dy*0.25*(visc(2:Nxb+1,2:Nyb+1)+visc(1:Nxb,2:Nyb+1)&
                                  +visc(1:Nxb,3:Nyb+2)+visc(2:Nxb+1,3:Nyb+2))

      tyym = inRe*(uP-uS)/dy*0.25*(visc(2:Nxb+1,2:Nyb+1)+visc(1:Nxb,2:Nyb+1)&
                                  +visc(1:Nxb,1:Nyb)+visc(2:Nxb+1,1:Nyb))

      !D1 = (inRe/dx)*(((uE-uP)/dx)-((uP-uW)/dx)) + (inRe/dy)*(((uN-uP)/dy)-((uP-uS)/dy))

      !D1 = (inRe/dx)*((uE-uP)/dx)&
      !    -(inRe/dx)*((uP-uW)/dx)&
      !    +(inRe/dy)*((uN-uP)/dy)&
      !    -(inRe/dy)*((uP-uS)/dy)

      D1 = (txxp-txxm)/dx + &
           (tyyp-tyym)/dy

      D1 = D1*(rho1x(2:Nxb+1,2:Nyb+1)+rho2x(2:Nxb+1,2:Nyb+1))

      !deallocate(uP,uN,uS,uE,uW)

end subroutine Diffusive_U_VD

!! DIFFUSIVE V !!
subroutine Diffusive_V_VD(vt,visc,rho1y,rho2y,dx,dy,inRe,D2)

#include "Solver.h"

      implicit none

      real,dimension(Nxb+2,Nyb+2), intent(in) :: vt
      real,dimension(Nxb+2,Nyb+2), intent(in) :: visc,rho1y,rho2y

      real, intent(in) :: dx
      real, intent(in) :: dy

      real, intent(in) :: inRe

      !real, allocatable, dimension(:,:) :: vP,vE,vW,vN,vS
      real,dimension(Nxb,Nyb) :: vP,vE,vW,vN,vS,txxp,txxm,tyyp,tyym

      real, dimension(Nxb,Nyb), intent(out) :: D2

      !allocate(vP(Nxb,Nyb))
      !allocate(vE(Nxb,Nyb))
      !allocate(vW(Nxb,Nyb)) 
      !allocate(vN(Nxb,Nyb))
      !allocate(vS(Nxb,Nyb))

      vP = vt(2:Nxb+1,2:Nyb+1)
      vE = vt(3:Nxb+2,2:Nyb+1)
      vW = vt(1:Nxb,2:Nyb+1)
      vN = vt(2:Nxb+1,3:Nyb+2)
      vS = vt(2:Nxb+1,1:Nyb)

      txxp = (inRe*(vE-vP)/dx)*0.25*(visc(2:Nxb+1,2:Nyb+1)+visc(2:Nxb+1,1:Nyb)&
                                    +visc(3:Nxb+2,2:Nyb+1)+visc(3:Nxb+2,1:Nyb))

      txxm = (inRe*(vP-vW)/dx)*0.25*(visc(2:Nxb+1,2:Nyb+1)+visc(2:Nxb+1,1:Nyb)&
                                    +visc(1:Nxb,2:Nyb+1)+visc(1:Nxb,1:Nyb))

      tyyp = (inRe*(vN-vP)/dy)*visc(2:Nxb+1,2:Nyb+1)
      tyym = (inRe*(vP-vS)/dy)*visc(2:Nxb+1,1:Nyb)

      !D2 = (inRe/dx)*(((vE-vP)/dx)-((vP-vW)/dx)) + (inRe/dy)*(((vN-vP)/dy)-((vP-vS)/dy))

      !D2 = (inRe/dx)*((vE-vP)/dx)&
      !    -(inRe/dx)*((vP-vW)/dx)&
      !    +(inRe/dy)*((vN-vP)/dy)&
      !    -(inRe/dy)*((vP-vS)/dy)

      D2 = (txxp-txxm)/dx + &
           (tyyp-tyym)/dy

      D2 = D2*(rho1y(2:Nxb+1,2:Nyb+1)+rho2y(2:Nxb+1,2:Nyb+1))
  
      !deallocate(vP,vE,vW,vN,vS)

end subroutine Diffusive_V_VD

