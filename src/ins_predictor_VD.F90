subroutine ins_predictor_VD(tstep,p_counter,p,u,v,ut,vt,visc,rho1x,rho1y,rho2x,rho2y,s,s2,sigp,sigx,sigy)

       use Poisson_interface, ONLY: Poisson_solver_VC            
       use Grid_data
       use Driver_data
       use MPI_data
       use IncompNS_data
       use MPI_interface, ONLY: MPI_applyBC, MPI_CollectResiduals, MPI_physicalBC_vel, MPI_applyBC_RMA
       use IncompNS_interface, ONLY: ins_rescaleVel
       use IBM_interface, ONLY: IBM_ApplyForcing

#include "Solver.h"

       implicit none

       !-Arguments
       integer, intent(in) :: tstep
       integer, intent(out) :: p_counter

       !-Local Variables

       !real, allocatable, dimension(:,:) :: ut,vt,u_old,v_old
       !real, allocatable, dimension(:,:) :: C1,G1,D1,C2,G2,D2,p_RHS

       real, dimension(Nxb+2,Nyb+2) :: u_old,v_old
       real, dimension(Nxb,Nyb) :: C1,G1,D1,C2,G2,D2,p_RHS
       real :: u_res1, v_res1, maxdiv, mindiv, umax, umin, vmax, vmin
       integer :: i,j
       real, intent(inout), dimension(:,:) :: u,v,visc,rho1x,rho1y,rho2x,rho2y,p,s,s2,sigp,sigx,sigy,ut,vt
       real, dimension(Nxb+2,Nyb+2) :: rhox,rhoy

       ins_v_res = 0
       ins_u_res = 0

       v_res1 = 0
       u_res1 = 0

       ! Rescaling velocities after first step
       !if(tstep>1) call ins_rescaleVel(u,v)

       u_old = u
       v_old = v

       ! Predictor Step

       call Convective_U_VD(u,v,gr_dx,gr_dy,C1)
       call Diffusive_U_VD(u,visc,rho1x,rho2x,gr_dx,gr_dy,ins_inRe,D1)

       G1 = C1 + D1 + ins_gravX

       if (tstep == 0) then

!              ut(2:Nxb+1,2:Nyb+1)=u(2:Nxb+1,2:Nyb+1)+(dr_dt/1)*(G1)
!              ins_G1_old = G1
       else

!              ut(2:Nxb+1,2:Nyb+1)=u(2:Nxb+1,2:Nyb+1)+(dr_dt/2)*(3*ins_G1_old-G1)
!              ins_G1_old = G1
       endif


       call Convective_V_VD(u,v,gr_dx,gr_dy,C2)
       call Diffusive_V_VD(v,visc,rho1y,rho2y,gr_dx,gr_dy,ins_inRe,D2)

       G2 = C2 + D2 + ins_gravY

       if (tstep == 0) then

!              vt(2:Nxb+1,2:Nyb+1)=v(2:Nxb+1,2:Nyb+1)+(dr_dt/1)*(G2)
!              ins_G2_old = G2
       else

!              vt(2:Nxb+1,2:Nyb+1)=v(2:Nxb+1,2:Nyb+1)+(dr_dt/2)*(3*ins_G2_old-G2)
!              ins_G2_old = G2
       endif

       ! Boundary Conditions

#ifdef MPI_DIST
       !call MPI_applyBC(ut)
       !call MPI_applyBC(vt)
#endif

#ifdef MPI_SHRD
       !call MPI_BARRIER(shared_comm,ierr)
       !call MPI_applyBC_shared(USTR_VAR,FACEX)
       !call MPI_applyBC_shared(USTR_VAR,FACEY)
#endif 

#ifdef MPI_RMA
       !call MPI_applyBC_RMA(ut)
       !call MPI_applyBC_RMA(vt)
#endif

!       call MPI_physicalBC_vel(ut,vt)

#ifdef IBM
       ! Immersed Boundary - Predictor BC

       call IBM_ApplyForcing(ut,vt,s,s2)

#endif
       ! Poisson Solver

       p_RHS = -((1/(gr_dy*dr_dt))*(vt(2:Nxb+1,2:Nyb+1)-vt(2:Nxb+1,1:Nyb)))&
               -((1/(gr_dx*dr_dt))*(ut(2:Nxb+1,2:Nyb+1)-ut(1:Nxb,2:Nyb+1)))&
               -sigp(2:Nxb+1,2:Nyb+1)

       rhox = 1./(rho1x+rho2x)
       rhoy = 1./(rho1y+rho2y)

       call Poisson_solver_VC(p_RHS,p,rhox,rhoy,ins_p_res,p_counter,PRES_VAR)

       ! Corrector Step

       u(2:Nxb+1,2:Nyb+1) = ut(2:Nxb+1,2:Nyb+1) - (dr_dt/gr_dx)*(rho1x(2:Nxb+1,2:Nyb+1)+rho2x(2:Nxb+1,2:Nyb+1))*&
                                                  (p(3:Nxb+2,2:Nyb+1)-p(2:Nxb+1,2:Nyb+1))&
                                                  +dr_dt*sigx(2:Nxb+1,2:Nyb+1)

       v(2:Nxb+1,2:Nyb+1) = vt(2:Nxb+1,2:Nyb+1) - (dr_dt/gr_dy)*(rho1y(2:Nxb+1,2:Nyb+1)+rho2y(2:Nxb+1,2:Nyb+1))*&
                                                  (p(2:Nxb+1,3:Nyb+2)-p(2:Nxb+1,2:Nyb+1))&
                                                  +dr_dt*sigy(2:Nxb+1,2:Nyb+1)

       ! Boundary Conditions

#ifdef MPI_DIST
       !call MPI_applyBC(u)
       !call MPI_applyBC(v)
#endif

#ifdef MPI_SHRD
       !call MPI_BARRIER(shared_comm,ierr)
       !call MPI_applyBC_shared(VELC_VAR,FACEX)
       !call MPI_applyBC_shared(VELC_VAR,FACEY)
#endif 

#ifdef MPI_RMA
       !call MPI_applyBC_RMA(u)
       !call MPI_applyBC_RMA(v)
#endif

 !      call MPI_physicalBC_vel(u,v)

       ! Divergence

       maxdiv = -10.**(10.)
       mindiv = 10.**(10.)

       maxdiv = max(maxdiv,maxval(((1/(gr_dy))*(v(2:Nxb+1,2:Nyb+1)-v(2:Nxb+1,1:Nyb)))&
                                  +((1/(gr_dx))*(u(2:Nxb+1,2:Nyb+1)-u(1:Nxb,2:Nyb+1)))))

       mindiv = min(mindiv,minval(((1/(gr_dy))*(v(2:Nxb+1,2:Nyb+1)-v(2:Nxb+1,1:Nyb)))&
                                  +((1/(gr_dx))*(u(2:Nxb+1,2:Nyb+1)-u(1:Nxb,2:Nyb+1)))))

       umax = maxval(u)
       umin = minval(u)

       vmax = maxval(v)
       vmin = minval(v)

       call MPI_CollectResiduals(maxdiv,ins_maxdiv,MAX_DATA)
       call MPI_CollectResiduals(mindiv,ins_mindiv,MIN_DATA)

       call MPI_CollectResiduals(umax,ins_umaxmin(1),MAX_DATA)
       call MPI_CollectResiduals(umin,ins_umaxmin(2),MIN_DATA)

       call MPI_CollectResiduals(vmax,ins_vmaxmin(1),MAX_DATA)
       call MPI_CollectResiduals(vmin,ins_vmaxmin(2),MIN_DATA)

       ! Residuals

       do i=1,Nyb+2
          ins_u_res = ins_u_res + sum((u(:,i)-u_old(:,i))**2)
       enddo

       call MPI_CollectResiduals(ins_u_res,u_res1,SUM_DATA)
       ins_u_res = sqrt(u_res1/((nblockx*nblocky)*(Nxb+2)*(Nyb+2)))

       do i=1,Nyb+1
          ins_v_res = ins_v_res + sum((v(:,i)-v_old(:,i))**2)
       enddo

       call MPI_CollectResiduals(ins_v_res,v_res1,SUM_DATA)
       ins_v_res = sqrt(v_res1/((nblockx*nblocky)*(Nxb+2)*(Nyb+2)))


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

