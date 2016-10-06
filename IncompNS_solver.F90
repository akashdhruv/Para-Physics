subroutine IncompNS_solver(tstep,p_counter)

       use Poisson_interface, ONLY: Poisson_solver            
       use Grid_data
       use physicaldata
       use Driver_data
       use MPI_data
       use IncompNS_data
       use MPI_interface, ONLY: MPI_applyBC, MPI_CollectResiduals, MPI_physicalBC_vel
       use IncompNS_interface, ONLY: ins_rescaleVel

#include "Solver.h"

       implicit none

       integer, intent(in) :: tstep

       real, dimension(Nxb+2,Nyb+2) :: ut
       real, dimension(Nxb+2,Nyb+2) :: vt

       real, dimension(Nxb+2,Nyb+2) :: u_old
       real, dimension(Nxb+2,Nyb+2) :: v_old

       real, dimension(Nxb,Nyb)   :: C1
       real, dimension(Nxb,Nyb)   :: G1
       real, dimension(Nxb,Nyb)   :: D1

       real, dimension(Nxb,Nyb)   :: C2
       real, dimension(Nxb,Nyb)   :: G2
       real, dimension(Nxb,Nyb)   :: D2

       real :: u_res1, v_res1, maxdiv, mindiv

       real, dimension(Nxb,Nyb) :: p_RHS
       integer :: i,j
       integer, intent(out) :: p_counter

       real, pointer, dimension(:,:) :: u, v, p, s

       p => ph_center(PRES_VAR,:,:)
       u => ph_facex(VELC_VAR,:,:)
       v => ph_facey(VELC_VAR,:,:)
       s => ph_center(DFUN_VAR,:,:)

       ins_v_res = 0
       ins_u_res = 0

       v_res1 = 0
       u_res1 = 0

       ! Rescaling velocities after first step
       !if(tstep>1) call ins_rescaleVel(u,v)

       u_old = u
       v_old = v

       ! Predictor Step

       call Convective_U(u,v,gr_dx,gr_dy,C1)
       call Diffusive_U(u,gr_dx,gr_dy,ins_inRe,D1)

       G1 = C1 + D1

       if (tstep == 0) then

              ut(2:Nxb+1,2:Nyb+1)=u(2:Nxb+1,2:Nyb+1)+(dr_dt/1)*(G1)
              ins_G1_old = G1
       else

              ut(2:Nxb+1,2:Nyb+1)=u(2:Nxb+1,2:Nyb+1)+(dr_dt/2)*(3*ins_G1_old-G1)
              ins_G1_old = G1
       endif


       call Convective_V(u,v,gr_dx,gr_dy,C2)
       call Diffusive_V(v,gr_dx,gr_dy,ins_inRe,D2)

       G2 = C2 + D2

       if (tstep == 0) then

              vt(2:Nxb+1,2:Nyb+1)=v(2:Nxb+1,2:Nyb+1)+(dr_dt/1)*(G2)
              ins_G2_old = G2
       else

              vt(2:Nxb+1,2:Nyb+1)=v(2:Nxb+1,2:Nyb+1)+(dr_dt/2)*(3*ins_G2_old-G2)
              ins_G2_old = G2
       endif

       ! Boundary Conditions

       call MPI_applyBC(ut)
       call MPI_applyBC(vt)
       call MPI_physicalBC_vel(ut,vt)

       do j=1,Nyb+2
         do i=1,Nxb+2

            if(s(i,j) .ge. 0.) then
                ut(i,j) = 0.0
                vt(i,j) = 0.0
            end if

         end do
       end do

       ! Poisson Solver

       p_RHS = -((1/(gr_dy*dr_dt))*(vt(2:Nxb+1,2:Nyb+1)-vt(2:Nxb+1,1:Nyb)))&
               -((1/(gr_dx*dr_dt))*(ut(2:Nxb+1,2:Nyb+1)-ut(1:Nxb,2:Nyb+1)))

       call Poisson_solver(p_RHS,p,ins_p_res,p_counter,PRES_VAR)

       ! Corrector Step

       u(2:Nxb+1,2:Nyb+1) = ut(2:Nxb+1,2:Nyb+1) - (dr_dt/gr_dx)*(p(3:Nxb+2,2:Nyb+1)-p(2:Nxb+1,2:Nyb+1))
       v(2:Nxb+1,2:Nyb+1) = vt(2:Nxb+1,2:Nyb+1) - (dr_dt/gr_dy)*(p(2:Nxb+1,3:Nyb+2)-p(2:Nxb+1,2:Nyb+1))

       ! Boundary Conditions

       call MPI_applyBC(u)
       call MPI_applyBC(v)
       call MPI_physicalBC_vel(u,v)

       do j=1,Nyb+2
         do i=1,Nxb+2

            if(s(i,j) .ge. 0.) then
                u(i,j) = 0.0
                v(i,j) = 0.0
            end if

         end do
       end do

       ! Divergence

       maxdiv = -10.**(10.)
       mindiv = 10.**(10.)

       maxdiv = max(maxdiv,maxval(-((1/(gr_dy))*(v(2:Nxb+1,2:Nyb+1)-v(2:Nxb+1,1:Nyb)))&
                                  -((1/(gr_dx))*(u(2:Nxb+1,2:Nyb+1)-u(1:Nxb,2:Nyb+1)))))

       mindiv = min(mindiv,minval(-((1/(gr_dy))*(v(2:Nxb+1,2:Nyb+1)-v(2:Nxb+1,1:Nyb)))&
                                  -((1/(gr_dx))*(u(2:Nxb+1,2:Nyb+1)-u(1:Nxb,2:Nyb+1)))))


       call MPI_CollectResiduals(maxdiv,ins_maxdiv,3)
       call MPI_CollectResiduals(mindiv,ins_mindiv,2)

       ! Residuals

       do i=1,Nyb+2
          ins_u_res = ins_u_res + sum((u(:,i)-u_old(:,i))**2)
       enddo

       call MPI_CollectResiduals(ins_u_res,u_res1,1)
       ins_u_res = sqrt(u_res1/((nblockx*nblocky)*(Nxb+2)*(Nyb+2)))

       do i=1,Nyb+1
          ins_v_res = ins_v_res + sum((v(:,i)-v_old(:,i))**2)
       enddo

       call MPI_CollectResiduals(ins_v_res,v_res1,1)
       ins_v_res = sqrt(v_res1/((nblockx*nblocky)*(Nxb+2)*(Nyb+2)))


       nullify(u)
       nullify(v)
       nullify(p)
       nullify(s)

end subroutine IncompNS_solver


!! CONVECTIVE U !!
subroutine Convective_U(ut,vt,dx,dy,C1)

#include "Solver.h"
       
      implicit none

      real,dimension(Nxb+2,Nyb+2), intent(in) :: ut
      real,dimension(Nxb+2,Nyb+2), intent(in) :: vt

      real,intent(in) :: dx
      real,intent(in) :: dy

      real, dimension(Nxb,Nyb) :: ue
      real, dimension(Nxb,Nyb) :: uw
      real, dimension(Nxb,Nyb) :: us
      real, dimension(Nxb,Nyb) :: un
      real, dimension(Nxb,Nyb) :: vs
      real, dimension(Nxb,Nyb) :: vn
      real, dimension(Nxb,Nyb), intent(out) :: C1

      ue = (ut(2:Nxb+1,2:Nyb+1)+ut(3:Nxb+2,2:Nyb+1))/2
      uw = (ut(2:Nxb+1,2:Nyb+1)+ut(1:Nxb,2:Nyb+1))/2
      us = (ut(2:Nxb+1,2:Nyb+1)+ut(2:Nxb+1,1:Nyb))/2
      un = (ut(2:Nxb+1,2:Nyb+1)+ut(2:Nxb+1,3:Nyb+2))/2
      vs = (vt(2:Nxb+1,1:Nyb)+vt(3:Nxb+2,1:Nyb))/2
      vn = (vt(2:Nxb+1,2:Nyb+1)+vt(3:Nxb+2,2:Nyb+1))/2

      C1 = -((ue**2)-(uw**2))/dx - ((un*vn)-(us*vs))/dy

end subroutine Convective_U

!! CONVECTIVE V !!
subroutine Convective_V(ut,vt,dx,dy,C2)

#include "Solver.h"

      implicit none

      real,dimension(Nxb+2,Nyb+2), intent(in) :: ut
      real,dimension(Nxb+2,Nyb+2), intent(in) :: vt

      real, intent(in) :: dx
      real, intent(in) :: dy      
      real, dimension(Nxb,Nyb) :: vn, vs, ve, vw, ue, uw
      real, dimension(Nxb,Nyb), intent(out) :: C2

      vs = (vt(2:Nxb+1,2:Nyb+1)+vt(2:Nxb+1,1:Nyb))/2
      vn = (vt(2:Nxb+1,2:Nyb+1)+vt(2:Nxb+1,3:Nyb+2))/2
      ve = (vt(2:Nxb+1,2:Nyb+1)+vt(3:Nxb+2,2:Nyb+1))/2
      vw = (vt(2:Nxb+1,2:Nyb+1)+vt(1:Nxb,2:Nyb+1))/2
      ue = (ut(2:Nxb+1,2:Nyb+1)+ut(2:Nxb+1,3:Nyb+2))/2
      uw = (ut(1:Nxb,2:Nyb+1)+ut(1:Nxb,3:Nyb+2))/2

      C2 = -((ue*ve)-(uw*vw))/dx - ((vn**2)-(vs**2))/dy

end subroutine Convective_V

!! DIFFUSIVE U !!
subroutine Diffusive_U(ut,dx,dy,inRe,D1)

#include "Solver.h"

      implicit none

      real,dimension(Nxb+2,Nyb+2), intent(in) :: ut

      real, intent(in) :: dx
      real, intent(in) :: dy

      real, intent(in) :: inRe

      real, dimension(Nxb,Nyb) :: uP
      real, dimension(Nxb,Nyb) :: uN
      real, dimension(Nxb,Nyb) :: uS
      real, dimension(Nxb,Nyb) :: uE
      real, dimension(Nxb,Nyb) :: uW

      real, dimension(Nxb,Nyb), intent(out) :: D1

      uP = ut(2:Nxb+1,2:Nyb+1)
      uE = ut(3:Nxb+2,2:Nyb+1)
      uW = ut(1:Nxb,2:Nyb+1)
      uN = ut(2:Nxb+1,3:Nyb+2)
      uS = ut(2:Nxb+1,1:Nyb)

      !D1 = (inRe/dx)*(((uE-uP)/dx)-((uP-uW)/dx)) + (inRe/dy)*(((uN-uP)/dy)-((uP-uS)/dy))
      D1 = (inRe/dx)*((uE-uP)/dx)&
          -(inRe/dx)*((uP-uW)/dx)&
          +(inRe/dy)*((uN-uP)/dy)&
          -(inRe/dy)*((uP-uS)/dy)

end subroutine Diffusive_U

!! DIFFUSIVE V !!
subroutine Diffusive_V(vt,dx,dy,inRe,D2)

#include "Solver.h"

      implicit none

      real,dimension(Nxb+2,Nyb+2), intent(in) :: vt

      real, intent(in) :: dx
      real, intent(in) :: dy

      real, intent(in) :: inRe

      real, dimension(Nxb,Nyb) :: vP,vE,vW,vN,vS

      real, dimension(Nxb,Nyb), intent(out) :: D2

      vP = vt(2:Nxb+1,2:Nyb+1)
      vE = vt(3:Nxb+2,2:Nyb+1)
      vW = vt(1:Nxb,2:Nyb+1)
      vN = vt(2:Nxb+1,3:Nyb+2)
      vS = vt(2:Nxb+1,1:Nyb)

      !D2 = (inRe/dx)*(((vE-vP)/dx)-((vP-vW)/dx)) + (inRe/dy)*(((vN-vP)/dy)-((vP-vS)/dy))
      D2 = (inRe/dx)*((vE-vP)/dx)&
          -(inRe/dx)*((vP-vW)/dx)&
          +(inRe/dy)*((vN-vP)/dy)&
          -(inRe/dy)*((vP-vS)/dy)

end subroutine Diffusive_V

