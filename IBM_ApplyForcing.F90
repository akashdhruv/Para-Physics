subroutine IBM_ApplyForcing(ut,vt,s,s2)

#include "Solver.h"

    use Grid_data
    use IBM_data

    implicit none

    real, dimension(:,:),intent(inout) :: ut,vt,s,s2

    real :: xcell_up,xcell_vp,xcell_um,xcell_vm,xcell_u,xcell_v
    real :: ycell_up,ycell_vp,ycell_um,ycell_vm,ycell_u,ycell_v
    real :: theta_uxp,theta_vxp,theta_uxm,theta_vxm
    real :: theta_uyp,theta_vyp,theta_uym,theta_vym

    integer :: i,j


    do j=2,Nyb+1

      ! if(j==2) then
      ! ycell_um = gr_y(1,j-1)-0.5*gr_dy
      ! ycell_up = (gr_y(1,j)+gr_y(1,j+1))*0.5
      ! ycell_vm = gr_y(1,j-1)
      ! ycell_vp = gr_y(1,j+1)

      ! else if(j==Nyb+1) then
      ! ycell_um = (gr_y(1,j-1)+gr_y(1,j-2))*0.5
      ! ycell_up = gr_y(1,j)+0.5*gr_dy
      ! ycell_vm = gr_y(1,j-1)
      ! ycell_vp = gr_y(1,j) + gr_dy

      ! else
      ! ycell_um = (gr_y(1,j-1)+gr_y(1,j-2))*0.5
      ! ycell_up = (gr_y(1,j)+gr_y(1,j+1))*0.5
      ! ycell_vm = gr_y(1,j-1)
      ! ycell_vp = gr_y(1,j+1)

      ! end if

      ! ycell_u = (gr_y(1,j)+gr_y(1,j-1))*0.5
      ! ycell_v = gr_y(1,j)

       do i=2,Nxb+1


         ! if(i==2) then
         ! xcell_um = gr_x(i-1,1)
         ! xcell_up = gr_x(i+1,1)
         ! xcell_vm = gr_x(i-1,1)-0.5*gr_dx
         ! xcell_vp = (gr_x(i,1)+gr_x(i+1,1))*0.5

         ! else if(j==Nxb+1) then
         ! xcell_um = gr_x(i-1,1)
         ! xcell_up = gr_x(i,1) + gr_dx
         ! xcell_vm = (gr_x(i-1,1)+gr_x(i-2,1))*0.5
         ! xcell_vp = gr_x(i,1)+0.5*gr_dx

         ! else
         ! xcell_um = gr_x(i-1,1)
         ! xcell_up = gr_x(i+1,1)
         ! xcell_vm = (gr_x(i-1,1)+gr_x(i-2,1))*0.5
         ! xcell_vp = (gr_x(i,1)+gr_x(i+1,1))*0.5

         ! end if

         ! xcell_u = gr_x(1,i)
         ! xcell_v = 0.5*(gr_x(1,i)+gr_x(1,i-1))

         ! theta_uxm = atan2(ycell_u-ibm_y0,xcell_um-ibm_x0)
         ! theta_uxp = atan2(ycell_u-ibm_y0,xcell_up-ibm_x0)
         ! theta_uym = atan2(ycell_um-ibm_y0,xcell_u-ibm_x0)
         ! theta_uyp = atan2(ycell_up-ibm_y0,xcell_u-ibm_x0)

         ! theta_vxm = atan2(ycell_v-ibm_y0,xcell_vm-ibm_x0)
         ! theta_vxp = atan2(ycell_v-ibm_y0,xcell_vp-ibm_x0)
         ! theta_vym = atan2(ycell_vm-ibm_y0,xcell_v-ibm_x0)
         ! theta_vyp = atan2(ycell_vp-ibm_y0,xcell_v-ibm_x0)

          if(s(i,j) > 0. .and. s(i-1,j) <= 0.) ut(i-1,j) = 0.0

          if(s(i,j) > 0. .and. s(i+1,j) <= 0.) ut(i+1,j) = 0.0

          if(s(i,j) > 0. .and. s(i,j-1) <= 0.) ut(i,j-1) = 0.0

          if(s(i,j) > 0. .and. s(i,j+1) <= 0.) ut(i,j+1) = 0.0

          if(s2(i,j) > 0. .and. s2(i-1,j) <= 0.) vt(i-1,j) = 0.0

          if(s2(i,j) > 0. .and. s2(i+1,j) <= 0.) vt(i+1,j) = 0.0

          if(s2(i,j) > 0. .and. s2(i,j-1) <= 0.) vt(i,j-1) = 0.0

          if(s2(i,j) > 0. .and. s2(i,j+1) <= 0.) vt(i,j+1) = 0.0

       end do
    end do

end subroutine IBM_ApplyForcing
