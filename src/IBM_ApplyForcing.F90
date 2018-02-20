subroutine IBM_ApplyForcing(ut,vt,s,s2)

#include "Solver.h"

    use Grid_data
    use IBM_data
    use IncompNS_data

    implicit none

    real, dimension(:,:),intent(inout) :: ut,vt,s,s2

    real :: xcell_up,xcell_vp,xcell_um,xcell_vm,xcell_u,xcell_v
    real :: ycell_up,ycell_vp,ycell_um,ycell_vm,ycell_u,ycell_v
    real :: theta_uxp,theta_vxp,theta_uxm,theta_vxm
    real :: theta_uyp,theta_vyp,theta_uym,theta_vym

    integer :: i,j


    do j=1,Nyb+1
      do i=1,Nxb+1

#ifdef BACKWARD_FACING_STEP
        if(s(i,j) .le. ibm_xr .and. s(i+1,j) .ge. ibm_xr .and. s2(i,j) .le. ibm_yr) then
           ut(i+1,j) = 0.0
        end if

        if(s2(i,j) .le. ibm_yr .and. s2(i,j+1) .ge. ibm_yr .and. s(i,j) .le. ibm_xr) then
           vt(i,j+1) = 0.0
        end if
#endif

#ifdef FORWARD_FACING_STEP
        if(s(i,j) .le. ibm_xl .and. s(i+1,j) .ge. ibm_xl .and. s2(i,j) .le. ibm_yr) then
           ut(i+1,j) = 0.0
        end if

        if(s2(i,j) .le. ibm_yr .and. s2(i,j+1) .ge. ibm_yr .and. s(i,j) .ge. ibm_xl) then
           vt(i,j+1) = 0.0
        end if
#endif

#ifdef HOME_HEATING_SYSTEM

        if(s(i,j) .le. ibm_xl .and. s(i+1,j) .ge. ibm_xl .and. &
           s2(i,j) .le. ibm_yr .and. s2(i,j) .ge. ibm_yl) then
           ut(i+1,j) = 0.0
           vt(i+1,j) = -vt(i,j)
        end if

        if(s2(i,j) .le. ibm_yr .and. s2(i,j+1) .ge. ibm_yr .and. &
           s(i,j) .le. ibm_xr .and. s(i,j) .ge. ibm_xl) then
           vt(i,j+1) = 0.0
           ut(i,j+1) = -ut(i,j)
        end if

        if(s2(i,j) .le. ibm_yl .and. s2(i,j+1) .ge. ibm_yl .and. &
           s(i,j) .le. ins_dnIn1 .and. s(i,j) .ge. ibm_xl) then
           vt(i,j+1) = 0.0
           ut(i,j+1) = -ut(i,j)
        end if

        if(s2(i,j) .le. ibm_yl .and. s2(i,j+1) .ge. ibm_yl .and. &
           s(i,j) .le. ins_dnIn2 .and. s(i,j) .ge. ins_dnIn1) then
           vt(i,j+1) = -1.0
           ut(i,j+1) = -ut(i,j)          
        end if

#endif

      end do
    end do

    !do j=2,Nyb+1
    !   do i=2,Nxb+1

    !      if(s(i,j) >= 0. .and. s(i-1,j) <0.) ut(i-1,j) = 0.0

    !      if(s(i,j) >= 0. .and. s(i+1,j) <0.) ut(i+1,j) = 0.0

    !      if(s(i,j) >= 0. .and. s(i,j-1) <0.) ut(i,j-1) = 0.0

    !      if(s(i,j) >= 0. .and. s(i,j+1) <0.) ut(i,j+1) = 0.0

    !      if(s2(i,j) >= 0. .and. s2(i-1,j) <0.) vt(i-1,j) = 0.0

    !      if(s2(i,j) >= 0. .and. s2(i+1,j) <0.) vt(i+1,j) = 0.0

    !      if(s2(i,j) >= 0. .and. s2(i,j-1) <0.) vt(i,j-1) = 0.0

    !      if(s2(i,j) >= 0. .and. s2(i,j+1) <0.) vt(i,j+1) = 0.0

    !   end do
    !end do

end subroutine IBM_ApplyForcing
