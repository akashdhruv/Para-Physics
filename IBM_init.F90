subroutine IBM_init()

#include "Solver.h"

   use physicaldata
   use Grid_data
   use IBM_data
   use Multiphase_interface, only: mph_FillVars

   implicit none

   real,pointer,dimension(:,:) :: s,pf,thco,cprs
   real :: x0,y0,r,xcell,ycell

   integer :: i,j


   s => ph_center(DFUN_VAR,:,:)

   s = 0.0

   x0 = -0.75
   y0 =  0.0
   r  =  0.1

   do j=1,Nyb+2

    if(j==1) then
        ycell = gr_y(1,j) - 0.5*gr_dy

    else if(j==Nyb+2) then
        ycell = gr_y(1,Nyb+1) + 0.5*gr_dy

    else
        ycell = 0.5*(gr_y(1,j) + gr_y(1,j-1))

    end if

    do i=1,Nxb+2

     if(i==1) then
          xcell = gr_x(i,1) - 0.5*gr_dx

     else if(i==Nxb+2) then
          xcell = gr_x(Nxb+1,1) + 0.5*gr_dx

     else
          xcell = 0.5*(gr_x(i,1) + gr_x(i-1,1))

     end if

     s(i,j) = r - sqrt((xcell-x0)**2+(ycell-y0)**2)

     end do

   end do

   ibm_thco2 = 1.0
   ibm_rho2 = 1.0
   ibm_cp2 = 1.0*ibm_rho2
   ibm_vis2 = 1.0

   ibm_thco1 = 1.2
   ibm_rho1 = 1.0
   ibm_cp1 = 1.0*ibm_rho1
   ibm_vis1 = 1.0

   nullify(s)

end subroutine IBM_init
