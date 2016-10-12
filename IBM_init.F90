subroutine IBM_init()

#include "Solver.h"

   use physicaldata
   use Grid_data
   use IBM_data
   use Multiphase_interface, only: mph_FillVars

   implicit none

   real,pointer,dimension(:,:) :: s,pf,thco,cprs,s2
   real :: x0,y0,r,xcell,ycell
   real :: x1,y1,r1

   integer :: i,j


   s => ph_facex(IBMF_VAR,:,:)
   s2 => ph_facey(IBMF_VAR,:,:)

   s = 0.0
   s2 = 0.0

   ibm_x0 = -0.4
   !ibm_x0 = 0.0
   ibm_y0 =  0.0
   ibm_r0 =  0.15

   do j=1,Nyb+2

    if(j==1) then
        ycell = gr_y(1,j) - 0.5*gr_dy

    else if(j==Nyb+2) then
        ycell = gr_y(1,Nyb+1) + 0.5*gr_dy

    else
        ycell = 0.5*(gr_y(1,j) + gr_y(1,j-1))

    end if

    do i=1,Nxb+2

     if(i==Nxb+2) then
          xcell = gr_x(Nxb+1,1) + gr_dx

     else
          xcell = gr_x(i,1)

     end if

     s(i,j) = ibm_r0 - sqrt((xcell-ibm_x0)**2+(ycell-ibm_y0)**2)

     end do

   end do

   do j=1,Nyb+2

    if(j==Nyb+2) then
        ycell = gr_y(1,Nyb+1) + gr_dy

    else
        ycell = gr_y(1,j)

    end if

    do i=1,Nxb+2

    if(i==1) then
          xcell = gr_x(i,1) - 0.5*gr_dx

    else if(i==Nxb+2) then
          xcell = gr_x(Nxb+1,1) + 0.5*gr_dx

    else
          xcell = 0.5*(gr_x(i,1) + gr_x(i-1,1))

    end if

    s2(i,j) = ibm_r0 - sqrt((xcell-ibm_x0)**2+(ycell-ibm_y0)**2)

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

   ibm_omega = 1.0

   nullify(s)
   nullify(s2)

end subroutine IBM_init
