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
   pf => ph_center(PFUN_VAR,:,:)
   thco => ph_center(THCO_VAR,:,:)
   cprs => ph_center(CPRS_VAR,:,:)

   s = 0.0
   pf = 0.0
   thco = 0.0
   cprs = 0.0

   x0 = 0.0
   y0 =  0.0
   r = 0.2

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

   ibm_rho1 = 0.597
   ibm_rho2 = 958.4

   ibm_thco1 = 0.025
   ibm_thco2 = 0.679

   ibm_cp1 = 2030.0*ibm_rho1
   ibm_cp2 = 4216.0*ibm_rho2

   call mph_FillVars(s,pf,thco,cprs,ibm_thco1,ibm_thco2,ibm_cp1,ibm_cp2)

   nullify(s)
   nullify(pf)
   nullify(thco)
   nullify(cprs)

end subroutine IBM_init
