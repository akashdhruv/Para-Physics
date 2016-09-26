subroutine Multiphase_init()

   use Multiphase_data
   use physicaldata
   use Grid_data

#include "Solver.h"

   implicit none

   real,pointer,dimension(:,:) :: sf,pf,th,cp
   real :: x0,y0,r,xcell,ycell

   integer :: i,j


   sf => ph_center(DFUN_VAR,:,:)
   pf => ph_center(PFUN_VAR,:,:)
   th => ph_center(THCO_VAR,:,:)
   cp => ph_center(CPRS_VAR,:,:)


   sf = 0.
   pf = 0.
   th = 0.
   cp = 0.

   !________________The Missing Data Simulation___________!

    x0 = 0.d0
    y0 = 0.d0
    r = 0.d0

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

          sf(i,j) = r - sqrt((xcell-x0)**2+(ycell-y0)**2)

      end do

    end do


   !___________________________End________________________!

   mph_rho1 = 0.597
   mph_rho2 = 958.4

   mph_thco1 = 0.025
   mph_thco2 = 0.679

   mph_cp1 = 2030.0*mph_rho1
   mph_cp2 = 4216.0*mph_rho2

   nullify(sf)
   nullify(pf)
   nullify(th)
   nullify(cp)

end subroutine Multiphase_init
