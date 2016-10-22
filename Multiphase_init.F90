subroutine Multiphase_init()

   use Multiphase_data
   use physicaldata
   use Grid_data
   use MPI_interface , only: MPI_CollectResiduals

#include "Solver.h"

   implicit none

   real,pointer,dimension(:,:,:) :: solnData,facexData,faceyData
   real :: xcell,ycell
   real :: min_s,max_s,all_min_s,all_max_s

   integer :: i,j

   solnData => ph_center
   facexData => ph_facex
   faceyData => ph_facey

   solnData(DFUN_VAR,:,:) = 0.0
   solnData(PFUN_VAR,:,:) = 0.0
   solnData(THCO_VAR,:,:) = 0.0
   solnData(CPRS_VAR,:,:) = 0.0
   solnData(VISC_VAR,:,:) = 1.0
   facexData(RHOF_VAR,:,:) = 1.0
   faceyData(ALPH_VAR,:,:) = 1.0

   !___Vorticity Test__!
   !mph_x0 = -0.2
   !mph_y0 =  0.0
   !mph_r0 =  0.05

   !___Multiphase Test__!
   mph_x0 = 0.0
   mph_y0 = 0.0
   mph_r0 = 0.125

   !___Conjugate Heat Test_!
   !mph_x0 = -5.0
   !mph_y0 =  0.0
   !mph_r0 = 0.5

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

              solnData(DFUN_VAR,i,j) = mph_r0 - sqrt((xcell-mph_x0)**2+(ycell-mph_y0)**2)

      end do

    end do

   !mph_rho1 = 0.597
   !mph_rho2 = 958.4

   !mph_thco1 = 0.025
   !mph_thco2 = 0.679

   !mph_cp1 = 2030.0*mph_rho1
   !mph_cp2 = 4216.0*mph_rho2

   mph_thco2 = 1.0
   mph_rho2 = 1.0
   mph_cp2 = 1.0*mph_rho2
   mph_vis2 = 1.0

   mph_thco1 = 1.2
   mph_rho1 = 1.0
   mph_cp1 = 1.0*mph_rho1
   mph_vis1 = 1.0
   mph_beta = 0.0002

   nullify(solnData)
   nullify(facexData)
   nullify(faceyData)

end subroutine Multiphase_init
