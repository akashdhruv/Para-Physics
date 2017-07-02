subroutine Multiphase_init()

   use Multiphase_data
   use physicaldata, only: localCENTER,localFACEX,localFACEY
   use Grid_data
   use MPI_interface , only: MPI_CollectResiduals
   use MPI_data, only: blockCount

#include "Solver.h"

   implicit none

   real :: xcell,ycell
   real :: min_s,max_s,all_min_s,all_max_s

   integer :: i,j

   real,pointer,dimension(:,:,:,:) :: solnData,facexData,faceyData

   solnData  => localCENTER
   facexData => localFACEX
   faceyData => localFACEY  

   mph_thco2 = 1.0
   mph_rho2  = 1.0
   mph_cp2   = 1.0*mph_rho2
   mph_vis2  = 1.0

   mph_thco1 = 1.0
   mph_rho1  = 1.0
   mph_cp1   = 1.0*mph_rho1
   mph_vis1  = 1.0

   mph_beta  = 0.1
   mph_sten  = 1.0

   solnData(:,:,blockCount,DFUN_VAR)  = 0.0
   solnData(:,:,blockCount,PFUN_VAR)  = 0.0
   solnData(:,:,blockCount,NRMX_VAR)  = 0.0
   solnData(:,:,blockCount,NRMY_VAR)  = 0.0
   solnData(:,:,blockCount,VISC_VAR)  = 1.0
   solnData(:,:,blockCount,THCO_VAR)  = 1.0
   solnData(:,:,blockCount,CPRS_VAR)  = 1.0
   solnData(:,:,blockCount,SMHV_VAR)  = 0.0
   solnData(:,:,blockCount,SMRH_VAR)  = 0.0
   solnData(:,:,blockCount,SIGP_VAR)  = 0.0
   solnData(:,:,blockCount,CURV_VAR)  = 0.0
   solnData(:,:,blockCount,MDOT_VAR)  = 0.0

   facexData(:,:,blockCount,RH1F_VAR) = 0.0
   faceyData(:,:,blockCount,RH1F_VAR) = 0.0
   facexData(:,:,blockCount,RH2F_VAR) = 1.0
   faceyData(:,:,blockCount,RH2F_VAR) = 1.0
   facexData(:,:,blockCount,AL1F_VAR) = 0.0
   faceyData(:,:,blockCount,AL1F_VAR) = 0.0
   facexData(:,:,blockCount,AL2F_VAR) = 1.0
   faceyData(:,:,blockCount,AL2F_VAR) = 1.0
   facexData(:,:,blockCount,SIGM_VAR) = 0.0
   faceyData(:,:,blockCount,SIGM_VAR) = 0.0
   facexData(:,:,blockCount,VELI_VAR) = 0.0
   faceyData(:,:,blockCount,VELI_VAR) = 0.0


   !___Vorticity Test__!
   !mph_x0 = -0.2
   !mph_y0 =  0.0
   !mph_r0 =  0.05

   !___Multiphase Test__!
   mph_x0 = 0.0
   mph_y0 = 0.0
   mph_r0 = 0.5

   !___Conjugate Heat Test_!
   !mph_x0 = -5.0
   !mph_y0 =  0.0
   !mph_r0 = 0.5

   do j=1,Nyb+2

      if(j==1) then     
          ycell = gr_y(1,j,blockCount) - 0.5*gr_dy

      else if(j==Nyb+2) then
          ycell = gr_y(1,Nyb+1,blockCount) + 0.5*gr_dy

      else
          ycell = 0.5*(gr_y(1,j,blockCount) + gr_y(1,j-1,blockCount))

      end if

      do i=1,Nxb+2

         if(i==1) then     
              xcell = gr_x(i,1,blockCount) - 0.5*gr_dx

         else if(i==Nxb+2) then
              xcell = gr_x(Nxb+1,1,blockCount) + 0.5*gr_dx

         else
              xcell = 0.5*(gr_x(i,1,blockCount) + gr_x(i-1,1,blockCount))

         end if

              solnData(i,j,blockCount,DFUN_VAR) = mph_r0 - sqrt((xcell-mph_x0)**2+(ycell-mph_y0)**2)

      end do

    end do

    !mph_redistance

    nullify(solnData,facexData,faceyData)

end subroutine Multiphase_init
