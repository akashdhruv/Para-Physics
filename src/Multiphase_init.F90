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

   integer :: i,j,blk

   real,pointer,dimension(:,:,:,:) :: solnData,facexData,faceyData

   solnData  => localCENTER
   facexData => localFACEX
   faceyData => localFACEY  

   mph_thco2 = 1.0           ! Fluid thermal conductivity
   mph_rho2  = 1.0           ! Fluid density
   mph_cp2   = 1.0*mph_rho2  ! Fluid specific heat
   mph_vis2  = 1.0           ! Fluid viscosity

   mph_thco1 = 1.0           ! Vapor/Gas thermal conductivity
   mph_rho1  = 1.0           ! Vapor/Gas density
   mph_cp1   = 1.0*mph_rho1  ! Vapor/Gas specific heat
   mph_vis1  = 1.0           ! Vapor/Gas viscosity

   mph_beta  = 1.0           ! Thermal expansion coefficient
   mph_sten  = 1.0           ! inverse Weber number 

   solnData(:,:,DFUN_VAR,:)  = 0.0
   solnData(:,:,PFUN_VAR,:)  = 0.0
   solnData(:,:,NRMX_VAR,:)  = 0.0
   solnData(:,:,NRMY_VAR,:)  = 0.0
   solnData(:,:,VISC_VAR,:)  = 1.0
   solnData(:,:,THCO_VAR,:)  = 1.0
   solnData(:,:,CPRS_VAR,:)  = 1.0
   solnData(:,:,SMHV_VAR,:)  = 0.0
   solnData(:,:,SMRH_VAR,:)  = 0.0
   solnData(:,:,SIGP_VAR,:)  = 0.0
   solnData(:,:,CURV_VAR,:)  = 0.0
   solnData(:,:,MDOT_VAR,:)  = 0.0

   facexData(:,:,RH1F_VAR,:) = 0.0
   faceyData(:,:,RH1F_VAR,:) = 0.0
   facexData(:,:,RH2F_VAR,:) = 1.0
   faceyData(:,:,RH2F_VAR,:) = 1.0
   facexData(:,:,AL1F_VAR,:) = 0.0
   faceyData(:,:,AL1F_VAR,:) = 0.0
   facexData(:,:,AL2F_VAR,:) = 1.0
   faceyData(:,:,AL2F_VAR,:) = 1.0
   facexData(:,:,SIGM_VAR,:) = 0.0
   faceyData(:,:,SIGM_VAR,:) = 0.0
   facexData(:,:,VELI_VAR,:) = 0.0
   faceyData(:,:,VELI_VAR,:) = 0.0


   ! Initial vapor phase definition

   !___Vorticity Test__!
   !mph_x0 = -0.2
   !mph_y0 =  0.0
   !mph_r0 =  0.05

   !___Multiphase Test__!
   !mph_x0 = 0.0
   !mph_y0 = 0.0
   !mph_r0 = 0.5

   !___Conjugate Heat Test_!
   mph_x0 =  3.0
   mph_y0 =  0.0
   mph_r0 =  0.5

   ! Distance function calculation

   !do blk=1,blockCount

   !do j=1,Nyb+2

   !   if(j==1) then     
   !       ycell = gr_y(1,j,blk) - 0.5*gr_dy

   !   else if(j==Nyb+2) then
   !       ycell = gr_y(1,Nyb+1,blk) + 0.5*gr_dy

   !   else
   !       ycell = 0.5*(gr_y(1,j,blk) + gr_y(1,j-1,blk))

   !   end if

   !   do i=1,Nxb+2

   !      if(i==1) then     
   !           xcell = gr_x(i,1,blk) - 0.5*gr_dx

   !      else if(i==Nxb+2) then
   !           xcell = gr_x(Nxb+1,1,blk) + 0.5*gr_dx

   !      else
   !           xcell = 0.5*(gr_x(i,1,blk) + gr_x(i-1,1,blk))

   !      end if

   !           solnData(i,j,DFUN_VAR,blk) = mph_r0 - sqrt((xcell-mph_x0)**2+(ycell-mph_y0)**2)

   !   end do

   ! end do

   ! end do

    !mph_redistance

    nullify(solnData,facexData,faceyData)

end subroutine Multiphase_init
