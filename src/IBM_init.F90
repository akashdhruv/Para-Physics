subroutine IBM_init()

#include "Solver.h"

   use physicaldata, only: localFACEX,localFACEY,localCENTER
   use Grid_data
   use IBM_data
   use Multiphase_interface, only: mph_FillVars
   use MPI_data, only: blockCount,solver_comm,ierr
   use MPI_interface, ONLY: MPI_applyBC, MPI_CollectResiduals, MPI_physicalBC_vort

   implicit none

   real :: x0,y0,r,xcell,ycell
   real :: x1,y1,r1

   integer :: i,j,blk

   real,pointer,dimension(:,:,:,:) :: solnData,facexData,faceyData

   facexData => localFACEX
   faceyData => localFACEY  
   solnData  => localCENTER

   ibm_thco2 = 1.0          ! Fluid thermal conductivity
   ibm_rho2 = 1.0           ! Fluid density
   ibm_cp2 = 1.0*ibm_rho2   ! Fluid specific heat
   ibm_vis2 = 1.0           ! Fluid viscosity

   ibm_thco1 = 1.0          ! Solid thermal conductivity
   ibm_rho1 = 1.0           ! Solid density
   ibm_cp1 = 1.0*ibm_rho1   ! Solid specific heat
   ibm_vis1 = 1.0           ! Solid viscosity

   ibm_omega = 1.0          ! Rotational velocity of immersed solid cylinder

   ! Distance function to store distances from immersed boundary
   facexData(:,:,IBMF_VAR,:) = 0.0 
   faceyData(:,:,IBMF_VAR,:) = 0.0
   solnData(:,:,DFUN_VAR,:)  = 0.0
   solnData(:,:,PFUN_VAR,:)  = 0.0 

   ! parameters for cylinder's location
   ibm_x0 = 3.0
   ibm_y0 = 0.0
   ibm_r0 = 0.5

   ! parameters for square block

#ifdef BACKWARD_FACING_STEP
   ! backward facing step
   ibm_xl = -5.0
   ibm_xr =  2.5
   ibm_yl = -2.5
   ibm_yr = -1.0
#endif

#ifdef FORWARD_FACING_STEP
   ! forward facing step
   ibm_xl =  7.5
   ibm_xr = 15.0
   ibm_yl = -2.5
   ibm_yr = -1.0
#endif

   ! Calculate distance function
   do blk=1,blockCount

   do j=1,Nyb+2

      if(j==1) then     
          ycell = gr_y(1,j,blk) - 0.5*gr_dy

      else if(j==Nyb+2) then
          ycell = gr_y(1,Nyb+1,blk) + 0.5*gr_dy

      else
          ycell = 0.5*(gr_y(1,j,blk) + gr_y(1,j-1,blk))

      end if

      do i=1,Nxb+2

         if(i==1) then     
              xcell = gr_x(i,1,blk) - 0.5*gr_dx

         else if(i==Nxb+2) then
              xcell = gr_x(Nxb+1,1,blk) + 0.5*gr_dx

         else
              xcell = 0.5*(gr_x(i,1,blk) + gr_x(i-1,1,blk))

         end if

              solnData(i,j,DFUN_VAR,blk) = xcell
              solnData(i,j,PFUN_VAR,blk) = ycell
                
      end do
    end do 
   end do

   call MPI_BARRIER(solver_comm,ierr)

   nullify(facexData,faceyData,solnData)

end subroutine IBM_init
