subroutine IBM_init()

#include "Solver.h"

   use physicaldata, only: localFACEX,localFACEY
   use Grid_data
   use IBM_data
   use Multiphase_interface, only: mph_FillVars
   use MPI_data, only: blockCount

   implicit none

   real :: x0,y0,r,xcell,ycell
   real :: x1,y1,r1

   integer :: i,j,blk

   real,pointer,dimension(:,:,:,:) :: facexData,faceyData

   facexData => localFACEX
   faceyData => localFACEY  

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

   ibm_x0 = 0.0             ! parameters for cylinder's location
   ibm_y0 = 0.5
   ibm_r0 = 0.1

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
     if(i==Nxb+2) then
          xcell = gr_x(Nxb+1,1,blk) + gr_dx
     else
          xcell = gr_x(i,1,blk)
     end if
     facexData(i,j,IBMF_VAR,blk) = ibm_r0 - sqrt((xcell-ibm_x0)**2+(ycell-ibm_y0)**2)
     end do
   end do

   do j=1,Nyb+2
    if(j==Nyb+2) then
        ycell = gr_y(1,Nyb+1,blk) + gr_dy
    else
        ycell = gr_y(1,j,blk)
    end if

    do i=1,Nxb+2
    if(i==1) then
          xcell = gr_x(i,1,blk) - 0.5*gr_dx
    else if(i==Nxb+2) then
          xcell = gr_x(Nxb+1,1,blk) + 0.5*gr_dx
    else
          xcell = 0.5*(gr_x(i,1,blk) + gr_x(i-1,1,blk))
    end if
    faceyData(i,j,IBMF_VAR,blk) = ibm_r0 - sqrt((xcell-ibm_x0)**2+(ycell-ibm_y0)**2)
    end do
   end do

  end do

   nullify(facexData,faceyData)

end subroutine IBM_init
