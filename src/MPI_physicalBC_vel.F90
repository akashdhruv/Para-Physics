subroutine MPI_physicalBC_vel(u_ex,v_ex,x,y)

#include "Solver.h"

       use MPI_data
       use Driver_data, only: dr_dt
       use Grid_data, only: gr_dx,gr_dy
       use MPI_interface, only: MPI_periodicBC
       use IBM_data, only: ibm_yr
       use IncompNS_data

       implicit none

       real, dimension(:,:,:), intent(inout) :: u_ex, v_ex
       real, dimension(:,:,:), intent(in) :: x,y
       integer :: status(MPI_STATUS_SIZE),blk
       integer :: i

   do blk=1,blockCount

       if ( xLC(blk) == 0) then

           if(ins_xl_bnd .eq. NOSLIP) then
                v_ex(1,2:Nyb+1,blk)=-v_ex(2,2:Nyb+1,blk)
                u_ex(1,2:Nyb+1,blk)=0

           else if(ins_xl_bnd .eq. SLIP) then
                v_ex(1,2:Nyb+1,blk)=v_ex(2,2:Nyb+1,blk)
                u_ex(1,2:Nyb+1,blk)=0
             
           else if(ins_xl_bnd .eq. INFLOW) then
                v_ex(1,2:Nyb+1,blk)=-v_ex(2,2:Nyb+1,blk)
                u_ex(1,2:Nyb+1,blk)=1.0

           end if
       end if

       if ( xLC(blk) == nblockx-1) then

           if(ins_xr_bnd .eq. NOSLIP) then
                v_ex(Nxb+2,2:Nyb+1,blk)=-v_ex(Nxb+1,2:Nyb+1,blk)
                u_ex(Nxb+1,2:Nyb+1,blk)=0
                u_ex(Nxb+2,2:Nyb+1,blk)=-u_ex(Nxb,2:Nyb+1,blk)

           else if(ins_xr_bnd .eq. SLIP) then
                v_ex(Nxb+2,2:Nyb+1,blk)=v_ex(Nxb+1,2:Nyb+1,blk)
                u_ex(Nxb+1,2:Nyb+1,blk)=0
                u_ex(Nxb+2,2:Nyb+1,blk)=-u_ex(Nxb,2:Nyb+1,blk)

           else if(ins_xr_bnd .eq. NEUMANN) then
                v_ex(Nxb+2,2:Nyb+1,blk)=v_ex(Nxb+1,2:Nyb+1,blk)
                u_ex(Nxb+2,2:Nyb+1,blk)=u_ex(Nxb,2:Nyb+1,blk)

           else if(ins_xr_bnd .eq. OUTFLOW) then
                u_ex(Nxb+1,2:Nyb+1,blk) = u_old(4,2:Nyb+1,blk) - &
                                          ins_convvel(2,1)*(dr_dt/gr_dx)* &
                                         (u_old(4,2:Nyb+1,blk)- &
                                          u_old(3,2:Nyb+1,blk))

                v_ex(Nxb+2,2:Nyb+1,blk)=v_ex(Nxb+1,2:Nyb+1,blk)
                u_ex(Nxb+2,2:Nyb+1,blk)=u_ex(Nxb,2:Nyb+1,blk)

           end if
       end if


       if ( yLC(blk) == 0) then

          if(ins_yl_bnd .eq. NOSLIP) then
                v_ex(1:Nxb+2,1,blk)=0
                u_ex(1:Nxb+2,1,blk)=-u_ex(1:Nxb+2,2,blk)

          else if(ins_yl_bnd .eq. SLIP) then
                v_ex(1:Nxb+2,1,blk)=0
                u_ex(1:Nxb+2,1,blk)=u_ex(1:Nxb+2,2,blk)
          
          end if
       end if

       if ( yLC(blk) == nblocky-1) then

          if(ins_yr_bnd .eq. NOSLIP) then
                v_ex(1:Nxb+2,Nyb+1,blk)=0
                v_ex(1:Nxb+2,Nyb+2,blk)=-v_ex(1:Nxb+2,Nyb,blk)
                u_ex(1:Nxb+2,Nyb+2,blk)=-u_ex(1:Nxb+2,Nyb+1,blk)

          else if(ins_yr_bnd .eq. SLIP) then
                v_ex(1:Nxb+2,Nyb+1,blk)=0
                v_ex(1:Nxb+2,Nyb+2,blk)=-v_ex(1:Nxb+2,Nyb,blk)
                u_ex(1:Nxb+2,Nyb+2,blk)=u_ex(1:Nxb+2,Nyb+1,blk)

          else if(ins_yr_bnd .eq. MOVLID) then
                v_ex(1:Nxb+2,Nyb+1,blk)=0
                v_ex(1:Nxb+2,Nyb+2,blk)=-v_ex(1:Nxb+2,Nyb,blk)
                u_ex(1:Nxb+2,Nyb+2,blk)=2-u_ex(1:Nxb+2,Nyb+1,blk)

          else if(ins_yr_bnd .eq. OUTFLOW) then
                v_ex(1:Nxb+2,Nyb+1,blk) = v_old(1:Nxb+2,4,blk) - &
                                          ins_convvel(2,2)*(dr_dt/gr_dy)* &
                                         (v_old(1:Nxb+2,4,blk)- &
                                          v_old(1:Nxb+2,3,blk))

                v_ex(1:Nxb+2,Nyb+2,blk)=v_ex(1:Nxb+2,Nyb,blk)
                u_ex(1:Nxb+2,Nyb+2,blk)=u_ex(1:Nxb+2,Nyb+1,blk)

          end if
       end if

   end do

#ifdef BACKWARD_FACING_STEP
           do blk=1,blockCount
           do i=2,Nyb+1
            if(xLC(blk) == 0 .and. y(1,i,blk) .le. ibm_yr) u_ex(1,i,blk)=0.0
           end do
           end do
#endif

#ifdef FORWARD_FACING_STEP
           do blk=1,blockCount
           do i=2,Nyb+1
              if(xLC(blk) == nblockx-1 .and. y(Nxb+1,i,blk) .le. ibm_yr) u_ex(Nxb+1,i,blk) = 0.0
           end do
           end do
#endif

end subroutine MPI_physicalBC_vel
