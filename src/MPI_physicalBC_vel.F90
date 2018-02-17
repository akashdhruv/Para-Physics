subroutine MPI_physicalBC_vel(u_ex,v_ex,y)

#include "Solver.h"

       use MPI_data
       use Driver_data, only: dr_dt
       use Grid_data, only: gr_dx,gr_dy
       use MPI_interface, only: MPI_periodicBC
       use IBM_data, only: ibm_yr
       use IncompNS_data, only: u_old, ins_convvel

       implicit none

       real, dimension(:,:,:), intent(inout) :: u_ex, v_ex
       real, dimension(:,:,:), intent(out) :: y
       integer :: status(MPI_STATUS_SIZE),blk
       integer :: i

#ifdef LID_DRIVEN_FLOW      
   do blk=1,blockCount
       if ( xLC(blk) == 0) then

           v_ex(1,:,blk)=-v_ex(2,:,blk)
           u_ex(1,:,blk)=0

       end if

       if ( xLC(blk) == nblockx-1) then

           v_ex(Nxb+2,:,blk)=-v_ex(Nxb+1,:,blk)
           u_ex(Nxb+1,:,blk)=0
           u_ex(Nxb+2,:,blk)=0

       end if


       if ( yLC(blk) == 0) then

           v_ex(:,1,blk)=0
           u_ex(:,1,blk)=-u_ex(:,2,blk)

       end if

       if ( yLC(blk) == nblocky-1) then

           v_ex(:,Nyb+2,blk)=0
           v_ex(:,Nyb+1,blk)=0
           u_ex(:,Nyb+2,blk)=2-u_ex(:,Nyb+1,blk)

       end if
   end do
#endif

#ifdef CHANNEL_FLOW
   do blk=1,blockCount
       if ( xLC(blk) == 0) then

           v_ex(1,2:Nyb+1,blk)=-v_ex(2,2:Nyb+1,blk)

#ifdef BACKWARD_FACING_STEP
           ! backward facing step
           do i=2,Nyb+1
            if(y(1,i,blk) .le. ibm_yr) then
               u_ex(1,i,blk)=0.0
            else
               u_ex(1,i,blk)=1.0
            end if
           end do
#endif

#ifdef FORWARD_FACING_STEP
           ! forward facing step
           u_ex(1,2:Nyb+1,blk) = 1.0
#endif

       end if

       if ( xLC(blk) == nblockx-1) then

#ifdef BACKWARD_FACING_STEP
           ! backward facing step
           u_ex(Nxb+1,2:Nyb+1,blk) = u_old(2,2:Nyb+1,blk) - &
                                     ins_convvel*(dr_dt/gr_dx)* &
                                    (u_old(2,2:Nyb+1,blk)- &
                                     u_old(1,2:Nyb+1,blk))

#endif

#ifdef FORWARD_FACING_STEP
           ! forward facing step
           do i=2,Nyb+1
              if(y(Nxb+1,i,blk) .le. ibm_yr) then
                 u_ex(Nxb+1,i,blk) = 0.0
              else
                 u_ex(Nxb+1,i,blk) = u_old(2,i,blk) - &
                                     ins_convvel*(dr_dt/gr_dx)* &
                                    (u_old(2,i,blk)- &
                                     u_old(1,i,blk))
              end if
           end do
#endif

           v_ex(Nxb+2,2:Nyb+1,blk) = v_ex(Nxb+1,2:Nyb+1,blk)
           u_ex(Nxb+2,2:Nyb+1,blk) = u_ex(Nxb+1,2:Nyb+1,blk)

       end if


       if ( yLC(blk) == 0) then

           v_ex(:,1,blk)=0.0
           u_ex(:,1,blk)=-u_ex(:,2,blk)

       end if

       if ( yLC(blk) == nblocky-1) then

           v_ex(:,Nyb+2,blk)=-v_ex(:,Nyb,blk)
           v_ex(:,Nyb+1,blk)=0.0
           u_ex(:,Nyb+2,blk)=u_ex(:,Nyb+1,blk)

       end if
   end do
#endif

end subroutine MPI_physicalBC_vel
