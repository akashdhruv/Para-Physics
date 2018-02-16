subroutine MPI_physicalBC_vel(u_ex,v_ex,y)

#include "Solver.h"

       use MPI_data
       use Driver_data, only: dr_dt
       use Grid_data, only: gr_dx,gr_dy
       use MPI_interface, only: MPI_periodicBC
       use IBM_data, only: ibm_yr

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

           v_ex(1,:,blk)=-v_ex(2,:,blk)

           do i=1,Nyb+2
            if(y(1,i,blk) .le. ibm_yr) then
               u_ex(1,i,blk)=0.0
            else
               u_ex(1,i,blk)=1.0
            end if
           end do

       end if

       if ( xLC(blk) == nblockx-1) then

           v_ex(Nxb+2,:,blk) = v_ex(Nxb+1,:,blk)
           u_ex(Nxb+2,:,blk) = u_ex(Nxb+1,:,blk)

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
