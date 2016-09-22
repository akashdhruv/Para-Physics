subroutine MPI_physicalBC_vel(u_ex,v_ex)

#include "Solver.h"

       use MPI_data
       use Driver_data, only: dr_dt
       use Grid_data, only: gr_dx
       use MPI_interface, only: MPI_periodicBC

       implicit none

       include "mpif.h"

       real, dimension(Nxb+2,Nyb+2), intent(inout) :: u_ex, v_ex
       integer :: status(MPI_STATUS_SIZE)


#ifdef LID_DRIVEN_FLOW
       
       if ( x_id == 0) then

           v_ex(1,:)=-v_ex(2,:)
           u_ex(1,:)=0

       end if

       if ( x_id == nblockx-1) then

           v_ex(Nxb+2,:)=-v_ex(Nxb+1,:)
           u_ex(Nxb+1,:)=0
           u_ex(Nxb+2,:)=0

       end if


       if ( y_id == 0) then

           v_ex(:,1)=0
           u_ex(:,1)=-u_ex(:,2)

       end if

       if ( y_id == nblocky-1) then

           v_ex(:,Nyb+2)=0
           v_ex(:,Nyb+1)=0
           u_ex(:,Nyb+2)=2-u_ex(:,Nyb+1)

       end if

#endif

#ifdef CHANNEL_FLOW

       if ( x_id == 0) then

           v_ex(1,:)=v_ex(2,:)
           u_ex(1,:)=1.0

       end if

       if ( x_id == nblockx-1) then

           v_ex(Nxb+2,:)=v_ex(Nxb+1,:)
           u_ex(Nxb+1,:)=u_ex(Nxb+1,:)-(dr_dt/gr_dx)*(u_ex(Nxb+1,:)-u_ex(Nxb,:))
           u_ex(Nxb+2,:)=u_ex(Nxb+1,:)

       end if


       if ( y_id == 0) then

           v_ex(:,1)=0
           u_ex(:,1)=-u_ex(:,2)

       end if

       if ( y_id == nblocky-1) then

           v_ex(:,Nyb+2)=0
           v_ex(:,Nyb+1)=0
           u_ex(:,Nyb+2)=-u_ex(:,Nyb+1)

       end if

#endif

#ifdef COUETTE_FLOW

       call MPI_periodicBC(u_ex,v_ex,1)

       if ( x_id == nblockx-1) then

           u_ex(Nxb+2,:)=u_ex(Nxb+1,:)

       end if


       if ( y_id == 0) then

           v_ex(:,1)=0
           u_ex(:,1)=-2-u_ex(:,2)

       end if

       if ( y_id == nblocky-1) then

           v_ex(:,Nyb+2)=0
           v_ex(:,Nyb+1)=0
           u_ex(:,Nyb+2)=2-u_ex(:,Nyb+1)

       end if

#endif

       call MPI_BARRIER(solver_comm,ierr)

end subroutine MPI_physicalBC_vel
