subroutine MPI_physicalBC_pres(p_ex)

#include "Solver.h"

       use MPI_data

       implicit none

       include "mpif.h"

       real, dimension(:,:), intent(inout) :: p_ex
       integer :: status(MPI_STATUS_SIZE)
       logical :: mask

       mask = .true.
    

#ifdef LID_DRIVEN_FLOW

       if ( x_id == 0) then

           p_ex(1,:)=p_ex(2,:)

       end if

       if ( x_id == nblockx-1) then

           p_ex(Nxb+2,:)=p_ex(Nxb+1,:)

       end if


       if ( y_id == 0) then

           p_ex(:,1)=p_ex(:,2)

       end if

       if ( y_id == nblocky-1) then

           p_ex(:,Nyb+2)=p_ex(:,Nyb+1)

       end if

       !call MPI_BARRIER(solver_comm,ierr)
   
       mask = .false.

#endif


#ifdef CHANNEL_FLOW

       if ( x_id == 0) then

           p_ex(1,:)=p_ex(2,:)

       end if

       if ( x_id == nblockx-1) then

           p_ex(Nxb+2,:)=-p_ex(Nxb+1,:)

       end if


       if ( y_id == 0) then

           p_ex(:,1)=p_ex(:,2)

       end if

       if ( y_id == nblocky-1) then

           p_ex(:,Nyb+2)=p_ex(:,Nyb+1)

       end if

       !call MPI_BARRIER(solver_comm,ierr)

       mask = .false.
   
#endif

#ifdef MPH_FLOW

       if ( x_id == 0) then

           p_ex(1,:)=-p_ex(2,:)

       end if

       if ( x_id == nblockx-1) then

           p_ex(Nxb+2,:)=-p_ex(Nxb+1,:)

       end if


       if ( y_id == 0) then

           p_ex(:,1)=-p_ex(:,2)

       end if

       if ( y_id == nblocky-1) then

           p_ex(:,Nyb+2)=-p_ex(:,Nyb+1)

       end if

       !call MPI_BARRIER(solver_comm,ierr)
 
       mask = .false.

#endif

      if (mask) then

       if ( x_id == 0) then

           p_ex(1,:)=p_ex(2,:)

       end if

       if ( x_id == nblockx-1) then

           p_ex(Nxb+2,:)=p_ex(Nxb+1,:)

       end if


       if ( y_id == 0) then

           p_ex(:,1)=p_ex(:,2)

       end if

       if ( y_id == nblocky-1) then

           p_ex(:,Nyb+2)=p_ex(:,Nyb+1)

       end if

       !call MPI_BARRIER(solver_comm,ierr)

       mask = .false.

      end if

end subroutine MPI_physicalBC_pres
