subroutine MPI_physicalBC_vort(d_ex)

#include "Solver.h"

       use MPI_data
       use physicaldata
       use Grid_data

       implicit none

       include "mpif.h"

       real, dimension(:,:), intent(inout) :: d_ex
       integer :: status(MPI_STATUS_SIZE)
       logical :: mask

       real, pointer, dimension(:,:) :: u,v

       u => ph_facex(VELC_VAR,:,:)
       v => ph_facey(VELC_VAR,:,:)

       mask = .true.
    
       if ( x_id == 0) then

           d_ex(1,:)=d_ex(:,2)

       end if

       if ( x_id == nblockx-1) then

           d_ex(Nxb+2,:)=d_ex(Nxb+1,:)

       end if


       if ( y_id == 0) then

           d_ex(:,1)=d_ex(:,2)

       end if

       if ( y_id == nblocky-1) then

           d_ex(:,Nyb+2)=d_ex(:,Nyb+1)

       end if

       call MPI_BARRIER(solver_comm,ierr)
   
       mask = .false.

       nullify(u,v)

end subroutine MPI_physicalBC_vort
