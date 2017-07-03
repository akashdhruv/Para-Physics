subroutine MPI_periodicBC(u_ex,v_ex,aDIM)

#include "Solver.h"

     use MPI_data

     implicit none
     
     real, dimension(:,:,:) , intent(inout) :: u_ex,v_ex
     integer, intent(in) :: aDIM

     integer :: status(MPI_STATUS_SIZE)


end subroutine MPI_periodicBC
