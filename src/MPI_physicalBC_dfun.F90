subroutine MPI_physicalBC_dfun(d_ex)

#include "Solver.h"

       use MPI_data

       implicit none

       real, dimension(:,:,:), intent(inout) :: d_ex
       integer :: status(MPI_STATUS_SIZE),blk
       logical :: mask

       mask = .true.

       do blk=1,blockCount
    
       if (xLC(blk) == 0) d_ex(1,:,blk)=d_ex(2,:,blk)
       if (xLC(blk) == nblockx-1) d_ex(Nxb+2,:,blk)=d_ex(Nxb+1,:,blk)

       if (yLC(blk) == 0) d_ex(:,1,blk)=d_ex(:,2,blk)
       if (yLC(blk) == nblocky-1) d_ex(:,Nyb+2,blk)=d_ex(:,Nyb+1,blk)

       end do

       mask = .false.

end subroutine MPI_physicalBC_dfun
