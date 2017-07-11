subroutine MPI_physicalBC_temp(t_ex)

#include "Solver.h"

       use MPI_data

       implicit none

       real, dimension(:,:,:), intent(inout) :: t_ex
       integer :: status(MPI_STATUS_SIZE),blk

       do blk=1,blockCount

       if (xLC(blk) == 0) t_ex(1,:,blk) = t_ex(2,:,blk)
       if (xLC(blk) == nblockx-1) t_ex(Nxb+2,:,blk) = t_ex(Nxb+1,:,blk)

       !if (yLC(blk) == 0) t_ex(:,1,blk) = 2.0-t_ex(:,2,blk)
       !if (yLC(blk) == nblocky-1) t_ex(:,Nyb+2,blk) = 2.0-t_ex(:,Nyb+1,blk)

       if (yLC(blk) == 0) t_ex(:,1,blk) = t_ex(:,2,blk)
       if (yLC(blk) == nblocky-1) t_ex(:,Nyb+2,blk) = t_ex(:,Nyb+1,blk)

       end do 

end subroutine MPI_physicalBC_temp
