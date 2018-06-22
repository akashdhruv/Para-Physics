subroutine MPI_physicalBC_pres(p_ex,x,y)

#include "Solver.h"

       use MPI_data
       use Grid_data
       use IncompNS_data

       implicit none

       real, dimension(:,:,:), intent(inout) :: p_ex
       real, dimension(:,:,:), intent(in) :: x,y
       integer :: status(MPI_STATUS_SIZE),blk,i,j
       logical :: mask

       mask = .true.
    
       do blk=1,blockCount

       if (xLC(blk) == 0) p_ex(1,:,blk)=p_ex(2,:,blk)
       if (xLC(blk) == nblockx-1) p_ex(Nxb+2,:,blk)=p_ex(Nxb+1,:,blk)

       if (yLC(blk) == 0) p_ex(:,1,blk)=p_ex(:,2,blk)
       if (yLC(blk) == nblocky-1) p_ex(:,Nyb+2,blk)=p_ex(:,Nyb+1,blk)

       end do   

       mask = .false.

end subroutine MPI_physicalBC_pres
