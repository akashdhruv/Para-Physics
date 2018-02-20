subroutine MPI_physicalBC_temp(t_ex,x,y)

#include "Solver.h"

       use MPI_data
       use Grid_data
       use IncompNS_data

       implicit none

       real, dimension(:,:,:), intent(inout) :: t_ex
       real, dimension(:,:,:), intent(in) :: x,y
       integer :: status(MPI_STATUS_SIZE),blk,i,j

       do blk=1,blockCount

       if (xLC(blk) == 0) t_ex(1,:,blk) = t_ex(2,:,blk)
       if (xLC(blk) == nblockx-1) t_ex(Nxb+2,:,blk) = t_ex(Nxb+1,:,blk)

       if (yLC(blk) == 0) t_ex(:,1,blk) = t_ex(:,2,blk)
       if (yLC(blk) == nblocky-1) t_ex(:,Nyb+2,blk) = t_ex(:,Nyb+1,blk)

       end do 

#ifdef HOME_HEATING_SYSTEM
       do blk=1,blockCount
          if(yLC(blk) == nblocky-1) then 
             do i=1,Nxb+2
                if(x(i,Nyb+2,blk) .ge. ins_upIn1 .and. x(i,Nyb+2,blk) .le. ins_upIn2) t_ex(i,Nyb+2,blk) = 2.0 - t_ex(i,Nyb+1,blk)
             end do
          end if
       end do
#endif

end subroutine MPI_physicalBC_temp
