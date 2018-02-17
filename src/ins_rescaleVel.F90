subroutine ins_rescaleVel(u)

#include "Solver.h"

    use IncompNS_data, only: ins_Qin,ins_Qout,ins_Qinout
    use Grid_data
    use MPI_data

    implicit none

    real,dimension(:,:,:),intent(inout) :: u
   
    real :: tol=1e-10
    integer :: blk

    if(ins_Qout == 0) then 
      ins_Qinout = 1.0
    else
      ins_Qinout = ins_Qin/ins_Qout
    end if
  
    do blk=1,blockCount
       if ( xLC(blk) == nblockx-1) then
          u(Nxb+1,2:Nyb+1,blk) = u(Nxb+1,2:Nyb+1,blk)*ins_Qinout
       end if
    end do 

end subroutine
