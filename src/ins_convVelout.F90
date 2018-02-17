subroutine ins_convVelout(u)

#include "Solver.h"

    use IncompNS_data, only: ins_convvel,u_old
    use MPI_interface, only: MPI_CollectResiduals
    use Grid_data
    use MPI_data

    implicit none

    real,dimension(:,:,:),intent(in) :: u
    real :: convvel
    real :: avg_OP
    integer :: blk,j
    real :: localCounter,globalCounter

    convvel = 0.0
    avg_OP  = gr_dy/(gr_Ly*nblocky)
    localCounter = 0
    globalCounter = 0

    do blk=1,blockCount
       if ( xLC(blk) == nblockx-1) then
          do j=2,Nyb+1
             convvel = convvel + u(Nxb+1,j,blk)*avg_OP
             u_old(1,j,blk) = u(Nxb  ,j,blk)
             u_old(2,j,blk) = u(Nxb+1,j,blk)
             localCounter = localCounter + 1
          end do
       end if
   end do 

  
   call MPI_CollectResiduals(convvel,ins_convvel,SUM_DATA)
   call MPI_CollectResiduals(localCounter,globalCounter,SUM_DATA)
   
end subroutine ins_convVelout
