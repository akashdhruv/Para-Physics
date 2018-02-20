subroutine ins_convVelout(u,v,x,y)

#include "Solver.h"

    use IncompNS_data
    use MPI_interface, only: MPI_CollectResiduals
    use Grid_data
    use MPI_data

    implicit none

    real,dimension(:,:,:),intent(in) :: u,v,x,y
    real :: convvel(2,2)
    real :: avg_OPy, avg_OPx
    integer :: blk,j

    convvel = 0.0

    avg_OPy  = gr_dy/(gr_Ly*nblocky)
    avg_OPx  = gr_dx/(gr_Lx*nblockx)

#ifdef HOME_HEATING_SYSTEM

    avg_OPy  = gr_dy/(ins_dnEx2-ins_dnEx1)
    avg_OPx  = gr_dx/(ins_upEx2-ins_upEx1)

    do blk=1,blockCount

       if ( xLC(blk) == nblockx-1) then
          do j=2,Nyb+1
             if(y(Nxb+1,j,blk) .ge. ins_dnEx1 .and. y(Nxb+1,j,blk) .le. ins_dnEx2) then
                convvel(2,1) = convvel(2,1) + u(Nxb+1,j,blk)*avg_OPy
                u_old(3,j,blk) = u(Nxb  ,j,blk)
                u_old(4,j,blk) = u(Nxb+1,j,blk)
             end if
          end do
       end if

       if ( yLC(blk) == nblocky-1) then
          do j=1,Nxb+2
             if(x(j,Nyb+1,blk) .ge. ins_upEx1 .and. x(j,Nyb+1,blk) .le. ins_upEx2) then
                convvel(2,2) = convvel(2,2) + v(j,Nyb+1,blk)*avg_OPx
                v_old(j,3,blk) = v(j,Nyb  ,blk)
                v_old(j,4,blk) = v(j,Nyb+1,blk)
             end if
          end do
       end if

    end do
#else
    do blk=1,blockCount

       if ( xLC(blk) == 0) then
          do j=2,Nyb+1
             convvel(1,1) = convvel(1,1) + u(1,j,blk)*avg_OPy
             u_old(1,j,blk) = u(1,j,blk)
             u_old(2,j,blk) = u(2,j,blk)
          end do
       end if

       if ( xLC(blk) == nblockx-1) then
          do j=2,Nyb+1
             convvel(2,1) = convvel(2,1) + u(Nxb+1,j,blk)*avg_OPy
             u_old(3,j,blk) = u(Nxb  ,j,blk)
             u_old(4,j,blk) = u(Nxb+1,j,blk)
          end do
       end if

       if ( yLC(blk) == 0) then
          do j=2,Nxb+1
             convvel(1,2) = convvel(1,2) + v(j,1,blk)*avg_OPx
             v_old(j,1,blk) = v(j,1,blk)
             v_old(j,2,blk) = v(j,2,blk)
          end do
       end if

       if ( yLC(blk) == nblocky-1) then
          do j=2,Nxb+1
             convvel(2,2) = convvel(2,2) + v(j,Nyb+1,blk)*avg_OPx
             v_old(j,3,blk) = v(j,Nyb  ,blk)
             v_old(j,4,blk) = v(j,Nyb+1,blk)
          end do
       end if

   end do 
#endif

   call MPI_ALLREDUCE(convvel,ins_convvel,4,MPI_REAL,MPI_SUM,solver_comm,ierr)
 
end subroutine ins_convVelout
