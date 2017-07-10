subroutine Grid_init()
           
     use Grid_data
     use physicaldata
     use MPI_data

#include "Solver.h"

     implicit none

     integer :: I,blk
     real :: pi=4.0*atan(1.0)

     gr_Lx = (D_xmax)-(D_xmin) ! Domain length in X DIR
     gr_Ly = (D_ymax)-(D_ymin) ! Domain length in Y DIR

     gr_Lx = gr_Lx/nblockx     ! Block length in X DIR
     gr_Ly = gr_Ly/nblocky     ! Block length in Y DIR

     gr_dx = gr_Lx/Nxb         ! dx
     gr_dy = gr_Ly/Nyb         ! dy
    
     allocate(gr_x(Nxb+1,Nyb+1,blockCount)) ! Array to store X location per block
     allocate(gr_y(Nxb+1,Nyb+1,blockCount)) ! Array to store Y location per block

     ! Calculate X and Y coordinates per block
     do blk=1,blockCount

       do i=1,Nyb+1
        gr_x(:,i,blk)=D_xmin+xLC(blk)*gr_Lx+gr_dx*(/(I,I=0,Nxb)/)
       enddo

       do i=1,Nxb+1
        gr_y(i,:,blk)=D_ymin+yLC(blk)*gr_Ly+gr_dy*(/(I,I=0,Nyb)/)
       enddo

    end do

end subroutine Grid_init
