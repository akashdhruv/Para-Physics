subroutine Grid_init()
           
     use Grid_data
     use physicaldata
     use MPI_data

#include "Solver.h"

     implicit none

     integer :: I
     real :: pi=4.0*atan(1.0)

     gr_Lx = (D_xmax)-(D_xmin)
     gr_Ly = (D_ymax)-(D_ymin)

     gr_Lx = gr_Lx/nblockx
     gr_Ly = gr_Ly/nblocky

     gr_dx = gr_Lx/Nxb
     gr_dy = gr_Ly/Nyb
    
     allocate(gr_x(Nxb+1,Nyb+1,blockCount))
     allocate(gr_y(Nxb+1,Nyb+1,blockCount))

     do i=1,Nyb+1
        gr_x(:,i,blockCount)=D_xmin+mod(myid,nblockx)*gr_Lx+gr_dx*(/(I,I=0,Nxb)/)
     enddo

     do i=1,Nxb+1
        gr_y(i,:,blockCount)=D_ymin+(myid/nblockx)*gr_Ly+gr_dy*(/(I,I=0,Nyb)/)
    enddo

end subroutine Grid_init
