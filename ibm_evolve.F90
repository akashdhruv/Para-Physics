subroutine ibm_evolve(x0,y0,r0,x1,y1,r1)

#include "Solver.h"

   use physicaldata
   use Grid_data

   implicit none
  
   real, intent(in) :: x0,y0,r0,x1,y1,r1
   real, pointer, dimension(:,:) :: s,s2
   real :: ycell,xcell
   integer :: i,j

   s => ph_center(IBM1_VAR,:,:)
   s2 => ph_center(IBM2_VAR,:,:)

   do j=1,Nyb+2

    if(j==1) then
        ycell = gr_y(1,j) - 0.5*gr_dy

    else if(j==Nyb+2) then
        ycell = gr_y(1,Nyb+1) + 0.5*gr_dy

    else
        ycell = 0.5*(gr_y(1,j) + gr_y(1,j-1))

    end if

    do i=1,Nxb+2

     if(i==1) then
          xcell = gr_x(i,1) - 0.5*gr_dx

     else if(i==Nxb+2) then
          xcell = gr_x(Nxb+1,1) + 0.5*gr_dx

     else
          xcell = 0.5*(gr_x(i,1) + gr_x(i-1,1))

     end if

     s(i,j) = r0 - sqrt((xcell-x0)**2+(ycell-y0)**2)

#if NBOD == 2
     s2(i,j) = r1 - sqrt((xcell-x1)**2+(ycell-y1)**2)
#endif

     end do
   end do

   nullify(s)
   nullify(s2)
   
end subroutine
