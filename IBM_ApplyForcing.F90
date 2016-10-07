subroutine IBM_ApplyForcing(ut,vt,s)

#include "Solver.h"

    implicit none

    real, dimension(:,:),intent(inout) :: ut,vt,s

    integer :: i,j

    
    do j=2,Nyb+1
       do i=2,Nxb+1

          if(s(i,j) >= 0. .and. s(i-1,j) < 0.) then

             if(abs(s(i,j))/(abs(s(i,j))+abs(s(i-1,j))) .le. 0.5) then
             ut(i,j) = 0.0

             else
             ut(i-1,j) = 0.0

             end if
          end if

          if(s(i,j) >= 0. .and. s(i+1,j) < 0.) then 

             if(abs(s(i,j))/(abs(s(i,j))+abs(s(i+1,j))) .le. 0.5) then
             ut(i,j) = 0.0

             else
             ut(i+1,j) = 0.0

             end if
          end if

          if(s(i,j) >= 0. .and. s(i,j-1) < 0.) then

             if(abs(s(i,j))/(abs(s(i,j))+abs(s(i,j-1))) .le. 0.5) then
             vt(i,j) = 0.0
            
             else
             vt(i,j-1) = 0.0
       
             end if
          end if

          if(s(i,j) >= 0. .and. s(i,j+1) < 0.) then

             if(abs(s(i,j))/(abs(s(i,j))+abs(s(i,j+1))) .le. 0.5) then
             vt(i,j) = 0.0
            
             else
             vt(i,j+1) = 0.0
       
             end if
          end if

       end do
    end do

end subroutine IBM_ApplyForcing
