subroutine IBM_ApplyForcing(ut,vt,s,s2,omega1,omega2)

#include "Solver.h"

    use Grid_data

    implicit none

    real, dimension(:,:),intent(inout) :: ut,vt,s,s2
    real, intent(in) :: omega1,omega2

    real :: xcell,ycell,x0,y0,r0,x1,y1,r1,theta

    integer :: i,j

    x0 = -0.40
    x1 = -0.10

    y0 = 0.0
    y1 = 0.0

    r0 = 0.15
    r1 = 0.10
    
    do j=2,Nyb+1

       ycell = 0.5*(gr_y(1,j) + gr_y(1,j-1))

       do i=2,Nxb+1

          xcell = 0.5*(gr_x(i,1) + gr_x(i-1,1)) 

          theta = atan2((ycell-y0),(xcell-x0))

          if(s(i,j) >= 0. .and. s(i-1,j) < 0.) then

             if(abs(s(i,j))/(abs(s(i,j))+abs(s(i-1,j))) .le. 0.5) then
             ut(i,j) = omega1*cos(theta)

             else
             ut(i-1,j) = omega1*cos(theta)

             end if
          end if

          if(s(i,j) >= 0. .and. s(i+1,j) < 0.) then 

             if(abs(s(i,j))/(abs(s(i,j))+abs(s(i+1,j))) .le. 0.5) then
             ut(i,j) = omega1*cos(theta)

             else
             ut(i+1,j) = omega1*cos(theta)

             end if
          end if

          if(s(i,j) >= 0. .and. s(i,j-1) < 0.) then

             if(abs(s(i,j))/(abs(s(i,j))+abs(s(i,j-1))) .le. 0.5) then
             vt(i,j) = omega1*sin(theta)
            
             else
             vt(i,j-1) = omega1*sin(theta)
       
             end if
          end if

          if(s(i,j) >= 0. .and. s(i,j+1) < 0.) then

             if(abs(s(i,j))/(abs(s(i,j))+abs(s(i,j+1))) .le. 0.5) then
             vt(i,j) = omega1*sin(theta)
            
             else
             vt(i,j+1) = omega1*sin(theta)
       
             end if
          end if

#if NBOD == 2

          theta = atan2((ycell-y1),(xcell-x1))

          if(s2(i,j) >= 0. .and. s2(i-1,j) < 0.) then

             if(abs(s2(i,j))/(abs(s2(i,j))+abs(s2(i-1,j))) .le. 0.5) then
             ut(i,j) = omega2*cos(theta)

             else
             ut(i-1,j) = omega2*cos(theta)

             end if
          end if

          if(s2(i,j) >= 0. .and. s2(i+1,j) < 0.) then 

             if(abs(s2(i,j))/(abs(s2(i,j))+abs(s2(i+1,j))) .le. 0.5) then
             ut(i,j) = omega2*cos(theta)

             else
             ut(i+1,j) = omega2*cos(theta)

             end if
          end if

          if(s2(i,j) >= 0. .and. s2(i,j-1) < 0.) then

             if(abs(s2(i,j))/(abs(s2(i,j))+abs(s2(i,j-1))) .le. 0.5) then
             vt(i,j) = omega2*sin(theta)
            
             else
             vt(i,j-1) = omega2*sin(theta)
       
             end if
          end if

          if(s2(i,j) >= 0. .and. s2(i,j+1) < 0.) then

             if(abs(s2(i,j))/(abs(s2(i,j))+abs(s2(i,j+1))) .le. 0.5) then
             vt(i,j) = omega2*sin(theta)
            
             else
             vt(i,j+1) = omega2*sin(theta)
       
             end if
          end if
#endif

       end do
    end do

end subroutine IBM_ApplyForcing
