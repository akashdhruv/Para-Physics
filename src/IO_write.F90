subroutine IO_write(x,y,uu,vv,pp,tt,ww,rr,blk,blockOffset)

#include "Solver.h"

       implicit none
       real, dimension(Nxb+1, Nyb+1), intent(in) :: x,y,uu,vv,pp,tt,ww,rr
       integer, intent(in) :: blk,blockOffset
       character(len=10) :: f1,f2,f3,f4,f7,f8,f9,f0
       integer :: id

       id = blk+blockOffset-1
     
       write (f1, '( "X", I4.4, ".dat" )' )id
       write (f2, '( "Y", I4.4, ".dat" )' )id
       write (f3, '( "U", I4.4, ".dat" )' )id
       write (f4, '( "V", I4.4, ".dat" )' )id
       write (f7, '( "P", I4.4, ".dat" )' )id
       write (f8, '( "T", I4.4, ".dat" )' )id
       write (f9, '( "W", I4.4, ".dat" )' )id
       write (f0, '( "R", I4.4, ".dat" )' )id

       open(unit = 1, file = f1)
       open(unit = 2, file = f2)
       open(unit = 3, file = f3)
       open(unit = 4, file = f4)
       open(unit = 7, file = f7)
       open(unit = 8, file = f8)
       open(unit = 9, file = f9)
       open(unit = 0, file = f0)

       write(1,*),x
       write(2,*),y
       write(3,*),uu
       write(4,*),vv
       write(7,*),pp
       write(8,*),tt
       write(9,*),ww
       write(0,*),rr

       close(1)
       close(2)
       close(3)
       close(4)
       close(7)
       close(8)
       close(9)
       close(0)

end subroutine IO_write
