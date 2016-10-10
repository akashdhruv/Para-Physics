module IBM_interface


   interface
      subroutine IBM_init()
      end subroutine IBM_init
   end interface

   interface
      subroutine IBM_ApplyForcing(ut,vt,s,s2,omega1,omega2)
      implicit none
      real,dimension(:,:),intent(inout) :: ut,vt,s,s2
      real, intent(in) :: omega1,omega2
      end subroutine
   end interface


   interface
      subroutine IBM_solver(tstep)
      implicit none
      integer, intent(in) :: tstep
      end subroutine
   end interface

   interface 
      subroutine ibm_evolve(x0,y0,r0,x1,y1,r1)
      implicit none
      real, intent(in) :: x0,y0,r0,x1,y1,r1
      end subroutine
   end interface

end module IBM_interface
