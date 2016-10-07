module IBM_interface


   interface
      subroutine IBM_init()
      end subroutine IBM_init
   end interface

   interface
      subroutine IBM_ApplyForcing(ut,vt,s)
      implicit none
      real,dimension(:,:),intent(inout) :: ut,vt,s
      end subroutine
   end interface

end module IBM_interface
