module morton_interface

   interface 
     subroutine morton_sort(blockCount,myid,procs,blockID,blockLC)
     implicit none
     integer,dimension(:),intent(inout) :: blockID,blockLC
     integer,intent(in) :: blockCount
     integer,intent(in) :: myid
     integer,intent(in) :: procs
     end subroutine
   end interface

end module morton_interface

