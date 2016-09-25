module morton_interface

   interface 
     subroutine morton_sort(blockID,blockCount,myid,procs,blockLC)
     implicit none
     integer,dimension(:),intent(inout) :: blockID
     integer,dimension(:,:),intent(inout) :: blockLC
     integer,intent(in) :: blockCount
     integer,intent(in) :: myid
     integer,intent(in) :: procs
     end subroutine
   end interface

end module morton_interface

