subroutine Grid_finalize()

   use Grid_data

   implicit none
   
   deallocate(gr_x)
   deallocate(gr_y)

end subroutine Grid_finalize 
