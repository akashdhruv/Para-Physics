subroutine Grid_finalize()

   use Grid_data
   use physicaldata

   implicit none
   
   deallocate(gr_x)
   deallocate(gr_y)
   deallocate(ph_center)
   deallocate(ph_facex)
   deallocate(ph_facey)

end subroutine Grid_finalize 
