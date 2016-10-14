subroutine ibm_evolve(x0,y0,r0,x1,y1,r1)

#include "Solver.h"

   use physicaldata
   use Grid_data

   implicit none
  
   real, intent(in) :: x0,y0,r0,x1,y1,r1
   real, pointer, dimension(:,:,:) :: facexData,faceyData
   real :: ycell,xcell
   integer :: i,j

   facexData => ph_facex
   faceyData => ph_facey

   nullify(facexData)
   nullify(faceyData)
   
end subroutine
