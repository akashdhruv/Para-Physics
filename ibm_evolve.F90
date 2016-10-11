subroutine ibm_evolve(x0,y0,r0,x1,y1,r1)

#include "Solver.h"

   use physicaldata
   use Grid_data

   implicit none
  
   real, intent(in) :: x0,y0,r0,x1,y1,r1
   real, pointer, dimension(:,:) :: s,s2
   real :: ycell,xcell
   integer :: i,j

   s => ph_facex(IBMF_VAR,:,:)
   s2 => ph_facey(IBMF_VAR,:,:)

   nullify(s)
   nullify(s2)
   
end subroutine
