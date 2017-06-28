subroutine HeatAD_init()

#include "Solver.h"

   use HeatAD_data
   use IncompNS_data
   use physicaldata, only: localCENTER

   implicit none

   real,pointer,dimension(:,:,:) :: solnData

   integer :: j,i
   real :: solnX
   real :: solnY
   real :: ycell

   solnData  => localCENTER

   ht_Pr  = 1.0
   ht_Nu  = 0.332*(ht_Pr**0.33)/(ins_inRe**0.5)
   ht_src = 0.0
   ht_St  = 1.0

   solnData(TEMP_VAR,:,:) = 0.0
   solnData(TOLD_VAR,:,:) = 0.0

   nullify(solnData)

end subroutine HeatAD_init
