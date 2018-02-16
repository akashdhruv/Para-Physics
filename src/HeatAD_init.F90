subroutine HeatAD_init()

#include "Solver.h"

   use HeatAD_data
   use IncompNS_data
   use physicaldata, only: localCENTER
   use MPI_data, only: blockCount

   implicit none

   real,pointer,dimension(:,:,:,:) :: solnData

   integer :: j,i
   real :: solnX
   real :: solnY
   real :: ycell

   solnData  => localCENTER

   ht_Pr  =  1.0                                 ! Prandtl number
   ht_Nu  =  0.332*(ht_Pr**0.33)/(ins_inRe**0.5) ! Nusselt number
   ht_src =  0.0                                 ! Heat source
   ht_St  =  1.0                                 ! Stefan number

   ! Initialize variables relevant to Heat AD unit
   solnData(:,:,TEMP_VAR,:) = 0.0               
   solnData(:,:,TOLD_VAR,:) = 0.0

   nullify(solnData)

end subroutine HeatAD_init
