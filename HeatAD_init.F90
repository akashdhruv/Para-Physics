subroutine HeatAD_init()

#include "Solver.h"

   use HeatAD_data
   use IncompNS_data
   use physicaldata
   use Grid_data
   use Multiphase_data, only: mph_thco2,mph_cp2,mph_vis2,mph_rho2
   use IBM_data,        only: ibm_thco2,ibm_cp2,ibm_rho2,ibm_vis2
   use Driver_data, only: dr_dt
   use MPI_interface, only: MPI_applyBC, MPI_physicalBC_temp

   implicit none

   real,pointer,dimension(:,:,:) :: solnData
   integer :: j,i
   real :: solnX
   real :: solnY
   real :: ycell

   ht_Pr  = 7.2495
   ht_Nu  = 0.332*(ht_Pr**0.33)/(ins_inRe**0.5)
   ht_src = 0.0
   ht_St  = 1.0

   solnData => ph_center

   solnData(TEMP_VAR,:,:) = 0.0
   solnData(TOLD_VAR,:,:) = 0.0

   nullify(solnData)

end subroutine HeatAD_init
