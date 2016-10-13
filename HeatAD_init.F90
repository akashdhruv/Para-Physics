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

   real,pointer,dimension(:,:) :: T,s
   integer :: j,i
   real :: solnX
   real :: solnY
   real :: ycell

   ht_Pr = 0.7
   ht_Nu = 0.332*(ht_Pr**0.33)/(ins_inRe**0.5)
   ht_src = 100.0

   T => ph_center(TEMP_VAR,:,:)
   s => ph_center(DFUN_VAR,:,:)

   T = 313.0

   !do j=1,Nyb+2
   !  do i=1,Nxb+2

   !     if(s(i,j) .ge. 0.0) then

   !          T(i,j) = 400.00

   !      end if

   !  end do
   !end do

   nullify(T)
   nullify(s)

end subroutine HeatAD_init
