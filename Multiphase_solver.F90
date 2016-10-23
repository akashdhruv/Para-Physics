subroutine Multiphase_solver(tstep,solnX)

#include "Solver.h"

    use Driver_data,  only: dr_dt
    use Multiphase_data, only: mph_thco1, mph_cp1, mph_thco2, mph_cp2,mph_beta
    use physicaldata
    use Grid_data
    use Multiphase_interface , only: mph_FillVars
    use MPI_interface, only: MPI_applyBC,MPI_physicalBC_dfun

    implicit none

    integer, intent(in) :: tstep
    

    real,intent(out) :: solnX
    real :: ycell
    real, pointer, dimension(:,:,:) :: solnData,facexData,faceyData
    integer :: j,i

    solnData => ph_center
    facexData => ph_facex
    faceyData => ph_facey

    call mph_FillVars(solnData(DFUN_VAR,:,:),solnData(PFUN_VAR,:,:),&
                      solnData(THCO_VAR,:,:),solnData(CPRS_VAR,:,:),&
                      solnData(VISC_VAR,:,:),&
                      facexData(RHOF_VAR,:,:),faceyData(RHOF_VAR,:,:),&
                      facexData(ALPH_VAR,:,:),faceyData(ALPH_VAR,:,:),&
                      solnData(TEMP_VAR,:,:),solnData(TOLD_VAR,:,:),&
                      mph_beta)


    nullify(solnData)  
    nullify(facexData)
    nullify(faceyData)
 
end subroutine Multiphase_solver

