subroutine Multiphase_solver(tstep,solnX,jump_flag)

#define MPH_DEBUG
#include "Solver.h"

    use Driver_data,  only: dr_dt
    use Multiphase_data, only: mph_thco1, mph_cp1, mph_thco2, mph_cp2,mph_beta
    use physicaldata
    use Grid_data
    use Multiphase_interface , only: mph_FillVars, mph_PressureJumps,&
                                     mph_getInterfaceVelocity,&
                                     mph_FillVars_ibm    

    use MPI_interface, only: MPI_applyBC,MPI_physicalBC_dfun

    implicit none

    integer, intent(in) :: tstep
    

    real,intent(out) :: solnX
    real :: ycell
    real, pointer, dimension(:,:,:) :: solnData,facexData,faceyData
    integer :: j,i
    logical,intent(in) :: jump_flag

    solnData => ph_center
    facexData => ph_facex
    faceyData => ph_facey


if (jump_flag .eqv. .FALSE.) then

#ifdef IBM

     call mph_FillVars_ibm(solnData(DFUN_VAR,:,:),solnData(PFUN_VAR,:,:),&
                           solnData(THCO_VAR,:,:),solnData(CPRS_VAR,:,:),&
                           solnData(VISC_VAR,:,:),&
                           facexData(RH2F_VAR,:,:),faceyData(RH2F_VAR,:,:),&
                           facexData(AL2F_VAR,:,:),faceyData(AL2F_VAR,:,:),&
                           solnData(TEMP_VAR,:,:),solnData(TOLD_VAR,:,:),&
                           mph_beta)

#else

#ifdef MPH_DEBUG

     call mph_FillVars_ibm(solnData(DFUN_VAR,:,:),solnData(PFUN_VAR,:,:),&
                           solnData(THCO_VAR,:,:),solnData(CPRS_VAR,:,:),&
                           solnData(VISC_VAR,:,:),&
                           facexData(RH2F_VAR,:,:),faceyData(RH2F_VAR,:,:),&
                           facexData(AL2F_VAR,:,:),faceyData(AL2F_VAR,:,:),&
                           solnData(TEMP_VAR,:,:),solnData(TOLD_VAR,:,:),&
                           mph_beta)


#else
    if(tstep > 0) then 

        call mph_getInterfaceVelocity(facexData(VELC_VAR,:,:),faceyData(VELC_VAR,:,:),&
                                      faceyData(VELI_VAR,:,:),faceyData(VELI_VAR,:,:),&
                                      solnData(SMRH_VAR,:,:),solnData(MDOT_VAR,:,:),&
                                      solnData(NRMX_VAR,:,:),solnData(NRMY_VAR,:,:))

    !   call mph_advect 
    !   call mph_redistance


    end if


    call mph_FillVars(solnData(DFUN_VAR,:,:),solnData(PFUN_VAR,:,:),&
                      solnData(CURV_VAR,:,:),&
                      solnData(THCO_VAR,:,:),solnData(CPRS_VAR,:,:),&
                      solnData(VISC_VAR,:,:),&
                      facexData(RH1F_VAR,:,:),faceyData(RH1F_VAR,:,:),&
                      facexData(RH2F_VAR,:,:),faceyData(RH2F_VAR,:,:),&
                      facexData(AL1F_VAR,:,:),faceyData(AL1F_VAR,:,:),&
                      facexData(AL2F_VAR,:,:),faceyData(AL2F_VAR,:,:),&
                      solnData(NRMX_VAR,:,:),solnData(NRMY_VAR,:,:),&
                      solnData(SMHV_VAR,:,:),solnData(SMRH_VAR,:,:))

#endif
#endif

else if (jump_flag .eqv. .TRUE.) then
                      
#ifdef IBM

#else

#ifdef MPH_DEBUG

#else
    call mph_PressureJumps(solnData(DFUN_VAR,:,:),solnData(PFUN_VAR,:,:),&
                           solnData(CURV_VAR,:,:),&
                           facexData(RH1F_VAR,:,:),faceyData(RH1F_VAR,:,:),&
                           facexData(RH2F_VAR,:,:),faceyData(RH2F_VAR,:,:),&
                           solnData(SIGP_VAR,:,:),&
                           facexData(SIGM_VAR,:,:),faceyData(SIGM_VAR,:,:),&
                           solnData(MDOT_VAR,:,:))

#endif
#endif

end if

    nullify(solnData)  
    nullify(facexData)
    nullify(faceyData)
 
end subroutine Multiphase_solver

