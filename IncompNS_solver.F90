subroutine IncompNS_solver(tstep,p_counter)

#include "Solver.h"

    use IncompNS_interface, only: ins_momentum,ins_vorticity,ins_momentum_VD
    use physicaldata

    implicit none

    integer, intent(in) :: tstep
    integer, intent(out) :: p_counter

    real, pointer,dimension(:,:,:) :: solnData,facexData,faceyData

    solnData => ph_center
    facexData => ph_facex
    faceyData => ph_facey

#ifdef SINGLEPHASE
    call ins_momentum(tstep,p_counter,solnData(PRES_VAR,:,:),&
                      facexData(VELC_VAR,:,:),faceyData(VELC_VAR,:,:),&
                      facexData(IBMF_VAR,:,:),faceyData(IBMF_VAR,:,:))
#endif

#ifdef MULTIPHASE
    call ins_momentum_VD(tstep,p_counter,solnData(PRES_VAR,:,:),&
                         facexData(VELC_VAR,:,:),faceyData(VELC_VAR,:,:),&
                         solnData(VISC_VAR,:,:),&
                         facexData(RHOF_VAR,:,:),faceyData(RHOF_VAR,:,:),&                                       
                         facexData(IBMF_VAR,:,:),faceyData(IBMF_VAR,:,:))

#endif

    call ins_vorticity(tstep,solnData(VORT_VAR,:,:),&
                       facexData(VELC_VAR,:,:),faceyData(VELC_VAR,:,:),&
                       solnData(DFUN_VAR,:,:))


    nullify(solnData)
    nullify(facexData)
    nullify(faceyData)

end subroutine
