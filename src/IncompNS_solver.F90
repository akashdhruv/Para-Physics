subroutine IncompNS_solver(tstep,p_counter)

#include "Solver.h"

    use IncompNS_interface, only: ins_momentum,ins_vorticity,ins_momentum_VD
    use physicaldata, only: solnData,facexData,faceyData

    implicit none

    integer, intent(in) :: tstep
    integer, intent(out) :: p_counter

#ifdef SINGLEPHASE
    call ins_momentum(tstep,p_counter,solnData(PRES_VAR,:,:),&
                      facexData(VELC_VAR,:,:),faceyData(VELC_VAR,:,:),&
                      facexData(USTR_VAR,:,:),faceyData(USTR_VAR,:,:),&
                      facexData(IBMF_VAR,:,:),faceyData(IBMF_VAR,:,:))
#endif

#ifdef MULTIPHASE
    call ins_momentum_VD(tstep,p_counter,solnData(PRES_VAR,:,:),&
                         facexData(VELC_VAR,:,:),faceyData(VELC_VAR,:,:),&
                         facexData(USTR_VAR,:,:),faceyData(USTR_VAR,:,:),&
                         solnData(VISC_VAR,:,:),&
                         facexData(RH1F_VAR,:,:),faceyData(RH1F_VAR,:,:),&                                       
                         facexData(RH2F_VAR,:,:),faceyData(RH2F_VAR,:,:),&
                         facexData(IBMF_VAR,:,:),faceyData(IBMF_VAR,:,:),&
                         solnData(SIGP_VAR,:,:),&
                         facexData(SIGM_VAR,:,:),faceyData(SIGM_VAR,:,:))

#endif

    !call ins_vorticity(tstep,solnData(VORT_VAR,:,:),&
    !                   facexData(VELC_VAR,:,:),faceyData(VELC_VAR,:,:),&
    !                   solnData(DFUN_VAR,:,:))

end subroutine
