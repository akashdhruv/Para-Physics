subroutine IncompNS_solver(tstep,p_counter)

#include "Solver.h"

    use IncompNS_interface, only: ins_momentum,ins_vorticity,ins_momentum_VD
    use physicaldata, only: localCENTER,localFACEX,localFACEY
    use MPI_data, only: blockCount

    implicit none

    integer, intent(in) :: tstep
    integer, intent(out) :: p_counter

    real, pointer, dimension(:,:,:,:) :: solnData,facexData,faceyData

    solnData  => localCENTER
    facexData => localFACEX
    faceyData => localFACEY

#ifdef SINGLEPHASE
    call ins_momentum(tstep,p_counter,solnData(:,:,blockCount,PRES_VAR),&
                      facexData(:,:,blockCount,VELC_VAR),faceyData(:,:,blockCount,VELC_VAR),&
                      facexData(:,:,blockCount,USTR_VAR),faceyData(:,:,blockCount,USTR_VAR),&
                      facexData(:,:,blockCount,IBMF_VAR),faceyData(:,:,blockCount,IBMF_VAR))
#endif

#ifdef MULTIPHASE
    call ins_momentum_VD(tstep,p_counter,solnData(:,:,blockCount,PRES_VAR),&
                         facexData(:,:,blockCount,VELC_VAR),faceyData(:,:,blockCount,VELC_VAR),&
                         facexData(:,:,blockCount,USTR_VAR),faceyData(:,:,blockCount,USTR_VAR),&
                         solnData(:,:,blockCount,VISC_VAR),&
                         facexData(:,:,blockCount,RH1F_VAR),faceyData(:,:,blockCount,RH1F_VAR),&                                       
                         facexData(:,:,blockCount,RH2F_VAR),faceyData(:,:,blockCount,RH2F_VAR),&
                         facexData(:,:,blockCount,IBMF_VAR),faceyData(:,:,blockCount,IBMF_VAR),&
                         solnData(:,:,blockCount,SIGP_VAR),&
                         facexData(:,:,blockCount,SIGM_VAR),faceyData(:,:,blockCount,SIGM_VAR))

#endif

    !call ins_vorticity(tstep,solnData(:,:,blockCount,VORT_VAR),&
    !                   facexData(:,:,blockCount,VELC_VAR),faceyData(:,:,blockCount,VELC_VAR),&
    !                   solnData(:,:,blockCount,DFUN_VAR))


    nullify(solnData,facexData,faceyData)

end subroutine
