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
    call ins_momentum(tstep,p_counter,solnData(:,:,PRES_VAR,blockCount),&
                      facexData(:,:,VELC_VAR,blockCount),faceyData(:,:,VELC_VAR,blockCount),&
                      facexData(:,:,USTR_VAR,blockCount),faceyData(:,:,USTR_VAR,blockCount),&
                      facexData(:,:,IBMF_VAR,blockCount),faceyData(:,:,IBMF_VAR,blockCount))
#endif

#ifdef MULTIPHASE
    call ins_momentum_VD(tstep,p_counter,solnData(:,:,PRES_VAR,blockCount),&
                         facexData(:,:,VELC_VAR,blockCount),faceyData(:,:,VELC_VAR,blockCount),&
                         facexData(:,:,USTR_VAR,blockCount),faceyData(:,:,USTR_VAR,blockCount),&
                         solnData(:,:,VISC_VAR,blockCount),&
                         facexData(:,:,RH1F_VAR,blockCount),faceyData(:,:,RH1F_VAR,blockCount),&                                       
                         facexData(:,:,RH2F_VAR,blockCount),faceyData(:,:,RH2F_VAR,blockCount),&
                         facexData(:,:,IBMF_VAR,blockCount),faceyData(:,:,IBMF_VAR,blockCount),&
                         solnData(:,:,SIGP_VAR,blockCount),&
                         facexData(:,:,SIGM_VAR,blockCount),faceyData(:,:,SIGM_VAR,blockCount))

#endif

    !call ins_vorticity(tstep,solnData(:,:,VORT_VAR,blockCount),&
    !                   facexData(:,:,VELC_VAR,blockCount),faceyData(:,:,VELC_VAR,blockCount),&
    !                   solnData(:,:,DFUN_VAR,blockCount))


    nullify(solnData,facexData,faceyData)

end subroutine
