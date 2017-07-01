subroutine HeatAD_solver(tstep)

#define MPH_DEBUG
#include "Solver.h"

   use HeatAD_interface, only: heat_tempSolver,heat_tempSolver_ibm,heat_tempSolver_mph
   use physicaldata, only: localCENTER,localFACEX,localFACEY
 
   implicit none
   
   integer, intent(in) :: tstep

   real, pointer, dimension(:,:,:) :: solnData,facexData,faceyData

   solnData  => localCENTER
   facexData => localFACEX
   faceyData => localFACEY

   solnData(:,:,TOLD_VAR) = solnData(:,:,TEMP_VAR) 

#ifdef SINGLEPHASE

   call heat_tempSolver(tstep,solnData(:,:,TEMP_VAR),solnData(:,:,TOLD_VAR),&
                         solnData(:,:,MDOT_VAR),solnData(:,:,SMRH_VAR),&
                         facexData(:,:,VELC_VAR),faceyData(:,:,VELC_VAR),&
                         facexData(:,:,AL1F_VAR),faceyData(:,:,AL1F_VAR),&
                         facexData(:,:,AL2F_VAR),faceyData(:,:,AL2F_VAR),&
                         solnData(:,:,DFUN_VAR),solnData(:,:,PFUN_VAR),&
                         solnData(:,:,THCO_VAR),solnData(:,:,CPRS_VAR))

#endif

#ifdef MULTIPHASE

#ifdef IBM

   call heat_tempSolver_ibm(tstep,solnData(:,:,TEMP_VAR),solnData(:,:,TOLD_VAR),&
                         solnData(:,:,MDOT_VAR),solnData(:,:,SMRH_VAR),&
                         facexData(:,:,VELC_VAR),faceyData(:,:,VELC_VAR),&
                         facexData(:,:,AL1F_VAR),faceyData(:,:,AL1F_VAR),&
                         facexData(:,:,AL2F_VAR),faceyData(:,:,AL2F_VAR),&
                         solnData(:,:,DFUN_VAR),solnData(:,:,PFUN_VAR),&
                         solnData(:,:,THCO_VAR),solnData(:,:,CPRS_VAR))
#else

#ifdef MPH_DEBUG

   call heat_tempSolver(tstep,solnData(:,:,TEMP_VAR),solnData(:,:,TOLD_VAR),&
                         solnData(:,:,MDOT_VAR),solnData(:,:,SMRH_VAR),&
                         facexData(:,:,VELC_VAR),faceyData(:,:,VELC_VAR),&
                         facexData(:,:,AL1F_VAR),faceyData(:,:,AL1F_VAR),&
                         facexData(:,:,AL2F_VAR),faceyData(:,:,AL2F_VAR),&
                         solnData(:,:,DFUN_VAR),solnData(:,:,PFUN_VAR),&
                         solnData(:,:,THCO_VAR),solnData(:,:,CPRS_VAR))
#else

   call heat_tempSolver_mph(tstep,solnData(:,:,TEMP_VAR),solnData(:,:,TOLD_VAR),&
                         solnData(:,:,MDOT_VAR),solnData(:,:,SMRH_VAR),&
                         facexData(:,:,VELC_VAR),faceyData(:,:,VELC_VAR),&
                         facexData(:,:,AL1F_VAR),faceyData(:,:,AL1F_VAR),&
                         facexData(:,:,AL2F_VAR),faceyData(:,:,AL2F_VAR),&
                         solnData(:,:,DFUN_VAR),solnData(:,:,PFUN_VAR),&
                         solnData(:,:,THCO_VAR),solnData(:,:,CPRS_VAR))

#endif

#endif

#endif

   nullify(solnData,facexData,faceyData)

end subroutine
