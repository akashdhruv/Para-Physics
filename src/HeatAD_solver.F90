subroutine HeatAD_solver(tstep)

#define MPH_DEBUG
#include "Solver.h"

   use HeatAD_interface, only: heat_tempSolver,heat_tempSolver_ibm,heat_tempSolver_mph
   use physicaldata, only: localCENTER,localFACEX,localFACEY
   use MPI_data, only: blockCount
 
   implicit none
   
   integer, intent(in) :: tstep

   real, pointer, dimension(:,:,:,:) :: solnData,facexData,faceyData

   solnData  => localCENTER
   facexData => localFACEX
   faceyData => localFACEY

   solnData(:,:,blockCount,TOLD_VAR) = solnData(:,:,blockCount,TEMP_VAR) 

#ifdef SINGLEPHASE

   call heat_tempSolver(tstep,solnData(:,:,blockCount,TEMP_VAR),solnData(:,:,blockCount,TOLD_VAR),&
                         solnData(:,:,blockCount,MDOT_VAR),solnData(:,:,blockCount,SMRH_VAR),&
                         facexData(:,:,blockCount,VELC_VAR),faceyData(:,:,blockCount,VELC_VAR),&
                         facexData(:,:,blockCount,AL1F_VAR),faceyData(:,:,blockCount,AL1F_VAR),&
                         facexData(:,:,blockCount,AL2F_VAR),faceyData(:,:,blockCount,AL2F_VAR),&
                         solnData(:,:,blockCount,DFUN_VAR),solnData(:,:,blockCount,PFUN_VAR),&
                         solnData(:,:,blockCount,THCO_VAR),solnData(:,:,blockCount,CPRS_VAR))

#endif

#ifdef MULTIPHASE

#ifdef IBM

   call heat_tempSolver_ibm(tstep,solnData(:,:,blockCount,TEMP_VAR),solnData(:,:,blockCount,TOLD_VAR),&
                         solnData(:,:,blockCount,MDOT_VAR),solnData(:,:,blockCount,SMRH_VAR),&
                         facexData(:,:,blockCount,VELC_VAR),faceyData(:,:,blockCount,VELC_VAR),&
                         facexData(:,:,blockCount,AL1F_VAR),faceyData(:,:,blockCount,AL1F_VAR),&
                         facexData(:,:,blockCount,AL2F_VAR),faceyData(:,:,blockCount,AL2F_VAR),&
                         solnData(:,:,blockCount,DFUN_VAR),solnData(:,:,blockCount,PFUN_VAR),&
                         solnData(:,:,blockCount,THCO_VAR),solnData(:,:,blockCount,CPRS_VAR))
#else

#ifdef MPH_DEBUG

   call heat_tempSolver(tstep,solnData(:,:,blockCount,TEMP_VAR),solnData(:,:,blockCount,TOLD_VAR),&
                         solnData(:,:,blockCount,MDOT_VAR),solnData(:,:,blockCount,SMRH_VAR),&
                         facexData(:,:,blockCount,VELC_VAR),faceyData(:,:,blockCount,VELC_VAR),&
                         facexData(:,:,blockCount,AL1F_VAR),faceyData(:,:,blockCount,AL1F_VAR),&
                         facexData(:,:,blockCount,AL2F_VAR),faceyData(:,:,blockCount,AL2F_VAR),&
                         solnData(:,:,blockCount,DFUN_VAR),solnData(:,:,blockCount,PFUN_VAR),&
                         solnData(:,:,blockCount,THCO_VAR),solnData(:,:,blockCount,CPRS_VAR))
#else

   call heat_tempSolver_mph(tstep,solnData(:,:,blockCount,TEMP_VAR),solnData(:,:,blockCount,TOLD_VAR),&
                         solnData(:,:,blockCount,MDOT_VAR),solnData(:,:,blockCount,SMRH_VAR),&
                         facexData(:,:,blockCount,VELC_VAR),faceyData(:,:,blockCount,VELC_VAR),&
                         facexData(:,:,blockCount,AL1F_VAR),faceyData(:,:,blockCount,AL1F_VAR),&
                         facexData(:,:,blockCount,AL2F_VAR),faceyData(:,:,blockCount,AL2F_VAR),&
                         solnData(:,:,blockCount,DFUN_VAR),solnData(:,:,blockCount,PFUN_VAR),&
                         solnData(:,:,blockCount,THCO_VAR),solnData(:,:,blockCount,CPRS_VAR))
#endif

#endif

#endif

   nullify(solnData,facexData,faceyData)

end subroutine
