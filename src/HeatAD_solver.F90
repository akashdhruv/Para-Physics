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

   solnData(:,:,TOLD_VAR,blockCount) = solnData(:,:,TEMP_VAR,blockCount) 

#ifdef SINGLEPHASE

   call heat_tempSolver(tstep,solnData(:,:,TEMP_VAR,blockCount),solnData(:,:,TOLD_VAR,blockCount),&
                         solnData(:,:,MDOT_VAR,blockCount),solnData(:,:,SMRH_VAR,blockCount),&
                         facexData(:,:,VELC_VAR,blockCount),faceyData(:,:,VELC_VAR,blockCount),&
                         facexData(:,:,AL1F_VAR,blockCount),faceyData(:,:,AL1F_VAR,blockCount),&
                         facexData(:,:,AL2F_VAR,blockCount),faceyData(:,:,AL2F_VAR,blockCount),&
                         solnData(:,:,DFUN_VAR,blockCount),solnData(:,:,PFUN_VAR,blockCount),&
                         solnData(:,:,THCO_VAR,blockCount),solnData(:,:,CPRS_VAR,blockCount))

#endif

#ifdef MULTIPHASE

#ifdef IBM

   call heat_tempSolver_ibm(tstep,solnData(:,:,TEMP_VAR,blockCount),solnData(:,:,TOLD_VAR,blockCount),&
                         solnData(:,:,MDOT_VAR,blockCount),solnData(:,:,SMRH_VAR,blockCount),&
                         facexData(:,:,VELC_VAR,blockCount),faceyData(:,:,VELC_VAR,blockCount),&
                         facexData(:,:,AL1F_VAR,blockCount),faceyData(:,:,AL1F_VAR,blockCount),&
                         facexData(:,:,AL2F_VAR,blockCount),faceyData(:,:,AL2F_VAR,blockCount),&
                         solnData(:,:,DFUN_VAR,blockCount),solnData(:,:,PFUN_VAR,blockCount),&
                         solnData(:,:,THCO_VAR,blockCount),solnData(:,:,CPRS_VAR,blockCount))
#else

#ifdef MPH_DEBUG

   call heat_tempSolver(tstep,solnData(:,:,TEMP_VAR,blockCount),solnData(:,:,TOLD_VAR,blockCount),&
                         solnData(:,:,MDOT_VAR,blockCount),solnData(:,:,SMRH_VAR,blockCount),&
                         facexData(:,:,VELC_VAR,blockCount),faceyData(:,:,VELC_VAR,blockCount),&
                         facexData(:,:,AL1F_VAR,blockCount),faceyData(:,:,AL1F_VAR,blockCount),&
                         facexData(:,:,AL2F_VAR,blockCount),faceyData(:,:,AL2F_VAR,blockCount),&
                         solnData(:,:,DFUN_VAR,blockCount),solnData(:,:,PFUN_VAR,blockCount),&
                         solnData(:,:,THCO_VAR,blockCount),solnData(:,:,CPRS_VAR,blockCount))
#else

   call heat_tempSolver_mph(tstep,solnData(:,:,TEMP_VAR,blockCount),solnData(:,:,TOLD_VAR,blockCount),&
                         solnData(:,:,MDOT_VAR,blockCount),solnData(:,:,SMRH_VAR,blockCount),&
                         facexData(:,:,VELC_VAR,blockCount),faceyData(:,:,VELC_VAR,blockCount),&
                         facexData(:,:,AL1F_VAR,blockCount),faceyData(:,:,AL1F_VAR,blockCount),&
                         facexData(:,:,AL2F_VAR,blockCount),faceyData(:,:,AL2F_VAR,blockCount),&
                         solnData(:,:,DFUN_VAR,blockCount),solnData(:,:,PFUN_VAR,blockCount),&
                         solnData(:,:,THCO_VAR,blockCount),solnData(:,:,CPRS_VAR,blockCount))

#endif

#endif

#endif

   nullify(solnData,facexData,faceyData)

end subroutine
