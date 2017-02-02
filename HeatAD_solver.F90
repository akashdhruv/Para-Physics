subroutine HeatAD_solver(tstep)

#include "Solver.h"

   use HeatAD_interface, only: heat_tempSolver
   use physicaldata
 
   implicit none
   
   integer, intent(in) :: tstep
   real, pointer, dimension(:,:,:) :: solnData,facexData,faceyData

   solnData => ph_center
   facexData => ph_facex
   faceyData => ph_facey

   solnData(TOLD_VAR,:,:) = solnData(TEMP_VAR,:,:) 

   call heat_tempSolver(tstep,solnData(TEMP_VAR,:,:),solnData(TOLD_VAR,:,:),&
                         solnData(MDOT_VAR,:,:),solnData(SMRH_VAR,:,:),&
                         facexData(VELC_VAR,:,:),faceyData(VELC_VAR,:,:),&
                         facexData(AL1F_VAR,:,:),faceyData(AL1F_VAR,:,:),&
                         facexData(AL2F_VAR,:,:),faceyData(AL2F_VAR,:,:),&
                         solnData(DFUN_VAR,:,:),solnData(PFUN_VAR,:,:),&
                         solnData(THCO_VAR,:,:),solnData(CPRS_VAR,:,:))

   nullify(solnData)
   nullify(facexData)
   nullify(faceyData)

end subroutine
