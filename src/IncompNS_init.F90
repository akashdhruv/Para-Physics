subroutine IncompNS_init()

      use physicaldata, only: solnData,facexData,faceyData
      use IncompNS_data

#include "Solver.h"
   
      implicit none
      
      solnData(PRES_VAR,:,:) = 0.0
      solnData(VORT_VAR,:,:) = 0.0

      facexData(VELC_VAR,:,:) = 0.0
      faceyData(VELC_VAR,:,:) = 0.0
     
      ins_inRe  = 0.001

      ins_sigma = 0.1
      ins_cfl   = 0.1
      ins_timePoisson = 0.0

      ins_gravX =  0.0
      ins_gravY =  0.0

end subroutine IncompNS_init
