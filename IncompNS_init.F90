subroutine IncompNS_init()

      use physicaldata
      use IncompNS_data

#include "Solver.h"
   
      implicit none
      
      real,pointer,dimension(:,:,:) :: facexData,faceyData,solnData

      solnData => ph_center
      facexData => ph_facex
      faceyData => ph_facey

      solnData(PRES_VAR,:,:) = 0.0
      facexData(VELC_VAR,:,:) = 0.0
      faceyData(VELC_VAR,:,:) = 0.0

      nullify(solnData)
      nullify(facexData)
      nullify(faceyData)
      
      ins_inRe = 0.001
      !ins_sigma = 0.001
      !ins_cfl = 0.001
      ins_sigma = 0.01
      ins_cfl = 0.15

end subroutine IncompNS_init
