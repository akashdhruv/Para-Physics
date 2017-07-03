subroutine IncompNS_init()

      use physicaldata, only: localCENTER,localFACEX,localFACEY
      use IncompNS_data
      use MPI_data, only: blockCount

#include "Solver.h"
   
      implicit none
        
      real,pointer,dimension(:,:,:,:) :: solnData,facexData,faceyData

      solnData  => localCENTER
      facexData => localFACEX
      faceyData => localFACEY      

      solnData(:,:,:,PRES_VAR) = 0.0
      solnData(:,:,:,VORT_VAR) = 0.0

      facexData(:,:,:,VELC_VAR) = 0.0
      faceyData(:,:,:,VELC_VAR) = 0.0

      facexData(:,:,:,USTR_VAR) = 0.0
      faceyData(:,:,:,USTR_VAR) = 0.0

      facexData(:,:,:,GOLD_VAR) = 0.0
      faceyData(:,:,:,GOLD_VAR) = 0.0
     
      ins_inRe  = 0.001

      ins_sigma = 0.1
      ins_cfl   = 0.1
      ins_timePoisson = 0.0

      ins_gravX =  0.0
      ins_gravY =  0.0

      nullify(facexData,faceyData,solnData)

end subroutine IncompNS_init
