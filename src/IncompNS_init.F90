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

      solnData(:,:,blockCount,PRES_VAR) = 0.0
      solnData(:,:,blockCount,VORT_VAR) = 0.0

      facexData(:,:,blockCount,VELC_VAR) = 0.0
      faceyData(:,:,blockCount,VELC_VAR) = 0.0

      facexData(:,:,blockCount,USTR_VAR) = 0.0
      faceyData(:,:,blockCount,USTR_VAR) = 0.0
     
      ins_inRe  = 0.001

      ins_sigma = 0.1
      ins_cfl   = 0.1
      ins_timePoisson = 0.0

      ins_gravX =  0.0
      ins_gravY =  0.0

      nullify(facexData,faceyData,solnData)

end subroutine IncompNS_init
