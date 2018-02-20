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

      ! Initialize variables relevant to INS unit
      solnData(:,:,PRES_VAR,:) = 0.0
      solnData(:,:,VORT_VAR,:) = 0.0

      facexData(:,:,VELC_VAR,:) = 1.0
      faceyData(:,:,VELC_VAR,:) = 1.0

      facexData(:,:,USTR_VAR,:) = 0.0
      faceyData(:,:,USTR_VAR,:) = 0.0

      facexData(:,:,GOLD_VAR,:) = 0.0
      faceyData(:,:,GOLD_VAR,:) = 0.0
     
      ins_inRe  = 0.010        ! Reynolds number

      ins_sigma = 0.010        ! CFL for diffusive  dt
      ins_cfl   = 0.010        ! CFL for convective dt
      ins_timePoisson = 0.0    ! variable to store Poisson solver time

      ins_gravX =  0.0         ! 1/(Fr**2) in X
      ins_gravY =  0.0         ! 1/(Fr**2) in Y


      ins_xl_bnd = NOSLIP
      ins_xr_bnd = NOSLIP
      ins_yl_bnd = NOSLIP
      ins_yr_bnd = NOSLIP

      ins_upEx1 = -10.0
      ins_upEx2 = -07.0
      ins_dnEx1 = -10.0
      ins_dnEx2 = -07.0

      ins_upIn1 = 07.0
      ins_upIn2 = 10.0
      ins_dnIn1 = 07.0
      ins_dnIn2 = 10.0

      nullify(facexData,faceyData,solnData)

end subroutine IncompNS_init
