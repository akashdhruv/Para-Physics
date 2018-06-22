subroutine Driver_init()
 
#include "Solver.h"

      use Driver_data
      use Grid_data
      use IncompNS_data
      use Multiphase_data
      use HeatAD_data, only : ht_Pr
      use physicaldata, only: localCENTER,localFACEX,localFACEY
      use MPI_data, only: blockCount

      implicit none
      
      real :: dt_sig, dt_cfl, dt_mph, dx_min, dy_min, dt_temp,dt_s,dt_g
      real :: velcoeff
      real, parameter :: eps = 1e-12
      real, parameter :: pi  = acos(-1.0)

      real, pointer, dimension(:,:,:,:) :: solnData,facexData,faceyData

      solnData  => localCENTER
      facexData => localFACEX
      faceyData => localFACEY

      dr_t   = TIME_END

      dx_min = gr_dx
      dy_min = gr_dy

      dt_cfl = ins_cfl*min(gr_dx,gr_dy)/max(maxval(abs(facexData(:,:,VELC_VAR,:))),maxval(abs(faceyData(:,:,VELC_VAR,:))))

      dt_g = ins_sigma*sqrt((abs(ins_gravX)/dx_min)+(abs(ins_gravY)/dy_min))

      dt_sig  = ins_sigma*(min(dx_min,dy_min)**2)/ins_inRe
   
      dt_temp = ht_Pr*ins_sigma*(min(dx_min,dy_min)**2)/ins_inRe

      dr_dt = min(dt_cfl,dt_sig)

      dr_nt = NEND_END

      dr_tile = (Nyb)/NTHREADS
      
      nullify(solnData,facexData,faceyData)

end subroutine Driver_init
