subroutine Driver_init()
 
#include "Solver.h"

      use Driver_data
      use Grid_data
      use IncompNS_data
      use Multiphase_data

      implicit none
      
      real :: dt_sig, dt_cfl, dt_mph, dx_min, dy_min

      dr_t  = TIME_END


      dx_min = gr_dx
      dy_min = gr_dy

      dt_sig = ins_sigma*(min(dx_min,dy_min)**2)/ins_inRe
      dt_cfl = ins_cfl*min(dx_min,dy_min)

#ifdef MULTIPHASE
      dt_mph = min(0.05*((mph_cp1)/mph_thco1)*min(dx_min,dy_min)**2,&
                   0.05*((mph_cp2)/mph_thco2)*min(dx_min,dy_min)**2)
#endif

      dr_dt = min(dt_sig,dt_cfl)

#ifdef MULTIPHASE
      dr_dt = min(dr_dt,dt_mph)
      dr_dt = 0.000001
#endif

      dr_nt = dr_t/dr_dt
 
      !dr_tile = (Nyb)/(ceiling((Nyb)/40.)*NTHREADS)
      dr_tile = (Nyb)/NTHREADS
      
end subroutine Driver_init
