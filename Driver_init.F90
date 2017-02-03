subroutine Driver_init()
 
#include "Solver.h"

      use Driver_data
      use Grid_data
      use IncompNS_data
      use Multiphase_data
      use HeatAD_data, only : ht_Pr

      implicit none
      
      real :: dt_sig, dt_cfl, dt_mph, dx_min, dy_min, dt_temp

      dr_t  = TIME_END


      dx_min = gr_dx
      dy_min = gr_dy

      dt_sig = ins_sigma*(min(dx_min,dy_min)**2)/ins_inRe
      dt_cfl = ins_cfl*min(dx_min,dy_min)

#ifdef ENERGY
      dt_temp = dt_sig*ht_Pr
#endif

#ifdef MULTIPHASE
     
      !dt_mph = min(0.05*((mph_cp1)/mph_thco1)*min(dx_min,dy_min)**2,&
      !             0.05*((mph_cp2)/mph_thco2)*min(dx_min,dy_min)**2)

      !_________The Missing Data simulation__________!

      dt_mph = (ins_sigma)*(ht_Pr) / (ins_inRe*MAX( 1.0/(gr_dx*gr_dx), 1.0/(gr_dy*gr_dy)))

      !____________________End_______________________!
#endif

      dr_dt = min(dt_sig,dt_cfl)

#ifdef ENERGY
 
      dr_dt = min(dr_dt,dt_temp)

#endif

#ifdef MULTIPHASE

      !_______For testing in Stefan_Problem_____!
      !dr_dt = min(dr_dt,dt_mph)
      !dr_dt = 0.000001
      !___________________End___________________!

      dr_dt = min(dr_dt,dt_mph)


#endif

      dr_dt = 1.0d-6

      dr_nt = dr_t/dr_dt
      dr_nt = 0
 
      !dr_tile = (Nyb)/(ceiling((Nyb)/40.)*NTHREADS)
      dr_tile = (Nyb)/NTHREADS
      
end subroutine Driver_init
