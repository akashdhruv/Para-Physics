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


      velcoeff =  MAX( MAXVAL(ABS(facexData(:,:,blockCount,VELC_VAR))/gr_dx), &
                       MAXVAL(ABS(faceyData(:,:,blockCount,VELC_VAR))/gr_dy))

      dt_cfl = ins_cfl*min(gr_dx,gr_dy)

      dt_g = ins_sigma*sqrt((abs(ins_gravX)/dx_min)+(abs(ins_gravY)/dy_min))

#ifdef SINGLEPHASE

      dt_sig  = ins_sigma*(min(dx_min,dy_min)**2)/ins_inRe
   
#ifdef ENERGY
      dt_temp = ht_Pr*ins_sigma*(min(dx_min,dy_min)**2)/ins_inRe
#endif

      dr_dt = min(dt_cfl,dt_sig)

      !dr_dt = 1./(1./dt_cfl + 1./dt_sig + &
      !        sqrt((1./(dt_sig+dt_cfl))**2)+4*(dt_g**2))

#endif


#ifdef MULTIPHASE
     
      !dt_mph = min(0.05*((mph_cp1)/mph_thco1)*min(dx_min,dy_min)**2,&
      !             0.05*((mph_cp2)/mph_thco2)*min(dx_min,dy_min)**2)

      !_________The Missing Data simulation__________!

      !dt_mph = (ins_sigma)*(ht_Pr) / (ins_inRe*MAX( 1.0/(gr_dx*gr_dx), 1.0/(gr_dy*gr_dy)))

      !____________________End_______________________!

      !dt_s    = ins_sigma*sqrt((1+(mph_rho1/mph_rho2))/(pi*mph_sten))*0.5*min(dx_min**(3./2.),dy_min**(3./2.))
      !dt_sig  = (1./max(1.0,(mph_vis1/mph_vis2)/(mph_rho1/mph_rho2)))*ins_sigma*(min(dx_min,dy_min)**2)/ins_inRe

      dt_sig  = ins_sigma*(min(dx_min,dy_min)**2)/ins_inRe

#ifdef ENERGY
      !dt_temp = ht_Pr*(1./max(1.0,(mph_thco1/mph_cp1)/(mph_thco2/mph_cp2)))*ins_sigma*(min(dx_min,dy_min)**2)/ins_inRe
      dt_temp = ht_Pr*ins_sigma*(min(dx_min,dy_min)**2)/ins_inRe
#endif

      !dr_dt =  1./(1./dt_cfl + 1./min(dt_sig,dt_temp) + &
      ! sqrt((1./(min(dt_sig,dt_temp)+dt_cfl))**2 + &
      !     4*(dt_g**2)+4*((1./dt_s)**2)))

       dr_dt = min(dt_cfl,dt_sig)
#endif

      dr_nt = dr_t/dr_dt

      dr_tile = (Nyb)/NTHREADS
      
      nullify(solnData,facexData,faceyData)

end subroutine Driver_init
