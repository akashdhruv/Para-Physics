subroutine IO_display(u_res,v_res,w_res,p_res,T_res,p_counter,simtime,maxdiv,mindiv,umaxmin,vmaxmin)

#include "Solver.h"

         implicit none

         real, intent(in) :: u_res,v_res,p_res,simtime,T_res,maxdiv,mindiv,w_res
         integer, intent(in) :: p_counter
         real, intent(in),dimension(2) :: umaxmin,vmaxmin

         print *,"**************PARAMETER DISPLAY***************"
         print *,"Simulation Time         : ",simtime
         print *,"U Velocity Residual     : ",u_res
         print *,"V Velocity Residual     : ",v_res
         print *,"Vorticity Residual      : ",w_res
         print *,"Pressure Residual       : ",p_res
         print *,"Temperature Residual    : ",T_res
         print *,"Pressure Poisson Counter: ",p_counter
         print *,"MAX DIV: ",maxdiv," MIN DIV: ",mindiv
         print *,"MAX U  : ",umaxmin(1)," MIN U  :",umaxmin(2)
         print *,"MAX V  : ",vmaxmin(1)," MIN V  :",vmaxmin(2)
        
end subroutine IO_display
