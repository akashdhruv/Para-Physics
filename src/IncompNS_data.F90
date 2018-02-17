module IncompNS_data

#include "Solver.h"

    implicit none

    real, save :: ins_u_res,ins_v_res,ins_p_res,ins_w_res
    real, save :: ins_inRe, ins_sigma, ins_cfl
    real, save :: ins_maxdiv, ins_mindiv
    real, save, dimension(2) :: ins_umaxmin,ins_vmaxmin

    real, save :: ins_Qin,ins_Qout,ins_convvel, ins_Qinout
    real, save :: ins_maxU
    double precision, save :: ins_timePoisson
    real, save :: ins_gravX,ins_gravY

    real, save, dimension(2,Nyb+2,MAX_BLOCKS) :: u_old

end module IncompNS_data
