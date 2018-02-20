module IncompNS_data

#include "Solver.h"

    implicit none

    real, save :: ins_u_res,ins_v_res,ins_p_res,ins_w_res
    real, save :: ins_inRe, ins_sigma, ins_cfl
    real, save :: ins_maxdiv, ins_mindiv
    real, save, dimension(2) :: ins_umaxmin,ins_vmaxmin

    real, save :: ins_Qin,ins_Qout,ins_convvel(2,2), ins_Qinout
    real, save :: ins_maxU
    double precision, save :: ins_timePoisson
    real, save :: ins_gravX,ins_gravY

    real, save, dimension(4,Nyb+2,MAX_BLOCKS) :: u_old
    real, save, dimension(Nxb+2,4,MAX_BLOCKS) :: v_old

    real, save :: ins_upEx1, ins_upEx2
    real, save :: ins_dnEx1, ins_dnEx2

    real, save :: ins_upIn1, ins_upIn2
    real, save :: ins_dnIn1, ins_dnIn2

    integer, save :: ins_xl_bnd, ins_xr_bnd, ins_yl_bnd, ins_yr_bnd

end module IncompNS_data
