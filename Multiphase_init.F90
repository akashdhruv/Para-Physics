subroutine Multiphase_init()

   use Multiphase_data
   use physicaldata

#include "Solver.h"

   implicit none

   real,pointer,dimension(:,:) :: sf,pf,th,cp

   sf => ph_center(DFUN_VAR,:,:)
   pf => ph_center(PFUN_VAR,:,:)
   th => ph_center(THCO_VAR,:,:)
   cp => ph_center(CPRS_VAR,:,:)


   sf = 0.
   pf = 0.
   th = 0.
   cp = 0.

   mph_rho1 = 0.597
   mph_rho2 = 958.4

   mph_thco1 = 0.025
   mph_thco2 = 0.679

   mph_cp1 = 2030.0*mph_rho1
   mph_cp2 = 4216.0*mph_rho2

   nullify(sf)
   nullify(pf)
   nullify(th)
   nullify(cp)

end subroutine Multiphase_init
