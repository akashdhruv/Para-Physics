subroutine ins_computeQinout()

      use IncompNS_data, only: ins_Qin,ins_Qout
      use physicaldata
      use Grid_data
      use MPI_data


      implicit none

      include "mpif.h"

      real :: Qin1, Qout1

      Qin1     = 0.
      Qout1    = 0.
      ins_Qin  = 0.
      ins_Qout = 0.
 
      




end subroutine
