subroutine ins_rescaleVel(u,v)

#include "Solver.h"

    use IncompNS_data, only: ins_maxU
    use MPI_interface, only: MPI_CollectResiduals

    implicit none

    real,dimension(:,:),intent(inout) :: u,v
    real :: maxU,maxV,maxVel

    maxU = maxval(abs(u))
    maxV = maxval(abs(v))

    maxVel = max(maxU,maxV)

    call MPI_CollectResiduals(maxVel,ins_maxU,3) 

    u = u/ins_maxU
    v = v/ins_maxU

end subroutine
