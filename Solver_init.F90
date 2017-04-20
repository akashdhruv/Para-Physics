subroutine Solver_init

#include "Solver.h"

    use MPI_interface, only: MPIsolver_init
    use Grid_interface, only: Grid_init
    use IncompNS_interface, only: IncompNS_init
    use HeatAD_interface, only: HeatAD_init
    use Driver_interface, only: Driver_init
    use Multiphase_interface, only: Multiphase_init
    use IBM_interface, only: IBM_init

    implicit none

    call MPIsolver_init()

    call Grid_init()

    call IncompNS_init()

    call Multiphase_init()

    call HeatAD_init()

    call IBM_init()

    call Driver_init()

end subroutine Solver_init
