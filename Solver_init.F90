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

#ifdef INS
    call IncompNS_init()
#endif

#ifdef MULTIPHASE
    call Multiphase_init()
#endif

#ifdef IBM
    call IBM_init()
#endif

#ifdef ENERGY
    call HeatAD_init()
#endif

    call Driver_init()

end subroutine Solver_init
