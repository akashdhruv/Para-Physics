subroutine MPI_applyBC(ivar,datatype)

#include "Solver.h"

        use MPI_data
        use physicaldata, only: localCENTER,localFACEX,localFACEY,&
                                sharedCENTER,sharedFACEX,sharedFACEY
        use MPI_interface, only: MPI_applyBC_SHM, MPI_applyBC_RMA, MPI_applyBC_DIS
                               
        implicit none

        integer, intent(in) :: ivar,datatype


#ifdef MPI_DIS
    if(datatype == CENTER) call MPI_applyBC_DIS(localCENTER(:,:,:,ivar))
    if(datatype == FACEX)  call MPI_applyBC_DIS(localFACEX(:,:,:,ivar))
    if(datatype == FACEY)  call MPI_applyBC_DIS(localFACEY(:,:,:,ivar))
#endif

#ifdef MPI_SHM
    if(datatype == CENTER) call MPI_applyBC_SHM(localCENTER(:,:,:,ivar),sharedCENTER(:,:,ivar,:))
    if(datatype == FACEX)  call MPI_applyBC_SHM(localFACEX(:,:,:,ivar),sharedFACEX(:,:,ivar,:))
    if(datatype == FACEY)  call MPI_applyBC_SHM(localFACEY(:,:,:,ivar),sharedFACEY(:,:,ivar,:))
#endif
   
#ifdef MPI_RMA
    if(datatype == CENTER) call MPI_applyBC_RMA(localCENTER(:,:,:,ivar))
    if(datatype == FACEX)  call MPI_applyBC_RMA(localFACEX(:,:,:,ivar))
    if(datatype == FACEY)  call MPI_applyBC_RMA(localFACEY(:,:,:,ivar))
#endif

end subroutine MPI_applyBC
