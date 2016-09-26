module physicaldata

#include "Solver.h"

       implicit none

       real, save, allocatable :: ph_center(:,:,:), ph_facex(:,:,:), ph_facey(:,:,:)
       integer,save :: blockID(MAX_BLOCKS),blockLC(nblockx*nblocky,2)
       integer,save :: blockCount

       target :: ph_center, ph_facex, ph_facey

end module physicaldata
