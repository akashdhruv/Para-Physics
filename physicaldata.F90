module physicaldata

       implicit none

       real, save, allocatable :: ph_center(:,:,:), ph_facex(:,:,:), ph_facey(:,:,:)
       integer,save,allocatable :: blockID(:)
       integer,save :: blockCount

       target :: ph_center, ph_facex, ph_facey

end module physicaldata
