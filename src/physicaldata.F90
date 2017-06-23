module physicaldata

       implicit none

       real, save, pointer, dimension(:,:,:) :: solnData, facexData, faceyData
       real, save, pointer, dimension(:,:,:) :: SHD_solnData, SHD_facexData, SHD_faceyData

end module physicaldata
