module physicaldata

       implicit none

       real, save, pointer, dimension(:,:,:,:) :: localCENTER, localFACEX, localFACEY          ! Local data per block
       real, save, pointer, dimension(:,:,:,:) :: sharedCENTER, sharedFACEX, sharedFACEY       ! Shared data per node

       real, save, pointer, dimension(:) :: dataTARGET                                         ! Target data for RMA
       real, save, pointer, dimension(:,:) :: northORIGIN, southORIGIN, eastORIGIN, westORIGIN ! Origin data for RMA

end module physicaldata
