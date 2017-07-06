module physicaldata

       implicit none

       real, save, pointer, dimension(:,:,:,:) :: localCENTER, localFACEX, localFACEY
       real, save, pointer, dimension(:,:,:,:) :: sharedCENTER, sharedFACEX, sharedFACEY

       real, save, pointer, dimension(:) :: dataTARGET
       real, save, pointer, dimension(:,:) :: northORIGIN, southORIGIN, eastORIGIN, westORIGIN

end module physicaldata
