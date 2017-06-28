module physicaldata

       implicit none

       real, save, pointer, dimension(:,:,:) :: localCENTER, localFACEX, localFACEY

       real, save, pointer, dimension(:,:,:) :: northCENTER, southCENTER, westCENTER, eastCENTER
       real, save, pointer, dimension(:,:,:) :: northFACEX, southFACEX, westFACEX, eastFACEX
       real, save, pointer, dimension(:,:,:) :: northFACEY, southFACEY, westFACEY, eastFACEY

       real, save, pointer, dimension(:) :: dataTARGET
       real, save, pointer, dimension(:) :: northORIGIN, southORIGIN, eastORIGIN, westORIGIN

end module physicaldata
