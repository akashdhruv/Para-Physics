module physicaldata

       implicit none

       real, save, pointer, dimension(:,:,:) :: solnData, facexData, faceyData
       real, save, pointer, dimension(:,:,:) :: northCENTER, southCENTER, westCENTER, eastCENTER
       real, save, pointer, dimension(:,:,:) :: northFACEX, southFACEX, westFACEX, eastFACEX
       real, save, pointer, dimension(:,:,:) :: northFACEY, southFACEY, westFACEY, eastFACEY

end module physicaldata
