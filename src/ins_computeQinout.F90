subroutine ins_computeQinout(u,flg)

#include"Solver.h"

      use IncompNS_data, only: ins_Qin,ins_Qout
      use physicaldata
      use Grid_data
      use MPI_data


      implicit none
      
      real, dimension(:,:,:), intent(inout) :: u
      logical, intent(in) :: flg

      real :: Qaux
      integer :: blk,j

      real :: dxdz, dydz

      Qaux = 0.0
      dxdz = gr_dx
      dydz = gr_dy

      if(flg) then

        do blk=1,blockCount
           if ( xLC(blk) == 0) then
              do j=2,Nyb+1                
                 Qaux = Qaux + u(1,j,blk)*dydz
                 !Qaux = Qaux - ((u(2,j,blk)-u(1,j,blk))/gr_dx)*dydz
              end do 
           end if
        end do 

        call MPI_CollectResiduals(Qaux,ins_Qin,SUM_DATA)

      else

        do blk=1,blockCount           
           if ( xLC(blk) == nblockx-1) then
              do j=2,Nyb+1
                 Qaux = Qaux + u(Nxb+1,j,blk)*dydz
                 !Qaux = Qaux + ((u(Nxb+1,j,blk)-u(Nxb,j,blk))/gr_dy)*dydz
              end do
           end if
        end do

        call MPI_CollectResiduals(Qaux,ins_Qout,SUM_DATA)

      end if
      
end subroutine ins_computeQinout
