subroutine ins_computeQinout(u,v,x,y,flg)

#include"Solver.h"

      use IncompNS_data
      use physicaldata
      use Grid_data
      use MPI_data
      use IBM_data

      implicit none
      
      real, dimension(:,:,:), intent(in) :: u,v,x,y
      logical, intent(in) :: flg

      real :: Qaux
      integer :: blk,j,i

      real :: dxdz, dydz

      Qaux = 0.0
      dxdz = gr_dx
      dydz = gr_dy


#ifdef HOME_HEATING_SYSTEM
      if(flg) then
       
        do blk=1,blockCount

           if ( yLC(blk) == nblocky-1 ) then
              do i=1,Nxb+2
                if(x(i,Nyb+1,blk) .ge. ins_upIn1 .and. x(i,Nyb+1,blk) .le. ins_upIn2) Qaux = Qaux - v(i,Nyb+1,blk)*dxdz
              end do
           end if

           do j=1,Nyb+1
              do i=1,Nxb+2

                if(y(i,j,blk) .le. ibm_yl .and. y(i,j+1,blk) .ge. ibm_yl .and. &
                   x(i,j,blk) .ge. ins_dnIn1 .and. x(i,j,blk) .le. ins_dnIn2) Qaux = Qaux - v(i,j,blk)*dxdz

              end do
           end do


          call MPI_CollectResiduals(Qaux,ins_Qin,SUM_DATA)

          ins_Qin = abs(ins_Qin)

        end do

      else

        do blk=1,blockCount
           if ( yLC(blk) == nblocky-1 ) then
              do i=1,Nxb+2
                if(x(i,Nyb+1,blk) .ge. ins_upEx1 .and. x(i,Nyb+1,blk) .le. ins_upEx2) Qaux = Qaux + v(i,Nyb+1,blk)*dxdz
              end do
           end if

           if ( xLC(blk) == nblockx-1 ) then
              do j=2,Nyb+1
                if(y(Nxb+1,j,blk) .ge. ins_dnEx1 .and. y(Nxb+1,j,blk) .le. ins_dnEx2) Qaux = Qaux + u(Nxb+1,j,blk)*dydz
              end do                       
           end if 

           call MPI_CollectResiduals(Qaux,ins_Qout,SUM_DATA)
        end do
      end if

#else
      if(flg) then

        do blk=1,blockCount

           if ( xLC(blk) == 0 .and. ins_xl_bnd .ne. OUTFLOW .and. ins_xl_bnd .ne. NEUMANN) then
              do j=2,Nyb+1                
                 Qaux = Qaux + u(1,j,blk)*dydz
              end do 
           end if

           if ( xLC(blk) == nblockx-1 .and. ins_xr_bnd .ne. OUTFLOW .and. ins_xr_bnd .ne. NEUMANN) then
              do j=2,Nyb+1
                 Qaux = Qaux - u(Nxb+1,j,blk)*dydz
              end do
           end if

           if ( yLC(blk) == 0 .and. ins_yl_bnd .ne. OUTFLOW .and. ins_yl_bnd .ne. NEUMANN) then
              do i=2,Nxb+1
                 Qaux = Qaux + v(i,1,blk)*dxdz
              end do
           end if

           if ( yLC(blk) == nblocky-1 .and. ins_yr_bnd .ne. OUTFLOW .and. ins_yr_bnd .ne. NEUMANN) then
              do i=2,Nxb+1
                Qaux = Qaux - v(i,Nyb+1,blk)*dxdz
              end do
           end if

        end do 

        call MPI_CollectResiduals(Qaux,ins_Qin,SUM_DATA)

     else

        do blk=1,blockCount

           if ( xLC(blk) == 0 .and. (ins_xl_bnd .eq. OUTFLOW .or. ins_xl_bnd .eq. NEUMANN)) then
              do j=2,Nyb+1                
                 Qaux = Qaux - u(1,j,blk)*dydz
              end do 
           end if

           if ( xLC(blk) == nblockx-1 .and. (ins_xr_bnd .eq. OUTFLOW .or. ins_xr_bnd .eq. NEUMANN)) then
              do j=2,Nyb+1
                 Qaux = Qaux + u(Nxb+1,j,blk)*dydz
              end do
           end if

           if ( yLC(blk) == 0 .and. (ins_yl_bnd .eq. OUTFLOW .or. ins_yl_bnd .eq. NEUMANN)) then
              do i=2,Nxb+1
                 Qaux = Qaux - v(i,1,blk)*dxdz
              end do
           end if

           if ( yLC(blk) == nblocky-1 .and. (ins_yr_bnd .eq. OUTFLOW .or. ins_yr_bnd .eq. NEUMANN)) then
              do i=2,Nxb+1
                Qaux = Qaux + v(i,Nyb+1,blk)*dxdz
              end do
           end if

        end do 

        call MPI_CollectResiduals(Qaux,ins_Qout,SUM_DATA)

     end if

#endif

end subroutine ins_computeQinout
