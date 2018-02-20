subroutine ins_rescaleVel(u,v,x,y)

#include "Solver.h"

    use IncompNS_data
    use Grid_data
    use MPI_data

    implicit none

    real,dimension(:,:,:),intent(inout) :: u,v
    real,dimension(:,:,:),intent(in) :: x,y
   
    real :: tol=1e-10
    integer :: blk,i,j

    if(ins_Qout == 0) then 
      ins_Qinout = 1.0
    else
      ins_Qinout = ins_Qin/ins_Qout
    end if
  

#ifdef HOME_HEATING_SYSTEM
    do blk=1,blockCount

       if ( xLC(blk) == nblockx-1) then
          do j=2,Nyb+1
             if(y(Nxb+1,j,blk) .ge. ins_dnEx1 .and. y(Nxb+1,j,blk) .le.ins_dnEx2) then
                u(Nxb+1,j,blk) = u(Nxb+1,j,blk)*ins_Qinout
             end if
          end do
       end if

       if ( yLC(blk) == nblocky-1) then
          do j=1,Nxb+2
             if(x(j,Nyb+1,blk) .ge. ins_upEx1 .and. x(j,Nyb+1,blk) .le. ins_upEx2) then
                v(j,Nyb+1,blk) = v(j,Nyb+1,blk)*ins_Qinout
             end if
          end do
       end if

      end do

#else 
    do blk=1,blockCount

       if ( xLC(blk) == 0 .and. (ins_xl_bnd .eq. OUTFLOW .or. ins_xl_bnd .eq. NEUMANN)) then
          u(1,2:Nyb+1,blk) = u(1,2:Nyb+1,blk)*ins_Qinout
       end if

       if ( xLC(blk) == nblockx-1 .and. (ins_xr_bnd .eq. OUTFLOW .or. ins_xr_bnd .eq. NEUMANN)) then
          u(Nxb+1,2:Nyb+1,blk) = u(Nxb+1,2:Nyb+1,blk)*ins_Qinout
       end if

       if( yLC(blk) == 0 .and. (ins_yl_bnd .eq. OUTFLOW .or. ins_yl_bnd .eq. NEUMANN)) then
         v(2:Nxb+1,1,blk) = v(2:Nxb+1,1,blk)*ins_Qinout
       end if

       if( yLC(blk) == nblocky-1 .and. (ins_yr_bnd .eq. OUTFLOW .or. ins_yr_bnd .eq. NEUMANN)) then
         v(2:Nxb+1,Nyb+1,blk) = v(2:Nxb+1,Nyb+1,blk)*ins_Qinout
       end if

    end do
 
#endif

end subroutine
