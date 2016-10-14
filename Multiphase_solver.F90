subroutine Multiphase_solver(tstep,solnX)

#include "Solver.h"

    use Driver_data,  only: dr_dt
    use Multiphase_data, only: mph_thco1, mph_cp1, mph_thco2, mph_cp2
    use physicaldata
    use Grid_data
    use Multiphase_interface , only: mph_FillVars
    use MPI_interface, only: MPI_applyBC,MPI_physicalBC_dfun

    implicit none

    integer, intent(in) :: tstep
    

    real,intent(out) :: solnX
    real :: ycell
    real, pointer, dimension(:,:,:) :: solnData
    integer :: j,i

    ! _________Analaytical position of interface - Stefan Problem____!

    !solnX = 0.001-(2*0.066916063714766546815307639566313*&
    !                sqrt(((tstep)*dr_dt)*((mph_thco1)/(mph_cp1))))

    !______________________________End_______________________________!

    solnData => ph_center

    !_________Distance Function calculation - Stefan Problem_________!
 
    !do j=1,Nyb+2

    !  if(j==1) then     
    !      ycell = gr_y(1,j) - 0.5*gr_dy

    !  else if(j==Nyb+2) then
    !      ycell = gr_y(1,Nyb+1) + 0.5*gr_dy

    !  else
    !      ycell = 0.5*(gr_y(1,j) + gr_y(1,j-1))

    !  end if

    !  solnData(DFUN_VAR,:,j) = ycell - solnX

    !end do

    !______________________________End_______________________________!

    call mph_FillVars(solnData(DFUN_VAR,:,:),solnData(PFUN_VAR,:,:),&
                      solnData(THCO_VAR,:,:),solnData(CPRS_VAR,:,:),&
                      mph_thco1,mph_thco2,mph_cp1,mph_cp2)   

    nullify(solnData)  
 
end subroutine Multiphase_solver

