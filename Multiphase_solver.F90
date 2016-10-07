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
    real, pointer, dimension(:,:) :: s,pf,thco,cprs
    integer :: j,i

    ! _________Analaytical position of interface - Stefan Problem____!

    !solnX = 0.001-(2*0.066916063714766546815307639566313*&
    !                sqrt(((tstep)*dr_dt)*((mph_thco1)/(mph_cp1))))

    !______________________________End_______________________________!

    s => ph_center(DFUN_VAR,:,:)
    pf => ph_center(PFUN_VAR,:,:)
    thco => ph_center(THCO_VAR,:,:)
    cprs => ph_center(CPRS_VAR,:,:) 

    !_________Distance Function calculation - Stefan Problem_________!
 
    !do j=1,Nyb+2

    !  if(j==1) then     
    !      ycell = gr_y(1,j) - 0.5*gr_dy

    !  else if(j==Nyb+2) then
    !      ycell = gr_y(1,Nyb+1) + 0.5*gr_dy

    !  else
    !      ycell = 0.5*(gr_y(1,j) + gr_y(1,j-1))

    !  end if

    !  s(:,j) = ycell - solnX

    !end do

    !______________________________End_______________________________!

   call mph_FillVars(s,pf,thco,cprs,mph_thco1,mph_thco2,mph_cp1,mph_cp2)

   call MPI_applyBC(s)
   call MPI_applyBC(pf)
   call MPI_applyBC(thco)
   call MPI_applyBC(cprs)

   call MPI_physicalBC_dfun(s)
   call MPI_physicalBC_dfun(pf)
   call MPI_physicalBC_dfun(thco)
   call MPI_physicalBC_dfun(cprs)
    
   nullify(s)
   nullify(pf)
   nullify(thco)
   nullify(cprs)


end subroutine Multiphase_solver

