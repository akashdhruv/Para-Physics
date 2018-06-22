subroutine Poisson_test(tstep,p_counter,ext_mean,num_mean,error_min,error_max)

#include "Solver.h"

    use physicaldata, only: localCENTER,localFACEX,localFACEY
    use MPI_data, only: procs,blockCount,solver_comm,ierr
    use Grid_data, only: gr_x,gr_y,gr_dx,gr_dy
    use Poisson_interface, only: Poisson_solver, Poisson_solver_VC
    use MPI_interface, ONLY: MPI_CollectResiduals

    implicit none

    integer, intent(in) :: tstep
    integer, intent(out) :: p_counter
    real, intent(out) :: ext_mean,num_mean,error_min,error_max

    real, pointer, dimension(:,:,:,:) :: solnData,facexData,faceyData
    real :: xcell,ycell
    integer :: blk,i,j
    real, parameter :: pi = acos(-1.0)
    real :: ext_mean_aux, num_mean_aux, error_min_aux, error_max_aux

    solnData  => localCENTER
    facexData => localFACEX
    faceyData => localFACEY

    do blk=1,blockCount

        do j=1,Nyb+2

           if(j==1) then
           ycell = gr_y(1,j,blk) - 0.5*gr_dy

           else if(j==Nyb+2) then
           ycell = gr_y(1,Nyb+1,blk) + 0.5*gr_dy

           else
           ycell = 0.5*(gr_y(1,j,blk) + gr_y(1,j-1,blk))

           end if

           do i=1,Nxb+2

                if(i==1) then
                xcell = gr_x(i,1,blk) - 0.5*gr_dx
 
                else if(i==Nxb+2) then
                xcell = gr_x(Nxb+1,1,blk) + 0.5*gr_dx

                else
                xcell = 0.5*(gr_x(i,1,blk) + gr_x(i-1,1,blk))

                end if

                solnData(i,j,EXCT_VAR,blk) = cos(2*pi*xcell)*cos(2*pi*ycell)
                solnData(i,j,PRHS_VAR,blk) = 8*pi*pi*cos(2*pi*xcell)*cos(2*pi*ycell)

           end do
        end do
    end do

    nullify(facexData,faceyData,solnData)

    call Poisson_solver(PRHS_VAR,PRES_VAR,p_counter)

    solnData  => localCENTER

    solnData(:,:,EROR_VAR,:) = abs(solnData(:,:,EXCT_VAR,:)-solnData(:,:,PRES_VAR,:))

    error_min_aux = minval(solnData(:,:,EROR_VAR,:))
    error_max_aux = maxval(solnData(:,:,EROR_VAR,:))

    ext_mean_aux = 0.0
    num_mean_aux = 0.0

    do blk=1,blockCount

    ext_mean_aux = ext_mean_aux + sum(sum(solnData(:,:,EXCT_VAR,blk),1))/((Nxb+2)*(Nyb+2)*blockCount)
    num_mean_aux = num_mean_aux + sum(sum(solnData(:,:,PRES_VAR,blk),1))/((Nxb+2)*(Nyb+2)*blockCount)

    end do

    call MPI_CollectResiduals(error_min_aux,error_min,MIN_DATA)
    call MPI_CollectResiduals(error_max_aux,error_max,MAX_DATA)
    
    call MPI_CollectResiduals(ext_mean_aux,ext_mean,SUM_DATA)
    call MPI_CollectResiduals(num_mean_aux,num_mean,SUM_DATA)

    ext_mean = ext_mean/procs
    num_mean = num_mean/procs
    
    nullify(solnData)

end subroutine
