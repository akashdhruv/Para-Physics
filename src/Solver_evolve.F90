subroutine Solver_evolve    

#include "Solver.h"

!#define SOLVER_DEBUG

    use IncompNS_interface, only: IncompNS_solver
    use HeatAD_interface, only: HeatAD_solver
    use Grid_data
    use Driver_data
    use physicaldata, only: localCENTER,localFACEX,localFACEY
    use IncompNS_data
    use HeatAD_data
    use MPI_data
    use IO_interface, only: IO_display, IO_write, IO_display_v2
    use Multiphase_interface, only:Multiphase_solver
    use IBM_interface, only: IBM_solver
    use Driver_interface, only:Driver_init
    use Poisson_interface, only: Poisson_test

    implicit none

    integer :: tstep,p_counter,blk
    real :: ext_mean, num_mean, error_min, error_max
    real, dimension(Nxb+1,Nyb+1) :: uu,vv,pp,tt,ww,rr
    real :: solnX
    real, pointer, dimension(:,:,:,:) :: solnData,facexData,faceyData 

    tstep     = 1
    dr_dt_old = gr_dx

    do while(tstep<=dr_nt) 

       !call Multiphase_solver(tstep,solnX,.FALSE.)
       !call HeatAD_solver(tstep)
       !call Multiphase_solver(tstep,solnX,.TRUE.)
       !call IncompNS_solver(tstep,p_counter)
       !call IBM_solver(tstep)
       call Poisson_test(tstep,p_counter,ext_mean,num_mean,error_min,error_max)

       !-------Head Up Display-------!

       if (mod(tstep,1) == 0 .and. myid == 0) then

       !call IO_display(ins_u_res,ins_v_res,ins_w_res,ins_p_res,ht_T_res,p_counter,tstep*dr_dt,ins_maxdiv,ins_mindiv,&
       !                ins_umaxmin,ins_vmaxmin)

       !print *,"Convective velocity: ",ins_convvel
       !print *,"ins_Qin: ",ins_Qin
       !print *,"ins_Qout: ",ins_Qout
       !print *,"Qinout: ",ins_Qinout

       print *,"Poisson Iterations: ",p_counter
       print *,"Poisson residual  : ",ins_p_res
       print *,"Analytical mean   : ",ext_mean
       print *,"Numerical  mean   : ",num_mean
       print *,"Min error: ",error_min,"Max error: ",error_max
       print '( "Poisson time      :   ",f20.10," seconds")',ins_timePoisson
       end if

       !if((ins_u_res .lt. 0.0000001) .and. (ins_u_res .ne. 0).and. (ins_v_res .lt. 0.0000001) .and. (ins_v_res .ne. 0) ) exit


       !----------Write data------!

       if (mod(tstep,1) == 0) then

       solnData  => localCENTER
       facexData => localFACEX
       faceyData => localFACEY

       do blk=1,blockCount

          uu = (facexData(1:Nxb+1,1:Nyb+1,VELC_VAR,blk)+facexData(1:Nxb+1,2:Nyb+2,VELC_VAR,blk))/2 

          vv = (faceyData(1:Nxb+1,1:Nyb+1,VELC_VAR,blk)+faceyData(2:Nxb+2,1:Nyb+1,VELC_VAR,blk))/2

          pp = ((solnData(1:Nxb+1,1:Nyb+1,PRES_VAR,blk)+solnData(2:Nxb+2,1:Nyb+1,PRES_VAR,blk))/2 &
               +(solnData(1:Nxb+1,2:Nyb+2,PRES_VAR,blk)+solnData(2:Nxb+2,2:Nyb+2,PRES_VAR,blk))/2)/2

          tt = ((solnData(1:Nxb+1,1:Nyb+1,EXCT_VAR,blk)+solnData(2:Nxb+2,1:Nyb+1,EXCT_VAR,blk))/2 &
               +(solnData(1:Nxb+1,2:Nyb+2,EXCT_VAR,blk)+solnData(2:Nxb+2,2:Nyb+2,EXCT_VAR,blk))/2)/2

          ww = -((solnData(1:Nxb+1,1:Nyb+1,PRHS_VAR,blk)+solnData(2:Nxb+2,1:Nyb+1,PRHS_VAR,blk))/2 &
                +(solnData(1:Nxb+1,2:Nyb+2,PRHS_VAR,blk)+solnData(2:Nxb+2,2:Nyb+2,PRHS_VAR,blk))/2)/2

          rr = ((solnData(1:Nxb+1,1:Nyb+1,EROR_VAR,blk)+solnData(2:Nxb+2,1:Nyb+1,EROR_VAR,blk))/2 &
               +(solnData(1:Nxb+1,2:Nyb+2,EROR_VAR,blk)+solnData(2:Nxb+2,2:Nyb+2,EROR_VAR,blk))/2)/2

          !rr = (1./(facexData(1:Nxb+1,1:Nyb+1,RH1F_VAR,blk)+facexData(1:Nxb+1,2:Nyb+2,RH2F_VAR,blk)) + &
          !      1./(faceyData(1:Nxb+1,1:Nyb+1,RH1F_VAR,blk)+faceyData(2:Nxb+2,1:Nyb+1,RH2F_VAR,blk)))/2

        call IO_write(gr_x(:,:,blk),gr_y(:,:,blk),uu,vv,pp,tt,ww,rr,blk,blockOffset)

        end do

        nullify(solnData,facexData,faceyData)

        end if

        tstep     = tstep +1
        dr_dt_old = dr_dt

    end do

end subroutine Solver_evolve
