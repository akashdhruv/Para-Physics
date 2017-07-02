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

    implicit none

    integer :: tstep,p_counter

    real, dimension(Nxb+1,Nyb+1) :: uu,vv,pp,tt,ww,rr

    real :: solnX
  
    real, pointer, dimension(:,:,:,:) :: solnData,facexData,faceyData 

    tstep     = 0
    dr_dt_old = gr_dx

    do while(tstep<=dr_nt) 

#ifdef MULTIPHASE
       call Multiphase_solver(tstep,solnX,.FALSE.)
#endif

#ifdef ENERGY
       call HeatAD_solver(tstep)
#endif

#ifdef MULTIPHASE
       call Multiphase_solver(tstep,solnX,.TRUE.)
#endif

#ifdef INS
       call IncompNS_solver(tstep,p_counter)
#endif

#ifdef IBM
       call IBM_solver(tstep)
#endif

#ifdef ONLY_POISSON
! call FOR Poisson_analytical GOES HERE
#endif

       if (mod(tstep,20) == 0 .and. myid == 0) then
          call IO_display(ins_u_res,ins_v_res,ins_w_res,ins_p_res,ht_T_res,p_counter,tstep*dr_dt,ins_maxdiv,ins_mindiv,&
                          ins_umaxmin,ins_vmaxmin)
       endif

       if((ins_u_res .lt. 0.0000001) .and. (ins_u_res .ne. 0).and. (ins_v_res .lt. 0.0000001) .and. (ins_v_res .ne. 0) ) exit

#ifdef SOLVER_DEBUG

       if(mod(tstep,10000) == 0) then

          solnData  => localCENTER
          facexData => localFACEX
          faceyData => localFACEY     

          uu = (facexData(1:Nxb+1,1:Nyb+1,blockCount,VELC_VAR)+facexData(1:Nxb+1,2:Nyb+2,blockCount,VELC_VAR))/2 

          vv = (faceyData(1:Nxb+1,1:Nyb+1,blockCount,VELC_VAR)+faceyData(2:Nxb+2,1:Nyb+1,blockCount,VELC_VAR))/2

          pp = ((solnData(1:Nxb+1,1:Nyb+1,blockCount,PRES_VAR)+solnData(2:Nxb+2,1:Nyb+1,blockCount,PRES_VAR))/2 &
               +(solnData(1:Nxb+1,2:Nyb+2,blockCount,PRES_VAR)+solnData(2:Nxb+2,2:Nyb+2,blockCount,PRES_VAR))/2)/2

          tt = ((solnData(1:Nxb+1,1:Nyb+1,blockCount,TEMP_VAR)+solnData(2:Nxb+2,1:Nyb+1,blockCount,TEMP_VAR))/2 &
               +(solnData(1:Nxb+1,2:Nyb+2,blockCount,TEMP_VAR)+solnData(2:Nxb+2,2:Nyb+2,blockCount,TEMP_VAR))/2)/2

          ww = ((solnData(1:Nxb+1,1:Nyb+1,blockCount,VISC_VAR)+solnData(2:Nxb+2,1:Nyb+1,blockCount,VISC_VAR))/2 &
               +(solnData(1:Nxb+1,2:Nyb+2,blockCount,VISC_VAR)+solnData(2:Nxb+2,2:Nyb+2,blockCount,VISC_VAR))/2)/2

          rr = ((facexData(1:Nxb+1,1:Nyb+1,blockCount,RH1F_VAR)+facexData(1:Nxb+1,2:Nyb+2,blockCount,RH2F_VAR)) + &
                (faceyData(1:Nxb+1,1:Nyb+1,blockCount,RH1F_VAR)+faceyData(2:Nxb+2,1:Nyb+1,blockCount,RH2F_VAR)))/2

         call IO_write(gr_x,gr_y,uu,vv,pp,tt,ww,rr,myid)

         nullify(solnData,facexData,faceyData)

        end if
#endif
        
        tstep     = tstep +1
        dr_dt_old = dr_dt

    end do

    solnData  => localCENTER
    facexData => localFACEX
    faceyData => localFACEY

    uu = (facexData(1:Nxb+1,1:Nyb+1,blockCount,VELC_VAR)+facexData(1:Nxb+1,2:Nyb+2,blockCount,VELC_VAR))/2 

    vv = (faceyData(1:Nxb+1,1:Nyb+1,blockCount,VELC_VAR)+faceyData(2:Nxb+2,1:Nyb+1,blockCount,VELC_VAR))/2

    pp = ((solnData(1:Nxb+1,1:Nyb+1,blockCount,PRES_VAR)+solnData(2:Nxb+2,1:Nyb+1,blockCount,PRES_VAR))/2 &
         +(solnData(1:Nxb+1,2:Nyb+2,blockCount,PRES_VAR)+solnData(2:Nxb+2,2:Nyb+2,blockCount,PRES_VAR))/2)/2

    tt = ((solnData(1:Nxb+1,1:Nyb+1,blockCount,TEMP_VAR)+solnData(2:Nxb+2,1:Nyb+1,blockCount,TEMP_VAR))/2 &
         +(solnData(1:Nxb+1,2:Nyb+2,blockCount,TEMP_VAR)+solnData(2:Nxb+2,2:Nyb+2,blockCount,TEMP_VAR))/2)/2

    ww = ((solnData(1:Nxb+1,1:Nyb+1,blockCount,VISC_VAR)+solnData(2:Nxb+2,1:Nyb+1,blockCount,VISC_VAR))/2 &
         +(solnData(1:Nxb+1,2:Nyb+2,blockCount,VISC_VAR)+solnData(2:Nxb+2,2:Nyb+2,blockCount,VISC_VAR))/2)/2

    rr = ((facexData(1:Nxb+1,1:Nyb+1,blockCount,RH1F_VAR)+facexData(1:Nxb+1,2:Nyb+2,blockCount,RH2F_VAR)) + &
          (faceyData(1:Nxb+1,1:Nyb+1,blockCount,RH1F_VAR)+faceyData(2:Nxb+2,1:Nyb+1,blockCount,RH2F_VAR)))/2

    call IO_write(gr_x,gr_y,uu,vv,pp,tt,ww,rr,myid)

    nullify(solnData,facexData,faceyData)

end subroutine Solver_evolve
