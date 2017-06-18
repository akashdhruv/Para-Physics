subroutine Solver_evolve    

#include "Solver.h"

!#define SOLVER_DEBUG

    use IncompNS_interface, only: IncompNS_solver
    use HeatAD_interface, only: HeatAD_solver
    use Grid_data
    use Driver_data
    use physicaldata
    use IncompNS_data
    use HeatAD_data
    use MPI_data
    use IO_interface, only: IO_display, IO_write, IO_display_v2
    use Multiphase_interface, only:Multiphase_solver
    use IBM_interface, only: IBM_solver
    use Driver_interface, only:Driver_init

    implicit none

    integer :: tstep,p_counter

    !real, allocatable,  dimension(:,:) :: uu,vv,pp,tt,df,pf,th,cp,ww
    real, dimension(Nxb+1,Nyb+1) :: uu,vv,pp,tt,ww,rx,ry,vs,ax,ay

    real, pointer, dimension(:,:,:) :: facexData,faceyData
    real, pointer, dimension(:,:,:) :: solnData

    real :: solnX
 
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

       !if((ins_u_res .lt. 0.0000001) .and. (ins_u_res .ne. 0).and. (ins_v_res .lt. 0.0000001) .and. (ins_v_res .ne. 0) ) exit

#ifdef SOLVER_DEBUG

       if(mod(tstep,10) == 0) then
   
          facexData => ph_facex
          faceyData => ph_facey
          solnData => ph_center

          uu = (facexData(VELC_VAR,1:Nxb+1,1:Nyb+1)+facexData(VELC_VAR,1:Nxb+1,2:Nyb+2))/2 

          vv = (faceyData(VELC_VAR,1:Nxb+1,1:Nyb+1)+faceyData(VELC_VAR,2:Nxb+2,1:Nyb+1))/2

          pp = ((solnData(PRES_VAR,1:Nxb+1,1:Nyb+1)+solnData(PRES_VAR,2:Nxb+2,1:Nyb+1))/2 &
               +(solnData(PRES_VAR,1:Nxb+1,2:Nyb+2)+solnData(PRES_VAR,2:Nxb+2,2:Nyb+2))/2)/2

          tt = ((solnData(TEMP_VAR,1:Nxb+1,1:Nyb+1)+solnData(TEMP_VAR,2:Nxb+2,1:Nyb+1))/2 &
               +(solnData(TEMP_VAR,1:Nxb+1,2:Nyb+2)+solnData(TEMP_VAR,2:Nxb+2,2:Nyb+2))/2)/2

          nullify(solnData)
          nullify(facexData)
          nullify(faceyData)

          call IO_write(gr_x,gr_y,uu,vv,pp,tt,myid)

        end if
#endif
        
        tstep     = tstep +1
        dr_dt_old = dr_dt

    end do

    facexData => ph_facex
    faceyData => ph_facey
    solnData => ph_center

    uu = (facexData(VELC_VAR,1:Nxb+1,1:Nyb+1)+facexData(VELC_VAR,1:Nxb+1,2:Nyb+2))/2

    vv = (faceyData(VELC_VAR,1:Nxb+1,1:Nyb+1)+faceyData(VELC_VAR,2:Nxb+2,1:Nyb+1))/2

    pp = ((solnData(PRES_VAR,1:Nxb+1,1:Nyb+1)+solnData(PRES_VAR,2:Nxb+2,1:Nyb+1))/2 + &
          (solnData(PRES_VAR,1:Nxb+1,2:Nyb+2)+solnData(PRES_VAR,2:Nxb+2,2:Nyb+2))/2)/2

    tt = ((solnData(TEMP_VAR,1:Nxb+1,1:Nyb+1)+solnData(TEMP_VAR,2:Nxb+2,1:Nyb+1))/2 + &
          (solnData(TEMP_VAR,1:Nxb+1,2:Nyb+2)+solnData(TEMP_VAR,2:Nxb+2,2:Nyb+2))/2)/2

    ww = ((solnData(VORT_VAR,1:Nxb+1,1:Nyb+1)+solnData(VORT_VAR,2:Nxb+2,1:Nyb+1))/2 + &
          (solnData(VORT_VAR,1:Nxb+1,2:Nyb+2)+solnData(VORT_VAR,2:Nxb+2,2:Nyb+2))/2)/2

    nullify(facexData)
    nullify(faceyData)
    nullify(solnData)

    call IO_write(gr_x,gr_y,uu,vv,pp,tt,myid)

end subroutine Solver_evolve
