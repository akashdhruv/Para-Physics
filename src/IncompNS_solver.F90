subroutine IncompNS_solver(tstep,p_counter)

#include "Solver.h"

    use IncompNS_interface, only: ins_predictor,ins_vorticity,ins_predictor_VD,&
                                  ins_convVelout,ins_rescaleVel,ins_computeQinout

    use physicaldata, only: localCENTER,localFACEX,localFACEY
    use MPI_data, only: blockCount,solver_comm,ierr
    use IncompNS_data, only: ins_v_res,ins_u_res,ins_maxdiv,ins_mindiv,ins_umaxmin,ins_vmaxmin,ins_w_res
    use Grid_data, only: gr_dx,gr_dy
    use Driver_data, only: dr_dt
    use Poisson_interface, only: Poisson_solver, Poisson_solver_VC
    use MPI_interface, ONLY: MPI_applyBC, MPI_CollectResiduals,MPI_physicalBC_vort,MPI_physicalBC_vel
    use IBM_interface, ONLY: IBM_ApplyForcing
            
    implicit none

    integer, intent(in) :: tstep
    integer, intent(out) :: p_counter

    real, pointer, dimension(:,:,:,:) :: solnData,facexData,faceyData
    integer :: blk,i,j
    real :: v_res1,u_res1,maxdiv,mindiv,umax,vmax,umin,vmin,w_res1
    real, dimension(Nxb+2,Nyb+2) :: rhox,rhoy

    solnData  => localCENTER
    facexData => localFACEX
    faceyData => localFACEY

    ins_v_res = 0
    ins_u_res = 0
    ins_w_res = 0

    v_res1 = 0
    u_res1 = 0

    facexData(:,:,UOLD_VAR,:) = facexData(:,:,VELC_VAR,:)
    faceyData(:,:,UOLD_VAR,:) = faceyData(:,:,VELC_VAR,:)
    solnData(:,:,VORO_VAR,:)  = solnData(:,:,VORT_VAR,:)

    ! Get Qin and convective velocity

    call ins_computeQinout(facexData(:,:,VELC_VAR,:),.true.)

    call ins_convVelout(facexData(:,:,VELC_VAR,:))

    ! Predictor Step

#ifdef SINGLEPHASE
    do blk=1,blockCount
       call ins_predictor(tstep,facexData(:,:,VELC_VAR,blk),faceyData(:,:,VELC_VAR,blk),&
                                facexData(:,:,USTR_VAR,blk),faceyData(:,:,USTR_VAR,blk),&
                                facexData(:,:,GOLD_VAR,blk),faceyData(:,:,GOLD_VAR,blk))
    end do
#else
    do blk=1,blockCount
       call ins_predictor_VD(tstep,facexData(:,:,VELC_VAR,blk),faceyData(:,:,VELC_VAR,blk),&
                                   facexData(:,:,USTR_VAR,blk),faceyData(:,:,USTR_VAR,blk),&
                                   facexData(:,:,GOLD_VAR,blk),faceyData(:,:,GOLD_VAR,blk),&
                                   solnData(:,:,VISC_VAR,blk),&
                                   facexData(:,:,RH1F_VAR,blk),faceyData(:,:,RH1F_VAR,blk),&
                                   facexData(:,:,RH2F_VAR,blk),faceyData(:,:,RH2F_VAR,blk))
    end do
#endif

    ! Immersed Boundary Calculation

#ifdef IBM
    do blk=1,blockCount
       !call IBM_ApplyForcing(facexData(:,:,USTR_VAR,blk),faceyData(:,:,USTR_VAR,blk),&
       !                      facexData(:,:,IBMF_VAR,blk),faceyData(:,:,IBMF_VAR,blk))

       call IBM_ApplyForcing(facexData(:,:,USTR_VAR,blk),faceyData(:,:,USTR_VAR,blk),&
                             solnData(:,:,DFUN_VAR,blk),solnData(:,:,PFUN_VAR,blk))
    end do
#endif    

    ! MPI and Domain Predictor BCs

    call MPI_BARRIER(solver_comm,ierr)
    call MPI_applyBC(USTR_VAR,FACEX)
    call MPI_applyBC(USTR_VAR,FACEY)
    call MPI_physicalBC_vel(facexData(:,:,USTR_VAR,:),faceyData(:,:,USTR_VAR,:),solnData(:,:,PFUN_VAR,:))

    ! Compute Qout
  
    call ins_computeQinout(facexData(:,:,VELC_VAR,:),.false.)

    ! Rescale velocity

    call ins_rescaleVel(facexData(:,:,VELC_VAR,:))

    ! Poisson RHS

#ifdef SINGLEPHASE
    do blk=1,blockCount
       solnData(2:Nxb+1,2:Nyb+1,PRHS_VAR,blk) = &
       -((1/(gr_dy*dr_dt))*(faceyData(2:Nxb+1,2:Nyb+1,USTR_VAR,blk)-faceyData(2:Nxb+1,1:Nyb,USTR_VAR,blk)))&
       -((1/(gr_dx*dr_dt))*(facexData(2:Nxb+1,2:Nyb+1,USTR_VAR,blk)-facexData(1:Nxb,2:Nyb+1,USTR_VAR,blk)))
    end do

#else
    do blk=1,blockCount
       solnData(2:Nxb+1,2:Nyb+1,PRHS_VAR,blk) = &
       -((1/(gr_dy*dr_dt))*(faceyData(2:Nxb+1,2:Nyb+1,USTR_VAR,blk)-faceyData(2:Nxb+1,1:Nyb,USTR_VAR,blk)))&
       -((1/(gr_dx*dr_dt))*(facexData(2:Nxb+1,2:Nyb+1,USTR_VAR,blk)-facexData(1:Nxb,2:Nyb+1,USTR_VAR,blk)))&
       -solnData(2:Nxb+1,2:Nyb+1,SIGP_VAR,blk)
    end do

#endif

    nullify(solnData,facexData,faceyData)

   ! Poisson Solve

#ifdef SINGLEPHASE
    call Poisson_solver(PRHS_VAR,PRES_VAR,p_counter)
#else
    call Poisson_solver_VC(PRHS_VAR,PRES_VAR,p_counter,RH1F_VAR,RH2F_VAR)
#endif

    solnData  => localCENTER
    facexData => localFACEX
    faceyData => localFACEY

   ! Corrector Step

#ifdef SINGLEPHASE
    do blk=1,blockCount
    facexData(2:Nxb+1,2:Nyb+1,VELC_VAR,blk) = &
    facexData(2:Nxb+1,2:Nyb+1,USTR_VAR,blk) - (dr_dt/gr_dx)*(solnData(3:Nxb+2,2:Nyb+1,PRES_VAR,blk)-&
                                                             solnData(2:Nxb+1,2:Nyb+1,PRES_VAR,blk))


    faceyData(2:Nxb+1,2:Nyb+1,VELC_VAR,blk) = &
    faceyData(2:Nxb+1,2:Nyb+1,USTR_VAR,blk) - (dr_dt/gr_dy)*(solnData(2:Nxb+1,3:Nyb+2,PRES_VAR,blk)-&
                                                             solnData(2:Nxb+1,2:Nyb+1,PRES_VAR,blk))
    end do
#else

    do blk=1,blockCount
    facexData(2:Nxb+1,2:Nyb+1,VELC_VAR,blk) = &
    facexData(2:Nxb+1,2:Nyb+1,USTR_VAR,blk) - (dr_dt/gr_dx)*&
                                              (facexData(2:Nxb+1,2:Nyb+1,RH1F_VAR,blk)+&
                                               facexData(2:Nxb+1,2:Nyb+1,RH2F_VAR,blk))*&
                                              (solnData(3:Nxb+2,2:Nyb+1,PRES_VAR,blk)-&
                                               solnData(2:Nxb+1,2:Nyb+1,PRES_VAR,blk))&
                                               +dr_dt*&
                                               facexData(2:Nxb+1,2:Nyb+1,SIGM_VAR,blk)


    faceyData(2:Nxb+1,2:Nyb+1,VELC_VAR,blk) = &
    faceyData(2:Nxb+1,2:Nyb+1,USTR_VAR,blk) - (dr_dt/gr_dy)*&
                                              (faceyData(2:Nxb+1,2:Nyb+1,RH1F_VAR,blk)+&
                                               faceyData(2:Nxb+1,2:Nyb+1,RH2F_VAR,blk))*&                                               
                                              (solnData(2:Nxb+1,3:Nyb+2,PRES_VAR,blk)-&
                                               solnData(2:Nxb+1,2:Nyb+1,PRES_VAR,blk))&
                                               +dr_dt*&
                                                faceyData(2:Nxb+1,2:Nyb+1,SIGM_VAR,blk)
    end do
#endif

    ! MPI and Domain Corrector BCs

    call MPI_BARRIER(solver_comm,ierr)
    call MPI_applyBC(VELC_VAR,FACEX)
    call MPI_applyBC(VELC_VAR,FACEY)
    call MPI_physicalBC_vel(facexData(:,:,VELC_VAR,:),faceyData(:,:,VELC_VAR,:),solnData(:,:,PFUN_VAR,:))

    ! Divergence

    maxdiv = -10.**(10.)
    mindiv = 10.**(10.)

    do blk=1,blockCount
    maxdiv = max(maxdiv,maxval(((1/(gr_dy))*(faceyData(2:Nxb+1,2:Nyb+1,VELC_VAR,blk)-faceyData(2:Nxb+1,1:Nyb,VELC_VAR,blk)))&
                              +((1/(gr_dx))*(facexData(2:Nxb+1,2:Nyb+1,VELC_VAR,blk)-facexData(1:Nxb,2:Nyb+1,VELC_VAR,blk)))))

    mindiv = min(mindiv,minval(((1/(gr_dy))*(faceyData(2:Nxb+1,2:Nyb+1,VELC_VAR,blk)-faceyData(2:Nxb+1,1:Nyb,VELC_VAR,blk)))&
                              +((1/(gr_dx))*(facexData(2:Nxb+1,2:Nyb+1,VELC_VAR,blk)-facexData(1:Nxb,2:Nyb+1,VELC_VAR,blk)))))

    end do

    ! Max-Min Velocity

    umax = maxval(facexData(2:Nxb+1,2:Nyb+1,VELC_VAR,:))
    umin = minval(facexData(2:Nxb+1,2:Nyb+1,VELC_VAR,:))

    vmax = maxval(faceyData(2:Nxb+1,2:Nyb+1,VELC_VAR,:))
    vmin = minval(faceyData(2:Nxb+1,2:Nyb+1,VELC_VAR,:))

    call MPI_CollectResiduals(maxdiv,ins_maxdiv,MAX_DATA)
    call MPI_CollectResiduals(mindiv,ins_mindiv,MIN_DATA)

    call MPI_CollectResiduals(umax,ins_umaxmin(1),MAX_DATA)
    call MPI_CollectResiduals(umin,ins_umaxmin(2),MIN_DATA)

    call MPI_CollectResiduals(vmax,ins_vmaxmin(1),MAX_DATA)
    call MPI_CollectResiduals(vmin,ins_vmaxmin(2),MIN_DATA)

    ! Residuals

    do blk=1,blockCount
    do i=1,Nyb+2
          ins_u_res = ins_u_res + sum((facexData(:,i,VELC_VAR,blk)-facexData(:,i,UOLD_VAR,blk))**2)
          ins_v_res = ins_v_res + sum((faceyData(:,i,VELC_VAR,blk)-faceyData(:,i,UOLD_VAR,blk))**2)
    enddo
    enddo

    call MPI_CollectResiduals(ins_u_res,u_res1,SUM_DATA)
    ins_u_res = sqrt(u_res1/((nblockx*nblocky)*(Nxb+2)*(Nyb+2)))

    call MPI_CollectResiduals(ins_v_res,v_res1,SUM_DATA)
    ins_v_res = sqrt(v_res1/((nblockx*nblocky)*(Nxb+2)*(Nyb+2)))

    ! Vorticity Calculation

    do blk=1,blockCount
  
    solnData(2:Nxb+1,2:Nyb+1,VORT_VAR,blk) = (((faceyData(2:Nxb+1,1:Nyb,VELC_VAR,blk)-faceyData(1:Nxb,1:Nyb,VELC_VAR,blk))/gr_dx &
                                            -  (facexData(1:Nxb,2:Nyb+1,VELC_VAR,blk)-facexData(1:Nxb,1:Nyb,VELC_VAR,blk))/gr_dy)&

                                            + ((faceyData(3:Nxb+2,1:Nyb,VELC_VAR,blk)-faceyData(2:Nxb+1,1:Nyb,VELC_VAR,blk))/gr_dx & 
                                            -  (facexData(2:Nxb+1,2:Nyb+1,VELC_VAR,blk)-facexData(2:Nxb+1,1:Nyb,VELC_VAR,blk))/gr_dy)&

                                            + ((faceyData(2:Nxb+1,2:Nyb+1,VELC_VAR,blk)-faceyData(1:Nxb,2:Nyb+1,VELC_VAR,blk))/gr_dx &
                                            -  (facexData(1:Nxb,3:Nyb+2,VELC_VAR,blk)-facexData(1:Nxb,2:Nyb+1,VELC_VAR,blk))/gr_dy) &

                                      + ((faceyData(3:Nxb+2,2:Nyb+1,VELC_VAR,blk)-faceyData(2:Nxb+1,2:Nyb+1,VELC_VAR,blk))/gr_dx)&
                                      -  (facexData(2:Nxb+1,3:Nyb+2,VELC_VAR,blk)-facexData(2:Nxb+1,2:Nyb+1,VELC_VAR,blk))/gr_dy)/4

    end do

    ! Vorticity BC

    call MPI_BARRIER(solver_comm,ierr)
    call MPI_applyBC(VORT_VAR,CENTER)
    call MPI_physicalBC_vort(solnData(:,:,VORT_VAR,:))

    !do blk=1,blockCount
    !call ins_vorticity(tstep,solnData(:,:,VORT_VAR,blk),solnData(:,:,VORO_VAR,blk),&
    !                   facexData(:,:,VELC_VAR,blk),faceyData(:,:,VELC_VAR,blk),&
    !                   solnData(:,:,DFUN_VAR,blk))
    !end do

    !! Vorticity BC

    !call MPI_BARRIER(solver_comm,ierr)
    !call MPI_applyBC(VORT_VAR,CENTER)
    !call MPI_physicalBC_vort(solnData(:,:,VORT_VAR,:))

    !do blk=1,blockCount
    !do i=1,Nyb+2
    !      ins_w_res = ins_w_res + sum((solnData(:,i,VORT_VAR,blk)-solnData(:,i,VORO_VAR,blk))**2)
    !enddo
    !enddo

    !call MPI_CollectResiduals(ins_w_res,w_res1,SUM_DATA)
    !ins_w_res = sqrt(w_res1/((nblockx*nblocky)*(Nxb+2)*(Nyb+2)))

    nullify(solnData,facexData,faceyData)

end subroutine
