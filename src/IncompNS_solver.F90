subroutine IncompNS_solver(tstep,p_counter)

#include "Solver.h"

    use IncompNS_interface, only: ins_predictor,ins_vorticity,ins_predictor_VD
    use physicaldata, only: localCENTER,localFACEX,localFACEY
    use MPI_data, only: blockCount,shared_comm,ierr
    use IncompNS_data, only: ins_v_res,ins_u_res,ins_maxdiv,ins_mindiv,ins_umaxmin,ins_vmaxmin,ins_w_res
    use Grid_data, only: gr_dx,gr_dy
    use Driver_data, only: dr_dt
    use MPI_interface, ONLY: MPI_applyBC, MPI_CollectResiduals,MPI_physicalBC_vort,MPI_physicalBC_vel,&
                             MPI_applyBC_shared,MPI_applyBC_RMA

    implicit none

    integer, intent(in) :: tstep
    integer, intent(out) :: p_counter

    real, pointer, dimension(:,:,:,:) :: solnData,facexData,faceyData
    integer :: blk,i,j
    real :: v_res1,u_res1,maxdiv,mindiv,umax,vmax,umin,vmin,w_res1

    solnData  => localCENTER
    facexData => localFACEX
    faceyData => localFACEY

    ins_v_res = 0
    ins_u_res = 0
    ins_w_res = 0

    v_res1 = 0
    u_res1 = 0

    facexData(:,:,:,UOLD_VAR) = facexData(:,:,:,VELC_VAR)
    faceyData(:,:,:,UOLD_VAR) = faceyData(:,:,:,VELC_VAR)
    solnData(:,:,:,VORO_VAR)  = solnData(:,:,:,VORT_VAR)

    do blk=1,blockCount
       call ins_predictor(tstep,facexData(:,:,blk,VELC_VAR),faceyData(:,:,blk,VELC_VAR),&
                                facexData(:,:,blk,USTR_VAR),faceyData(:,:,blk,USTR_VAR),&
                                facexData(:,:,blk,GOLD_VAR),faceyData(:,:,blk,GOLD_VAR))
    end do

#ifdef MPI_DIST
    call MPI_applyBC(USTR_VAR,FACEX)
    call MPI_applyBC(USTR_VAR,FACEY)
#endif

#ifdef MPI_SHRD
    call MPI_BARRIER(shared_comm,ierr)
    call MPI_applyBC_shared(USTR_VAR,FACEX)
    call MPI_applyBC_shared(USTR_VAR,FACEY)
#endif

#ifdef MPI_RMA
    call MPI_applyBC_RMA(USTR_VAR,FACEX)
    call MPI_applyBC_RMA(USTR_VAR,FACEY)
#endif

    call MPI_physicalBC_vel(facexData(:,:,:,USTR_VAR),faceyData(:,:,:,USTR_VAR))

#ifdef IBM
    do blk=1,blockCount
       call IBM_ApplyForcing(facexData(:,:,blk,USTR_VAR),faceyData(:,:,blk,USTR_VAR),&
                             facexData(:,:,blk,IBMF_VAR),faceyData(:,:,blk,IBMF_VAR))
    end do
#endif    

    do blk=1,blockCount
       solnData(2:Nxb+1,2:Nyb+1,blk,PRHS_VAR) = &

       -((1/(gr_dy*dr_dt))*(faceyData(2:Nxb+1,2:Nyb+1,blk,USTR_VAR)-faceyData(2:Nxb+1,1:Nyb,blk,USTR_VAR)))&
       -((1/(gr_dx*dr_dt))*(facexData(2:Nxb+1,2:Nyb+1,blk,USTR_VAR)-facexData(1:Nxb,2:Nyb+1,blk,USTR_VAR)))
    end do

    nullify(solnData,facexData,faceyData)

    call Poisson_solver(PRHS_VAR,PRES_VAR,p_counter)

    solnData  => localCENTER
    facexData => localFACEX
    faceyData => localFACEY

    do blk=1,blockCount
    facexData(2:Nxb+1,2:Nyb+1,blk,VELC_VAR) = &
    facexData(2:Nxb+1,2:Nyb+1,blk,USTR_VAR) - (dr_dt/gr_dx)*(solnData(3:Nxb+2,2:Nyb+1,blk,PRES_VAR)-&
                                                             solnData(2:Nxb+1,2:Nyb+1,blk,PRES_VAR))


    faceyData(2:Nxb+1,2:Nyb+1,blk,VELC_VAR) = &
    faceyData(2:Nxb+1,2:Nyb+1,blk,USTR_VAR) - (dr_dt/gr_dy)*(solnData(2:Nxb+1,3:Nyb+2,blk,PRES_VAR)-&
                                                             solnData(2:Nxb+1,2:Nyb+1,blk,PRES_VAR))
    end do

#ifdef MPI_DIST
    call MPI_applyBC(VELC_VAR,FACEX)
    call MPI_applyBC(VELC_VAR,FACEY)
#endif

#ifdef MPI_SHRD
    call MPI_BARRIER(shared_comm,ierr)
    call MPI_applyBC_shared(VELC_VAR,FACEX)
    call MPI_applyBC_shared(VELC_VAR,FACEY)
#endif

#ifdef MPI_RMA
    call MPI_applyBC_RMA(VELC_VAR,FACEX)
    call MPI_applyBC_RMA(VELC_VAR,FACEY)
#endif

    call MPI_physicalBC_vel(facexData(:,:,:,VELC_VAR),faceyData(:,:,:,VELC_VAR))

    ! Divergence

    maxdiv = -10.**(10.)
    mindiv = 10.**(10.)

    do blk=1,blockCount
    maxdiv = max(maxdiv,maxval(((1/(gr_dy))*(faceyData(2:Nxb+1,2:Nyb+1,blk,VELC_VAR)-faceyData(2:Nxb+1,1:Nyb,blk,VELC_VAR)))&
                              +((1/(gr_dx))*(facexData(2:Nxb+1,2:Nyb+1,blk,VELC_VAR)-facexData(1:Nxb,2:Nyb+1,blk,VELC_VAR)))))

    mindiv = min(mindiv,minval(((1/(gr_dy))*(faceyData(2:Nxb+1,2:Nyb+1,blk,VELC_VAR)-faceyData(2:Nxb+1,1:Nyb,blk,VELC_VAR)))&
                              +((1/(gr_dx))*(facexData(2:Nxb+1,2:Nyb+1,blk,VELC_VAR)-facexData(1:Nxb,2:Nyb+1,blk,VELC_VAR)))))

    end do

    umax = maxval(facexData(:,:,:,VELC_VAR))
    umin = minval(facexData(:,:,:,VELC_VAR))

    vmax = maxval(faceyData(:,:,:,VELC_VAR))
    vmin = minval(faceyData(:,:,:,VELC_VAR))

    call MPI_CollectResiduals(maxdiv,ins_maxdiv,MAX_DATA)
    call MPI_CollectResiduals(mindiv,ins_mindiv,MIN_DATA)

    call MPI_CollectResiduals(umax,ins_umaxmin(1),MAX_DATA)
    call MPI_CollectResiduals(umin,ins_umaxmin(2),MIN_DATA)

    call MPI_CollectResiduals(vmax,ins_vmaxmin(1),MAX_DATA)
    call MPI_CollectResiduals(vmin,ins_vmaxmin(2),MIN_DATA)

    ! Residuals

    do blk=1,blockCount
    do i=1,Nyb+2
          ins_u_res = ins_u_res + sum((facexData(:,i,blk,VELC_VAR)-facexData(:,i,blk,UOLD_VAR))**2)
          ins_v_res = ins_v_res + sum((faceyData(:,i,blk,VELC_VAR)-faceyData(:,i,blk,UOLD_VAR))**2)
    enddo
    enddo

    call MPI_CollectResiduals(ins_u_res,u_res1,SUM_DATA)
    ins_u_res = sqrt(u_res1/((nblockx*nblocky)*(Nxb+2)*(Nyb+2)))

    call MPI_CollectResiduals(ins_v_res,v_res1,SUM_DATA)
    ins_v_res = sqrt(v_res1/((nblockx*nblocky)*(Nxb+2)*(Nyb+2)))

    !do blk=1,blockCount
    !call ins_vorticity(tstep,solnData(:,:,blk,VORT_VAR),solnData(:,:,blk,VORO_VAR),&
    !                   facexData(:,:,blk,VELC_VAR),faceyData(:,:,blk,VELC_VAR),&
    !                   solnData(:,:,blk,DFUN_VAR))
    !end do

#ifdef MPI_DIST
    !call MPI_applyBC(VORT_VAR,CENTER)
#endif

#ifdef MPI_SHRD
    !call MPI_applyBC_shared(VORT_VAR,CENTER)
#endif

    !call MPI_physicalBC_vort(solnData(:,:,:,VORT_VAR))

    !do blk=1,blockCount
    !do i=1,Nyb+2
    !      ins_w_res = ins_w_res + sum((solnData(:,i,blk,VORT_VAR)-solnData(:,i,blk,VORO_VAR))**2)
    !enddo
    !enddo

    !call MPI_CollectResiduals(ins_w_res,w_res1,SUM_DATA)
    !ins_w_res = sqrt(w_res1/((nblockx*nblocky)*(Nxb+2)*(Nyb+2)))

    nullify(solnData,facexData,faceyData)

end subroutine
