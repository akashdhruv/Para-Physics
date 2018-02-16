subroutine Multiphase_solver(tstep,solnX,jump_flag)

#define MPH_DEBUG
#include "Solver.h"

    use Driver_data,  only: dr_dt
    use Multiphase_data, only: mph_thco1, mph_cp1, mph_thco2, mph_cp2,mph_beta
    use physicaldata, only: localCENTER,localFACEX,localFACEY
    use Grid_data
    use Multiphase_interface , only: mph_FillVars, mph_PressureJumps,&
                                     mph_getInterfaceVelocity,&
                                     mph_FillVars_ibm    

    use MPI_interface, only: MPI_applyBC,MPI_physicalBC_dfun,MPI_physicalBC_vel

    use MPI_data, only: solver_comm,ierr,blockCount
    use HeatAD_data, only: ht_St

    implicit none

    integer, intent(in) :: tstep
    
    real,intent(out) :: solnX
    real :: ycell
    integer :: j,i,blk
    logical,intent(in) :: jump_flag

    real, pointer, dimension(:,:,:,:) :: solnData,facexData,faceyData

    solnData  => localCENTER
    facexData => localFACEX
    faceyData => localFACEY

    solnX = 0.0

  if (jump_flag .eqv. .FALSE.) then

    do blk=1,blockCount
    call mph_FillVars_ibm(solnData(:,:,DFUN_VAR,blk),solnData(:,:,PFUN_VAR,blk),&
                          solnData(:,:,THCO_VAR,blk),solnData(:,:,CPRS_VAR,blk),&
                          solnData(:,:,VISC_VAR,blk),&
                          facexData(:,:,RH2F_VAR,blk),faceyData(:,:,RH2F_VAR,blk),&
                          facexData(:,:,AL2F_VAR,blk),faceyData(:,:,AL2F_VAR,blk),&
                          solnData(:,:,TEMP_VAR,blk),solnData(:,:,TOLD_VAR,blk),&
                          mph_beta,ht_St)
    end do

    call MPI_BARRIER(solver_comm,ierr)
    call MPI_applyBC(VISC_VAR,CENTER)
    call MPI_applyBC(RH2F_VAR,FACEX)
    call MPI_applyBC(RH2F_VAR,FACEY)
    call MPI_applyBC(AL2F_VAR,FACEX)
    call MPI_applyBC(AL2F_VAR,FACEY)    
    call MPI_physicalBC_dfun(solnData(:,:,VISC_VAR,:))
    call MPI_physicalBC_dfun(facexData(:,:,RH2F_VAR,:))
    call MPI_physicalBC_dfun(faceyData(:,:,RH2F_VAR,:))
    call MPI_physicalBC_dfun(facexData(:,:,AL2F_VAR,:))
    call MPI_physicalBC_dfun(faceyData(:,:,AL2F_VAR,:))

    if(tstep > 0) then 

        do blk=1,blockCount
        call mph_getInterfaceVelocity(facexData(:,:,VELC_VAR,blk),faceyData(:,:,VELC_VAR,blk),&
                                      facexData(:,:,VELI_VAR,blk),faceyData(:,:,VELI_VAR,blk),&
                                      solnData(:,:,SMRH_VAR,blk),solnData(:,:,MDOT_VAR,blk),&
                                      solnData(:,:,NRMX_VAR,blk),solnData(:,:,NRMY_VAR,blk))
        end do

    !   call mph_advect 
    !   call mph_redistance

      call MPI_BARRIER(solver_comm,ierr)
      call MPI_applyBC(VELI_VAR,FACEX)
      call MPI_applyBC(VELI_VAR,FACEY)
      call MPI_physicalBC_vel(facexData(:,:,VELI_VAR,:),faceyData(:,:,VELI_VAR,:),solnData(:,:,PFUN_VAR,:))

    end if

 else if (jump_flag .eqv. .TRUE.) then
                      
#ifdef MPH_DEBUG

#else
    do blk=1,blockCount
    call mph_PressureJumps(solnData(:,:,DFUN_VAR,blk),solnData(:,:,PFUN_VAR,blk),&
                           solnData(:,:,CURV_VAR,blk),&
                           facexData(:,:,RH1F_VAR,blk),faceyData(:,:,RH1F_VAR,blk),&
                           facexData(:,:,RH2F_VAR,blk),faceyData(:,:,RH2F_VAR,blk),&
                           solnData(:,:,SIGP_VAR,blk),&
                           facexData(:,:,SIGM_VAR,blk),faceyData(:,:,SIGM_VAR,blk),&
                           solnData(:,:,MDOT_VAR,blk))
    end do

    call MPI_BARRIER(solver_comm,ierr)
    call MPI_applyBC(RH1F_VAR,FACEX)
    call MPI_applyBC(RH1F_VAR,FACEY)
    call MPI_applyBC(RH2F_VAR,FACEX)
    call MPI_applyBC(RH2F_VAR,FACEY)
    call MPI_physicalBC_dfun(facexData(:,:,RH1F_VAR,:))
    call MPI_physicalBC_dfun(facexData(:,:,RH2F_VAR,:))
    call MPI_physicalBC_dfun(faceyData(:,:,RH1F_VAR,:))
    call MPI_physicalBC_dfun(faceyData(:,:,RH2F_VAR,:))

#endif

 end if

   nullify(solnData,facexData,faceyData)

end subroutine Multiphase_solver
