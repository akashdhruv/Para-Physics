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

    use MPI_interface, only: MPI_applyBC_shared,MPI_applyBC,&
                             MPI_physicalBC_dfun,MPI_physicalBC_vel,&
                             MPI_applyBC_RMA

    use MPI_data, only: shared_comm,ierr,blockCount

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
    call mph_FillVars_ibm(solnData(:,:,blk,DFUN_VAR),solnData(:,:,blk,PFUN_VAR),&
                           solnData(:,:,blk,THCO_VAR),solnData(:,:,blk,CPRS_VAR),&
                           solnData(:,:,blk,VISC_VAR),&
                           facexData(:,:,blk,RH2F_VAR),faceyData(:,:,blk,RH2F_VAR),&
                           facexData(:,:,blk,AL2F_VAR),faceyData(:,:,blk,AL2F_VAR),&
                           solnData(:,:,blk,TEMP_VAR),solnData(:,:,blk,TOLD_VAR),&
                           mph_beta)
    end do

#ifdef MPI_DIST
    call MPI_applyBC(VISC_VAR,CENTER)
    call MPI_applyBC(RH2F_VAR,FACEX)
    call MPI_applyBC(RH2F_VAR,FACEY)
    call MPI_applyBC(AL2F_VAR,FACEX)
    call MPI_applyBC(AL2F_VAR,FACEY)    
#endif

#ifdef MPI_SHRD
    call MPI_BARRIER(shared_comm,ierr)
    call MPI_applyBC_shared(VISC_VAR,CENTER)
    call MPI_applyBC_shared(RH2F_VAR,FACEX)
    call MPI_applyBC_shared(RH2F_VAR,FACEY)
    call MPI_applyBC_shared(AL2F_VAR,FACEX)
    call MPI_applyBC_shared(AL2F_VAR,FACEY)
#endif

#ifdef MPI_RMA
    call MPI_applyBC_RMA(VISC_VAR,CENTER)
    call MPI_applyBC_RMA(RH2F_VAR,FACEX)
    call MPI_applyBC_RMA(RH2F_VAR,FACEY)
    call MPI_applyBC_RMA(AL2F_VAR,FACEX)
    call MPI_applyBC_RMA(AL2F_VAR,FACEY)   
#endif

    call MPI_physicalBC_dfun(solnData(:,:,:,VISC_VAR))
    call MPI_physicalBC_dfun(facexData(:,:,:,RH2F_VAR))
    call MPI_physicalBC_dfun(faceyData(:,:,:,RH2F_VAR))
    call MPI_physicalBC_dfun(facexData(:,:,:,AL2F_VAR))
    call MPI_physicalBC_dfun(faceyData(:,:,:,AL2F_VAR))

    if(tstep > 0) then 

        do blk=1,blockCount
        call mph_getInterfaceVelocity(facexData(:,:,blk,VELC_VAR),faceyData(:,:,blk,VELC_VAR),&
                                      facexData(:,:,blk,VELI_VAR),faceyData(:,:,blk,VELI_VAR),&
                                      solnData(:,:,blk,SMRH_VAR),solnData(:,:,blk,MDOT_VAR),&
                                      solnData(:,:,blk,NRMX_VAR),solnData(:,:,blk,NRMY_VAR))
        end do

    !   call mph_advect 
    !   call mph_redistance

#ifdef MPI_DIST
      call MPI_applyBC(VELI_VAR,FACEX)
      call MPI_applyBC(VELI_VAR,FACEY)
#endif

#ifdef MPI_SHRD
      call MPI_BARRIER(shared_comm,ierr)
      call MPI_applyBC_shared(VELI_VAR,FACEX)
      call MPI_applyBC_shared(VELI_VAR,FACEY)
#endif

#ifdef MPI_RMA
      call MPI_applyBC_RMA(VELI_VAR,FACEX)
      call MPI_applyBC_RMA(VELI_VAR,FACEY)
#endif

      call MPI_physicalBC_vel(facexData(:,:,:,VELI_VAR),faceyData(:,:,:,VELI_VAR))


    end if

else if (jump_flag .eqv. .TRUE.) then
                      
#ifdef MPH_DEBUG

#else
    do blk=1,blockCount
    call mph_PressureJumps(solnData(:,:,blk,DFUN_VAR),solnData(:,:,blk,PFUN_VAR),&
                           solnData(:,:,blk,CURV_VAR),&
                           facexData(:,:,blk,RH1F_VAR),faceyData(:,:,blk,RH1F_VAR),&
                           facexData(:,:,blk,RH2F_VAR),faceyData(:,:,blk,RH2F_VAR),&
                           solnData(:,:,blk,SIGP_VAR),&
                           facexData(:,:,blk,SIGM_VAR),faceyData(:,:,blk,SIGM_VAR),&
                           solnData(:,:,blk,MDOT_VAR))
    end do

#ifdef MPI_DIST
    call MPI_applyBC(RH1F_VAR,FACEX)
    call MPI_applyBC(RH1F_VAR,FACEY)
    call MPI_applyBC(RH2F_VAR,FACEX)
    call MPI_applyBC(RH2F_VAR,FACEY)
#endif

#ifdef MPI_SHRD
    call MPI_BARRIER(shared_comm,ierr)
    call MPI_applyBC_shared(RH1F_VAR,FACEX)
    call MPI_applyBC_shared(RH1F_VAR,FACEY)
    call MPI_applyBC_shared(RH2F_VAR,FACEX)
    call MPI_applyBC_shared(RH2F_VAR,FACEY)
#endif

#ifdef MPI_RMA
    call MPI_applyBC_RMA(RH1F_VAR,FACEX)
    call MPI_applyBC_RMA(RH1F_VAR,FACEY)
    call MPI_applyBC_RMA(RH2F_VAR,FACEX)
    call MPI_applyBC_RMA(RH2F_VAR,FACEY)
#endif

    call MPI_physicalBC_dfun(facexData(:,:,:,RH1F_VAR))
    call MPI_physicalBC_dfun(facexData(:,:,:,RH2F_VAR))
    call MPI_physicalBC_dfun(faceyData(:,:,:,RH1F_VAR))
    call MPI_physicalBC_dfun(faceyData(:,:,:,RH2F_VAR))

#endif

end if

   nullify(solnData,facexData,faceyData)

end subroutine Multiphase_solver

