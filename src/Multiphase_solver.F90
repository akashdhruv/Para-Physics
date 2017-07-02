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
    integer :: j,i
    logical,intent(in) :: jump_flag

    real, pointer, dimension(:,:,:,:) :: solnData,facexData,faceyData

    solnData  => localCENTER
    facexData => localFACEX
    faceyData => localFACEY

    solnX = 0.0

if (jump_flag .eqv. .FALSE.) then

    call mph_FillVars_ibm(solnData(:,:,blockCount,DFUN_VAR),solnData(:,:,blockCount,PFUN_VAR),&
                           solnData(:,:,blockCount,THCO_VAR),solnData(:,:,blockCount,CPRS_VAR),&
                           solnData(:,:,blockCount,VISC_VAR),&
                           facexData(:,:,blockCount,RH2F_VAR),faceyData(:,:,blockCount,RH2F_VAR),&
                           facexData(:,:,blockCount,AL2F_VAR),faceyData(:,:,blockCount,AL2F_VAR),&
                           solnData(:,:,blockCount,TEMP_VAR),solnData(:,:,blockCount,TOLD_VAR),&
                           mph_beta)

#ifdef MPI_DIST
    call MPI_applyBC(solnData(:,:,blockCount,VISC_VAR))
    call MPI_applyBC(facexData(:,:,blockCount,RH2F_VAR))
    call MPI_applyBC(faceyData(:,:,blockCount,RH2F_VAR))
    call MPI_applyBC(facexData(:,:,blockCount,AL2F_VAR))
    call MPI_applyBC(faceyData(:,:,blockCount,AL2F_VAR))    
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
    call MPI_applyBC_RMA(solnData(:,:,blockCount,VISC_VAR))
    call MPI_applyBC_RMA(facexData(:,:,blockCount,RH2F_VAR))
    call MPI_applyBC_RMA(faceyData(:,:,blockCount,RH2F_VAR))
    call MPI_applyBC_RMA(facexData(:,:,blockCount,AL2F_VAR))
    call MPI_applyBC_RMA(faceyData(:,:,blockCount,AL2F_VAR))    
#endif

    call MPI_physicalBC_dfun(solnData(:,:,blockCount,VISC_VAR))
    call MPI_physicalBC_dfun(facexData(:,:,blockCount,RH2F_VAR))
    call MPI_physicalBC_dfun(faceyData(:,:,blockCount,RH2F_VAR))
    call MPI_physicalBC_dfun(facexData(:,:,blockCount,AL2F_VAR))
    call MPI_physicalBC_dfun(faceyData(:,:,blockCount,AL2F_VAR))


    if(tstep > 0) then 

        call mph_getInterfaceVelocity(facexData(:,:,blockCount,VELC_VAR),faceyData(:,:,blockCount,VELC_VAR),&
                                      facexData(:,:,blockCount,VELI_VAR),faceyData(:,:,blockCount,VELI_VAR),&
                                      solnData(:,:,blockCount,SMRH_VAR),solnData(:,:,blockCount,MDOT_VAR),&
                                      solnData(:,:,blockCount,NRMX_VAR),solnData(:,:,blockCount,NRMY_VAR))

    !   call mph_advect 
    !   call mph_redistance

#ifdef MPI_DIST
      call MPI_applyBC(facexData(:,:,blockCount,VELI_VAR))
      call MPI_applyBC(faceyData(:,:,blockCount,VELI_VAR))
#endif

#ifdef MPI_SHRD
      call MPI_BARRIER(shared_comm,ierr)
      call MPI_applyBC_shared(VELI_VAR,FACEX)
      call MPI_applyBC_shared(VELI_VAR,FACEY)
#endif

#ifdef MPI_RMA
      call MPI_applyBC_RMA(facexData(:,:,blockCount,VELI_VAR))
      call MPI_applyBC_RMA(faceyData(:,:,blockCount,VELI_VAR))
#endif

      call MPI_physicalBC_vel(facexData(:,:,blockCount,VELI_VAR),faceyData(:,:,blockCount,VELI_VAR))


    end if

else if (jump_flag .eqv. .TRUE.) then
                      
#ifdef MPH_DEBUG

#else
    call mph_PressureJumps(solnData(:,:,blockCount,DFUN_VAR),solnData(:,:,blockCount,PFUN_VAR),&
                           solnData(:,:,blockCount,CURV_VAR),&
                           facexData(:,:,blockCount,RH1F_VAR),faceyData(:,:,blockCount,RH1F_VAR),&
                           facexData(:,:,blockCount,RH2F_VAR),faceyData(:,:,blockCount,RH2F_VAR),&
                           solnData(:,:,blockCount,SIGP_VAR),&
                           facexData(:,:,blockCount,SIGM_VAR),faceyData(:,:,blockCount,SIGM_VAR),&
                           solnData(:,:,blockCount,MDOT_VAR))

#ifdef MPI_DIST
    call MPI_applyBC(facexData(:,:,blockCount,RH1F_VAR))
    call MPI_applyBC(facexData(:,:,blockCount,RH2F_VAR))
    call MPI_applyBC(faceyData(:,:,blockCount,RH1F_VAR))
    call MPI_applyBC(faceyData(:,:,blockCount,RH2F_VAR))
#endif

#ifdef MPI_SHRD
    call MPI_BARRIER(shared_comm,ierr)
    call MPI_applyBC_shared(RH1F_VAR,FACEX)
    call MPI_applyBC_shared(RH1F_VAR,FACEY)
    call MPI_applyBC_shared(RH2F_VAR,FACEX)
    call MPI_applyBC_shared(RH2F_VAR,FACEY)
#endif

#ifdef MPI_RMA
    call MPI_applyBC_RMA(facexData(:,:,blockCount,RH1F_VAR))
    call MPI_applyBC_RMA(facexData(:,:,blockCount,RH2F_VAR))
    call MPI_applyBC_RMA(faceyData(:,:,blockCount,RH1F_VAR))
    call MPI_applyBC_RMA(faceyData(:,:,blockCount,RH2F_VAR))
#endif

    call MPI_physicalBC_dfun(facexData(:,:,blockCount,RH1F_VAR))
    call MPI_physicalBC_dfun(facexData(:,:,blockCount,RH2F_VAR))
    call MPI_physicalBC_dfun(faceyData(:,:,blockCount,RH1F_VAR))
    call MPI_physicalBC_dfun(faceyData(:,:,blockCount,RH2F_VAR))

#endif

end if

   nullify(solnData,facexData,faceyData)

end subroutine Multiphase_solver

