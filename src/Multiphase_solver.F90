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

    call mph_FillVars_ibm(solnData(:,:,DFUN_VAR,blockCount),solnData(:,:,PFUN_VAR,blockCount),&
                           solnData(:,:,THCO_VAR,blockCount),solnData(:,:,CPRS_VAR,blockCount),&
                           solnData(:,:,VISC_VAR,blockCount),&
                           facexData(:,:,RH2F_VAR,blockCount),faceyData(:,:,RH2F_VAR,blockCount),&
                           facexData(:,:,AL2F_VAR,blockCount),faceyData(:,:,AL2F_VAR,blockCount),&
                           solnData(:,:,TEMP_VAR,blockCount),solnData(:,:,TOLD_VAR,blockCount),&
                           mph_beta)

#ifdef MPI_DIST
    call MPI_applyBC(solnData(:,:,VISC_VAR,blockCount))
    call MPI_applyBC(facexData(:,:,RH2F_VAR,blockCount))
    call MPI_applyBC(faceyData(:,:,RH2F_VAR,blockCount))
    call MPI_applyBC(facexData(:,:,AL2F_VAR,blockCount))
    call MPI_applyBC(faceyData(:,:,AL2F_VAR,blockCount))    
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
    call MPI_applyBC_RMA(solnData(:,:,VISC_VAR,blockCount))
    call MPI_applyBC_RMA(facexData(:,:,RH2F_VAR,blockCount))
    call MPI_applyBC_RMA(faceyData(:,:,RH2F_VAR,blockCount))
    call MPI_applyBC_RMA(facexData(:,:,AL2F_VAR,blockCount))
    call MPI_applyBC_RMA(faceyData(:,:,AL2F_VAR,blockCount))    
#endif

    call MPI_physicalBC_dfun(solnData(:,:,VISC_VAR,blockCount))
    call MPI_physicalBC_dfun(facexData(:,:,RH2F_VAR,blockCount))
    call MPI_physicalBC_dfun(faceyData(:,:,RH2F_VAR,blockCount))
    call MPI_physicalBC_dfun(facexData(:,:,AL2F_VAR,blockCount))
    call MPI_physicalBC_dfun(faceyData(:,:,AL2F_VAR,blockCount))


    if(tstep > 0) then 

        call mph_getInterfaceVelocity(facexData(:,:,VELC_VAR,blockCount),faceyData(:,:,VELC_VAR,blockCount),&
                                      facexData(:,:,VELI_VAR,blockCount),faceyData(:,:,VELI_VAR,blockCount),&
                                      solnData(:,:,SMRH_VAR,blockCount),solnData(:,:,MDOT_VAR,blockCount),&
                                      solnData(:,:,NRMX_VAR,blockCount),solnData(:,:,NRMY_VAR,blockCount))

    !   call mph_advect 
    !   call mph_redistance

#ifdef MPI_DIST
      call MPI_applyBC(facexData(:,:,VELI_VAR,blockCount))
      call MPI_applyBC(faceyData(:,:,VELI_VAR,blockCount))
#endif

#ifdef MPI_SHRD
      call MPI_BARRIER(shared_comm,ierr)
      call MPI_applyBC_shared(VELI_VAR,FACEX)
      call MPI_applyBC_shared(VELI_VAR,FACEY)
#endif

#ifdef MPI_RMA
      call MPI_applyBC_RMA(facexData(:,:,VELI_VAR,blockCount))
      call MPI_applyBC_RMA(faceyData(:,:,VELI_VAR,blockCount))
#endif

      call MPI_physicalBC_vel(facexData(:,:,VELI_VAR,blockCount),faceyData(:,:,VELI_VAR,blockCount))


    end if

else if (jump_flag .eqv. .TRUE.) then
                      
#ifdef MPH_DEBUG

#else
    call mph_PressureJumps(solnData(:,:,DFUN_VAR,blockCount),solnData(:,:,PFUN_VAR,blockCount),&
                           solnData(:,:,CURV_VAR,blockCount),&
                           facexData(:,:,RH1F_VAR,blockCount),faceyData(:,:,RH1F_VAR,blockCount),&
                           facexData(:,:,RH2F_VAR,blockCount),faceyData(:,:,RH2F_VAR,blockCount),&
                           solnData(:,:,SIGP_VAR,blockCount),&
                           facexData(:,:,SIGM_VAR,blockCount),faceyData(:,:,SIGM_VAR,blockCount),&
                           solnData(:,:,MDOT_VAR,blockCount))

#ifdef MPI_DIST
    call MPI_applyBC(facexData(:,:,RH1F_VAR,blockCount))
    call MPI_applyBC(facexData(:,:,RH2F_VAR,blockCount))
    call MPI_applyBC(faceyData(:,:,RH1F_VAR,blockCount))
    call MPI_applyBC(faceyData(:,:,RH2F_VAR,blockCount))
#endif

#ifdef MPI_SHRD
    call MPI_BARRIER(shared_comm,ierr)
    call MPI_applyBC_shared(RH1F_VAR,FACEX)
    call MPI_applyBC_shared(RH1F_VAR,FACEY)
    call MPI_applyBC_shared(RH2F_VAR,FACEX)
    call MPI_applyBC_shared(RH2F_VAR,FACEY)
#endif

#ifdef MPI_RMA
    call MPI_applyBC_RMA(facexData(:,:,RH1F_VAR,blockCount))
    call MPI_applyBC_RMA(facexData(:,:,RH2F_VAR,blockCount))
    call MPI_applyBC_RMA(faceyData(:,:,RH1F_VAR,blockCount))
    call MPI_applyBC_RMA(faceyData(:,:,RH2F_VAR,blockCount))
#endif

    call MPI_physicalBC_dfun(facexData(:,:,RH1F_VAR,blockCount))
    call MPI_physicalBC_dfun(facexData(:,:,RH2F_VAR,blockCount))
    call MPI_physicalBC_dfun(faceyData(:,:,RH1F_VAR,blockCount))
    call MPI_physicalBC_dfun(faceyData(:,:,RH2F_VAR,blockCount))

#endif

end if

   nullify(solnData,facexData,faceyData)

end subroutine Multiphase_solver

