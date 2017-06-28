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

    use MPI_data, only: shared_comm,ierr

    implicit none

    integer, intent(in) :: tstep
    
    real,intent(out) :: solnX
    real :: ycell
    integer :: j,i
    logical,intent(in) :: jump_flag

    real, pointer, dimension(:,:,:) :: solnData,facexData,faceyData

    solnData  => localCENTER
    facexData => localFACEX
    faceyData => localFACEY

    solnX = 0.0

if (jump_flag .eqv. .FALSE.) then

    call mph_FillVars_ibm(solnData(DFUN_VAR,:,:),solnData(PFUN_VAR,:,:),&
                           solnData(THCO_VAR,:,:),solnData(CPRS_VAR,:,:),&
                           solnData(VISC_VAR,:,:),&
                           facexData(RH2F_VAR,:,:),faceyData(RH2F_VAR,:,:),&
                           facexData(AL2F_VAR,:,:),faceyData(AL2F_VAR,:,:),&
                           solnData(TEMP_VAR,:,:),solnData(TOLD_VAR,:,:),&
                           mph_beta)

#ifdef MPI_DIST
    call MPI_applyBC(solnData(VISC_VAR,:,:))
    call MPI_applyBC(facexData(RH2F_VAR,:,:))
    call MPI_applyBC(faceyData(RH2F_VAR,:,:))
    call MPI_applyBC(facexData(AL2F_VAR,:,:))
    call MPI_applyBC(faceyData(AL2F_VAR,:,:))    
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
    call MPI_applyBC_RMA(solnData(VISC_VAR,:,:))
    call MPI_applyBC_RMA(facexData(RH2F_VAR,:,:))
    call MPI_applyBC_RMA(faceyData(RH2F_VAR,:,:))
    call MPI_applyBC_RMA(facexData(AL2F_VAR,:,:))
    call MPI_applyBC_RMA(faceyData(AL2F_VAR,:,:))    
#endif

    call MPI_physicalBC_dfun(solnData(VISC_VAR,:,:))
    call MPI_physicalBC_dfun(facexData(RH2F_VAR,:,:))
    call MPI_physicalBC_dfun(faceyData(RH2F_VAR,:,:))
    call MPI_physicalBC_dfun(facexData(AL2F_VAR,:,:))
    call MPI_physicalBC_dfun(faceyData(AL2F_VAR,:,:))


    if(tstep > 0) then 

        call mph_getInterfaceVelocity(facexData(VELC_VAR,:,:),faceyData(VELC_VAR,:,:),&
                                      facexData(VELI_VAR,:,:),faceyData(VELI_VAR,:,:),&
                                      solnData(SMRH_VAR,:,:),solnData(MDOT_VAR,:,:),&
                                      solnData(NRMX_VAR,:,:),solnData(NRMY_VAR,:,:))

    !   call mph_advect 
    !   call mph_redistance

#ifdef MPI_DIST
      call MPI_applyBC(facexData(VELI_VAR,:,:))
      call MPI_applyBC(faceyData(VELI_VAR,:,:))
#endif

#ifdef MPI_SHRD
      call MPI_BARRIER(shared_comm,ierr)
      call MPI_applyBC_shared(VELI_VAR,FACEX)
      call MPI_applyBC_shared(VELI_VAR,FACEY)
#endif

#ifdef MPI_RMA
      call MPI_applyBC_RMA(facexData(VELI_VAR,:,:))
      call MPI_applyBC_RMA(faceyData(VELI_VAR,:,:))
#endif

      call MPI_physicalBC_vel(facexData(VELI_VAR,:,:),faceyData(VELI_VAR,:,:))


    end if

else if (jump_flag .eqv. .TRUE.) then
                      
#ifdef MPH_DEBUG

#else
    call mph_PressureJumps(solnData(DFUN_VAR,:,:),solnData(PFUN_VAR,:,:),&
                           solnData(CURV_VAR,:,:),&
                           facexData(RH1F_VAR,:,:),faceyData(RH1F_VAR,:,:),&
                           facexData(RH2F_VAR,:,:),faceyData(RH2F_VAR,:,:),&
                           solnData(SIGP_VAR,:,:),&
                           facexData(SIGM_VAR,:,:),faceyData(SIGM_VAR,:,:),&
                           solnData(MDOT_VAR,:,:))

#ifdef MPI_DIST
    call MPI_applyBC(facexData(RH1F_VAR,:,:))
    call MPI_applyBC(facexData(RH2F_VAR,:,:))
    call MPI_applyBC(faceyData(RH1F_VAR,:,:))
    call MPI_applyBC(faceyData(RH2F_VAR,:,:))
#endif

#ifdef MPI_SHRD
    call MPI_BARRIER(shared_comm,ierr)
    call MPI_applyBC_shared(RH1F_VAR,FACEX)
    call MPI_applyBC_shared(RH1F_VAR,FACEY)
    call MPI_applyBC_shared(RH2F_VAR,FACEX)
    call MPI_applyBC_shared(RH2F_VAR,FACEY)
#endif

#ifdef MPI_RMA
    call MPI_applyBC_RMA(facexData(RH1F_VAR,:,:))
    call MPI_applyBC_RMA(facexData(RH2F_VAR,:,:))
    call MPI_applyBC_RMA(faceyData(RH1F_VAR,:,:))
    call MPI_applyBC_RMA(faceyData(RH2F_VAR,:,:))
#endif

    call MPI_physicalBC_dfun(facexData(RH1F_VAR,:,:))
    call MPI_physicalBC_dfun(facexData(RH2F_VAR,:,:))
    call MPI_physicalBC_dfun(faceyData(RH1F_VAR,:,:))
    call MPI_physicalBC_dfun(faceyData(RH2F_VAR,:,:))

#endif

end if

   nullify(solnData,facexData,faceyData)

end subroutine Multiphase_solver

