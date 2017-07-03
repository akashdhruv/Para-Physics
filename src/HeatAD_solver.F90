subroutine HeatAD_solver(tstep)

#define MPH_DEBUG
#include "Solver.h"

   use HeatAD_interface, only: heat_tempSolver,heat_tempSolver_ibm,heat_tempSolver_mph
   use physicaldata, only: localCENTER,localFACEX,localFACEY
   use MPI_data, only: blockCount,shared_comm,ierr
   use HeatAD_data, only: ht_T_res
   use MPI_interface, ONLY: MPI_applyBC, MPI_CollectResiduals,MPI_physicalBC_temp,&
                             MPI_applyBC_shared,MPI_applyBC_RMA
 
   implicit none
   
   integer, intent(in) :: tstep

   real, pointer, dimension(:,:,:,:) :: solnData,facexData,faceyData
   integer :: blk,j
   real :: T_res1

   ht_T_res = 0.0
   T_res1 = 0.0

   solnData  => localCENTER
   facexData => localFACEX
   faceyData => localFACEY

   solnData(:,:,:,TOLD_VAR) = solnData(:,:,:,TEMP_VAR) 

#ifdef SINGLEPHASE
   do blk=1,blockCount
   call heat_tempSolver(tstep,solnData(:,:,blk,TEMP_VAR),solnData(:,:,blk,TOLD_VAR),&
                         solnData(:,:,blk,MDOT_VAR),solnData(:,:,blk,SMRH_VAR),&
                         facexData(:,:,blk,VELC_VAR),faceyData(:,:,blk,VELC_VAR),&
                         facexData(:,:,blk,AL1F_VAR),faceyData(:,:,blk,AL1F_VAR),&
                         facexData(:,:,blk,AL2F_VAR),faceyData(:,:,blk,AL2F_VAR),&
                         solnData(:,:,blk,DFUN_VAR),solnData(:,:,blk,PFUN_VAR),&
                         solnData(:,:,blk,THCO_VAR),solnData(:,:,blk,CPRS_VAR))
   end do
#endif

#ifdef MULTIPHASE
#ifdef IBM
   do blk=1,blockCount
   call heat_tempSolver_ibm(tstep,solnData(:,:,blk,TEMP_VAR),solnData(:,:,blk,TOLD_VAR),&
                         solnData(:,:,blk,MDOT_VAR),solnData(:,:,blk,SMRH_VAR),&
                         facexData(:,:,blk,VELC_VAR),faceyData(:,:,blk,VELC_VAR),&
                         facexData(:,:,blk,AL1F_VAR),faceyData(:,:,blk,AL1F_VAR),&
                         facexData(:,:,blk,AL2F_VAR),faceyData(:,:,blk,AL2F_VAR),&
                         solnData(:,:,blk,DFUN_VAR),solnData(:,:,blk,PFUN_VAR),&
                         solnData(:,:,blk,THCO_VAR),solnData(:,:,blk,CPRS_VAR))
   end do
#else

#ifdef MPH_DEBUG
   do blk=1,blockCount
   call heat_tempSolver(tstep,solnData(:,:,blk,TEMP_VAR),solnData(:,:,blk,TOLD_VAR),&
                         solnData(:,:,blk,MDOT_VAR),solnData(:,:,blk,SMRH_VAR),&
                         facexData(:,:,blk,VELC_VAR),faceyData(:,:,blk,VELC_VAR),&
                         facexData(:,:,blk,AL1F_VAR),faceyData(:,:,blk,AL1F_VAR),&
                         facexData(:,:,blk,AL2F_VAR),faceyData(:,:,blk,AL2F_VAR),&
                         solnData(:,:,blk,DFUN_VAR),solnData(:,:,blk,PFUN_VAR),&
                         solnData(:,:,blk,THCO_VAR),solnData(:,:,blk,CPRS_VAR))
   end do
#else
   do blk=1,blockCount
   call heat_tempSolver_mph(tstep,solnData(:,:,blk,TEMP_VAR),solnData(:,:,blk,TOLD_VAR),&
                         solnData(:,:,blk,MDOT_VAR),solnData(:,:,blk,SMRH_VAR),&
                         facexData(:,:,blk,VELC_VAR),faceyData(:,:,blk,VELC_VAR),&
                         facexData(:,:,blk,AL1F_VAR),faceyData(:,:,blk,AL1F_VAR),&
                         facexData(:,:,blk,AL2F_VAR),faceyData(:,:,blk,AL2F_VAR),&
                         solnData(:,:,blk,DFUN_VAR),solnData(:,:,blk,PFUN_VAR),&
                         solnData(:,:,blk,THCO_VAR),solnData(:,:,blk,CPRS_VAR))
   end do
#endif
#endif
#endif

#ifdef MPI_DIST
   call MPI_applyBC(TEMP_VAR,CENTER)
#endif

#ifdef MPI_SHRD
   call MPI_BARRIER(shared_comm,ierr)
   call MPI_applyBC_shared(TEMP_VAR,CENTER)
#endif

#ifdef MPI_RMA
   call MPI_applyBC_RMA(solnData(:,:,:,TEMP_VAR))
#endif

   call MPI_physicalBC_temp(solnData(:,:,:,TEMP_VAR))

   do blk =1,blockCount
     do j=1,Nyb+2
       ht_T_res = ht_T_res + sum((solnData(:,j,blk,TEMP_VAR)-solnData(:,j,blk,TOLD_VAR))**2)
     end do
   end do

   call MPI_CollectResiduals(ht_T_res,T_res1,SUM_DATA)

   ht_T_res = sqrt(T_res1/((Nxb+2)*(Nyb+2)*(nblockx*nblocky)))

   nullify(solnData,facexData,faceyData)

end subroutine
