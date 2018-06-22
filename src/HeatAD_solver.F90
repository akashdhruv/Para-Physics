subroutine HeatAD_solver(tstep)

#include "Solver.h"

   use HeatAD_interface, only: heat_tempSolver,heat_tempSolver_ibm,heat_tempSolver_mph
   use physicaldata, only: localCENTER,localFACEX,localFACEY
   use MPI_data, only: blockCount,solver_comm,ierr
   use HeatAD_data, only: ht_T_res
   use MPI_interface, ONLY: MPI_applyBC, MPI_CollectResiduals, MPI_physicalBC_temp
 
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

   solnData(:,:,TOLD_VAR,:) = solnData(:,:,TEMP_VAR,:) 

   do blk=1,blockCount
   call heat_tempSolver_ibm(tstep,solnData(:,:,TEMP_VAR,blk),solnData(:,:,TOLD_VAR,blk),&
                         solnData(:,:,MDOT_VAR,blk),solnData(:,:,SMRH_VAR,blk),&
                         facexData(:,:,VELC_VAR,blk),faceyData(:,:,VELC_VAR,blk),&
                         facexData(:,:,AL1F_VAR,blk),faceyData(:,:,AL1F_VAR,blk),&
                         facexData(:,:,AL2F_VAR,blk),faceyData(:,:,AL2F_VAR,blk),&
                         solnData(:,:,DFUN_VAR,blk),solnData(:,:,PFUN_VAR,blk),&
                         solnData(:,:,THCO_VAR,blk),solnData(:,:,CPRS_VAR,blk))   
   end do

   !do blk=1,blockCount
   !call heat_tempSolver(tstep,solnData(:,:,TEMP_VAR,blk),solnData(:,:,TOLD_VAR,blk),&
   !                      solnData(:,:,MDOT_VAR,blk),solnData(:,:,SMRH_VAR,blk),&
   !                      facexData(:,:,VELC_VAR,blk),faceyData(:,:,VELC_VAR,blk),&
   !                      facexData(:,:,AL1F_VAR,blk),faceyData(:,:,AL1F_VAR,blk),&
   !                      facexData(:,:,AL2F_VAR,blk),faceyData(:,:,AL2F_VAR,blk),&
   !                      solnData(:,:,DFUN_VAR,blk),solnData(:,:,PFUN_VAR,blk),&
   !                      solnData(:,:,THCO_VAR,blk),solnData(:,:,CPRS_VAR,blk))
   !end do

   !do blk=1,blockCount
   !call heat_tempSolver_mph(tstep,solnData(:,:,TEMP_VAR,blk),solnData(:,:,TOLD_VAR,blk),&
   !                      solnData(:,:,MDOT_VAR,blk),solnData(:,:,SMRH_VAR,blk),&
   !                      facexData(:,:,VELC_VAR,blk),faceyData(:,:,VELC_VAR,blk),&
   !                      facexData(:,:,AL1F_VAR,blk),faceyData(:,:,AL1F_VAR,blk),&
   !                      facexData(:,:,AL2F_VAR,blk),faceyData(:,:,AL2F_VAR,blk),&
   !                      solnData(:,:,DFUN_VAR,blk),solnData(:,:,PFUN_VAR,blk),&
   !                      solnData(:,:,THCO_VAR,blk),solnData(:,:,CPRS_VAR,blk))
   !end do

   call MPI_BARRIER(solver_comm,ierr)
   call MPI_applyBC(TEMP_VAR,CENTER)
   call MPI_physicalBC_temp(solnData(:,:,TEMP_VAR,:),solnData(:,:,DFUN_VAR,:),solnData(:,:,PFUN_VAR,:))

   do blk =1,blockCount
     do j=1,Nyb+2
       ht_T_res = ht_T_res + sum((solnData(:,j,TEMP_VAR,blk)-solnData(:,j,TOLD_VAR,blk))**2)
     end do
   end do

   call MPI_CollectResiduals(ht_T_res,T_res1,SUM_DATA)

   ht_T_res = sqrt(T_res1/((Nxb+2)*(Nyb+2)*(nblockx*nblocky)))

   nullify(solnData,facexData,faceyData)

end subroutine
