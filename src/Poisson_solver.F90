subroutine Poisson_solver(rvar,ivar,ps_counter)

  use Grid_data
  use MPI_data
  use MPI_interface, ONLY: MPI_applyBC, MPI_CollectResiduals, MPI_physicalBC_pres, MPI_applyBC_shared, MPI_applyBC_RMA
  use Driver_data, ONLY: dr_tile
  use IncompNS_data, ONLY: ins_timePoisson, ins_p_res
  use physicaldata, ONLY: localCENTER

#include "Solver.h"
                
  implicit none
 
  integer, intent(in)  :: rvar,ivar
  integer, intent(out) :: ps_counter

  real, dimension(Nxb+2,Nyb+2,blockCount) :: ps_old
  real, pointer, dimension(:,:,:) :: ps,ps_RHS

  real :: ps_res1
  integer :: i,j,thread_id,ii,jj,blk
  double precision :: poisson_start, poisson_finish

  ps     => localCENTER(:,:,:,ivar)
  ps_RHS => localCENTER(:,:,:,rvar)

  ps_old = 0
  ps_counter = 0

  poisson_start = MPI_Wtime()

  do while(ps_counter<MaxIt)

     ins_p_res = 0  
     ps_old = ps

     do blk=1,blockCount

#ifdef POISSON_SOLVER_JACOBI
     ps(2:Nxb+1,2:Nyb+1,blk)= ((ps_old(2:Nxb+1,3:Nyb+2,blk)/(gr_dy*gr_dy))&
                              +(ps_old(2:Nxb+1,1:Nyb,blk)/(gr_dy*gr_dy))&
                              +(ps_old(3:Nxb+2,2:Nyb+1,blk)/(gr_dx*gr_dx))&
                              +(ps_old(1:Nxb,2:Nyb+1,blk)/(gr_dx*gr_dx))&
                              +ps_RHS(2:Nxb+1,2:Nyb+1,blk)&
                              *(1/((1/(gr_dx*gr_dx))&
                                  +(1/(gr_dy*gr_dy))&
                                  +(1/(gr_dx*gr_dx))&
                                  +(1/(gr_dy*gr_dy))))
#endif

#ifdef POISSON_SOLVER_GS
     do j=2,Nyb+1
        do i=2,Nxb+1

           ps(i,j,blk)=((ps_old(i,j+1,blk)/(gr_dy*gr_dy))+(ps(i,j-1,blk)/(gr_dy*gr_dy))&
                       +(ps_old(i+1,j,blk)/(gr_dx*gr_dx))+(ps(i-1,j,blk)/(gr_dx*gr_dx))&
                       +ps_RHS(i,j,blk))&
                       *(1/((1/(gr_dx*gr_dx))+(1/(gr_dy*gr_dy))+&
                            (1/(gr_dx*gr_dx))+(1/(gr_dy*gr_dy))))
     
        end do
     end do
#endif

#ifdef POISSON_SOLVER_GS_SKEW
     do j=4,Nxb+Nyb+2
        do i=max(2,j-Nyb-1),min(Nxb+1,j-2)
           ps(i,j-i,blk)=((ps_old(i,j-i+1,blk)/(gr_dy*gr_dy))+(ps(i,j-i-1,blk)/(gr_dy*gr_dy))&
                         +(ps_old(i+1,j-i,blk)/(gr_dx*gr_dx))+(ps(i-1,j-i,blk)/(gr_dx*gr_dx))&
                         +ps_RHS(i,j-i,blk))&
                         *(1/((1/(gr_dx*gr_dx))+(1/(gr_dy*gr_dy))+&
                              (1/(gr_dx*gr_dx))+(1/(gr_dy*gr_dy))))     
        end do
     end do
#endif

#ifdef POISSON_SOLVER_GS_SOR
     do j=2,Nyb+1
        do i=2,Nxb+1

           ps(i,j,blk)=((ps_old(i,j+1,blk)/(gr_dy*gr_dy))+(ps(i,j-1,blk)/(gr_dy*gr_dy))&
                       +(ps_old(i+1,j,blk)/(gr_dx*gr_dx))+(ps(i-1,j,blk)/(gr_dx*gr_dx))&
                       +ps_RHS(i,j,blk))&
                       *(1/((1/(gr_dx*gr_dx))+(1/(gr_dy*gr_dy))+&
                            (1/(gr_dx*gr_dx))+(1/(gr_dy*gr_dy))))*omega + (1-omega)*ps_old(i,j,blk)
                  
        end do
     end do
#endif

    end do

    ! Pressure BC

#ifdef MPI_DIST
    call MPI_applyBC(ivar,CENTER)
#endif

#ifdef MPI_SHRD
    call MPI_BARRIER(shared_comm,ierr)
    call MPI_applyBC_shared(ivar,CENTER)
#endif

#ifdef MPI_RMA
    call MPI_applyBC_RMA(ps)
#endif

    if(ivar == PRES_VAR) call MPI_physicalBC_pres(ps)

    ps_counter = ps_counter + 1

    do blk=1,blockCount
       ins_p_res = ins_p_res + sum(sum((ps(:,:,blk)-ps_old(:,:,blk))**2,1))
    end do

    call MPI_CollectResiduals(ins_p_res,ps_res1,SUM_DATA)

    ins_p_res = sqrt(ps_res1/((Nxb+2)*(Nyb+2)*(nblockx*nblocky)))

    if( (ins_p_res < 0.000001 ) .and. (ins_p_res .ne. 0) ) exit

  end do

  nullify(ps,ps_RHS)

  poisson_finish = MPI_Wtime()

  ins_timePoisson = ins_timePoisson + (poisson_finish - poisson_start) 

end subroutine Poisson_solver
