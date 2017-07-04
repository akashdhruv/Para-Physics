subroutine Poisson_solver_VC(rvar,ivar,ps_counter,cvar,dvar)

  !$ use omp_lib
  use Grid_data
  use MPI_data
  use MPI_interface, ONLY: MPI_applyBC, MPI_CollectResiduals, MPI_physicalBC_pres
  use Driver_data, ONLY: dr_tile
  use IncompNS_data, ONLY: ins_timePoisson, ins_p_res
  use physicaldata, only: localCENTER,localFACEX,localFACEY

#include "Solver.h"
                
  implicit none

  integer, intent(in) :: rvar,ivar,cvar
  integer, optional, intent(in) :: dvar
  integer, intent(out) :: ps_counter


  real, pointer, dimension(:,:,:) :: ps, ps_RHS
  real, dimension(Nxb+2,Nyb+2,blockCount) :: ps_old,ps_rx,ps_ry
      
  real :: ps_res1
  integer :: i,j,thread_id,ii,jj,blk
  double precision :: poisson_start, poisson_finish

  ps     => localCENTER(:,:,:,ivar)
  ps_RHS => localCENTER(:,:,:,rvar)

  ps_old = 0
  ps_counter = 0

  if(present(dvar)) then
     ps_rx = 1./(localFACEX(:,:,:,cvar)+localFACEX(:,:,:,dvar))
     ps_ry = 1./(localFACEY(:,:,:,cvar)+localFACEY(:,:,:,dvar))
  else
     ps_rx = 1./(localFACEX(:,:,:,cvar))
     ps_ry = 1./(localFACEY(:,:,:,cvar))

  end if

  poisson_start = MPI_Wtime()

  do while(ps_counter<MaxIt)

     ins_p_res = 0  
     ps_old = ps

#ifdef POISSON_SOLVER_GS
    do blk=1,blockCount
     do j=2,Nyb+1
        do i=2,Nxb+1

           ps(i,j,blk)=((ps_old(i,j+1,blk)/(ps_ry(i,j,blk)*gr_dy*gr_dy))+(ps(i,j-1,blk)/(ps_ry(i,j-1,blk)*gr_dy*gr_dy))&
                       +(ps_old(i+1,j,blk)/(ps_rx(i,j,blk)*gr_dx*gr_dx))+(ps(i-1,j,blk)/(ps_rx(i-1,j,blk)*gr_dx*gr_dx))&
                       +ps_RHS(i,j,blk))&
                       *(1/((1/(ps_rx(i,j,blk)*gr_dx*gr_dx))+(1/(ps_ry(i,j,blk)*gr_dy*gr_dy))+&
                            (1/(ps_rx(i-1,j,blk)*gr_dx*gr_dx))+(1/(ps_ry(i,j-1,blk)*gr_dy*gr_dy))))

        end do
     end do
    end do
#endif

#ifdef POISSON_SOLVER_GS_SKEW
    do blk=1,blockCount
     do j=4,Nxb+Nyb+2 
        do i=max(2,j-Nyb-1),min(Nxb+1,j-2)

           ps(i,j-i,blk)=((ps_old(i,j-i+1,blk)/(ps_ry(i,j-i,blk)*gr_dy*gr_dy))+(ps(i,j-i-1,blk)/(ps_ry(i,j-i-1,blk)*gr_dy*gr_dy))&
                         +(ps_old(i+1,j-i,blk)/(ps_rx(i,j-i,blk)*gr_dx*gr_dx))+(ps(i-1,j-i,blk)/(ps_rx(i-1,j-i,blk)*gr_dx*gr_dx))&
                          +ps_RHS(i,j-i,blk))&
                          *(1/((1/(ps_rx(i,j-i,blk)*gr_dx*gr_dx))+(1/(ps_ry(i,j-i,blk)*gr_dy*gr_dy))+&
                               (1/(ps_rx(i-1,j-i,blk)*gr_dx*gr_dx))+(1/(ps_ry(i,j-i-1,blk)*gr_dy*gr_dy))))

        end do
     end do
   end do
#endif

#ifdef POISSON_SOLVER_GS_SOR
    do blk=1,blockCount
     do j=2,Nyb+1
        do i=2,Nxb+1

           ps(i,j,blk)=((ps_old(i,j+1,blk)/(ps_ry(i,j,blk)*gr_dy*gr_dy))+(ps(i,j-1,blk)/(ps_ry(i,j-1,blk)*gr_dy*gr_dy))&
                   +(ps_old(i+1,j,blk)/(ps_rx(i,j,blk)*gr_dx*gr_dx))+(ps(i-1,j,blk)/(ps_rx(i-1,j,blk)*gr_dx*gr_dx))&
                    +ps_RHS(i,j,blk))&
                    *(1/((1/(ps_rx(i,j,blk)*gr_dx*gr_dx))+(1/(ps_ry(i,j,blk)*gr_dy*gr_dy))+&
                         (1/(ps_rx(i-1,j,blk)*gr_dx*gr_dx))+(1/(ps_ry(i,j-1,blk)*gr_dy*gr_dy))))*omega + (1-omega)*ps_old(i,j,blk)

        end do
     end do
    end do
#endif

     call MPI_applyBC(ivar,CENTER)

     if(ivar == PRES_VAR) call MPI_physicalBC_pres(ps)
 
     ps_counter = ps_counter + 1

     ins_p_res = ins_p_res + sum(sum((ps-ps_old)**2,1))

     call MPI_CollectResiduals(ins_p_res,ps_res1,SUM_DATA)

     ins_p_res = sqrt(ps_res1/((Nxb+2)*(Nyb+2)*(nblockx*nblocky)))

     if((ins_p_res < 0.000001 ) .and. (ins_p_res .ne. 0) ) exit

  end do

  poisson_finish = MPI_Wtime()

  ins_timePoisson = ins_timePoisson + (poisson_finish - poisson_start)

end subroutine Poisson_solver_VC
