subroutine Poisson_solver(ps_RHS,ps,ps_res,ps_counter,ps_quant)

  !$ use omp_lib
  use Grid_data
  use MPI_data
  use MPI_interface, ONLY: MPI_applyBC, MPI_CollectResiduals, MPI_physicalBC_pres, MPI_applyBC_shared, MPI_applyBC_RMA
  use Driver_data, ONLY: dr_tile
  use IncompNS_data, ONLY: ins_timePoisson

#include "Solver.h"
                
  implicit none

  real, dimension(:,:), intent(in) :: ps_RHS

  real, dimension(:,:), intent(inout) :: ps

  integer, intent(in) :: ps_quant

  real, dimension(Nxb+2,Nyb+2) :: ps_old,ps_new

  real, intent(out) :: ps_res

  real, dimension(:,:), allocatable :: p_priv
        
  real :: ps_res1

  integer, intent(out) :: ps_counter
  integer :: i,j,thread_id,ii,jj

  double precision :: poisson_start, poisson_finish

  thread_id = 0
  ps_old = 0
  ps_counter = 0

  poisson_start = MPI_Wtime()

  !!DIR$ OFFLOAD BEGIN TARGET(mic) in(ps_old,gr_dy,gr_dx,ps_RHS,i,j,thread_id,ps_res1,ps_quant,dr_tile,ii,jj) inout(ps,ps_res,ps_counter)

  !$OMP PARALLEL PRIVATE(i,j,thread_id,ii,jj) DEFAULT(NONE) NUM_THREADS(NTHREADS) &
  !$OMP SHARED(ps_old,gr_dy,gr_dx,ps_RHS,ps,ps_res,ps_counter,ps_res1,ps_quant,dr_tile,shared_comm,ierr)

#if NTHREADS > 1
  thread_id = OMP_GET_THREAD_NUM()
#endif

  do while(ps_counter<MaxIt)

     if (thread_id == 0) then

     ps_res = 0  
     ps_old = ps

     end if
     
     !$OMP BARRIER

#ifdef POISSON_SOLVER_JACOBI
     ps(2:Nxb+1,2:Nyb+1)= ((ps_old(2:Nxb+1,3:Nyb+2)/(gr_dy*gr_dy))&
                          +(ps_old(2:Nxb+1,1:Nyb)/(gr_dy*gr_dy))&
                          +(ps_old(3:Nxb+2,2:Nyb+1)/(gr_dx*gr_dx))&
                          +(ps_old(1:Nxb,2:Nyb+1)/(gr_dx*gr_dx))&
                          +ps_RHS)&
                          *(1/((1/(gr_dx*gr_dx))&
                          +(1/(gr_dy*gr_dy))&
                          +(1/(gr_dx*gr_dx))&
                          +(1/(gr_dy*gr_dy))))
#endif

#ifdef POISSON_SOLVER_GS
     !$OMP DO COLLAPSE(2) SCHEDULE(STATIC) 

     !___Double Tile____
     !do jj=2,Nyb+1,dr_tile
     !do ii=2,Nxb+1,dr_tile
     !do j=jj,jj+dr_tile-1
     !   do i=ii,ii+dr_tile-1

     !___Single Tile____
     !do jj=2,Nyb+1,dr_tile
     ! do j=jj,jj+dr_tile-1
     !  do i=2,Nxb+1

     !____No Tile______
     do j=2,Nyb+1
        do i=2,Nxb+1

           ps(i,j)=((ps_old(i,j+1)/(gr_dy*gr_dy))+(ps(i,j-1)/(gr_dy*gr_dy))&
                  +(ps_old(i+1,j)/(gr_dx*gr_dx))+(ps(i-1,j)/(gr_dx*gr_dx))&
                  +ps_RHS(i-1,j-1))&
                  *(1/((1/(gr_dx*gr_dx))+(1/(gr_dy*gr_dy))+&
                   (1/(gr_dx*gr_dx))+(1/(gr_dy*gr_dy))))
     
        end do
     end do
     !___No Tile______



     !   end do
     !  end do
     !end do
     !___Single Tile___

     !   end do
     !  end do
     !end do
     !end do
     !__Double Tile___

     !$OMP END DO
#endif

#ifdef POISSON_SOLVER_GS_SKEW
     !$OMP DO SCHEDULE(STATIC)
     do j=4,Nxb+Nyb+2
        do i=max(2,j-Nyb-1),min(Nxb+1,j-2)
           ps(i,j-i)=((ps_old(i,j-i+1)/(gr_dy*gr_dy))+(ps(i,j-i-1)/(gr_dy*gr_dy))&
                     +(ps_old(i+1,j-i)/(gr_dx*gr_dx))+(ps(i-1,j-i)/(gr_dx*gr_dx))&
                     +ps_RHS(i-1,j-i-1))&
                     *(1/((1/(gr_dx*gr_dx))+(1/(gr_dy*gr_dy))+&
                          (1/(gr_dx*gr_dx))+(1/(gr_dy*gr_dy))))     
        end do
     end do
     !$OMP END DO
#endif

#ifdef POISSON_SOLVER_GS_SOR
     !$OMP DO COLLAPSE(2) SCHEDULE(STATIC)

     do j=2,Nyb+1
        do i=2,Nxb+1

           ps(i,j)=((ps_old(i,j+1)/(gr_dy*gr_dy))+(ps(i,j-1)/(gr_dy*gr_dy))&
                   +(ps_old(i+1,j)/(gr_dx*gr_dx))+(ps(i-1,j)/(gr_dx*gr_dx))&
                   +ps_RHS(i-1,j-1))&
                   *(1/((1/(gr_dx*gr_dx))+(1/(gr_dy*gr_dy))+&
                        (1/(gr_dx*gr_dx))+(1/(gr_dy*gr_dy))))*omega + (1-omega)*ps_old(i,j)
                  
        end do
     end do

     !$OMP END DO
#endif

     ! Pressure BC

     if (thread_id == 0) then

#ifdef MPI_DIST
     call MPI_applyBC(ps)
#endif

#ifdef MPI_SHRD
     call MPI_BARRIER(shared_comm,ierr)
     call MPI_applyBC_shared(ps_quant,CENTER)
#endif

#ifdef MPI_RMA
     call MPI_applyBC_RMA(ps)
#endif

     if(ps_quant == PRES_VAR) call MPI_physicalBC_pres(ps)

     ps_counter = ps_counter + 1

     ps_res = ps_res + sum(sum((ps-ps_old)**2,1))

     call MPI_CollectResiduals(ps_res,ps_res1,SUM_DATA)

     ps_res = sqrt(ps_res1/((Nxb+2)*(Nyb+2)*(nblockx*nblocky)))

     end if
   
     !$OMP BARRIER

     if( (ps_res < 0.000001 ) .and. (ps_res .ne. 0) ) exit

  end do

  !$OMP END PARALLEL

  !!DIR$ END OFFLOAD

  poisson_finish = MPI_Wtime()

  ins_timePoisson = ins_timePoisson + (poisson_finish - poisson_start) 

end subroutine Poisson_solver
