subroutine ins_vorticity(tstep)

#include "Solver.h"

    use Grid_data
    use physicaldata
    use Driver_data
    use MPI_data
    use IncompNS_data
    use MPI_interface, ONLY: MPI_applyBC, MPI_CollectResiduals,MPI_physicalBC_vort

    implicit none
    integer, intent(in) :: tstep
    real, pointer, dimension(:,:) :: w,u,v
    real :: u_conv,v_conv,u_plus,u_mins,v_plus,v_mins
    real :: wx_plus,wx_mins,wy_plus,wy_mins
    real, allocatable,dimension(:,:) :: w_old
    integer :: i,j
    real :: w_res1

    allocate(w_old(Nxb+2,Nyb+2))

    w => ph_center(OMGA_VAR,:,:)
    u => ph_facex(VELC_VAR,:,:)
    v => ph_facey(VELC_VAR,:,:)

    w_old = w

    !$OMP PARALLEL DEFAULT(NONE) PRIVATE(i,j,u_conv,v_conv,u_plus,u_mins,&
    !$OMP v_plus,v_mins,wx_plus,wx_mins,wy_plus,wy_mins) NUM_THREADS(NTHREADS) &
    !$OMP SHARED(w,w_old,dr_dt,gr_dy,gr_dx,ins_inRe,u,v,dr_tile)

    !$OMP DO COLLAPSE(2) SCHEDULE(STATIC)

    do j=2,Nyb+1
      do i=2,Nxb+1

     u_conv = (u(i,j)+u(i-1,j))/2.
     v_conv = (v(i,j)+v(i,j-1))/2.

     u_plus = max(u_conv, 0.)
     u_mins = min(u_conv, 0.)

     v_plus = max(v_conv, 0.)
     v_mins = min(v_conv, 0.)

     wx_plus = (w_old(i+1,j)-w_old(i,j))/gr_dx
     wx_mins = (w_old(i,j)-w_old(i-1,j))/gr_dx

     wy_plus = (w_old(i,j+1)-w_old(i,j))/gr_dy
     wy_mins = (w_old(i,j)-w_old(i,j-1))/gr_dy

     w(i,j) = w_old(i,j)+((dr_dt*ins_inRe)/(gr_dx*gr_dx))*(w_old(i+1,j)+w_old(i-1,j)-2*w_old(i,j))&
                        +((dr_dt*ins_inRe)/(gr_dy*gr_dy))*(w_old(i,j+1)+w_old(i,j-1)-2*w_old(i,j))&
                        -((dr_dt))*(u_plus*wx_mins + u_mins*wx_plus)&
                        -((dr_dt))*(v_plus*wy_mins + v_mins*wy_plus)

       end do
     end do

     !$OMP END DO
     !$OMP END PARALLEL

     call MPI_applyBC(w)
     call MPI_physicalBC_vort(w)
   
     do i=1,Nyb+2
          ins_w_res = ins_w_res + sum((w(:,i)-w_old(:,i))**2)
     enddo

     call MPI_CollectResiduals(ins_w_res,w_res1,1)
     ins_w_res = sqrt(w_res1/((nblockx*nblocky)*(Nxb+2)*(Nyb+2))) 

     nullify(u,v,w)
     deallocate(w_old)

end subroutine ins_vorticity
