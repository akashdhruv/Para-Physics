subroutine ins_vorticity(tstep)

#include "Solver.h"

    use Grid_data
    use physicaldata
    use Driver_data
    use MPI_data
    use IncompNS_data
    use MPI_interface, ONLY: MPI_applyBC, MPI_CollectResiduals,MPI_physicalBC_vort,MPI_physicalBC_dfun

    implicit none
    integer, intent(in) :: tstep
    real, pointer, dimension(:,:) :: w,u,v,s
    real :: u_conv,v_conv,u_plus,u_mins,v_plus,v_mins
    real :: wx_plus,wx_mins,wy_plus,wy_mins
    real, allocatable,dimension(:,:) :: w_old
    integer :: i,j
    real :: w_res1
    real :: w_sat,th,tol,wij

    allocate(w_old(Nxb+2,Nyb+2))

    w => ph_center(VORT_VAR,:,:)
    u => ph_facex(VELC_VAR,:,:)
    v => ph_facey(VELC_VAR,:,:)
    s => ph_center(DFUN_VAR,:,:)

    w_old = w
    
    tol = 0.01

    !$OMP PARALLEL DEFAULT(NONE) PRIVATE(w_sat,th,wij,i,j,u_conv,v_conv,u_plus,u_mins,&
    !$OMP v_plus,v_mins,wx_plus,wx_mins,wy_plus,wy_mins) NUM_THREADS(NTHREADS) &
    !$OMP SHARED(s,tol,w,w_old,dr_dt,gr_dy,gr_dx,ins_inRe,u,v,dr_tile)

    w_sat = 1.0

    !$OMP DO COLLAPSE(2) SCHEDULE(STATIC)

    do j=2,Nyb+1
      do i=2,Nxb+1

     u_conv = (u(i,j)+u(i-1,j))/2.
     v_conv = (v(i,j)+v(i,j-1))/2.

     u_plus = max(u_conv, 0.)
     u_mins = min(u_conv, 0.)

     v_plus = max(v_conv, 0.)
     v_mins = min(v_conv, 0.)

     wx_plus = w_old(i+1,j)
     wx_mins = w_old(i-1,j)

     wy_plus = w_old(i,j+1)
     wy_mins = w_old(i,j-1)

     wij = w_old(i,j)


     ! Case 1 !
     if(s(i,j)*s(i+1,j).le.0.d0) then

       if (abs(s(i,j))/(abs(s(i,j))+abs(s(i+1,j))) .gt. tol) then

       th = abs(s(i,j))/(abs(s(i,j))+abs(s(i+1,j)))
       wx_plus = (w_sat-w_old(i,j))/th + wij

       else

       th = abs(s(i-1,j))/(abs(s(i-1,j))+abs(s(i+1,j)))
       wx_plus = (w_sat-w_old(i-1,j))/th + w_old(i-1,j)
       
       end if
     end if
     ! End of Case 1 !


     ! Case 2 !
     if(s(i,j)*s(i-1,j).le.0.d0) then

       if (abs(s(i,j))/(abs(s(i,j))+abs(s(i-1,j))) .gt. tol) then

       th = abs(s(i,j))/(abs(s(i,j))+abs(s(i-1,j)))
       wx_mins = (w_sat-w_old(i,j))/th + wij

       else

       th = abs(s(i+1,j))/(abs(s(i+1,j))+abs(s(i-1,j)))
       wx_mins = (w_sat-w_old(i+1,j))/th + w_old(i+1,j)
       
       end if
     end if
     ! End of Case 2 !


    ! Case 3 !
    if(s(i,j)*s(i,j+1).le.0.d0) then

      if (abs(s(i,j))/(abs(s(i,j))+abs(s(i,j+1))) .gt. tol) then

      th = abs(s(i,j))/(abs(s(i,j))+abs(s(i,j+1)))
      wy_plus = (w_sat-w_old(i,j))/th + wij

      else

      th = abs(s(i,j-1))/(abs(s(i,j-1))+abs(s(i,j+1)))
      wy_plus = (w_sat-w_old(i,j-1))/th + w_old(i,j-1)
      
      end if
    end if
    ! End of Case 3 !

    ! Case 4 !
    if(s(i,j)*s(i,j-1).le.0.d0) then

      if (abs(s(i,j))/(abs(s(i,j))+abs(s(i,j-1))) .gt. tol) then

      th = abs(s(i,j))/(abs(s(i,j))+abs(s(i,j-1)))
      wy_mins = (w_sat-w_old(i,j))/th + wij

      else

      th = abs(s(i,j+1))/(abs(s(i,j+1))+abs(s(i,j-1)))
      wy_mins = (w_sat-w_old(i,j+1))/th + w_old(i,j+1)

      end if
    end if
    ! End of Case 4 ! 

    w(i,j) = w_old(i,j)+((dr_dt*ins_inRe)/(gr_dx*gr_dx))*(wx_plus+wx_mins-2*wij)&
                       +((dr_dt*ins_inRe)/(gr_dy*gr_dy))*(wy_plus+wy_mins-2*wij)&
                       -((dr_dt))*(u_plus*(wij-wx_mins)/gr_dx +u_mins*(wx_plus-wij)/gr_dx)&
                       -((dr_dt))*(v_plus*(wij-wy_mins)/gr_dy +v_mins*(wy_plus-wij)/gr_dy)
 
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
