subroutine heat_tempSolver_ibm(tstep,T,T_old,mdot,smrh,u,v,a1x,a1y,a2x,a2y,s,pf,thco,cp)

#include "Solver.h"

   !$ use omp_lib
   use Grid_data
   use IncompNS_data
   use HeatAD_data
   use Driver_data
   use MPI_interface, only: MPI_applyBC, MPI_physicalBC_temp, MPI_CollectResiduals
   use Multiphase_data, only: mph_cp2,mph_thco2,mph_max_s,mph_min_s
   use IBM_data, only: ibm_cp1,ibm_thco1

   implicit none
      
   integer, intent(in) :: tstep
   real, intent(inout), dimension(:,:) :: T,T_old,mdot,smrh,u,v,a1x,a1y,s,pf,thco,cp,a2x,a2y
   integer :: i,j,ii,jj

   real :: u_plus, u_mins, v_plus, v_mins, u_conv, v_conv
   real :: Tx_plus, Tx_mins, Ty_plus, Ty_mins
   real :: Tij, Tipj, Timj, Tijp, Tijm
   real :: Txx, Tyy, th, dxp, dxm, dyp, dym
   real :: alphax_plus, alphax_mins, alphay_plus, alphay_mins, alpha_interface
   real :: E_source

  real :: tol

  real :: T_res1
  real :: Tsat

  ht_T_res = 0.0
  T_res1 = 0.0

  tol = 0.01
  Tsat = 1.0

  !$OMP PARALLEL DEFAULT(NONE) PRIVATE(i,j,u_conv,v_conv,u_plus,u_mins,&
  !$OMP v_plus,v_mins,Tx_plus,Tx_mins,Ty_plus,Ty_mins,ii,jj,th,Tipj,Timj,Txx,Tyy,&
  !$OMP Tij,Tijp,Tijm,alphax_plus,alphay_plus,alphax_mins,alphay_mins,alpha_interface)&
  !$OMP NUM_THREADS(NTHREADS) &
  !$OMP SHARED(T,Tsat,a1x,a1y,a2x,a2y,T_old,dr_dt,gr_dy,gr_dx,ht_Pr,ins_inRe,u,v,dr_tile,s,tol,thco,cp,ht_Nu,ibm_cp1,ibm_thco1)

  !$OMP DO COLLAPSE(2) SCHEDULE(STATIC)
  !do jj=2,Nyb+1,dr_tile
  !do ii=2,Nxb+1,dr_tile
  !do j=jj,jj+dr_tile-1
     !do i=ii,ii+dr_tile-1
  do j=2,Nyb+1
     do i=2,Nxb+1

     u_conv = (u(i,j)+u(i-1,j))/2.
     v_conv = (v(i,j)+v(i,j-1))/2.

     u_plus = max(u_conv, 0.)
     u_mins = min(u_conv, 0.)

     v_plus = max(v_conv, 0.)
     v_mins = min(v_conv, 0.)

     Tx_plus = T_old(i+1,j)
     Tx_mins = T_old(i-1,j)

     Ty_plus = T_old(i,j+1)
     Ty_mins = T_old(i,j-1)

     Tij = T_old(i,j)

     ! Case 1 !
     if(s(i,j)*s(i+1,j).le.0.d0) then

       if (abs(s(i,j))/(abs(s(i,j))+abs(s(i+1,j))) .gt. tol) then

        th = abs(s(i,j))/(abs(s(i,j))+abs(s(i+1,j)))
 
       else

        th = tol

       end if

       Tx_plus = (Tsat-Tij)/th + Tij

     end if
     ! End of Case 1 !

     ! Case 2 !
     if(s(i,j)*s(i-1,j).le.0.d0) then

       if (abs(s(i,j))/(abs(s(i,j))+abs(s(i-1,j))) .gt. tol) then

        th = abs(s(i,j))/(abs(s(i,j))+abs(s(i-1,j)))

       else

        th = tol

       end if

       Tx_mins = (Tsat-Tij)/th + Tij

     end if
     ! End of Case 2 !

    ! Case 3 !
    if(s(i,j)*s(i,j+1).le.0.d0) then

      if (abs(s(i,j))/(abs(s(i,j))+abs(s(i,j+1))) .gt. tol) then

       th = abs(s(i,j))/(abs(s(i,j))+abs(s(i,j+1)))

      else

       th = tol

      end if

      Ty_plus = (Tsat-Tij)/th + Tij

    end if
    ! End of Case 3 !

    ! Case 4 !
    if(s(i,j)*s(i,j-1).le.0.d0) then

      if (abs(s(i,j))/(abs(s(i,j))+abs(s(i,j-1))) .gt. tol) then

        th = abs(s(i,j))/(abs(s(i,j))+abs(s(i,j-1)))

      else

        th = tol

      end if

      Ty_mins = (Tsat-Tij)/th + Tij

    end if
    ! End of Case 4 !

    if(s(i,j) .ge. 0.0) then
    T(i,j) = T_old(i,j)+((dr_dt*ins_inRe*(ibm_thco1/ibm_cp1))/(ht_Pr*gr_dx*gr_dx))*(Tx_plus+Tx_mins-2*Tij)&
                       +((dr_dt*ins_inRe*(ibm_thco1/ibm_cp1))/(ht_Pr*gr_dy*gr_dy))*(Ty_plus+Ty_mins-2*Tij)&
                       -((dr_dt))*(u_plus*(Tij-Tx_mins)/gr_dx + u_mins*(Tx_plus-Tij)/gr_dx)&
                       -((dr_dt))*(v_plus*(Tij-Ty_mins)/gr_dy + v_mins*(Ty_plus-Tij)/gr_dy)

    else

    T(i,j) = T_old(i,j)+dr_dt*((a1x(i,j)+a2x(i,j))*(ins_inRe/ht_Pr)*(Tx_plus-Tij)/gr_dx - &
                               (a1x(i-1,j)+a2x(i-1,j))*(ins_inRe/ht_Pr)*(Tij-Tx_mins)/gr_dx)/gr_dx&
                       +dr_dt*((a1y(i,j)+a2y(i,j))*(ins_inRe/ht_Pr)*(Ty_plus-Tij)/gr_dy - &
                               (a1y(i,j-1)+a2y(i,j-1))*(ins_inRe/ht_Pr)*(Tij-Ty_mins)/gr_dy)/gr_dy&
                       -((dr_dt))*(u_plus*(Tij-Tx_mins)/gr_dx + u_mins*(Tx_plus-Tij)/gr_dx)&
                       -((dr_dt))*(v_plus*(Tij-Ty_mins)/gr_dy + v_mins*(Ty_plus-Tij)/gr_dy)


    end if

    end do
  end do
  !end do
  !end do

  !$OMP END DO
  !$OMP END PARALLEL

  call MPI_applyBC(T)
  call MPI_physicalBC_temp(T)

  do i=1,Nxb+2
          ht_T_res = ht_T_res + sum((T(i,:)-T_old(i,:))**2)
  enddo

  call MPI_CollectResiduals(ht_T_res,T_res1,SUM_DATA)

  ht_T_res = sqrt(T_res1/((Nxb+2)*(Nyb+2)*(nblockx*nblocky)))

end subroutine heat_tempSolver_ibm