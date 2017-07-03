subroutine heat_tempSolver_mph(tstep,T,T_old,mdot,smrh,u,v,a1x,a1y,a2x,a2y,s,pf,thco,cp)

#define MPH_DEBUG
#include "Solver.h"

  !$ use omp_lib
  use Grid_data
  use IncompNS_data
  use HeatAD_data
  use Driver_data
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

  real :: Tsat

  tol  = 0.01
  Tsat = 0.0

  !$OMP PARALLEL DEFAULT(NONE) PRIVATE(i,j,u_conv,v_conv,u_plus,u_mins,&
  !$OMP v_plus,v_mins,Tx_plus,Tx_mins,Ty_plus,Ty_mins,ii,jj,th,Tipj,Timj,Txx,Tyy,&
  !$OMP Tij,Tijp,Tijm,alphax_plus,alphay_plus,alphax_mins,alphay_mins,alpha_interface) &
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

#ifdef MPH_DEBUG

#else
     ! Case 1 !
     if(s(i,j)*s(i+1,j).le.0.d0) then

       th = abs(s(i,j))/(abs(s(i,j))+abs(s(i+1,j)))

       if (th .lt. tol) th = tol

       Tx_plus = (Tsat-Tij)/th + Tij

     end if
     ! End of Case 1 !


     ! Case 2 !
     if(s(i,j)*s(i-1,j).le.0.d0) then

       th = abs(s(i,j))/(abs(s(i,j))+abs(s(i-1,j)))

       if (th .lt. tol) th = tol

       Tx_mins = (Tsat-Tij)/th + Tij

     end if
     ! End of Case 2 !


    ! Case 3 !
    if(s(i,j)*s(i,j+1).le.0.d0) then

      th = abs(s(i,j))/(abs(s(i,j))+abs(s(i,j+1)))

      if (th .lt. tol) th = tol

      Ty_plus = (Tsat-Tij)/th + Tij

    end if
    ! End of Case 3 !

    ! Case 4 !
    if(s(i,j)*s(i,j-1).le.0.d0) then

      th = abs(s(i,j))/(abs(s(i,j))+abs(s(i,j-1)))

      if (th .lt. tol) th = tol

      Ty_mins = (Tsat-T_old(i,j))/th + Tij

    end if
    ! End of Case 4 !
#endif

    alphax_plus = (thco(i,j)/cp(i,j))*(ins_inRe/ht_Pr)
    alphax_mins = (thco(i,j)/cp(i,j))*(ins_inRe/ht_Pr)
    alphay_plus = (thco(i,j)/cp(i,j))*(ins_inRe/ht_Pr)
    alphay_mins = (thco(i,j)/cp(i,j))*(ins_inRe/ht_Pr)

    Txx = (alphax_plus*(Tx_plus-Tij)/gr_dx - alphax_mins*(Tij-Tx_mins)/gr_dx)/gr_dx
    Tyy = (alphay_plus*(Ty_plus-Tij)/gr_dy - alphay_mins*(Tij-Ty_mins)/gr_dy)/gr_dy

    T(i,j) = T_old(i,j) + dr_dt*((-(u_plus*(Tij-Tx_mins)/gr_dx+u_mins*(Tx_plus-Tij)/gr_dx)&
                                    -(v_plus*(Tij-Ty_mins)/gr_dy+v_mins*(Ty_plus-Tij)/gr_dy))&
                                    +(Txx +Tyy))

    end do
  end do
  !end do
  !end do

  !$OMP END DO
  !$OMP END PARALLEL

end subroutine heat_tempSolver_mph
