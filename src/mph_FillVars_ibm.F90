subroutine mph_FillVars_ibm(s,pf,thco,cprs,visc,rhox,rhoy,alpx,alpy,T,T_old,beta)

#include "Solver.h"

    implicit none
    real,intent(inout),dimension(Nxb+2,Nyb+2) :: s,pf,thco,cprs,visc,rhox,rhoy,alpx,alpy
    real,intent(in),dimension(Nxb+2,Nyb+2) :: T,T_old
    real,intent(in) :: beta

    real, parameter :: eps = 1E-10
    integer :: i,j

    do j=2,Nyb+1
     do i=2,Nxb+1

         visc(i,j) = visc(i,j)*(1+beta*(T(i,j)-T_old(i,j)))

         rhox(i,j) = rhox(i,j)*(1+beta*((T(i,j)+T(i+1,j))*0.5-&
                                         (T_old(i,j)+T_old(i+1,j))*0.5))

         rhoy(i,j) = rhoy(i,j)*(1+beta*((T(i,j)+T(i,j+1))*0.5-&
                                        (T_old(i,j)+T_old(i,j+1))*0.5))

         alpx(i,j) = alpx(i,j)*(1+beta*((T(i,j)+T(i+1,j))*0.5-&
                                         (T_old(i,j)+T_old(i+1,j))*0.5))

         alpy(i,j) = alpy(i,j)*(1+beta*((T(i,j)+T(i,j+1))*0.5-&
                                        (T_old(i,j)+T_old(i,j+1))*0.5))


     end do
    end do

end subroutine mph_FillVars_ibm
