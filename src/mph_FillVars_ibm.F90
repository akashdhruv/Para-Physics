subroutine mph_FillVars_ibm(s,pf,thco,cprs,visc,rhox,rhoy,alpx,alpy,T,T_old,beta)

#include "Solver.h"

    use MPI_interface, ONLY: MPI_applyBC,MPI_physicalBC_dfun, MPI_applyBC_shared
    use physicaldata, only: SHD_solnData, SHD_facexData, SHD_faceyData

    implicit none
    real,intent(inout),dimension(Nxb+2,Nyb+2) :: s,pf,thco,cprs,visc,rhox,rhoy,alpx,alpy
    real,intent(in),dimension(Nxb+2,Nyb+2) :: T,T_old
    real,intent(in) :: beta

    real, parameter :: eps = 1E-10
    integer :: i,j

    do j=2,Nyb+1
     do i=2,Nxb+1

         !visc(i,j) = visc(i,j)*(((T(i,j))/(T_old(i,j)))**0.7)

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

#ifdef MPI_DIST
    call MPI_applyBC(visc)
    call MPI_applyBC(rhox)
    call MPI_applyBC(rhoy)
    call MPI_applyBC(alpx)
    call MPI_applyBC(alpy)    
#endif

#ifdef MPI_SHRD
    call MPI_applyBC_shared(visc,SHD_solnData(VISC_VAR,:,:))
    call MPI_applyBC_shared(rhox,SHD_facexData(RH2F_VAR,:,:))
    call MPI_applyBC_shared(rhoy,SHD_faceyData(RH2F_VAR,:,:))
    call MPI_applyBC_shared(alpx,SHD_facexData(AL2F_VAR,:,:))
    call MPI_applyBC_shared(alpy,SHD_faceyData(AL2F_VAR,:,:))
#endif

    call MPI_physicalBC_dfun(visc)
    call MPI_physicalBC_dfun(rhox)
    call MPI_physicalBC_dfun(rhoy)
    call MPI_physicalBC_dfun(alpx)
    call MPI_physicalBC_dfun(alpy)


end subroutine mph_FillVars_ibm
