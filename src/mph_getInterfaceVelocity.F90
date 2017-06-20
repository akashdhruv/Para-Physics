subroutine mph_getInterfaceVelocity(u,v,u_int,v_int,smrh,mdot,nrmx,nrmy)


#include "Solver.h"

       use MPI_interface, ONLY: MPI_applyBC,MPI_physicalBC_vel

       implicit none
      
       real, dimension(:,:), intent(in) :: u, v, smrh, mdot, nrmx, nrmy
       real, dimension(:,:), intent(out) :: u_int,v_int


       integer :: i,j
       real :: rhox,rhoy,mdotx,mdoty,normx,normy
  
        !++++++++++  U-COMPONENT  ++++++++++

        do j = 2,Nyb+1
        do i = 2,Nxb+1

             rhox  = (smrh(i-1,j)+smrh(i,j))/2.d0
             mdotx = (mdot(i-1,j)+mdot(i,j))/2.d0
             normx = (nrmy(i-1,j)+nrmx(i,j))/2.d0

             u_int(i,j)=u(i,j) + mdotx*normx*rhox

        end do
        end do

      !++++++++++  V-COMPONENT  ++++++++++

        do j = 2,Nyb+1
        do i = 2,Nxb+1

             rhoy  = (smrh(i,j-1)+smrh(i,j))/2.d0
             mdoty = (mdot(i,j-1)+mdot(i,j))/2.d0
             normy = (nrmx(i,j-1)+nrmy(i,j))/2.d0

             v_int(i,j)=v(i,j) + mdoty*normy*rhoy

        enddo
        enddo

      call MPI_applyBC(u_int)
      call MPI_applyBC(v_int)

      call MPI_physicalBC_vel(u_int,v_int)

end subroutine
