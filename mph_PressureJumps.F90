subroutine mph_PressureJumps(s,pf,crv,rho1x,rho1y,rho2x,rho2y,w,sigx,sigy,mdot)

#include "Solver.h"

        use Multiphase_data, ONLY: mph_rho1,mph_rho2,mph_sten,mph_crmx,mph_crmn
        use Grid_data, ONLY: gr_dx,gr_dy

        implicit none

        real, dimension(:,:), intent(inout):: s,crv,rho1x,rho2x,rho1y, &
                                               rho2y,pf,w,sigx,sigy

        real, dimension(:,:), intent(in) :: mdot

        real :: th,aa,xijl,xijr, &
                a1,a2,cri,xid,xij,xidl,xidr,yid,yidr,yidl,yij,yijl,yijr,bb,mT

        integer :: i,j

        real, parameter :: eps = 1E-13

        integer :: iSmear

        iSmear  = 1

        mph_crmx = -1E10
        mph_crmn = 1E10
        sigx = 0.
        sigy = 0.
        w = 0.

        do j = 2,Nyb+1
           do i = 2,Nxb+1

              if(pf(i,j).eq.0..and.pf(i+1,j).eq.1.) then

                 !          = (+)            = (+)           = (-)
                 th = abs(s(i+1,j))/(abs(s(i+1,j))+abs(s(i,j)))

                 cri = crv(i+1,j)*(1.-th) + crv(i,j)*th

                 mT = mdot(i+1,j)*(1.-th) + mdot(i,j)*th

                 xijl = mph_sten*crv(i,j)
                 xijr = mph_sten*crv(i+1,j)
                 xidl = 0.
                 xidr = 0.

                 xij = xijl*th + xijr*(1.-th)
                 xid = xidl*th + xidr*(1.-th)

                 if (iSmear  .eq. 1) then

                 aa = th*(mph_rho1/mph_rho2) + (1.-th)*(mph_rho2/mph_rho2) 
                 rho1x(i+1,j) = rho1x(i+1,j)*(mph_rho1/mph_rho2)/aa
                 rho2x(i+1,j) = rho2x(i+1,j)*(mph_rho2/mph_rho2)/aa

                 bb = (1./aa) - 1.

                 w(i,j)   = w(i,j)   - xij/aa/gr_dx**2 - xid*th*(mph_rho1/mph_rho2)/aa/gr_dx      + (bb*mT*mT)/aa/gr_dx**2
                 w(i+1,j) = w(i+1,j) + xij/aa/gr_dx**2 - xid*(1.-th)*(mph_rho2/mph_rho2)/aa/gr_dx - (bb*mT*mT)/aa/gr_dx**2

                 sigx(i+1,j) = - xij/aa/gr_dx + (bb*mT*mT)/aa/gr_dx

                 else

                 end if

                 mph_crmx = max(abs(cri),mph_crmx)
                 mph_crmn = min(abs(cri),mph_crmn)

              end if

          end do
       end do

end subroutine mph_PressureJumps
