subroutine mph_PressureJumps(s,pf,crv,rho1x,rho1y,rho2x,rho2y,w,sigx,sigy,mdot)

#include "Solver.h"

        use Multiphase_data, ONLY: mph_rho1,mph_rho2,mph_sten,mph_crmx,mph_crmn
        use Grid_data, ONLY: gr_dx,gr_dy

        implicit none

        real, dimension(:,:), intent(inout):: crv,rho1x,rho2x,rho1y, &
                                               rho2y,pf,w,sigx,sigy

        real, dimension(:,:), intent(in) :: mdot,s

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

        do j = 1,Nyb+1
           do i = 1,Nxb+1

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

             if(pf(i,j).eq.1..and.pf(i+1,j).eq.0.) then

                 th = abs(s(i,j))/(abs(s(i,j))+abs(s(i+1,j)))

                 cri = crv(i,j)*(1.-th) + crv(i+1,j)*th
                 mT = mdot(i,j)*(1.-th) + mdot(i+1,j)*th

                 xijl = mph_sten*crv(i,j)
                 xijr = mph_sten*crv(i+1,j)
                 xidl = 0.
                 xidr = 0.

                 xij = xijl*(1.-th) + xijr*th
                 xid = xidl*(1.-th) + xidr*th

                 if (iSmear  .eq. 1) then

                 aa = th*(mph_rho1/mph_rho2) + (1.-th)*(mph_rho2/mph_rho2)
                 rho1x(i+1,j) = rho1x(i+1,j)*(mph_rho1/mph_rho2)/aa
                 rho2x(i+1,j) = rho2x(i+1,j)*(mph_rho2/mph_rho2)/aa

                 bb = (1./aa) - 1.

                 w(i,j)   = w(i,j)   + xij/aa/gr_dx**2 + xid*(1.-th)*(mph_rho2/mph_rho2)/aa/gr_dx - (bb*mT*mT)/aa/gr_dx**2
                 w(i+1,j) = w(i+1,j) - xij/aa/gr_dx**2 + xid*th*(mph_rho1/mph_rho2)/aa/gr_dx      + (bb*mT*mT)/aa/gr_dx**2

                 sigx(i+1,j) = xij/aa/gr_dx - (bb*mT*mT)/aa/gr_dx

                 else

                 end if

                 mph_crmx = max(abs(cri),mph_crmx)
                 mph_crmn = min(abs(cri),mph_crmn)

              end if

              if(pf(i,j).eq.0..and.pf(i,j+1).eq.1.) then

                 th = abs(s(i,j+1))/(abs(s(i,j+1))+abs(s(i,j)))

                 cri = crv(i,j+1)*(1.-th) + crv(i,j)*th
                 mT = mdot(i,j+1)*(1.-th) + mdot(i,j)*th

                 yijl = mph_sten*crv(i,j)
                 yijr = mph_sten*crv(i,j+1)
                 yidl = 0.
                 yidr = 0.

                 yij = yijl*th + yijr*(1.-th)
                 yid = yidl*th + yidr*(1.-th)

                 if (iSmear  .eq. 1) then

                 !- kpd - All Densities are relative to rho2...
                 aa = th*(mph_rho1/mph_rho2) + (1.-th)*(mph_rho2/mph_rho2)
                 rho1y(i,j+1) = rho1y(i,j+1)*(mph_rho1/mph_rho2)/aa
                 rho2y(i,j+1) = rho2y(i,j+1)*(mph_rho2/mph_rho2)/aa

                 bb = (1./aa) - 1.

                 w(i,j)   = w(i,j) - yij/aa/gr_dy**2 - yid*th*(mph_rho1/mph_rho2)/aa/gr_dy + (bb*mT*mT)/aa/gr_dy**2
                 w(i,j+1) = w(i,j+1)   + yij/aa/gr_dy**2 - yid*(1.-th)*(mph_rho2/mph_rho2)/aa/gr_dy- (bb*mT*mT)/aa/gr_dy**2

                 sigy(i,j+1) = - yij/aa/gr_dy + (bb*mT*mT)/aa/gr_dy

                 else

                 end if

                 mph_crmx = max(abs(cri),mph_crmx)
                 mph_crmn = min(abs(cri),mph_crmn)

              end if

              if(pf(i,j).eq.1..and.pf(i,j+1).eq.0.) then

                 th = abs(s(i,j))/(abs(s(i,j))+abs(s(i,j+1)))

                 cri = crv(i,j)*(1.-th) + crv(i,j+1)*th
                 mT = mdot(i,j)*(1.-th) + mdot(i,j+1)*th

                 yijl = mph_sten*crv(i,j)
                 yijr = mph_sten*crv(i,j+1)
                 yidl = 0.
                 yidr = 0.

                 yij = yijl*(1.-th) + yijr*th
                 yid = yidl*(1.-th) + yidr*th

                 if (iSmear  .eq. 1) then

                 aa = th*(mph_rho1/mph_rho2) + (1.-th)*(mph_rho2/mph_rho2)
                 rho1y(i,j+1) = rho1y(i,j+1)*(mph_rho1/mph_rho2)/aa
                 rho2y(i,j+1) = rho2y(i,j+1)*(mph_rho2/mph_rho2)/aa

                 bb = (1./aa) - 1.

                 w(i,j)   = w(i,j)   + yij/aa/gr_dy**2 + yid*(1.-th)*(mph_rho2/mph_rho2)/aa/gr_dy - (bb*mT*mT)/aa/gr_dy**2
                 w(i,j+1) = w(i,j+1) - yij/aa/gr_dy**2 + yid*th*(mph_rho1/mph_rho2)/aa/gr_dy      + (bb*mT*mT)/aa/gr_dy**2

                 sigy(i,j+1) = yij/aa/gr_dy - (bb*mT*mT)/aa/gr_dy

                 else

                 end if

                 mph_crmx = max(abs(cri),mph_crmx)
                 mph_crmn = min(abs(cri),mph_crmn)

              end if

          end do
       end do

end subroutine mph_PressureJumps
