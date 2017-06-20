module IncompNS_interface

       implicit none

       interface
           subroutine IncompNS_init()
           end subroutine IncompNS_init
       end interface

       interface
           subroutine IncompNS_solver(tstep,p_counter)
            implicit none
            integer, intent(in) :: tstep
            integer, intent(out) :: p_counter
           end subroutine IncompNS_solver
       end interface

       interface
           subroutine ins_momentum(tstep,p_counter,p,u,v,s,s2)
            implicit none
            integer, intent(in) :: tstep
            integer, intent(out) :: p_counter
            real, intent(inout), dimension(:,:) :: u, v, p, s, s2
           end subroutine ins_momentum
       end interface

       interface
           subroutine ins_momentum_VD(tstep,p_counter,p,u,v,visc,rho1x,rho1y,rho2x,rho2y,s,s2,sigp,sigx,sigy)
            implicit none
            integer, intent(in) :: tstep
            integer, intent(out) :: p_counter
            real, intent(inout), dimension(:,:) :: u, v, p, visc, rho1x, rho1y, rho2x, rho2y, s, s2, sigp, sigx, sigy
           end subroutine ins_momentum_VD
       end interface

       interface
           subroutine ins_vorticity(tstep,w,u,v,s)
            implicit none
            integer, intent(in) :: tstep
            real, intent(inout), dimension(:,:) :: w,u,v,s
           end subroutine ins_vorticity
       end interface

       interface
           subroutine ins_computeQinout()
           end subroutine
       end interface


       interface 
           subroutine ins_rescaleVel(u,v)
            implicit none
            real, dimension(:,:), intent(inout) :: u,v
           end subroutine
       end interface

       interface 
           subroutine ins_convVelout()
           end subroutine
       end interface

end module IncompNS_interface
