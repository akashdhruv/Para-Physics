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
           subroutine ins_predictor(tstep,u,v,ut,vt,g1_old,g2_old)
            implicit none
            integer, intent(in) :: tstep
            real, intent(inout), dimension(:,:) :: u, v
            real, intent(inout), dimension(:,:) :: ut,vt,g1_old,g2_old
           end subroutine ins_predictor
       end interface

       interface
           subroutine ins_predictor_VD(tstep,p_counter,p,u,v,ut,vt,visc,rho1x,rho1y,rho2x,rho2y,s,s2,sigp,sigx,sigy)
            implicit none
            integer, intent(in) :: tstep
            integer, intent(out) :: p_counter
            real, intent(inout), dimension(:,:) :: u, v, p, visc, rho1x, rho1y, rho2x, rho2y, s, s2, sigp, sigx, sigy, ut, vt
           end subroutine ins_predictor_VD
       end interface

       interface
           subroutine ins_vorticity(tstep,w,w_old,u,v,s)
            implicit none
            integer, intent(in) :: tstep
            real, intent(inout), dimension(:,:) :: w,u,v,s,w_old
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
