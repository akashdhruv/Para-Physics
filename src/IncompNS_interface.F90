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
           subroutine ins_predictor_VD(tstep,u,v,ut,vt,g1_old,g2_old,visc,rho1x,rho1y,rho2x,rho2y)
            implicit none
            integer, intent(in) :: tstep
            real, intent(inout), dimension(:,:) :: u, v, visc, rho1x, rho1y, rho2x, rho2y, &
                                                   ut, vt,g1_old,g2_old
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
           subroutine ins_computeQinout(u,v,x,y,flg)
           implicit none
           real, dimension(:,:,:), intent(in) :: u,v,x,y
           logical, intent(in) :: flg
           end subroutine ins_computeQinout
       end interface


       interface 
           subroutine ins_rescaleVel(u,v,x,y)
            implicit none
            real, dimension(:,:,:), intent(inout) :: u,v
            real, dimension(:,:,:), intent(in) :: x,y
           end subroutine ins_rescaleVel
       end interface

       interface 
           subroutine ins_convVelout(u,v,x,y)
           implicit none
           real, dimension(:,:,:), intent(in) :: u,v,x,y
           end subroutine ins_convVelout
       end interface

end module IncompNS_interface
