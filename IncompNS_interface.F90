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
