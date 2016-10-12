subroutine IncompNS_solver(tstep,p_counter)

    use IncompNS_interface, only: ins_momentum,ins_vorticity

    implicit none

    integer, intent(in) :: tstep
    integer, intent(out) :: p_counter

    call ins_momentum(tstep,p_counter)
    call ins_vorticity(tstep)

end subroutine IncompNS_solver
