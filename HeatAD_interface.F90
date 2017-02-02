module HeatAD_interface

     implicit none

     interface
         subroutine HeatAD_init()
         end subroutine HeatAD_init
     end interface


     interface
         subroutine HeatAD_solver(tstep)
         implicit none
         integer, intent(in) :: tstep
         end subroutine HeatAD_solver
     end interface

     interface
         subroutine heat_tempSolver(tstep,T,T_old,mdot,smrh,u,v,a1x,a1y,a2x,a2y,s,pf,thco,cp)
         implicit none
         integer, intent(in) :: tstep
         real, intent(inout), dimension(:,:) :: T,u,v,a1x,a1y,s,pf,thco,cp,a2x,a2y,T_old,mdot,smrh
         end subroutine heat_tempSolver
     end interface

end module HeatAD_interface
