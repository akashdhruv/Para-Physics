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
         subroutine heat_tempSolver(tstep,T,u,v,s,pf,thco,cp)
         implicit none
         integer, intent(in) :: tstep
         real, intent(inout), dimension(:,:) :: T,u,v,s,pf,thco,cp
         end subroutine heat_tempSolver
     end interface

end module HeatAD_interface
