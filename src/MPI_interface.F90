module MPI_interface
#include "Solver.h"

    implicit none

    interface
          subroutine MPIsolver_init()
          implicit none
          end subroutine MPIsolver_init
    end interface

    interface 
          subroutine MPIsolver_finalize(sim_Complete)
          implicit none
          logical, intent(in) :: sim_Complete
          end subroutine MPIsolver_finalize
    end interface

    interface
          subroutine MPI_applyBC(ivar,datatype)
          implicit none
          integer, intent(in) :: ivar,datatype
          end subroutine MPI_applyBC
    end interface 

    interface
          subroutine MPI_applyBC_DIS(local)
          implicit none
          real, intent(inout), dimension(:,:,:) :: local
          end subroutine MPI_applyBC_DIS
    end interface

    interface
          subroutine MPI_applyBC_SHM(local,sharedData)
          implicit none
          real, intent(inout), dimension(:,:,:) :: local
          real, intent(inout), dimension(:,:,:) :: sharedData
          end subroutine MPI_applyBC_SHM
    end interface

    interface
          subroutine MPI_applyBC_RMA(local)
          implicit none
          real, intent(inout), dimension(:,:,:) :: local
          end subroutine MPI_applyBC_RMA
    end interface

    interface
          subroutine MPI_applyBC_ORIG(local)
          implicit none
          real, dimension(:,:),intent(inout) :: local
          end subroutine MPI_applyBC_ORIG
    end interface

    interface 
       subroutine MPI_CollectResiduals(res,res1,collect_type)
       implicit none
       real, intent(inout) :: res,res1
       integer,intent(in) :: collect_type
       end subroutine MPI_CollectResiduals
    end interface

    interface 
        subroutine MPI_physicalBC_vel(u_ex,v_ex,x,y)
        implicit none
        real, dimension(:,:,:), intent(inout) :: u_ex, v_ex
        real, dimension(:,:,:), intent(out) :: x,y
        end subroutine MPI_physicalBC_vel
    end interface

    interface
        subroutine MPI_physicalBC_pres(p_ex,x,y)
        implicit none
        real, dimension(:,:,:), intent(inout) :: p_ex
        real, dimension(:,:,:), intent(in) :: x,y
        end subroutine MPI_physicalBC_pres
    end interface

   interface
        subroutine MPI_physicalBC_temp(t_ex,x,y)
        implicit none
        real, dimension(:,:,:), intent(inout) :: t_ex
        real, dimension(:,:,:), intent(in) :: x,y
        end subroutine MPI_physicalBC_temp
   end interface

   interface
        subroutine MPI_physicalBC_dfun(d_ex)
        implicit none
        real, dimension(:,:,:), intent(inout) :: d_ex
        end subroutine MPI_physicalBC_dfun
   end interface

   interface
        subroutine MPI_physicalBC_vort(d_ex)
        implicit none
        real, dimension(:,:,:), intent(inout) :: d_ex
        end subroutine MPI_physicalBC_vort
   end interface

   interface
        subroutine MPI_periodicBC(u_ex,v_ex,aDIM)
        implicit none
        real, dimension(:,:,:), intent(inout) :: u_ex,v_ex
        integer,intent(in) :: aDIM
        end subroutine MPI_periodicBC
   end interface

end module MPI_interface
