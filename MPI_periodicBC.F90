subroutine MPI_periodicBC(u_ex,v_ex,aDIM)

#include "Solver.h"

     use MPI_data

     implicit none
     
     include "mpif.h"

     real, dimension(:,:) , intent(inout) :: u_ex,v_ex
     integer, intent(in) :: aDIM

     integer :: status(MPI_STATUS_SIZE)
     
     select case (aDIM)

           case (1)

              if(x_id == x_procs - 1) then

                call MPI_SENDRECV(u_ex(Nxb,:),   Nyb+2, MPI_REAL, 0, 1,&
                                  u_ex(Nxb+1,:), Nyb+2, MPI_REAL, 0, 2, x_comm, status, ierr)

                call MPI_SENDRECV(v_ex(Nxb+1,:), Nyb+2, MPI_REAL, 0, 3,&
                                  v_ex(Nxb+2,:), Nyb+2, MPI_REAL, 0, 4, x_comm, status, ierr)

             else if (x_id == 0) then

              call MPI_SENDRECV(u_ex(2,:), Nyb+2, MPI_REAL, x_procs - 1, 2,&
                                u_ex(1,:), Nyb+2, MPI_REAL, x_procs - 1, 1,x_comm,status, ierr)

              call MPI_SENDRECV(v_ex(2,:), Nyb+2, MPI_REAL, x_procs - 1, 4,&
                                v_ex(1,:), Nyb+2, MPI_REAL, x_procs - 1, 3, x_comm, status, ierr)

             end if

           case (2)

              if(y_id == y_procs - 1) then

                call MPI_SENDRECV(v_ex(:,Nyb),   Nxb+2, MPI_REAL, 0, 1,&
                                  v_ex(:,Nyb+1), Nxb+2, MPI_REAL, 0, 2, y_comm, status, ierr)

                call MPI_SENDRECV(u_ex(:,Nyb+1), Nxb+2, MPI_REAL, 0, 3,&
                                  u_ex(:,Nyb+2), Nxb+2, MPI_REAL, 0, 4, y_comm, status, ierr)

             else if (y_id == 0) then

              call MPI_SENDRECV(v_ex(:,2), Nxb+2, MPI_REAL, y_procs - 1, 2,&
                                v_ex(:,1), Nxb+2, MPI_REAL, y_procs - 1, 1,y_comm,status, ierr)

              call MPI_SENDRECV(u_ex(:,2), Nxb+2, MPI_REAL, y_procs - 1, 4,&
                                u_ex(:,1), Nxb+2, MPI_REAL, y_procs - 1, 3,y_comm, status, ierr)

             end if

           case default
       
              print *,"MPI_periodicBC: The solver supports only two dimensions for now."
    
      end select

      !call MPI_BARRIER(solver_comm,ierr)

end subroutine MPI_periodicBC
