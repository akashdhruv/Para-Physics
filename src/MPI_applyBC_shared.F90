subroutine MPI_applyBC_shared(local,shared)

#include "Solver.h"

        use MPI_data

        implicit none

        include "mpif.h"

        real,intent(inout), dimension(:,:) :: local,shared

        integer :: send_req,recv_req

        if(x_id < x_procs - 1) then


                if(shared_part(myid+1+1) /= MPI_UNDEFINED) then

                        local(Nxb+2,:) = shared(2,(shared_id+1)*(Nyb+2)+1:(shared_id+1)*(Nyb+2)+(Nyb+2)+1)

                else

                        call MPI_IRECV(local(Nxb+2,:), Nyb+2, MPI_REAL, x_id+1, 1, x_comm, recv_req, ierr)

                end if

        end if

        if(shared_part(myid-1+1) == MPI_UNDEFINED .and. x_id>0) &
        call MPI_ISEND(local(2,:), Nyb+2, MPI_REAL, x_id-1, 1, x_comm, send_req,ierr)


        if(x_id > 0) then

                if(shared_part(myid-1+1) /= MPI_UNDEFINED) then

                        local(1,:) = shared(Nxb+1,(shared_id-1)*(Nyb+2)+1:(shared_id-1)*(Nyb+2)+(Nyb+2)+1)


                else

                        call MPI_IRECV(local(1,:), Nyb+2, MPI_REAL, x_id-1, 2, x_comm, recv_req, ierr)

                end if

        end if

        if(shared_part(myid+1+1) == MPI_UNDEFINED .and. x_id < x_procs-1) &
        call MPI_ISEND(local(Nxb+1,:), Nyb+2, MPI_REAL, x_id+1, 2, x_comm, send_req,ierr)

end subroutine MPI_applyBC_shared
