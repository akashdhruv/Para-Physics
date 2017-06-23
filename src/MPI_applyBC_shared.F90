subroutine MPI_applyBC_shared(local,shared)

#include "Solver.h"

        use MPI_data

        implicit none

        real,intent(inout), dimension(:,:) :: local,shared

        integer :: status(MPI_STATUS_SIZE)
        integer :: sendreq, recvreq

        !_______________________MPI BC for High X______________________________!
        if(x_id < x_procs - 1) then

                if(shared_part(myid+1+1) /= MPI_UNDEFINED) then

                        local(Nxb+2,:) = shared(2,(shared_id+1)*(Nyb+2)+1:(shared_id+1)*(Nyb+2)+(Nyb+2)+1)
                else
  
                        call MPI_RECV(local(Nxb+2,:), Nyb+2, MPI_REAL, x_id+1, 1, x_comm, status, ierr)

                end if
        end if

        if(shared_part(myid+1-1) == MPI_UNDEFINED .and. x_id>0) &
        call MPI_ISEND(local(2,:), Nyb+2, MPI_REAL, x_id-1, 1, x_comm, sendreq, ierr)

        !_______________________MPI BC for Low X______________________________!
        if(x_id > 0) then

                if(shared_part(myid+1-1) /= MPI_UNDEFINED) then

                        local(1,:) = shared(Nxb+1,(shared_id-1)*(Nyb+2)+1:(shared_id-1)*(Nyb+2)+(Nyb+2)+1)
                else
                        call MPI_RECV(local(1,:), Nyb+2, MPI_REAL, x_id-1, 2, x_comm, status, ierr)

                end if
        end if

        if(shared_part(myid+1+1) == MPI_UNDEFINED .and. x_id < x_procs-1) &
        call MPI_ISEND(local(Nxb+1,:), Nyb+2, MPI_REAL, x_id+1, 2, x_comm, sendreq, ierr)
                
        !_______________________MPI BC for High Y______________________________!
        if(y_id < y_procs - 1) then

                if(shared_part(myid+1+x_procs) /= MPI_UNDEFINED) then

                        local(:,Nyb+2) = shared(:,(shared_id+x_procs)*(Nyb+2)+2)
                else
                        call MPI_RECV(local(:,Nyb+2), Nxb+2, MPI_REAL, y_id+1, 3, y_comm, status, ierr)


                end if
        end if

        if(shared_part(myid+1-x_procs) == MPI_UNDEFINED .and. y_id>0) &
        call MPI_ISEND(local(:,2), Nxb+2, MPI_REAL, y_id-1, 3, y_comm, sendreq, ierr)

       !_______________________MPI BC for Low Y______________________________!
       if(y_id > 0) then

                if(shared_part(myid+1-x_procs) /= MPI_UNDEFINED) then

                        local(:,1) = shared(:,(shared_id-x_procs)*(Nyb+2)+Nyb+1)
                else
                        call MPI_RECV(local(:,1), Nxb+2, MPI_REAL, y_id-1, 4, y_comm, status, ierr)

                end if
        end if

        if(shared_part(myid+1+x_procs) == MPI_UNDEFINED .and. y_id < y_procs -1) &
        call MPI_ISEND(local(:,Nyb+1), Nxb+2, MPI_REAL, y_id+1, 4 ,y_comm, sendreq, ierr)
        
        call MPI_BARRIER(shared_comm,ierr)

end subroutine MPI_applyBC_shared
