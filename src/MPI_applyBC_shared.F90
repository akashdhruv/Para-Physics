subroutine MPI_applyBC_shared(local,shared,myid,procs,solver_comm,world_prt,shared_id,shared_procs,shared_comm,shared_prt,Nx,Ny)

        implicit none

        include "mpif.h"

        real,intent(inout), dimension(:,:) :: local,shared
        integer,intent(in), dimension(:) :: world_prt,shared_prt
        integer,intent(in) :: myid,procs,shared_id,shared_procs,Nx,Ny,solver_comm,shared_comm

        integer :: send_req,recv_req,ierr

        !_____Check if right neighbour is on the same node__________!

        if(shared_prt(myid+1+1) /= MPI_UNDEFINED) then

         if(myid < procs - 1) local(2,2) = shared(1,(shared_id+1)*Ny+1) ! Write data from shared memory to local chunk

        else

         if(myid < procs - 1) call MPI_IRECV(local(2,2), 1, MPI_REAL, myid+1, 1, solver_comm, recv_req, ierr) ! Recieve data from other node

        end if

        !______Send data to left node _______________!
        if(shared_prt(myid-1+1) == MPI_UNDEFINED .and. myid > 0) &
        call MPI_ISEND(local(1,1), 1, MPI_REAL, myid-1, 1,solver_comm,send_req, ierr)

end subroutine MPI_applyBC_shared
