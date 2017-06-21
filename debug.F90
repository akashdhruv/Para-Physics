module MPI_interface
     contains
        subroutine MPI_applyBC(local_data,shared_data,myid,procs,solver_comm,shared_id,shared_procs,shared_comm,Nx,Ny)

        implicit none

        include "mpif.h"

        real,intent(inout), dimension(:,:) :: local_data,shared_data
        integer,intent(in) :: myid,procs,shared_id,shared_procs,Nx,Ny,solver_comm,shared_comm

        integer :: send_req,recv_req,ierr

        if(myid/shared_procs == (myid+1)/shared_procs) then

                if(myid < procs - 1) local_data(2,2) = shared_data(1,(shared_id+1)*Ny+1)

        else

                if(myid < procs - 1) call MPI_IRECV(local_data(2,2), 1, MPI_REAL, myid+1, 1, solver_comm, recv_req, ierr)

                if(myid/shared_procs /= (myid-1)/shared_procs .and. myid > 0) &
                call MPI_ISEND(local_data(1,1), 1, MPI_REAL, myid-1, 1, solver_comm, send_req, ierr)
 
        end if

        end subroutine MPI_applyBC
end module

program debug


        use MPI_interface
        use, intrinsic :: ISO_C_BINDING, ONLY: C_PTR, C_F_POINTER

        implicit none

        include "mpif.h"

        integer :: solver_comm, shared_comm, myid, procs
        integer :: shared_id,shared_procs,ierr
        integer :: resultlen
        character*(MPI_MAX_PROCESSOR_NAME) :: nameproc
        integer(kind=MPI_ADDRESS_KIND) :: sze
        integer :: disp_unit

        type(C_PTR) :: baseptr,shareptr
        integer :: win,i,j

        real,pointer,dimension(:,:) :: local_data,shared_data

        integer :: Nx,Ny

        real :: A

        integer :: status(MPI_STATUS_SIZE), send_req, recv_req

        Nx = 2
        Ny = 2

        solver_comm = MPI_COMM_WORLD

        call MPI_INIT(ierr)
        call MPI_COMM_RANK(solver_comm, myid, ierr)
        call MPI_COMM_SIZE(solver_comm, procs, ierr)

        call MPI_GET_PROCESSOR_NAME(nameproc,resultlen,ierr)

        call MPI_COMM_SPLIT_TYPE(solver_comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, shared_comm, ierr)

        call MPI_COMM_RANK(shared_comm,shared_id,ierr)
        call MPI_COMM_SIZE(shared_comm,shared_procs,ierr)

        !print *,"Global rank: ",myid,"Shared rank: ",shared_id,"Proc: ",nameproc(1:resultlen)
 
        sze       = Nx*Ny*sizeof(A)
        disp_unit = sizeof(A)

        call MPI_WIN_ALLOCATE_SHARED(sze,disp_unit,MPI_INFO_NULL,shared_comm,baseptr,win,ierr)

        call MPI_WIN_SHARED_QUERY(win, shared_id, sze, disp_unit, baseptr,ierr)
        call MPI_BARRIER(shared_comm,ierr)
        call C_F_POINTER(baseptr, local_data,[Nx,Ny])

        call MPI_WIN_SHARED_QUERY(win, 0 ,sze,disp_unit,baseptr,ierr)
        call MPI_BARRIER(shared_comm,ierr)
        call C_F_POINTER(baseptr,shared_data,[Nx,Ny*shared_procs])


        local_data(1,1) = myid + 1

        call MPI_BARRIER(shared_comm,ierr)

        call MPI_applyBC(local_data,shared_data,myid,procs,solver_comm,shared_id,shared_procs,shared_comm,Nx,Ny)
       
        call MPI_BARRIER(solver_comm,ierr)
        call MPI_BARRIER(shared_comm,ierr)

        print *,"Global Rank: ",myid," Shared Rank: ",shared_id," data(1,1): ",local_data(1,1)," data(2,2): ",local_data(2,2)

        call MPI_BARRIER(solver_comm, ierr)
        call MPI_BARRIER(shared_comm, ierr)

         if (shared_id == 0) then
         print *,"Data on rank: ",shared_id
          do j=1,Ny*shared_procs
              print *,shared_data(:,j)
          end do
         end if

        !call MPI_BARRIER(solver_comm,ierr)

        ! if (shared_id == 1) then

        ! print *,"Data on rank: ",shared_id
        !  do j=1,Ny*shared_procs
        !      print *,shared_data(:,j)
        !  end do
        !end if

        ! call MPI_BARRIER(shared_comm,ierr)

        ! if (shared_id == 2) then
        ! print *,"Data on rank: ",shared_id
        !  do j=1,Ny*shared_procs
        !      print *,shared_data(:,j)
        !  end do
       !end if

       ! call MPI_BARRIER(shared_comm,ierr)

       ! if (shared_id == 3) then
       ! print *,"Data on rank: ",shared_id
       !  do j=1,Ny*shared_procs
       !      print *,shared_data(:,j)
       !  end do
       !end if

       call MPI_BARRIER(shared_comm,ierr)
       call MPI_WIN_FREE(win,ierr)
       call MPI_COMM_FREE(shared_comm,ierr)
       call MPI_FINALIZE(ierr)

end program debug


