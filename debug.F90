module MPI_interface
     contains
        subroutine MPI_applyBC(local,shared,myid,procs,solver_comm,world_prt,shared_id,shared_procs,shared_comm,shared_prt,Nx,Ny)

        implicit none

        include "mpif.h"

        real,intent(inout), dimension(:,:) :: local,shared
        integer,intent(in), dimension(:) :: world_prt,shared_prt
        integer,intent(in) :: myid,procs,shared_id,shared_procs,Nx,Ny,solver_comm,shared_comm

        integer :: send_req,recv_req,ierr

        if(shared_prt(myid+1) /= MPI_UNDEFINED) then

                if(myid < procs - 1) local(2,2) = shared(1,(shared_id+1)*Ny+1)

        else

                if(myid < procs - 1) call MPI_IRECV(local(2,2), 1, MPI_REAL, myid+1, 1, solver_comm, recv_req, ierr)

        end if

        if(shared_prt(myid-1) == MPI_UNDEFINED .and. myid > 0) &
        call MPI_ISEND(local(1,1), 1, MPI_REAL, myid-1, 1, solver_comm,send_req, ierr)

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

        integer, allocatable,dimension(:) :: world_part, shared_part

        integer :: Nx,Ny

        real :: A

        integer :: status(MPI_STATUS_SIZE), send_req, recv_req, world_grp, shared_grp

        Nx = 2
        Ny = 2

        solver_comm = MPI_COMM_WORLD

        call MPI_INIT(ierr)
        call MPI_COMM_RANK(solver_comm, myid, ierr)
        call MPI_COMM_SIZE(solver_comm, procs, ierr)

        call MPI_COMM_SPLIT_TYPE(solver_comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, shared_comm, ierr)

        call MPI_COMM_RANK(shared_comm,shared_id,ierr)
        call MPI_COMM_SIZE(shared_comm,shared_procs,ierr)

        allocate(world_part(procs))
        allocate(shared_part(procs))

        world_part = (/(I,I=0,procs-1)/)

        call MPI_COMM_GROUP(solver_comm, world_grp, ierr)
        call MPI_COMM_GROUP(shared_comm, shared_grp,ierr)
        call MPI_GROUP_TRANSLATE_RANKS(world_grp,procs,world_part,shared_grp,shared_part,ierr)

        !print *,"Global rank: ",myid,"Shared rank: ",shared_id,"Proc: ",nameproc(1:resultlen)
 
        !if(myid == 0) print *,world_part
        !if(myid == 0) print *,shared_part

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

        call MPI_applyBC(local_data,shared_data,myid,procs,solver_comm,world_part,&
                        shared_id,shared_procs,shared_comm,shared_part,Nx,Ny)
       
        call MPI_BARRIER(solver_comm,ierr)
        call MPI_BARRIER(shared_comm,ierr)

        print *,"Global Rank: ",myid," Shared Rank: ",shared_id," data(1,1): ",local_data(1,1)," data(2,2): ",local_data(2,2)

        call MPI_BARRIER(solver_comm, ierr)
        call MPI_BARRIER(shared_comm, ierr)

        ! if (shared_id == 0) then
        ! print *,"Data on rank: ",shared_id
        !  do j=1,Ny*shared_procs
        !      print *,shared_data(:,j)
        !  end do
        ! end if

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

       deallocate(world_part)
       deallocate(shared_part)

       call MPI_BARRIER(shared_comm,ierr)
       call MPI_WIN_FREE(win,ierr)
       call MPI_COMM_FREE(shared_comm,ierr)
       call MPI_FINALIZE(ierr)

end program debug


