module modify

     contains
     subroutine modify_data(local_data)

#include "Solver.h"     

     use MPI_data

     implicit none
     real,intent(inout),dimension(:,:) :: local_data

     local_data(2:Nxb+1,2:Nyb+1) = myid + 1

     end subroutine

end module modify

program MPI_SHMtest

#include "Solver.h"

        use modify
        use MPI_interface, only: MPI_applyBC_shared
        use MPI_data
        use, intrinsic :: ISO_C_BINDING, ONLY: C_PTR, C_F_POINTER

        implicit none

        !include "mpif.h"

        integer(kind=MPI_ADDRESS_KIND) :: sze
        !integer :: disp_unit
        type(C_PTR) :: baseptr,shareptr
        integer :: win,i,j

        real,save,pointer,dimension(:,:) :: local_data,shared_data
        real :: A

        !_____Defined Block Size________!

        ! Get block size from Solver.h !

        !_________Define Global Communication Environment___________!
        solver_comm = MPI_COMM_WORLD

        call MPI_INIT(ierr)
        call MPI_COMM_RANK(solver_comm, myid, ierr)
        call MPI_COMM_SIZE(solver_comm, procs, ierr)
        call MPI_COMM_GROUP(solver_comm, world_group, ierr)

        allocate(world_part(procs))
        world_part = (/(I,I=0,procs-1)/)

        !_________Define Communication Based On Grid______________!
        call MPI_COMM_SPLIT(solver_comm,myid/nblockx,myid,x_comm,ierr)
        call MPI_COMM_SPLIT(solver_comm,mod(myid,nblockx),myid,y_comm,ierr)

        call MPI_COMM_RANK(x_comm,x_id,ierr)   
        call MPI_COMM_SIZE(x_comm,x_procs,ierr)

        call MPI_COMM_RANK(y_comm,y_id,ierr)
        call MPI_COMM_size(y_comm,y_procs,ierr)

        !__________Define Shared Communication Environment__________!

        call MPI_COMM_SPLIT_TYPE(solver_comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, shared_comm, ierr)
        call MPI_COMM_RANK(shared_comm,shared_id,ierr)
        call MPI_COMM_SIZE(shared_comm,shared_procs,ierr)
        call MPI_COMM_GROUP(shared_comm, shared_group,ierr)

        allocate(shared_part(procs))
        call MPI_GROUP_TRANSLATE_RANKS(world_group,procs,world_part,shared_group,shared_part,ierr)

        !________Make on-node processes allocate their chunk of shared Memory____________________!
        sze       = (Nxb+2)*(Nyb+2)*sizeof(A)
        disp_unit = sizeof(A)

        call MPI_WIN_ALLOCATE_SHARED(sze,disp_unit,MPI_INFO_NULL,shared_comm,baseptr,win,ierr)

        !__________________Point to local chunk of the shared data_______________________________!
        call MPI_WIN_SHARED_QUERY(win, shared_id, sze, disp_unit, baseptr,ierr)
        call MPI_BARRIER(shared_comm,ierr)
        call C_F_POINTER(baseptr, local_data,[Nxb+2,Nyb+2])

        !_____________________Point to the enitre shared data_____________________________!
        call MPI_WIN_SHARED_QUERY(win, 0 ,sze,disp_unit,baseptr,ierr)
        call MPI_BARRIER(shared_comm,ierr)
        call C_F_POINTER(baseptr,shared_data,[Nxb+2,(Nyb+2)*shared_procs])


        !___Modify Local Data____!
        call modify_data(local_data)

        call MPI_BARRIER(shared_comm,ierr)
        call MPI_BARRIER(solver_comm,ierr)

        !_____Exchange Information_____!
        call MPI_applyBC_shared(local_data,shared_data)

        call MPI_BARRIER(solver_comm,ierr)
        call MPI_BARRIER(shared_comm,ierr)

        !_______Print Results_____________!
        !if(myid == 0) print *,"world_part:",world_part
        !if(myid == 0) print *,"shared_part: ",shared_part
       
        !print *,"x_rank : ",x_id," y_rank: ", y_id," comm_x: ",mod(x_id-1+x_procs,x_procs)," comm_y: ",mod(y_id-1+y_procs,y_procs)
 
        !do i=1,procs

         if (myid == 27) then
         print *,"Global rank:",myid," Shared rank:",shared_id
          do j=1,Nyb+2
              print *,local_data(:,j)
          end do
         end if

        call MPI_BARRIER(solver_comm,ierr)
        call MPI_BARRIER(shared_comm,ierr)

         if (myid == 35) then
         print *,"Global rank:",myid," Shared rank:",shared_id
          do j=1,Nyb+2
              print *,local_data(:,j)
          end do
         end if

        call MPI_BARRIER(solver_comm,ierr)
        call MPI_BARRIER(shared_comm,ierr)

         if (myid == 43) then
         print *,"Global rank:",myid," Shared rank:",shared_id
          do j=1,Nyb+2
              print *,local_data(:,j)
          end do
         end if

        call MPI_BARRIER(solver_comm,ierr)
        call MPI_BARRIER(shared_comm,ierr)

        !end do

     !___Deallocate memory and clean up_____!
       deallocate(world_part)
       deallocate(shared_part)

       call MPI_BARRIER(shared_comm,ierr)
       call MPI_WIN_FREE(win,ierr)
       call MPI_COMM_FREE(shared_comm,ierr)
       call MPI_FINALIZE(ierr)

end program MPI_SHMtest


