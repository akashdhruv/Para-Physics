subroutine MPI_applyBC_RMA(local)

#include "Solver.h"

        use MPI_data
        use physicaldata, only : dataTARGET,&
                                 eastORIGIN,westORIGIN,northORIGIN,southORIGIN

        implicit none

        real, intent(inout), dimension(:,:) :: local

        integer :: status(MPI_STATUS_SIZE)
        integer(kind=MPI_ADDRESS_KIND) :: target_disp = 0 
        integer :: index1=Nyb+2, index2=Nyb+2+Nyb+2, index3=Nyb+2+Nyb+2+Nxb+2

        dataTARGET(1:Nyb+2)               = local(2,:)
        dataTARGET(index1+1:index1+Nyb+2) = local(Nxb+1,:)
        dataTARGET(index2+1:index2+Nxb+2) = local(:,2)
        dataTARGET(index3+1:index3+Nxb+2) = local(:,Nyb+1)

#ifdef MPI_RMA_ACTIVE
        call MPI_WIN_FENCE(0,RMA_win,ierr)

        !_______________________MPI BC for High X______________________________!
        if(x_id < x_procs - 1) then
          call MPI_GET(eastORIGIN(1),Nyb+2,MPI_REAL,myid+1,target_disp,Nyb+2,MPI_REAL,RMA_win,ierr)
        end if

        !_______________________MPI BC for Low X______________________________!
        if(x_id > 0) then
          call MPI_GET(westORIGIN(1),Nyb+2,MPI_REAL,myid-1,target_disp+index1,Nyb+2,MPI_REAL,RMA_win,ierr)
        end if
              
        !_______________________MPI BC for High Y______________________________!
        if(y_id < y_procs - 1) then
          call MPI_GET(northORIGIN(1),Nxb+2,MPI_REAL,myid+x_procs,target_disp+index2,Nxb+2,MPI_REAL,RMA_win,ierr)
        end if

        !_______________________MPI BC for Low Y______________________________!
        if(y_id > 0) then
          call MPI_GET(southORIGIN(1),Nxb+2,MPI_REAL,myid-x_procs,target_disp+index3,Nxb+2,MPI_REAL,RMA_win,ierr)
        end if

        call MPI_WIN_FENCE(0,RMA_win,ierr)
#endif

#ifdef MPI_RMA_PASSIVE
        call MPI_BARRIER(solver_comm,ierr)

        !_______________________MPI BC for High X______________________________!
        if(x_id < x_procs - 1) then
          call MPI_GET(eastORIGIN(1),Nyb+2,MPI_REAL,myid+1,target_disp,Nyb+2,MPI_REAL,RMA_win,ierr)
        end if

        !_______________________MPI BC for Low X______________________________!
        if(x_id > 0) then
          call MPI_GET(westORIGIN(1),Nyb+2,MPI_REAL,myid-1,target_disp+index1,Nyb+2,MPI_REAL,RMA_win,ierr)
        end if
              
        !_______________________MPI BC for High Y______________________________!
        if(y_id < y_procs - 1) then
          call MPI_GET(northORIGIN(1),Nxb+2,MPI_REAL,myid+x_procs,target_disp+index2,Nxb+2,MPI_REAL,RMA_win,ierr)
        end if

        !_______________________MPI BC for Low Y______________________________!
        if(y_id > 0) then
          call MPI_GET(southORIGIN(1),Nxb+2,MPI_REAL,myid-x_procs,target_disp+index3,Nxb+2,MPI_REAL,RMA_win,ierr)
        end if

        call MPI_BARRIER(solver_comm,ierr)
#endif

        local(Nxb+2,:) = eastORIGIN
        local(1,:)     = westORIGIN
        local(:,Nyb+2) = northORIGIN
        local(:,1)     = southORIGIN

end subroutine MPI_applyBC_RMA
