subroutine MPI_applyBC_RMA(local)

#include "Solver.h"

        use MPI_data
        use physicaldata, only : eastTARGET,westTARGET,northTARGET,southTARGET,&
                                 eastORIGIN,westORIGIN,northORIGIN,southORIGIN

        implicit none

        real, intent(inout), dimension(:,:) :: local

        integer :: status(MPI_STATUS_SIZE)
        integer(kind=MPI_ADDRESS_KIND) :: target_disp = 0 

        eastTARGET  = local(2,:)
        westTARGET  = local(Nxb+1,:)
        northTARGET = local(:,2)
        southTARGET = local(:,Nyb+1)

        call MPI_WIN_FENCE(0,north_win,ierr)
        call MPI_WIN_FENCE(0,south_win,ierr)
        call MPI_WIN_FENCE(0,east_win,ierr)
        call MPI_WIN_FENCE(0,west_win,ierr)

        !_______________________MPI BC for High X______________________________!
        if(x_id < x_procs - 1) then
          call MPI_GET(eastORIGIN(1),Nyb+2,MPI_REAL,myid+1,target_disp,Nyb+2,MPI_REAL,east_win,ierr)
        end if

        !_______________________MPI BC for Low X______________________________!
        if(x_id > 0) then
          call MPI_GET(westORIGIN(1),Nyb+2,MPI_REAL,myid-1,target_disp,Nyb+2,MPI_REAL,west_win,ierr)
        end if
              
        !_______________________MPI BC for High Y______________________________!
        if(y_id < y_procs - 1) then
          call MPI_GET(northORIGIN(1),Nxb+2,MPI_REAL,myid+x_procs,target_disp,Nxb+2,MPI_REAL,north_win,ierr)
        end if

        !_______________________MPI BC for Low Y______________________________!
        if(y_id > 0) then
          call MPI_GET(southORIGIN(1),Nxb+2,MPI_REAL,myid-x_procs,target_disp,Nxb+2,MPI_REAL,south_win,ierr)
        end if

        call MPI_WIN_FENCE(0,north_win,ierr)
        call MPI_WIN_FENCE(0,south_win,ierr)
        call MPI_WIN_FENCE(0,east_win,ierr)
        call MPI_WIN_FENCE(0,west_win,ierr)

        local(Nxb+2,:) = eastORIGIN
        local(1,:)     = westORIGIN
        local(:,Nyb+2) = northORIGIN
        local(:,1)     = southORIGIN

end subroutine MPI_applyBC_RMA
