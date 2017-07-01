subroutine MPI_applyBC_shared(ivar,datatype)

#include "Solver.h"

        use MPI_data
        use physicaldata, only: localCENTER,localFACEX,localFACEY,&
                                northCENTER,southCENTER,eastCENTER,westCENTER,&
                                northFACEX,southFACEX,eastFACEX,westFACEX,&
                                northFACEY,southFACEY,eastFACEY,westFACEY

        implicit none

        integer, intent(in) :: ivar,datatype

        integer :: status(MPI_STATUS_SIZE)

        if(datatype == CENTER) then

                !call MPI_WIN_FENCE(0,center_win,ierr)

                !_______________________MPI BC for High X______________________________!
                if(x_id < x_procs - 1) then

                if(shared_part(myid+1+1) /= MPI_UNDEFINED) then

                        localCENTER(Nxb+2,:,ivar) = eastCENTER(2,:,ivar)
                else
  
                        call MPI_SENDRECV(localCENTER(Nxb+1,:,ivar), Nyb+2, MPI_REAL, x_id+1, 2,&
                                          localCENTER(Nxb+2,:,ivar), Nyb+2, MPI_REAL, x_id+1, 1, x_comm, status, ierr)

                end if
                end if

                !_______________________MPI BC for Low X______________________________!
                if(x_id > 0) then

                if(shared_part(myid+1-1) /= MPI_UNDEFINED) then

                        localCENTER(1,:,ivar) = westCENTER(Nxb+1,:,ivar)
                else
                        call MPI_SENDRECV(localCENTER(2,:,ivar), Nyb+2, MPI_REAL, x_id-1, 1,&
                                          localCENTER(1,:,ivar), Nyb+2, MPI_REAL, x_id-1, 2, x_comm, status, ierr)

                end if
                end if
              
                !_______________________MPI BC for High Y______________________________!
                if(y_id < y_procs - 1) then

                if(shared_part(myid+1+x_procs) /= MPI_UNDEFINED) then

                        localCENTER(:,Nyb+2,ivar) = northCENTER(:,2,ivar)
                else
                        call MPI_SENDRECV(localCENTER(:,Nyb+1,ivar), Nxb+2, MPI_REAL, y_id+1, 4,& 
                                          localCENTER(:,Nyb+2,ivar), Nxb+2, MPI_REAL, y_id+1, 3, y_comm, status, ierr)


                end if
                end if

                !_______________________MPI BC for Low Y______________________________!
                if(y_id > 0) then

                if(shared_part(myid+1-x_procs) /= MPI_UNDEFINED) then

                        localCENTER(:,1,ivar) = southCENTER(:,Nyb+1,ivar)
                else
                        call MPI_SENDRECV(localCENTER(:,2,ivar), Nxb+2, MPI_REAL, y_id-1, 3,&
                                          localCENTER(:,1,ivar), Nxb+2, MPI_REAL, y_id-1, 4, y_comm, status, ierr)

                end if
                end if

                !call MPI_WIN_FENCE(0,center_win,ierr)       
                return

        else if(datatype == FACEX) then

                !call MPI_WIN_FENCE(0,facex_win,ierr)

                !_______________________MPI BC for High X______________________________!
                if(x_id < x_procs - 1) then

                if(shared_part(myid+1+1) /= MPI_UNDEFINED) then

                        localFACEX(Nxb+2,:,ivar) = eastFACEX(2,:,ivar)
                else
  
                        call MPI_SENDRECV(localFACEX(Nxb+1,:,ivar), Nyb+2, MPI_REAL, x_id+1, 2,&
                                          localFACEX(Nxb+2,:,ivar), Nyb+2, MPI_REAL, x_id+1, 1, x_comm, status, ierr)

                end if
                end if

                !_______________________MPI BC for Low X______________________________!
                if(x_id > 0) then

                if(shared_part(myid+1-1) /= MPI_UNDEFINED) then

                        localFACEX(1,:,ivar) = westFACEX(Nxb+1,:,ivar)
                else
                        call MPI_SENDRECV(localFACEX(2,:,ivar), Nyb+2, MPI_REAL, x_id-1, 1,&
                                          localFACEX(1,:,ivar), Nyb+2, MPI_REAL, x_id-1, 2, x_comm, status, ierr)

                end if
                end if
              
                !_______________________MPI BC for High Y______________________________!
                if(y_id < y_procs - 1) then

                if(shared_part(myid+1+x_procs) /= MPI_UNDEFINED) then

                        localFACEX(:,Nyb+2,ivar) = northFACEX(:,2,ivar)
                else
                        call MPI_SENDRECV(localFACEX(:,Nyb+1,ivar), Nxb+2, MPI_REAL, y_id+1, 4,& 
                                          localFACEX(:,Nyb+2,ivar), Nxb+2, MPI_REAL, y_id+1, 3, y_comm, status, ierr)


                end if
                end if

                !_______________________MPI BC for Low Y______________________________!
                if(y_id > 0) then

                if(shared_part(myid+1-x_procs) /= MPI_UNDEFINED) then

                        localFACEX(:,1,ivar) = southFACEX(:,Nyb+1,ivar)
                else
                        call MPI_SENDRECV(localFACEX(:,2,ivar), Nxb+2, MPI_REAL, y_id-1, 3,&
                                          localFACEX(:,1,ivar), Nxb+2, MPI_REAL, y_id-1, 4, y_comm, status, ierr)

                end if
                end if

                !call MPI_WIN_FENCE(0,facex_win,ierr)
                return


        else if(datatype == FACEY) then

                !call MPI_WIN_FENCE(0,facey_win,ierr)

                !_______________________MPI BC for High X______________________________!
                if(x_id < x_procs - 1) then

                if(shared_part(myid+1+1) /= MPI_UNDEFINED) then

                        localFACEY(Nxb+2,:,ivar) = eastFACEY(2,:,ivar)
                else
  
                        call MPI_SENDRECV(localFACEY(Nxb+1,:,ivar), Nyb+2, MPI_REAL, x_id+1, 2,&
                                          localFACEY(Nxb+2,:,ivar), Nyb+2, MPI_REAL, x_id+1, 1, x_comm, status, ierr)

                end if
                end if

                !_______________________MPI BC for Low X______________________________!
                if(x_id > 0) then

                if(shared_part(myid+1-1) /= MPI_UNDEFINED) then

                        localFACEY(1,:,ivar) = westFACEY(Nxb+1,:,ivar)
                else
                        call MPI_SENDRECV(localFACEY(2,:,ivar), Nyb+2, MPI_REAL, x_id-1, 1,&
                                          localFACEY(1,:,ivar), Nyb+2, MPI_REAL, x_id-1, 2, x_comm, status, ierr)

                end if
                end if
              
                !_______________________MPI BC for High Y______________________________!
                if(y_id < y_procs - 1) then

                if(shared_part(myid+1+x_procs) /= MPI_UNDEFINED) then

                        localFACEY(:,Nyb+2,ivar) = northFACEY(:,2,ivar)
                else
                        call MPI_SENDRECV(localFACEY(:,Nyb+1,ivar), Nxb+2, MPI_REAL, y_id+1, 4,& 
                                          localFACEY(:,Nyb+2,ivar), Nxb+2, MPI_REAL, y_id+1, 3, y_comm, status, ierr)


                end if
                end if

                !_______________________MPI BC for Low Y______________________________!
                if(y_id > 0) then

                if(shared_part(myid+1-x_procs) /= MPI_UNDEFINED) then

                        localFACEY(:,1,ivar) = southFACEY(:,Nyb+1,ivar)
                else
                        call MPI_SENDRECV(localFACEY(:,2,ivar), Nxb+2, MPI_REAL, y_id-1, 3,&
                                          localFACEY(:,1,ivar), Nxb+2, MPI_REAL, y_id-1, 4, y_comm, status, ierr)

                end if
                end if

                !call MPI_WIN_FENCE(0,facey_win,ierr)       
                return

        end if

end subroutine MPI_applyBC_shared
