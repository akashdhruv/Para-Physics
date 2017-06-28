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

                        localCENTER(ivar,Nxb+2,:) = eastCENTER(ivar,2,:)
                else
  
                        call MPI_SENDRECV(localCENTER(ivar,Nxb+1,:), Nyb+2, MPI_REAL, x_id+1, 2,&
                                          localCENTER(ivar,Nxb+2,:), Nyb+2, MPI_REAL, x_id+1, 1, x_comm, status, ierr)

                end if
                end if

                !_______________________MPI BC for Low X______________________________!
                if(x_id > 0) then

                if(shared_part(myid+1-1) /= MPI_UNDEFINED) then

                        localCENTER(ivar,1,:) = westCENTER(ivar,Nxb+1,:)
                else
                        call MPI_SENDRECV(localCENTER(ivar,2,:), Nyb+2, MPI_REAL, x_id-1, 1,&
                                          localCENTER(ivar,1,:), Nyb+2, MPI_REAL, x_id-1, 2, x_comm, status, ierr)

                end if
                end if
              
                !_______________________MPI BC for High Y______________________________!
                if(y_id < y_procs - 1) then

                if(shared_part(myid+1+x_procs) /= MPI_UNDEFINED) then

                        localCENTER(ivar,:,Nyb+2) = northCENTER(ivar,:,2)
                else
                        call MPI_SENDRECV(localCENTER(ivar,:,Nyb+1), Nxb+2, MPI_REAL, y_id+1, 4,& 
                                          localCENTER(ivar,:,Nyb+2), Nxb+2, MPI_REAL, y_id+1, 3, y_comm, status, ierr)


                end if
                end if

                !_______________________MPI BC for Low Y______________________________!
                if(y_id > 0) then

                if(shared_part(myid+1-x_procs) /= MPI_UNDEFINED) then

                        localCENTER(ivar,:,1) = southCENTER(ivar,:,Nyb+1)
                else
                        call MPI_SENDRECV(localCENTER(ivar,:,2), Nxb+2, MPI_REAL, y_id-1, 3,&
                                          localCENTER(ivar,:,1), Nxb+2, MPI_REAL, y_id-1, 4, y_comm, status, ierr)

                end if
                end if

                !call MPI_WIN_FENCE(0,center_win,ierr)       
                return

        else if(datatype == FACEX) then

                !call MPI_WIN_FENCE(0,facex_win,ierr)

                !_______________________MPI BC for High X______________________________!
                if(x_id < x_procs - 1) then

                if(shared_part(myid+1+1) /= MPI_UNDEFINED) then

                        localFACEX(ivar,Nxb+2,:) = eastFACEX(ivar,2,:)
                else
  
                        call MPI_SENDRECV(localFACEX(ivar,Nxb+1,:), Nyb+2, MPI_REAL, x_id+1, 2,&
                                          localFACEX(ivar,Nxb+2,:), Nyb+2, MPI_REAL, x_id+1, 1, x_comm, status, ierr)

                end if
                end if

                !_______________________MPI BC for Low X______________________________!
                if(x_id > 0) then

                if(shared_part(myid+1-1) /= MPI_UNDEFINED) then

                        localFACEX(ivar,1,:) = westFACEX(ivar,Nxb+1,:)
                else
                        call MPI_SENDRECV(localFACEX(ivar,2,:), Nyb+2, MPI_REAL, x_id-1, 1,&
                                          localFACEX(ivar,1,:), Nyb+2, MPI_REAL, x_id-1, 2, x_comm, status, ierr)

                end if
                end if
              
                !_______________________MPI BC for High Y______________________________!
                if(y_id < y_procs - 1) then

                if(shared_part(myid+1+x_procs) /= MPI_UNDEFINED) then

                        localFACEX(ivar,:,Nyb+2) = northFACEX(ivar,:,2)
                else
                        call MPI_SENDRECV(localFACEX(ivar,:,Nyb+1), Nxb+2, MPI_REAL, y_id+1, 4,& 
                                          localFACEX(ivar,:,Nyb+2), Nxb+2, MPI_REAL, y_id+1, 3, y_comm, status, ierr)


                end if
                end if

                !_______________________MPI BC for Low Y______________________________!
                if(y_id > 0) then

                if(shared_part(myid+1-x_procs) /= MPI_UNDEFINED) then

                        localFACEX(ivar,:,1) = southFACEX(ivar,:,Nyb+1)
                else
                        call MPI_SENDRECV(localFACEX(ivar,:,2), Nxb+2, MPI_REAL, y_id-1, 3,&
                                          localFACEX(ivar,:,1), Nxb+2, MPI_REAL, y_id-1, 4, y_comm, status, ierr)

                end if
                end if

                !call MPI_WIN_FENCE(0,facex_win,ierr)
                return


        else if(datatype == FACEY) then

                !call MPI_WIN_FENCE(0,facey_win,ierr)

                !_______________________MPI BC for High X______________________________!
                if(x_id < x_procs - 1) then

                if(shared_part(myid+1+1) /= MPI_UNDEFINED) then

                        localFACEY(ivar,Nxb+2,:) = eastFACEY(ivar,2,:)
                else
  
                        call MPI_SENDRECV(localFACEY(ivar,Nxb+1,:), Nyb+2, MPI_REAL, x_id+1, 2,&
                                          localFACEY(ivar,Nxb+2,:), Nyb+2, MPI_REAL, x_id+1, 1, x_comm, status, ierr)

                end if
                end if

                !_______________________MPI BC for Low X______________________________!
                if(x_id > 0) then

                if(shared_part(myid+1-1) /= MPI_UNDEFINED) then

                        localFACEY(ivar,1,:) = westFACEY(ivar,Nxb+1,:)
                else
                        call MPI_SENDRECV(localFACEY(ivar,2,:), Nyb+2, MPI_REAL, x_id-1, 1,&
                                          localFACEY(ivar,1,:), Nyb+2, MPI_REAL, x_id-1, 2, x_comm, status, ierr)

                end if
                end if
              
                !_______________________MPI BC for High Y______________________________!
                if(y_id < y_procs - 1) then

                if(shared_part(myid+1+x_procs) /= MPI_UNDEFINED) then

                        localFACEY(ivar,:,Nyb+2) = northFACEY(ivar,:,2)
                else
                        call MPI_SENDRECV(localFACEY(ivar,:,Nyb+1), Nxb+2, MPI_REAL, y_id+1, 4,& 
                                          localFACEY(ivar,:,Nyb+2), Nxb+2, MPI_REAL, y_id+1, 3, y_comm, status, ierr)


                end if
                end if

                !_______________________MPI BC for Low Y______________________________!
                if(y_id > 0) then

                if(shared_part(myid+1-x_procs) /= MPI_UNDEFINED) then

                        localFACEY(ivar,:,1) = southFACEY(ivar,:,Nyb+1)
                else
                        call MPI_SENDRECV(localFACEY(ivar,:,2), Nxb+2, MPI_REAL, y_id-1, 3,&
                                          localFACEY(ivar,:,1), Nxb+2, MPI_REAL, y_id-1, 4, y_comm, status, ierr)

                end if
                end if

                !call MPI_WIN_FENCE(0,facey_win,ierr)       
                return

        end if

end subroutine MPI_applyBC_shared
