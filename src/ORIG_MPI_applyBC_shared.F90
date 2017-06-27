subroutine MPI_applyBC_shared(ivar,datatype)

#include "Solver.h"

        use MPI_data
        use physicaldata

        implicit none

        integer, intent(in) :: ivar,datatype

        integer :: status(MPI_STATUS_SIZE)

        if(datatype == CENTER) then

                !_______________________MPI BC for High X______________________________!
                if(x_id < x_procs - 1) then

                if(shared_part(myid+1+1) /= MPI_UNDEFINED) then

                        solnData(ivar,Nxb+2,:) = eastCENTER(ivar,2,:)
                else
  
                        call MPI_SENDRECV(solnData(ivar,Nxb+1,:), Nyb+2, MPI_REAL, x_id+1, 2,&
                                          solnData(ivar,Nxb+2,:), Nyb+2, MPI_REAL, x_id+1, 1, x_comm, status, ierr)

                end if
                end if

                !_______________________MPI BC for Low X______________________________!
                if(x_id > 0) then

                if(shared_part(myid+1-1) /= MPI_UNDEFINED) then

                        solnData(ivar,1,:) = westCENTER(ivar,Nxb+1,:)
                else
                        call MPI_SENDRECV(solnData(ivar,2,:), Nyb+2, MPI_REAL, x_id-1, 1,&
                                          solnData(ivar,1,:), Nyb+2, MPI_REAL, x_id-1, 2, x_comm, status, ierr)

                end if
                end if
              
                !_______________________MPI BC for High Y______________________________!
                if(y_id < y_procs - 1) then

                if(shared_part(myid+1+x_procs) /= MPI_UNDEFINED) then

                        solnData(ivar,:,Nyb+2) = northCENTER(ivar,:,2)
                else
                        call MPI_SENDRECV(solnData(ivar,:,Nyb+1), Nxb+2, MPI_REAL, y_id+1, 4,& 
                                          solnData(ivar,:,Nyb+2), Nxb+2, MPI_REAL, y_id+1, 3, y_comm, status, ierr)


                end if
                end if

                !_______________________MPI BC for Low Y______________________________!
                if(y_id > 0) then

                if(shared_part(myid+1-x_procs) /= MPI_UNDEFINED) then

                        solnData(ivar,:,1) = southCENTER(ivar,:,Nyb+1)
                else
                        call MPI_SENDRECV(solnData(ivar,:,2), Nxb+2, MPI_REAL, y_id-1, 3,&
                                          solnData(ivar,:,1), Nxb+2, MPI_REAL, y_id-1, 4, y_comm, status, ierr)

                end if
                end if
       
                call MPI_BARRIER(shared_comm,ierr)
                return

        else if(datatype == FACEX) then

                !_______________________MPI BC for High X______________________________!
                if(x_id < x_procs - 1) then

                if(shared_part(myid+1+1) /= MPI_UNDEFINED) then

                        facexData(ivar,Nxb+2,:) = eastFACEX(ivar,2,:)
                else
  
                        call MPI_SENDRECV(facexData(ivar,Nxb+1,:), Nyb+2, MPI_REAL, x_id+1, 2,&
                                          facexData(ivar,Nxb+2,:), Nyb+2, MPI_REAL, x_id+1, 1, x_comm, status, ierr)

                end if
                end if

                !_______________________MPI BC for Low X______________________________!
                if(x_id > 0) then

                if(shared_part(myid+1-1) /= MPI_UNDEFINED) then

                        facexData(ivar,1,:) = westFACEX(ivar,Nxb+1,:)
                else
                        call MPI_SENDRECV(facexData(ivar,2,:), Nyb+2, MPI_REAL, x_id-1, 1,&
                                          facexData(ivar,1,:), Nyb+2, MPI_REAL, x_id-1, 2, x_comm, status, ierr)

                end if
                end if
              
                !_______________________MPI BC for High Y______________________________!
                if(y_id < y_procs - 1) then

                if(shared_part(myid+1+x_procs) /= MPI_UNDEFINED) then

                        facexData(ivar,:,Nyb+2) = northFACEX(ivar,:,2)
                else
                        call MPI_SENDRECV(facexData(ivar,:,Nyb+1), Nxb+2, MPI_REAL, y_id+1, 4,& 
                                          facexData(ivar,:,Nyb+2), Nxb+2, MPI_REAL, y_id+1, 3, y_comm, status, ierr)


                end if
                end if

                !_______________________MPI BC for Low Y______________________________!
                if(y_id > 0) then

                if(shared_part(myid+1-x_procs) /= MPI_UNDEFINED) then

                        facexData(ivar,:,1) = southFACEX(ivar,:,Nyb+1)
                else
                        call MPI_SENDRECV(facexData(ivar,:,2), Nxb+2, MPI_REAL, y_id-1, 3,&
                                          facexData(ivar,:,1), Nxb+2, MPI_REAL, y_id-1, 4, y_comm, status, ierr)

                end if
                end if
       
                call MPI_BARRIER(shared_comm,ierr)
                return


        else if(datatype == FACEY) then

                !_______________________MPI BC for High X______________________________!
                if(x_id < x_procs - 1) then

                if(shared_part(myid+1+1) /= MPI_UNDEFINED) then

                        faceyData(ivar,Nxb+2,:) = eastFACEY(ivar,2,:)
                else
  
                        call MPI_SENDRECV(faceyData(ivar,Nxb+1,:), Nyb+2, MPI_REAL, x_id+1, 2,&
                                          faceyData(ivar,Nxb+2,:), Nyb+2, MPI_REAL, x_id+1, 1, x_comm, status, ierr)

                end if
                end if

                !_______________________MPI BC for Low X______________________________!
                if(x_id > 0) then

                if(shared_part(myid+1-1) /= MPI_UNDEFINED) then

                        faceyData(ivar,1,:) = westFACEY(ivar,Nxb+1,:)
                else
                        call MPI_SENDRECV(faceyData(ivar,2,:), Nyb+2, MPI_REAL, x_id-1, 1,&
                                          faceyData(ivar,1,:), Nyb+2, MPI_REAL, x_id-1, 2, x_comm, status, ierr)

                end if
                end if
              
                !_______________________MPI BC for High Y______________________________!
                if(y_id < y_procs - 1) then

                if(shared_part(myid+1+x_procs) /= MPI_UNDEFINED) then

                        faceyData(ivar,:,Nyb+2) = northFACEY(ivar,:,2)
                else
                        call MPI_SENDRECV(faceyData(ivar,:,Nyb+1), Nxb+2, MPI_REAL, y_id+1, 4,& 
                                          faceyData(ivar,:,Nyb+2), Nxb+2, MPI_REAL, y_id+1, 3, y_comm, status, ierr)


                end if
                end if

                !_______________________MPI BC for Low Y______________________________!
                if(y_id > 0) then

                if(shared_part(myid+1-x_procs) /= MPI_UNDEFINED) then

                        faceyData(ivar,:,1) = southFACEY(ivar,:,Nyb+1)
                else
                        call MPI_SENDRECV(faceyData(ivar,:,2), Nxb+2, MPI_REAL, y_id-1, 3,&
                                          faceyData(ivar,:,1), Nxb+2, MPI_REAL, y_id-1, 4, y_comm, status, ierr)

                end if
                end if
       
                call MPI_BARRIER(shared_comm,ierr)
                return

        end if

end subroutine MPI_applyBC_shared
