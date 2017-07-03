subroutine MPI_applyBC_shared(ivar,datatype)

#include "Solver.h"

        use MPI_data
        use physicaldata, only: localCENTER,localFACEX,localFACEY,&
                                sharedCENTER,sharedFACEX,sharedFACEY

        implicit none

        integer, intent(in) :: ivar,datatype

        integer :: status(MPI_STATUS_SIZE), blk

        if(datatype == CENTER) then

            do blk = 1,blockCount

                !_______________________MPI BC for High X______________________________!
                if(xLC(blk) < nblockx - 1) then

                  if(blockLC(blk+blockOffset) == blockLC(blk+blockOffset+1)) then

                    localCENTER(Nxb+2,:,blockID(blk+blockOffset),ivar) = &
                    localCENTER(2,:,blockID(blk+blockOffset+1),ivar)

                  else

                   if(shared_part(blockLC(blk+blockOffset+1)+1) /= MPI_UNDEFINED) then

                    localCENTER(Nxb+2,:,blockID(blk+blockOffset),ivar) = &
                    sharedCENTER(2,:,blockID(blk+blockOffset+1),ivar+shared_part(blockLC(blk+blockOffset+1)+1)*CENT_VAR)

                   else
  
                    call MPI_SENDRECV(localCENTER(Nxb+1,:,blockID(blk+blockOffset),ivar), Nyb+2, MPI_REAL, &
                                      blockLC(blk+blockOffset+1), blk+blockOffset,&
                                      localCENTER(Nxb+2,:,blockID(blk+blockOffset),ivar), Nyb+2, MPI_REAL, &
                                      blockLC(blk+blockOffset+1), blk+blockOffset+1, solver_comm, status, ierr)

                   end if
                  end if
                end if

                !_______________________MPI BC for Low X______________________________!
                if(xLC(blk) > 0) then

                   if(blockLC(blk+blockOffset) == blockLC(blk+blockOffset-1)) then

                     localCENTER(1,:,blockID(blk+blockOffset),ivar) = &
                     localCENTER(Nxb+1,:,blockID(blk+blockOffset-1),ivar)

                  else

                   if(shared_part(blockLC(blk+blockOffset-1)+1) /= MPI_UNDEFINED) then

                    localCENTER(1,:,blockID(blk+blockOffset),ivar) = &
                    sharedCENTER(Nxb+1,:,blockID(blk+blockOffset-1),ivar+shared_part(blockLC(blk+blockOffset-1)+1)*CENT_VAR)

                   else
  
                    call MPI_SENDRECV(localCENTER(2,:,blockID(blk+blockOffset),ivar), Nyb+2, MPI_REAL, &
                                      blockLC(blk+blockOffset-1), blk+blockOffset,&
                                      localCENTER(1,:,blockID(blk+blockOffset),ivar), Nyb+2, MPI_REAL, &
                                      blockLC(blk+blockOffset-1), blk+blockOffset-1, solver_comm, status, ierr)

                   end if
                  end if
                end if
            
                !_______________________MPI BC for High Y______________________________!
                if(yLC(blk) < nblocky - 1) then

                  if(blockLC(blk+blockOffset) == blockLC(blk+blockOffset+nblockx)) then

                    localCENTER(:,Nyb+2,blockID(blk+blockOffset),ivar) = &
                    localCENTER(:,2,blockID(blk+blockOffset+nblockx),ivar)

                  else

                   if(shared_part(blockLC(blk+blockOffset+nblockx)+1) /= MPI_UNDEFINED) then

                    localCENTER(:,Nyb+2,blockID(blk+blockOffset),ivar) = &
                    sharedCENTER(:,2,blockID(blk+blockOffset+nblockx),ivar+shared_part(blockLC(blk+blockOffset+nblockx)+1)*CENT_VAR)

                   else
  
                    call MPI_SENDRECV(localCENTER(:,Nyb+1,blockID(blk+blockOffset),ivar), Nxb+2, MPI_REAL, &
                                      blockLC(blk+blockOffset+nblockx), blk+blockOffset,&
                                      localCENTER(:,Nyb+2,blockID(blk+blockOffset),ivar), Nxb+2, MPI_REAL, &
                                      blockLC(blk+blockOffset+nblockx), blk+blockOffset+nblockx, solver_comm, status, ierr)

                   end if
                  end if
                end if

                !_______________________MPI BC for Low Y______________________________!
                if(yLC(blk) > 0) then

                  if(blockLC(blk+blockOffset) == blockLC(blk+blockOffset-nblockx)) then

                    localCENTER(:,1,blockID(blk+blockOffset),ivar) = &
                    localCENTER(:,Nyb+1,blockID(blk+blockOffset-nblockx),ivar)

                  else

                   if(shared_part(blockLC(blk+blockOffset-nblockx)+1) /= MPI_UNDEFINED) then

                    localCENTER(:,1,blockID(blk+blockOffset),ivar) = &
                    sharedCENTER(:,Nyb+1,blockID(blk+blockOffset-nblockx),ivar+shared_part(blockLC(blk+blockOffset-nblockx)+1)*CENT_VAR)

                   else
  
                    call MPI_SENDRECV(localCENTER(:,2,blockID(blk+blockOffset),ivar), Nxb+2, MPI_REAL, &
                                      blockLC(blk+blockOffset-nblockx), blk+blockOffset,&
                                      localCENTER(:,1,blockID(blk+blockOffset),ivar), Nxb+2, MPI_REAL, &
                                      blockLC(blk+blockOffset-nblockx), blk+blockOffset-nblockx, solver_comm, status, ierr)

                   end if
                  end if
                end if
           end do
           return

        else if(datatype == FACEX) then

            do blk = 1,blockCount

                !_______________________MPI BC for High X______________________________!
                if(xLC(blk) < nblockx - 1) then

                  if(blockLC(blk+blockOffset) == blockLC(blk+blockOffset+1)) then

                    localFACEX(Nxb+2,:,blockID(blk+blockOffset),ivar) = &
                    localFACEX(2,:,blockID(blk+blockOffset+1),ivar)

                  else

                   if(shared_part(blockLC(blk+blockOffset+1)+1) /= MPI_UNDEFINED) then

                    localFACEX(Nxb+2,:,blockID(blk+blockOffset),ivar) = &
                    sharedFACEX(2,:,blockID(blk+blockOffset+1),ivar+shared_part(blockLC(blk+blockOffset+1)+1)*FACE_VAR)

                   else
  
                    call MPI_SENDRECV(localFACEX(Nxb+1,:,blockID(blk+blockOffset),ivar), Nyb+2, MPI_REAL, &
                                      blockLC(blk+blockOffset+1), blk+blockOffset,&
                                      localFACEX(Nxb+2,:,blockID(blk+blockOffset),ivar), Nyb+2, MPI_REAL, &
                                      blockLC(blk+blockOffset+1), blk+blockOffset+1, solver_comm, status, ierr)

                   end if
                  end if
                end if

                !_______________________MPI BC for Low X______________________________!
                if(xLC(blk) > 0) then

                   if(blockLC(blk+blockOffset) == blockLC(blk+blockOffset-1)) then

                     localFACEX(1,:,blockID(blk+blockOffset),ivar) = &
                     localFACEX(Nxb+1,:,blockID(blk+blockOffset-1),ivar)

                  else

                   if(shared_part(blockLC(blk+blockOffset-1)+1) /= MPI_UNDEFINED) then

                     localFACEX(1,:,blockID(blk+blockOffset),ivar) = &
                     sharedFACEX(Nxb+1,:,blockID(blk+blockOffset-1),ivar+shared_part(blockLC(blk+blockOffset-1)+1)*FACE_VAR)

                   else
  
                    call MPI_SENDRECV(localFACEX(2,:,blockID(blk+blockOffset),ivar), Nyb+2, MPI_REAL, &
                                      blockLC(blk+blockOffset-1), blk+blockOffset,&
                                      localFACEX(1,:,blockID(blk+blockOffset),ivar), Nyb+2, MPI_REAL, &
                                      blockLC(blk+blockOffset-1), blk+blockOffset-1, solver_comm, status, ierr)

                   end if
                  end if
                end if
            
                !_______________________MPI BC for High Y______________________________!
                if(yLC(blk) < nblocky - 1) then

                   if(blockLC(blk+blockOffset) == blockLC(blk+blockOffset+nblockx)) then

                     localFACEX(:,Nyb+2,blockID(blk+blockOffset),ivar) = &
                     localFACEX(:,2,blockID(blk+blockOffset+nblockx),ivar)

                  else

                   if(shared_part(blockLC(blk+blockOffset+nblockx)+1) /= MPI_UNDEFINED) then

                     localFACEX(:,Nyb+2,blockID(blk+blockOffset),ivar) = &
                     sharedFACEX(:,2,blockID(blk+blockOffset+nblockx),ivar+shared_part(blockLC(blk+blockOffset+nblockx)+1)*FACE_VAR)

                   else
  
                    call MPI_SENDRECV(localFACEX(:,Nyb+1,blockID(blk+blockOffset),ivar), Nxb+2, MPI_REAL, &
                                      blockLC(blk+blockOffset+nblockx), blk+blockOffset,&
                                      localFACEX(:,Nyb+2,blockID(blk+blockOffset),ivar), Nxb+2, MPI_REAL, &
                                      blockLC(blk+blockOffset+nblockx), blk+blockOffset+nblockx, solver_comm, status, ierr)

                   end if
                  end if
                end if

                !_______________________MPI BC for Low Y______________________________!
                if(yLC(blk) > 0) then

                   if(blockLC(blk+blockOffset) == blockLC(blk+blockOffset-nblockx)) then

                     localFACEX(:,1,blockID(blk+blockOffset),ivar) = &
                     localFACEX(:,Nyb+1,blockID(blk+blockOffset-nblockx),ivar)

                  else

                   if(shared_part(blockLC(blk+blockOffset-nblockx)+1) /= MPI_UNDEFINED) then

                     localFACEX(:,1,blockID(blk+blockOffset),ivar) = &
                     sharedFACEX(:,Nyb+1,blockID(blk+blockOffset-nblockx),ivar+shared_part(blockLC(blk+blockOffset-nblockx)+1)*FACE_VAR)

                   else
  
                    call MPI_SENDRECV(localFACEX(:,2,blockID(blk+blockOffset),ivar), Nxb+2, MPI_REAL, &
                                      blockLC(blk+blockOffset-nblockx), blk+blockOffset,&
                                      localFACEX(:,1,blockID(blk+blockOffset),ivar), Nxb+2, MPI_REAL, &
                                      blockLC(blk+blockOffset-nblockx), blk+blockOffset-nblockx, solver_comm, status, ierr)

                   end if
                  end if
                end if
           end do
           return


        else if(datatype == FACEY) then

            do blk = 1,blockCount

                !_______________________MPI BC for High X______________________________!
                if(xLC(blk) < nblockx - 1) then

                   if(blockLC(blk+blockOffset) == blockLC(blk+blockOffset+1)) then

                     localFACEY(Nxb+2,:,blockID(blk+blockOffset),ivar) = &
                     localFACEY(2,:,blockID(blk+blockOffset+1),ivar)

                  else

                   if(shared_part(blockLC(blk+blockOffset+1)+1) /= MPI_UNDEFINED) then

                     localFACEY(Nxb+2,:,blockID(blk+blockOffset),ivar) = &
                     sharedFACEY(2,:,blockID(blk+blockOffset+1),ivar+shared_part(blockLC(blk+blockOffset+1)+1)*FACE_VAR)

                   else
  
                    call MPI_SENDRECV(localFACEY(Nxb+1,:,blockID(blk+blockOffset),ivar), Nyb+2, MPI_REAL, &
                                      blockLC(blk+blockOffset+1), blk+blockOffset,&
                                      localFACEY(Nxb+2,:,blockID(blk+blockOffset),ivar), Nyb+2, MPI_REAL, &
                                      blockLC(blk+blockOffset+1), blk+blockOffset+1, solver_comm, status, ierr)

                   end if
                  end if
                end if

                !_______________________MPI BC for Low X______________________________!
                if(xLC(blk) > 0) then

                   if(blockLC(blk+blockOffset) == blockLC(blk+blockOffset-1)) then

                     localFACEY(1,:,blockID(blk+blockOffset),ivar) = &
                     localFACEY(Nxb+1,:,blockID(blk+blockOffset-1),ivar)

                  else

                   if(shared_part(blockLC(blk+blockOffset-1)+1) /= MPI_UNDEFINED) then

                     localFACEY(1,:,blockID(blk+blockOffset),ivar) = &
                     sharedFACEY(Nxb+1,:,blockID(blk+blockOffset-1),ivar+shared_part(blockLC(blk+blockOffset-1)+1)*FACE_VAR)

                   else
  
                    call MPI_SENDRECV(localFACEY(2,:,blockID(blk+blockOffset),ivar), Nyb+2, MPI_REAL, &
                                      blockLC(blk+blockOffset-1), blk+blockOffset,&
                                      localFACEY(1,:,blockID(blk+blockOffset),ivar), Nyb+2, MPI_REAL, &
                                      blockLC(blk+blockOffset-1), blk+blockOffset-1, solver_comm, status, ierr)

                   end if
                  end if
                end if
            
                !_______________________MPI BC for High Y______________________________!
                if(yLC(blk) < nblocky - 1) then

                  if(blockLC(blk+blockOffset) == blockLC(blk+blockOffset+nblockx)) then

                    localFACEY(:,Nyb+2,blockID(blk+blockOffset),ivar) = &
                    localFACEY(:,2,blockID(blk+blockOffset+nblockx),ivar)

                  else

                   if(shared_part(blockLC(blk+blockOffset+nblockx)+1) /= MPI_UNDEFINED) then

                     localFACEY(:,Nyb+2,blockID(blk+blockOffset),ivar) = &
                     sharedFACEY(:,2,blockID(blk+blockOffset+nblockx),ivar+shared_part(blockLC(blk+blockOffset+nblockx)+1)*FACE_VAR)

                   else
  
                    call MPI_SENDRECV(localFACEY(:,Nyb+1,blockID(blk+blockOffset),ivar), Nxb+2, MPI_REAL, &
                                      blockLC(blk+blockOffset+nblockx), blk+blockOffset,&
                                      localFACEY(:,Nyb+2,blockID(blk+blockOffset),ivar), Nxb+2, MPI_REAL, &
                                      blockLC(blk+blockOffset+nblockx), blk+blockOffset+nblockx, solver_comm, status, ierr)

                   end if
                  end if
                end if

                !_______________________MPI BC for Low Y______________________________!
                if(yLC(blk) > 0) then

                  if(blockLC(blk+blockOffset) == blockLC(blk+blockOffset-nblockx)) then

                    localFACEY(:,1,blockID(blk+blockOffset),ivar) = &
                    localFACEY(:,Nyb+1,blockID(blk+blockOffset-nblockx),ivar)

                  else

                   if(shared_part(blockLC(blk+blockOffset-nblockx)+1) /= MPI_UNDEFINED) then

                     localFACEY(:,1,blockID(blk+blockOffset),ivar) = &
                     sharedFACEY(:,Nyb+1,blockID(blk+blockOffset-nblockx),ivar+shared_part(blockLC(blk+blockOffset-nblockx)+1)*FACE_VAR)

                   else
  
                    call MPI_SENDRECV(localFACEY(:,2,blockID(blk+blockOffset),ivar), Nxb+2, MPI_REAL, &
                                      blockLC(blk+blockOffset-nblockx), blk+blockOffset,&
                                      localFACEY(:,1,blockID(blk+blockOffset),ivar), Nxb+2, MPI_REAL, &
                                      blockLC(blk+blockOffset-nblockx), blk+blockOffset-nblockx, solver_comm, status, ierr)

                     end if
                  end if
                end if
           end do
           return

        end if

end subroutine MPI_applyBC_shared
