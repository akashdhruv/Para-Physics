subroutine MPI_applyBC_SHM(local,sharedDATA,ivar,totvar)

#include "Solver.h"

        use MPI_data

        implicit none

        real, intent(inout), dimension(:,:,:) :: local
        real, intent(inout), dimension(:,:,:,:) :: sharedData
        integer, intent(in) :: ivar,totvar

        integer :: status(MPI_STATUS_SIZE), blk

        do blk = 1,blockCount

           !_______________________MPI BC for High X______________________________!
           if(xLC(blk) < nblockx - 1) then

              if(blockLC(blk+blockOffset) == blockLC(blk+blockOffset+1)) then

                 local(Nxb+2,:,blockID(blk+blockOffset)) = &
                 local(2,:,blockID(blk+blockOffset+1))

              else

                 if(shared_part(blockLC(blk+blockOffset+1)+1) /= MPI_UNDEFINED) then

                    local(Nxb+2,:,blockID(blk+blockOffset)) = &
                    sharedDATA(2,:,blockID(blk+blockOffset+1),ivar+shared_part(blockLC(blk+blockOffset+1)+1)*totvar)

                 else
  
                    call MPI_SENDRECV(local(Nxb+1,:,blockID(blk+blockOffset)), Nyb+2, MPI_REAL, &
                                      blockLC(blk+blockOffset+1), blk+blockOffset,&
                                      local(Nxb+2,:,blockID(blk+blockOffset)), Nyb+2, MPI_REAL, &
                                      blockLC(blk+blockOffset+1), blk+blockOffset+1, solver_comm, status, ierr)

                 end if
              end if
            end if

            !_______________________MPI BC for Low X______________________________!
            if(xLC(blk) > 0) then

               if(blockLC(blk+blockOffset) == blockLC(blk+blockOffset-1)) then

                  local(1,:,blockID(blk+blockOffset)) = &
                  local(Nxb+1,:,blockID(blk+blockOffset-1))

               else

                   if(shared_part(blockLC(blk+blockOffset-1)+1) /= MPI_UNDEFINED) then

                    local(1,:,blockID(blk+blockOffset)) = &
                    sharedDATA(Nxb+1,:,blockID(blk+blockOffset-1),ivar+shared_part(blockLC(blk+blockOffset-1)+1)*totvar)

                   else
  
                    call MPI_SENDRECV(local(2,:,blockID(blk+blockOffset)), Nyb+2, MPI_REAL, &
                                      blockLC(blk+blockOffset-1), blk+blockOffset,&
                                      local(1,:,blockID(blk+blockOffset)), Nyb+2, MPI_REAL, &
                                      blockLC(blk+blockOffset-1), blk+blockOffset-1, solver_comm, status, ierr)

                   end if
               end if
             end if
            
             !_______________________MPI BC for High Y______________________________!
             if(yLC(blk) < nblocky - 1) then

                if(blockLC(blk+blockOffset) == blockLC(blk+blockOffset+nblockx)) then

                   local(:,Nyb+2,blockID(blk+blockOffset)) = &
                   local(:,2,blockID(blk+blockOffset+nblockx))

                else

                   if(shared_part(blockLC(blk+blockOffset+nblockx)+1) /= MPI_UNDEFINED) then

                    local(:,Nyb+2,blockID(blk+blockOffset)) = &
                    sharedDATA(:,2,blockID(blk+blockOffset+nblockx),ivar+shared_part(blockLC(blk+blockOffset+nblockx)+1)*totvar)

                   else
  
                    call MPI_SENDRECV(local(:,Nyb+1,blockID(blk+blockOffset)), Nxb+2, MPI_REAL, &
                                      blockLC(blk+blockOffset+nblockx), blk+blockOffset,&
                                      local(:,Nyb+2,blockID(blk+blockOffset)), Nxb+2, MPI_REAL, &
                                      blockLC(blk+blockOffset+nblockx), blk+blockOffset+nblockx, solver_comm, status, ierr)

                   end if
                end if
             end if

             !_______________________MPI BC for Low Y______________________________!
             if(yLC(blk) > 0) then

                if(blockLC(blk+blockOffset) == blockLC(blk+blockOffset-nblockx)) then

                   local(:,1,blockID(blk+blockOffset)) = &
                   local(:,Nyb+1,blockID(blk+blockOffset-nblockx))

                else

                   if(shared_part(blockLC(blk+blockOffset-nblockx)+1) /= MPI_UNDEFINED) then

                    local(:,1,blockID(blk+blockOffset)) = &
                    sharedDATA(:,Nyb+1,blockID(blk+blockOffset-nblockx),ivar+shared_part(blockLC(blk+blockOffset-nblockx)+1)*totvar)

                   else
  
                    call MPI_SENDRECV(local(:,2,blockID(blk+blockOffset)), Nxb+2, MPI_REAL, &
                                      blockLC(blk+blockOffset-nblockx), blk+blockOffset,&
                                      local(:,1,blockID(blk+blockOffset)), Nxb+2, MPI_REAL, &
                                      blockLC(blk+blockOffset-nblockx), blk+blockOffset-nblockx, solver_comm, status, ierr)

                   end if
                end if
             end if

           end do

end subroutine MPI_applyBC_SHM
