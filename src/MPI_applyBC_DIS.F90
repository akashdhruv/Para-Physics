subroutine MPI_applyBC_DIS(local)

#include "Solver.h"

        use MPI_data
                               
        implicit none

        real, intent(inout), dimension(:,:,:) :: local

        integer :: status(MPI_STATUS_SIZE), blk, req_count

!===============================SWEEP IN X DIRECTION===================================================
        req_count = 0

        do blk = 1,blockCount

        !_______________________MPI BC for High X______________________________!
           if(xLC(blk) < nblockx - 1) then
             if(blockLC(blk+blockOffset) == blockLC(blk+blockOffset+1)) then

                local(Nxb+2,:,blockID(blk+blockOffset)) = &
                local(2,:,blockID(blk+blockOffset+1))

             else
                req_count = req_count + 1

                call MPI_ISEND(local(Nxb+1,:,blockID(blk+blockOffset)), Nyb+2, MPI_REAL, &
                               blockLC(blk+blockOffset+1), blk+blockOffset, solver_comm, reqs(req_count), ierr)

             end if
           end if

        !_______________________MPI BC for Low X________________________________!
           if(xLC(blk) > 0) then
              if(blockLC(blk+blockOffset) == blockLC(blk+blockOffset-1)) then

                 local(1,:,blockID(blk+blockOffset)) = &
                 local(Nxb+1,:,blockID(blk+blockOffset-1))

              else
                 req_count = req_count + 1

                 call MPI_ISEND(local(2,:,blockID(blk+blockOffset)), Nyb+2, MPI_REAL, &
                                blockLC(blk+blockOffset-1), blk+blockOffset, solver_comm, reqs(req_count), ierr)

              end if
            end if

        end do

        do blk = 1,blockCount

        !_______________________MPI BC for High X______________________________!
           if((xLC(blk) < nblockx - 1) .and. (blockLC(blk+blockOffset) /= blockLC(blk+blockOffset+1))) &

           call MPI_RECV(local(Nxb+2,:,blockID(blk+blockOffset)), Nyb+2, MPI_REAL, &
                         blockLC(blk+blockOffset+1), blk+blockOffset+1, solver_comm, status, ierr)
 
        !_______________________MPI BC for Low X________________________________!
           if((xLC(blk) > 0) .and. (blockLC(blk+blockOffset) /= blockLC(blk+blockOffset-1))) &

           call MPI_RECV(local(1,:,blockID(blk+blockOffset)), Nyb+2, MPI_REAL, &
                         blockLC(blk+blockOffset-1), blk+blockOffset-1, solver_comm, status, ierr)
 
        end do


!===============================SWEEP IN Y DIRECTION===================================================
        req_count = 0

        do blk = 1,blockCount            

         !_______________________MPI BC for High Y______________________________!
            if(yLC(blk) < nblocky - 1) then
               if(blockLC(blk+blockOffset) == blockLC(blk+blockOffset+nblockx)) then

                  local(:,Nyb+2,blockID(blk+blockOffset)) = &
                  local(:,2,blockID(blk+blockOffset+nblockx))

               else
                  req_count = req_count + 1

                  call MPI_ISEND(local(:,Nyb+1,blockID(blk+blockOffset)), Nxb+2, MPI_REAL, &
                                 blockLC(blk+blockOffset+nblockx), blk+blockOffset, solver_comm, reqs(req_count), ierr)

               end if
             end if

         !_______________________MPI BC for Low Y________________________________!
             if(yLC(blk) > 0) then
                if(blockLC(blk+blockOffset) == blockLC(blk+blockOffset-nblockx)) then

                   local(:,1,blockID(blk+blockOffset)) = &
                   local(:,Nyb+1,blockID(blk+blockOffset-nblockx))

                else
                   req_count = req_count + 1
                                         
                   call MPI_ISEND(local(:,2,blockID(blk+blockOffset)), Nxb+2, MPI_REAL, &
                                  blockLC(blk+blockOffset-nblockx), blk+blockOffset, solver_comm, reqs(req_count), ierr)

                end if
              end if

        end do

        do blk = 1,blockCount            

         !_______________________MPI BC for High Y______________________________!
            if((yLC(blk) < nblocky - 1) .and. (blockLC(blk+blockOffset) /= blockLC(blk+blockOffset+nblockx))) &

            call MPI_RECV(local(:,Nyb+2,blockID(blk+blockOffset)), Nxb+2, MPI_REAL, &
                          blockLC(blk+blockOffset+nblockx), blk+blockOffset+nblockx, solver_comm, status, ierr)

         !_______________________MPI BC for Low Y________________________________!
             if((yLC(blk) > 0) .and. (blockLC(blk+blockOffset) /= blockLC(blk+blockOffset-nblockx))) &

             call MPI_RECV(local(:,1,blockID(blk+blockOffset)), Nxb+2, MPI_REAL, &
                           blockLC(blk+blockOffset-nblockx), blk+blockOffset-nblockx, solver_comm, status, ierr)

        end do

end subroutine MPI_applyBC_DIS
