subroutine MPI_applyBC_SHM(local,sharedDATA)

#include "Solver.h"

        use MPI_data

        implicit none

        real, intent(inout), dimension(:,:,:) :: local
        real, intent(inout), dimension(:,:,:) :: sharedData

        integer :: status(MPI_STATUS_SIZE), blk, req_count

!===============================SWEEP IN X DIRECTION===================================================
        req_count = 0

        do blk = 1,blockCount

           !_______________________MPI BC for High X______________________________!
           if(xLC(blk) < nblockx - 1) then

              if(blockLC(blk+blockOffset) == blockLC(blk+blockOffset+1)) then           ! Check if neighbouring block on same process

                 local(Nxb+2,:,blockID(blk+blockOffset)) = &
                 local(2,:,blockID(blk+blockOffset+1))


              else if(shared_part(blockLC(blk+blockOffset+1)+1) /= MPI_UNDEFINED) then  ! If not on same process then check if on the same node

                 local(Nxb+2,:,blockID(blk+blockOffset)) = &
                 sharedDATA(2,:,blockID(blk+blockOffset+1)+shared_part(blockLC(blk+blockOffset+1)+1)*blockCount)

              else                                                                      ! If on different node perform MPI exchange

                req_count = req_count + 1

                call MPI_ISEND(local(Nxb+1,:,blockID(blk+blockOffset)), Nyb+2, MPI_REAL, &
                               blockLC(blk+blockOffset+1), blk+blockOffset, solver_comm, reqs(req_count), ierr)

              end if
            end if

            !_______________________MPI BC for Low X______________________________!
            if(xLC(blk) > 0) then

               if(blockLC(blk+blockOffset) == blockLC(blk+blockOffset-1)) then          ! Check if neighbouring block on same process

                  local(1,:,blockID(blk+blockOffset)) = &
                  local(Nxb+1,:,blockID(blk+blockOffset-1))

               else if(shared_part(blockLC(blk+blockOffset-1)+1) /= MPI_UNDEFINED) then ! If not on same process then check if on the same node

                  local(1,:,blockID(blk+blockOffset)) = &
                  sharedDATA(Nxb+1,:,blockID(blk+blockOffset-1)+shared_part(blockLC(blk+blockOffset-1)+1)*blockCount)

               else                                                                     ! If on different node perform MPI exchange

                 req_count = req_count + 1

                 call MPI_ISEND(local(2,:,blockID(blk+blockOffset)), Nyb+2, MPI_REAL, &
                                blockLC(blk+blockOffset-1), blk+blockOffset, solver_comm, reqs(req_count), ierr)

               end if
             end if

        end do

        ! Receive information to complete MPI exchange
        do blk = 1,blockCount

        !_______________________MPI BC for High X______________________________!
           if((xLC(blk) < nblockx - 1) .and. &
              (blockLC(blk+blockOffset) /= blockLC(blk+blockOffset+1)) .and. &
              (shared_part(blockLC(blk+blockOffset+1)+1) == MPI_UNDEFINED)) &

           call MPI_RECV(local(Nxb+2,:,blockID(blk+blockOffset)), Nyb+2, MPI_REAL, &
                         blockLC(blk+blockOffset+1), blk+blockOffset+1, solver_comm, status, ierr)

        !_______________________MPI BC for Low X________________________________!
           if((xLC(blk) > 0) .and. &
              (blockLC(blk+blockOffset) /= blockLC(blk+blockOffset-1)) .and. &
              (shared_part(blockLC(blk+blockOffset-1)+1) == MPI_UNDEFINED)) &

           call MPI_RECV(local(1,:,blockID(blk+blockOffset)), Nyb+2, MPI_REAL, &
                         blockLC(blk+blockOffset-1), blk+blockOffset-1, solver_comm, status, ierr)

        end do

        call MPI_WAITALL(req_count,reqs(1:req_count),req_stat(:,1:req_count),ierr)

!===============================SWEEP IN Y DIRECTION===================================================
          req_count = 0

          do blk = 1,blockCount

             !_______________________MPI BC for High Y______________________________!
             if(yLC(blk) < nblocky - 1) then

                if(blockLC(blk+blockOffset) == blockLC(blk+blockOffset+nblockx)) then           ! Check if neighbouring block on same process

                   local(:,Nyb+2,blockID(blk+blockOffset)) = &
                   local(:,2,blockID(blk+blockOffset+nblockx))

                else if(shared_part(blockLC(blk+blockOffset+nblockx)+1) /= MPI_UNDEFINED) then  ! If not on same process then check if on the same node

                   local(:,Nyb+2,blockID(blk+blockOffset)) = &
                   sharedDATA(:,2,blockID(blk+blockOffset+nblockx)+shared_part(blockLC(blk+blockOffset+nblockx)+1)*blockCount)

                else                                                                            ! If on different node perform MPI exchange

                  req_count = req_count + 1

                  call MPI_ISEND(local(:,Nyb+1,blockID(blk+blockOffset)), Nxb+2, MPI_REAL, &
                                 blockLC(blk+blockOffset+nblockx), blk+blockOffset, solver_comm, reqs(req_count), ierr)  

                end if
             end if

             !_______________________MPI BC for Low Y______________________________!
             if(yLC(blk) > 0) then

                if(blockLC(blk+blockOffset) == blockLC(blk+blockOffset-nblockx)) then           ! Check if neighbouring block on same process

                   local(:,1,blockID(blk+blockOffset)) = &
                   local(:,Nyb+1,blockID(blk+blockOffset-nblockx))

                else if(shared_part(blockLC(blk+blockOffset-nblockx)+1) /= MPI_UNDEFINED) then  ! If not on same process then check if on the same node

                   local(:,1,blockID(blk+blockOffset)) = &
                   sharedDATA(:,Nyb+1,blockID(blk+blockOffset-nblockx)+shared_part(blockLC(blk+blockOffset-nblockx)+1)*blockCount)

                else                                                                            ! If on different node perform MPI exchange

                   req_count = req_count + 1

                   call MPI_ISEND(local(:,2,blockID(blk+blockOffset)), Nxb+2, MPI_REAL, &
                                  blockLC(blk+blockOffset-nblockx), blk+blockOffset, solver_comm, reqs(req_count), ierr)
  
                end if
             end if

          end do

         ! Receive information to complete MPI exchange
         do blk = 1,blockCount

         !_______________________MPI BC for High Y______________________________!
            if((yLC(blk) < nblocky - 1) .and. &
               (blockLC(blk+blockOffset) /= blockLC(blk+blockOffset+nblockx)) .and. &
               (shared_part(blockLC(blk+blockOffset+nblockx)+1) == MPI_UNDEFINED)) &

            call MPI_RECV(local(:,Nyb+2,blockID(blk+blockOffset)), Nxb+2, MPI_REAL, &
                          blockLC(blk+blockOffset+nblockx), blk+blockOffset+nblockx, solver_comm, status, ierr)

         !_______________________MPI BC for Low Y________________________________!
             if((yLC(blk) > 0) .and. &
                (blockLC(blk+blockOffset) /= blockLC(blk+blockOffset-nblockx)) .and. &
                (shared_part(blockLC(blk+blockOffset-nblockx)+1) == MPI_UNDEFINED)) &

             call MPI_RECV(local(:,1,blockID(blk+blockOffset)), Nxb+2, MPI_REAL, &
                           blockLC(blk+blockOffset-nblockx), blk+blockOffset-nblockx, solver_comm, status, ierr)

         end do

         call MPI_WAITALL(req_count,reqs(1:req_count),req_stat(:,1:req_count),ierr)

end subroutine MPI_applyBC_SHM
