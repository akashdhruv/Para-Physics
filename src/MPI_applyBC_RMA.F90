subroutine MPI_applyBC_RMA(local)

#include "Solver.h"

        use MPI_data
        use physicaldata, only : dataTARGET,&
                                 eastORIGIN,westORIGIN,northORIGIN,southORIGIN

        implicit none

        real, intent(inout), dimension(:,:,:) :: local

        integer :: status(MPI_STATUS_SIZE)
        integer(kind=MPI_ADDRESS_KIND) :: target_disp = 0 
        integer :: index1=Nyb+2, index2=Nyb+2+Nyb+2, index3=Nyb+2+Nyb+2+Nxb+2
        integer :: blk
        integer :: bSize

        bSize = Nxb+2+Nxb+2+Nyb+2+Nyb+2

        do blk=1,blockCount
        dataTARGET(1+(blk-1)*bsize:Nyb+2+(blk-1)*bSize)               = local(2,:,blk)
        dataTARGET(index1+1+(blk-1)*bsize:index1+Nyb+2+(blk-1)*bsize) = local(Nxb+1,:,blk)
        dataTARGET(index2+1+(blk-1)*bsize:index2+Nxb+2+(blk-1)*bsize) = local(:,2,blk)
        dataTARGET(index3+1+(blk-1)*bsize:index3+Nxb+2+(blk-1)*bsize) = local(:,Nyb+1,blk)
        end do

#ifdef MPI_RMA_ACTIVE
        call MPI_WIN_FENCE(0,RMA_win,ierr)

        do blk=1,blockCount

        !_______________________MPI BC for High X______________________________!
        if(xLC(blk) < nblockx - 1) then
           if(blockLC(blk+blockOffset) == blockLC(blk+blockOffset+1)) then
              eastORIGIN(:,blk) = local(2,:,blockID(blk+blockOffset+1))
           else
             call MPI_GET(eastORIGIN(1,blk),Nyb+2,MPI_REAL,blockLC(blk+blockOffset+1),&
                  target_disp+(blockID(blk+blockOffset+1)-1)*bSize,Nyb+2,MPI_REAL,RMA_win,ierr)
           end if
        end if

        !_______________________MPI BC for Low X______________________________!
        if(xLC(blk) > 0) then
           if(blockLC(blk+blockOffset) == blockLC(blk+blockOffset-1)) then
              westORIGIN(:,blk) = local(Nxb+1,:,blockID(blk+blockOffset-1))
           else
              call MPI_GET(westORIGIN(1,blk),Nyb+2,MPI_REAL,blockLC(blk+blockOffset-1),&
                   target_disp+index1+(blockID(blk+blockOffset-1)-1)*bSize,Nyb+2,MPI_REAL,RMA_win,ierr)
           end if
        end if         
      
        !_______________________MPI BC for High Y______________________________!
        if(yLC(blk) < nblocky - 1) then
           if(blockLC(blk+blockOffset) == blockLC(blk+blockOffset+nblockx)) then
              northORIGIN(:,blk) = local(:,2,blockID(blk+blockOffset+nblockx))
           else
              call MPI_GET(northORIGIN(1,blk),Nxb+2,MPI_REAL,blockLC(blk+blockOffset+nblockx),&
                   target_disp+index2+(blockID(blk+blockOffset+nblockx)-1)*bSize,Nxb+2,MPI_REAL,RMA_win,ierr)
           
           end if
        end if

        !_______________________MPI BC for Low Y______________________________!
        if(yLC(blk) > 0) then
           if(blockLC(blk+blockOffset) == blockLC(blk+blockOffset-nblockx)) then
              southORIGIN(:,blk) = local(:,Nyb+1,blockID(blk+blockOffset-nblockx))
           else
              call MPI_GET(southORIGIN(1,blk),Nxb+2,MPI_REAL,blockLC(blk+blockOffset-nblockx),&
                   target_disp+index3+(blockID(blk+blockOffset-nblockx)-1)*bSize,Nxb+2,MPI_REAL,RMA_win,ierr)
    
           end if
        end if

        end do

        call MPI_WIN_FENCE(0,RMA_win,ierr)
#endif

#ifdef MPI_RMA_PASSIVE
        call MPI_BARRIER(solver_comm,ierr)

        do blk=1,blockCount

        !_______________________MPI BC for High X______________________________!
        if(xLC(blk) < nblockx - 1) then
           if(blockLC(blk+blockOffset) == blockLC(blk+blockOffset+1)) then
              eastORIGIN(:,blk) = local(2,:,blockID(blk+blockOffset+1))
           else
             call MPI_GET(eastORIGIN(1,blk),Nyb+2,MPI_REAL,blockLC(blk+blockOffset+1),&
                  target_disp+(blockID(blk+blockOffset+1)-1)*bSize,Nyb+2,MPI_REAL,RMA_win,ierr)
           end if
        end if

        !_______________________MPI BC for Low X______________________________!
        if(xLC(blk) > 0) then
           if(blockLC(blk+blockOffset) == blockLC(blk+blockOffset-1)) then
              westORIGIN(:,blk) = local(Nxb+1,:,blockID(blk+blockOffset-1))
           else
              call MPI_GET(westORIGIN(1,blk),Nyb+2,MPI_REAL,blockLC(blk+blockOffset-1),&
                   target_disp+index1+(blockID(blk+blockOffset-1)-1)*bSize,Nyb+2,MPI_REAL,RMA_win,ierr)
           end if
        end if         
      
        !_______________________MPI BC for High Y______________________________!
        if(yLC(blk) < nblocky - 1) then
           if(blockLC(blk+blockOffset) == blockLC(blk+blockOffset+nblockx)) then
              northORIGIN(:,blk) = local(:,2,blockID(blk+blockOffset+nblockx))
           else
              call MPI_GET(northORIGIN(1,blk),Nxb+2,MPI_REAL,blockLC(blk+blockOffset+nblockx),&
                   target_disp+index2+(blockID(blk+blockOffset+nblockx)-1)*bSize,Nxb+2,MPI_REAL,RMA_win,ierr)
           
           end if
        end if

        !_______________________MPI BC for Low Y______________________________!
        if(yLC(blk) > 0) then
           if(blockLC(blk+blockOffset) == blockLC(blk+blockOffset-nblockx)) then
              southORIGIN(:,blk) = local(:,Nyb+1,blockID(blk+blockOffset-nblockx))
           else
              call MPI_GET(southORIGIN(1,blk),Nxb+2,MPI_REAL,blockLC(blk+blockOffset-nblockx),&
                   target_disp+index3+(blockID(blk+blockOffset-nblockx)-1)*bSize,Nxb+2,MPI_REAL,RMA_win,ierr)
    
           end if
        end if

        end do

        call MPI_BARRIER(solver_comm,ierr)
#endif

        do blk=1,blockCount
        local(Nxb+2,:,blk) = eastORIGIN(:,blk)
        local(1,:,blk)     = westORIGIN(:,blk)
        local(:,Nyb+2,blk) = northORIGIN(:,blk)
        local(:,1,blk)     = southORIGIN(:,blk)
        end do

end subroutine MPI_applyBC_RMA
