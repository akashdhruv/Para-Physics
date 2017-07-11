subroutine MPIsolver_init()

#include "Solver.h"

    !$ use omp_lib

    use MPI_data
    use morton_interface, only: morton_sort
    use physicaldata

    implicit none

    integer :: checkSumMPI
    integer :: status = 0
    integer :: i,blk
    real :: A

    !_________Define Global Communication Environment___________!
    solver_comm = MPI_COMM_WORLD

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(solver_comm, myid, ierr)         ! Get global rank
    call MPI_COMM_SIZE(solver_comm, procs, ierr)        ! Get global comm size
    call MPI_COMM_GROUP(solver_comm, world_group, ierr) ! Create global group

    allocate(world_part(procs))
    world_part = (/(I,I=0,procs-1)/) ! Create array of global id's on each process

    !_______________Split Domain Into Blocks For Cache Optimization_________!

    blockCount = ((nblockx*nblocky)/procs)   ! Divide total blocks among all the processes
    checkSumMPI = blockCount*procs           ! checksum for load balancing
   
    blockOffset = myid*blockCount            ! Block offset to identify global position
                                             ! of local blocks with respect to blocks 
                                             ! on other processes. Used during MPI boundary conditions

    if (checkSumMPI /= nblockx*nblocky) then ! Checking if block distribution is uniform for balancing
       
        if (myid == 0) &
        print *,"RUNTIME ERROR: The number of blocks should be greater than/equal to &
                 and exactly divisible by total number of MPI processes."

        call MPI_FINALIZE(ierr)
        call exit(status)

    else if(blockCount > MAX_BLOCKS) then    ! Checking if blockCount per process is less than the threshold
    
        if (myid == 0) &
        print *,"RUNTIME ERROR: The total number of blocks per process exceed the maximum limit. &
                 Increase the number of MPI jobs."

        call MPI_FINALIZE(ierr)
        call exit(status)

    end if

    allocate(blockID(nblockx*nblocky)) ! Array to store local block number of each global block
    allocate(blockLC(nblockx*nblocky)) ! Array to store processes id for each block
   
    do i=1,procs

       blockLC(1+(i-1)*blockCount:blockCount+(i-1)*blockCount) = i-1
       blockID(1+(i-1)*blockCount:blockCount+(i-1)*blockCount) = (/(I,I=1,blockCount)/)

    end do

    !_______________Apply Refinement Using AMR - still in debug_________!

    !call morton_sort(blockCount,myid,procs,blockID,blockLC) ! Morton sort for AMR blocks
 
    !_________________Define Communication Based On Grid________________!

    ! Splitting communication environment in X and Y direction.
    ! Used in previous version of para-physics, not relevant anymore.
    ! Required if want to use subroutine MPI_applyBC_ORIG with blockCount = 1

    call MPI_COMM_SPLIT(solver_comm,myid/nblockx,myid,x_comm,ierr)
    call MPI_COMM_SPLIT(solver_comm,mod(myid,nblockx),myid,y_comm,ierr)

    call MPI_COMM_RANK(x_comm,x_id,ierr)
    call MPI_COMM_SIZE(x_comm,x_procs,ierr)

    call MPI_COMM_RANK(y_comm,y_id,ierr)
    call MPI_COMM_size(y_comm,y_procs,ierr)

    ! Spliting block IDs in X and Y direction. 
    ! Analogous to splitting communication in the previous version of para-physics.
    ! Used in MPI_applyBC_DIS, MPI_applyBC_RMA, MPI_applyBC_SHM and
    ! MPI physical BCs. 

    ! Required for the latest version

    allocate(xLC(blockCount))
    allocate(yLC(blockCount))

    xLC = mod(((/(I,I=0,blockCount-1)/) + blockOffset),nblockx)
    yLC = ((/(I,I=0,blockCount-1)/) + blockOffset)/nblockx

    ! Allocating request and status array for non-blocking MPI communications.
    ! Used in MPI boundary conditions

    allocate(reqs(blockCount*2))
    allocate(req_stat(MPI_STATUS_SIZE,blockCount*2))

#ifdef MPI_DIS
    ! Distributed Memory Model

    ! Allocate physical data for distributed memory MPI communication
    allocate(localCENTER(Nxb+2,Nyb+2,CENT_VAR,blockCount))
    allocate(localFACEX(Nxb+2,Nyb+2,FACE_VAR,blockCount))
    allocate(localFACEY(Nxb+2,Nyb+2,FACE_VAR,blockCount))
#endif

#ifdef MPI_SHM
    ! Shared Memory Model

    !__________Define Shared Communication Environment__________!

    ! Splitting processes into shared groups based on node locality
    call MPI_COMM_SPLIT_TYPE(solver_comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, shared_comm, ierr)

    call MPI_COMM_RANK(shared_comm,shared_id,ierr)     ! Shared process id
    call MPI_COMM_SIZE(shared_comm,shared_procs,ierr)  ! Shared process comm size
    call MPI_COMM_GROUP(shared_comm,shared_group,ierr) ! Create shared group

    ! Translate global rank array to shared rank array
    ! This will be used to identify if the neighbours are on the same node or a different one
    allocate(shared_part(procs))
    call MPI_GROUP_TRANSLATE_RANKS(world_group,procs,world_part,shared_group,shared_part,ierr)

    ! Create info key for shared memory window allocation
    call MPI_INFO_CREATE(mpi_info_key,ierr)
    call MPI_INFO_SET(mpi_info_key,"alloc_shared_noncontig","false",ierr)

   !_________Make on-node processes allocate their chunk of shared Memory____________________!

    ! Size of windows for cell-center and face-center data
    center_size = blockCount*CENT_VAR*(Nxb+2)*(Nyb+2)*sizeof(A)
    facex_size  = blockCount*FACE_VAR*(Nxb+2)*(Nyb+2)*sizeof(A)
    facey_size  = blockCount*FACE_VAR*(Nxb+2)*(Nyb+2)*sizeof(A) 

    ! Memory displacement between each element of array
    disp_unit = sizeof(A)

    ! Allocate shared memory window - each processes will allocate their chunk which is shared between all
    ! the processes in shared communication environment
    call MPI_WIN_ALLOCATE_SHARED(center_size,disp_unit,mpi_info_key,shared_comm,center_ptr,center_win,ierr)
    call MPI_WIN_ALLOCATE_SHARED(facex_size,disp_unit,mpi_info_key,shared_comm,facex_ptr,facex_win,ierr)
    call MPI_WIN_ALLOCATE_SHARED(facey_size,disp_unit,mpi_info_key,shared_comm,facey_ptr,facey_win,ierr)

    !__________________Point to local chunk of the shared data_______________________________!
    call MPI_WIN_SHARED_QUERY(center_win, shared_id, center_size, disp_unit, center_ptr,ierr)
    call MPI_BARRIER(shared_comm,ierr)
    call C_F_POINTER(center_ptr, localCENTER,[Nxb+2,Nyb+2,CENT_VAR,blockCount])

    call MPI_WIN_SHARED_QUERY(facex_win, shared_id, facex_size, disp_unit, facex_ptr,ierr)
    call MPI_BARRIER(shared_comm,ierr)
    call C_F_POINTER(facex_ptr, localFACEX, [Nxb+2,Nyb+2,FACE_VAR,blockCount])

    call MPI_WIN_SHARED_QUERY(facey_win, shared_id, facey_size, disp_unit, facey_ptr,ierr)
    call MPI_BARRIER(shared_comm,ierr)
    call C_F_POINTER(facey_ptr, localFACEY, [Nxb+2,Nyb+2,FACE_VAR,blockCount])

    !____________________________Point to the shared data________________________________!
    call MPI_WIN_SHARED_QUERY(center_win, 0 ,center_size,disp_unit,center_ptr,ierr)
    call C_F_POINTER(center_ptr,sharedCENTER,[Nxb+2,Nyb+2,CENT_VAR,blockCount*shared_procs])
    call MPI_BARRIER(shared_comm,ierr)

    call MPI_WIN_SHARED_QUERY(facex_win, 0,facex_size,disp_unit,facex_ptr,ierr)
    call C_F_POINTER(facex_ptr,sharedFACEX,[Nxb+2,Nyb+2,FACE_VAR,blockCount*shared_procs])
    call MPI_BARRIER(shared_comm,ierr)

    call MPI_WIN_SHARED_QUERY(facey_win, 0,facey_size,disp_unit,facey_ptr,ierr)
    call C_F_POINTER(facey_ptr,sharedFACEY,[Nxb+2,Nyb+2,FACE_VAR,blockCount*shared_procs])
    call MPI_BARRIER(shared_comm,ierr)
#endif

#ifdef MPI_RMA
    ! Remote Memory Access Model

    !__________Define Shared Communication Environment__________!

    ! Splitting processes into shared groups based on node locality
    call MPI_COMM_SPLIT_TYPE(solver_comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, shared_comm, ierr)
    call MPI_COMM_RANK(shared_comm,shared_id,ierr)
    call MPI_COMM_SIZE(shared_comm,shared_procs,ierr)
    call MPI_COMM_GROUP(shared_comm, shared_group,ierr)

    ! Translate global rank array to shared rank array
    ! This will be used to identify if the neighbours are on the same node or a different one
    allocate(shared_part(procs))
    call MPI_GROUP_TRANSLATE_RANKS(world_group,procs,world_part,shared_group,shared_part,ierr)

    ! Create info key for RMA window allocation
    call MPI_INFO_CREATE(mpi_info_key,ierr)
    call MPI_INFO_SET(mpi_info_key,"no_locks","true",ierr)

    ! Allocate physical data
    allocate(localCENTER(Nxb+2,Nyb+2,CENT_VAR,blockCount))
    allocate(localFACEX(Nxb+2,Nyb+2,FACE_VAR,blockCount))
    allocate(localFACEY(Nxb+2,Nyb+2,FACE_VAR,blockCount))

    ! Allocate array for ORIGIN data per process
    allocate(eastORIGIN(Nyb+2,blockCount))
    allocate(westORIGIN(Nyb+2,blockCount))
    allocate(northORIGIN(Nxb+2,blockCount))
    allocate(southORIGIN(Nxb+2,blockCount))

    ! RMA window size
    RMA_size   = (Nyb+2+Nxb+2+Nyb+2+Nxb+2)*blockCount*sizeof(A)

    ! Memory displacement between each element of array
    disp_unit  = sizeof(A)

    ! Allocate RMA window for active and passive communications
#ifdef MPI_RMA_ACTIVE
    call MPI_WIN_ALLOCATE(RMA_size,disp_unit,mpi_info_key,solver_comm,RMA_ptr,RMA_win,ierr)
    call C_F_POINTER(RMA_ptr,dataTARGET,[Nyb+2+Nyb+2+Nxb+2+Nxb+2*blockCount])
#endif

#ifdef MPI_RMA_PASSIVE
    call MPI_WIN_ALLOCATE(RMA_size,disp_unit,MPI_INFO_NULL,solver_comm,RMA_ptr,RMA_win,ierr)
    call C_F_POINTER(RMA_ptr,dataTARGET,[Nyb+2+Nyb+2+Nxb+2+Nxb+2*blockCount])
#endif

#endif

    !call cpu_time(start)
    !start = omp_get_wtime()
    start = MPI_Wtime()

end subroutine MPIsolver_init 
