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
    call MPI_COMM_RANK(solver_comm, myid, ierr)
    call MPI_COMM_SIZE(solver_comm, procs, ierr)
    call MPI_COMM_GROUP(solver_comm, world_group, ierr)

    allocate(world_part(procs))
    world_part = (/(I,I=0,procs-1)/)

    !_______________Split Domain Into Blocks For Cache Optimization_________!

    blockCount = ((nblockx*nblocky)/procs)
    checkSumMPI = blockCount*procs
   
    blockOffset = myid*blockCount

    if (checkSumMPI /= nblockx*nblocky) then
       
        if (myid == 0) &
        print *,"RUNTIME ERROR: The number of blocks should be greater than/equal to &
                 and exactly divisible by total number of MPI processes."

        call MPI_FINALIZE(ierr)
        call exit(status)

    else if(blockCount > MAX_BLOCKS) then
    
        if (myid == 0) &
        print *,"RUNTIME ERROR: The total number of blocks per process exceed the maximum limit. &
                 Increase the number of MPI jobs."

        call MPI_FINALIZE(ierr)
        call exit(status)

    end if

    allocate(blockID(nblockx*nblocky))
    allocate(blockLC(nblockx*nblocky))
   
    do i=1,procs

       blockLC(1+(i-1)*blockCount:blockCount+(i-1)*blockCount) = i-1
       blockID(1+(i-1)*blockCount:blockCount+(i-1)*blockCount) = (/(I,I=1,blockCount)/)

    end do

    !_______________Apply Refinement Using AMR - still in debug_________!

    !call morton_sort(blockCount,myid,procs,blockID,blockLC)
 
    !_________Define Communication Based On Grid______________!
    call MPI_COMM_SPLIT(solver_comm,myid/nblockx,myid,x_comm,ierr)
    call MPI_COMM_SPLIT(solver_comm,mod(myid,nblockx),myid,y_comm,ierr)

    call MPI_COMM_RANK(x_comm,x_id,ierr)
    call MPI_COMM_SIZE(x_comm,x_procs,ierr)

    call MPI_COMM_RANK(y_comm,y_id,ierr)
    call MPI_COMM_size(y_comm,y_procs,ierr)

    allocate(xLC(blockCount))
    allocate(yLC(blockCount))
    allocate(reqs(blockCount*4*2))
    allocate(req_stat(blockCount*4*2*MPI_STATUS_SIZE))

    xLC = mod(((/(I,I=0,blockCount-1)/) + blockOffset),nblockx)
    yLC = ((/(I,I=0,blockCount-1)/) + blockOffset)/nblockx

#ifdef MPI_DIS
    allocate(localCENTER(Nxb+2,Nyb+2,blockCount,CENT_VAR))
    allocate(localFACEX(Nxb+2,Nyb+2,blockCount,FACE_VAR))
    allocate(localFACEY(Nxb+2,Nyb+2,blockCount,FACE_VAR))
#endif

#ifdef MPI_SHM
    !__________Define Shared Communication Environment__________!
    call MPI_COMM_SPLIT_TYPE(solver_comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, shared_comm, ierr)
    call MPI_COMM_RANK(shared_comm,shared_id,ierr)
    call MPI_COMM_SIZE(shared_comm,shared_procs,ierr)
    call MPI_COMM_GROUP(shared_comm, shared_group,ierr)

    allocate(shared_part(procs))
    call MPI_GROUP_TRANSLATE_RANKS(world_group,procs,world_part,shared_group,shared_part,ierr)

    call MPI_INFO_CREATE(mpi_info_key,ierr)
    call MPI_INFO_SET(mpi_info_key,"alloc_shared_noncontig","true",ierr)

   !_________Make on-node processes allocate their chunk of shared Memory____________________!
    center_size = blockCount*CENT_VAR*(Nxb+2)*(Nyb+2)*sizeof(A)
    facex_size  = blockCount*FACE_VAR*(Nxb+2)*(Nyb+2)*sizeof(A)
    facey_size  = blockCount*FACE_VAR*(Nxb+2)*(Nyb+2)*sizeof(A) 

    disp_unit = sizeof(A)

    call MPI_WIN_ALLOCATE_SHARED(center_size,disp_unit,MPI_INFO_NULL,shared_comm,center_ptr,center_win,ierr)
    call MPI_WIN_ALLOCATE_SHARED(facex_size,disp_unit,MPI_INFO_NULL,shared_comm,facex_ptr,facex_win,ierr)
    call MPI_WIN_ALLOCATE_SHARED(facey_size,disp_unit,MPI_INFO_NULL,shared_comm,facey_ptr,facey_win,ierr)

    !__________________Point to local chunk of the shared data_______________________________!
    call MPI_WIN_SHARED_QUERY(center_win, shared_id, center_size, disp_unit, center_ptr,ierr)
    call MPI_BARRIER(shared_comm,ierr)
    call C_F_POINTER(center_ptr, localCENTER,[Nxb+2,Nyb+2,blockCount,CENT_VAR])

    call MPI_WIN_SHARED_QUERY(facex_win, shared_id, facex_size, disp_unit, facex_ptr,ierr)
    call MPI_BARRIER(shared_comm,ierr)
    call C_F_POINTER(facex_ptr, localFACEX, [Nxb+2,Nyb+2,blockCount,FACE_VAR])

    call MPI_WIN_SHARED_QUERY(facey_win, shared_id, facey_size, disp_unit, facey_ptr,ierr)
    call MPI_BARRIER(shared_comm,ierr)
    call C_F_POINTER(facey_ptr, localFACEY, [Nxb+2,Nyb+2,blockCount,FACE_VAR])

    !____________________________Point to the shared data________________________________!
    call MPI_WIN_SHARED_QUERY(center_win, 0 ,center_size,disp_unit,center_ptr,ierr)
    call C_F_POINTER(center_ptr,sharedCENTER,[Nxb+2,Nyb+2,blockCount,CENT_VAR*shared_procs])
    call MPI_BARRIER(shared_comm,ierr)

    call MPI_WIN_SHARED_QUERY(facex_win, 0,facex_size,disp_unit,facex_ptr,ierr)
    call C_F_POINTER(facex_ptr,sharedFACEX,[Nxb+2,Nyb+2,blockCount,FACE_VAR*shared_procs])
    call MPI_BARRIER(shared_comm,ierr)

    call MPI_WIN_SHARED_QUERY(facey_win, 0,facey_size,disp_unit,facey_ptr,ierr)
    call C_F_POINTER(facey_ptr,sharedFACEY,[Nxb+2,Nyb+2,blockCount,FACE_VAR*shared_procs])
    call MPI_BARRIER(shared_comm,ierr)
#endif

#ifdef MPI_RMA
    !__________Define Shared Communication Environment__________!
    call MPI_COMM_SPLIT_TYPE(solver_comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, shared_comm, ierr)
    call MPI_COMM_RANK(shared_comm,shared_id,ierr)
    call MPI_COMM_SIZE(shared_comm,shared_procs,ierr)
    call MPI_COMM_GROUP(shared_comm, shared_group,ierr)

    allocate(shared_part(procs))
    call MPI_GROUP_TRANSLATE_RANKS(world_group,procs,world_part,shared_group,shared_part,ierr)

    call MPI_INFO_CREATE(mpi_info_key,ierr)
    call MPI_INFO_SET(mpi_info_key,"no_locks","true",ierr)

    allocate(localCENTER(Nxb+2,Nyb+2,blockCount,CENT_VAR))
    allocate(localFACEX(Nxb+2,Nyb+2,blockCount,FACE_VAR))
    allocate(localFACEY(Nxb+2,Nyb+2,blockCount,FACE_VAR))

    allocate(eastORIGIN(Nyb+2,blockCount))
    allocate(westORIGIN(Nyb+2,blockCount))
    allocate(northORIGIN(Nxb+2,blockCount))
    allocate(southORIGIN(Nxb+2,blockCount))

    RMA_size   = (Nyb+2+Nxb+2+Nyb+2+Nxb+2)*blockCount*sizeof(A)

    disp_unit  = sizeof(A)

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
