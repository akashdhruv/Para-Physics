subroutine MPIsolver_init()

#include "Solver.h"

    !$ use omp_lib

    use MPI_data
    use morton_interface, only: morton_sort
    use physicaldata

    implicit none

    integer :: checkSumMPI
    integer :: status = 0
    integer :: i
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
       
        call MPI_FINALIZE(ierr)
        if (myid == 0) &
        print *,"RUNTIME ERROR: The number of blocks should be greater than/equal to &
                 and exactly divisible by total number of MPI processes."
        call exit(status)

    end if

    allocate(blockID(nblockx*nblocky))
    allocate(blockLC(nblockx*nblocky))
   
    do i=0,procs-1

       blockLC(1+i*blockCount:blockCount+i*blockCount) = i
       blockID(1+i*blockCount:blockCount+i*blockCount) = (/(I,I=1,blockCount)/)

    end do

    !_______________Apply Refinement Using AMR - still in debug_________!

    call morton_sort(blockCount,myid,procs,blockID,blockLC)
 
    !_________Define Communication Based On Grid______________!
    call MPI_COMM_SPLIT(solver_comm,myid/nblockx,myid,x_comm,ierr)
    call MPI_COMM_SPLIT(solver_comm,mod(myid,nblockx),myid,y_comm,ierr)

    call MPI_COMM_RANK(x_comm,x_id,ierr)
    call MPI_COMM_SIZE(x_comm,x_procs,ierr)

    call MPI_COMM_RANK(y_comm,y_id,ierr)
    call MPI_COMM_size(y_comm,y_procs,ierr)

#ifdef MPI_DIST
    allocate(localCENTER(Nxb+2,Nyb+2,CENT_VAR,blockCount))
    allocate(localFACEX(Nxb+2,Nyb+2,FACE_VAR,blockCount))
    allocate(localFACEY(Nxb+2,Nyb+2,FACE_VAR,blockCount))
#endif

#ifdef MPI_SHRD
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
    call C_F_POINTER(center_ptr, localCENTER,[Nxb+2,Nyb+2,CENT_VAR,blockCount])

    call MPI_WIN_SHARED_QUERY(facex_win, shared_id, facex_size, disp_unit, facex_ptr,ierr)
    call MPI_BARRIER(shared_comm,ierr)
    call C_F_POINTER(facex_ptr, localFACEX, [Nxb+2,Nyb+2,FACE_VAR,blockCount])

    call MPI_WIN_SHARED_QUERY(facey_win, shared_id, facey_size, disp_unit, facey_ptr,ierr)
    call MPI_BARRIER(shared_comm,ierr)
    call C_F_POINTER(facey_ptr, localFACEY, [Nxb+2,Nyb+2,FACE_VAR,blockCount])

    !_____________________Point to the neighbour's data________________________________!
    if(x_id < x_procs-1 .and. shared_part(myid+1+1) /= MPI_UNDEFINED) then
    call MPI_WIN_SHARED_QUERY(center_win, shared_id+1 ,center_size,disp_unit,center_ptr,ierr)
    call C_F_POINTER(center_ptr,eastCENTER,[Nxb+2,Nyb+2,CENT_VAR,blockCount])

    call MPI_WIN_SHARED_QUERY(facex_win, shared_id+1,facex_size,disp_unit,facex_ptr,ierr)
    call C_F_POINTER(facex_ptr,eastFACEX,[Nxb+2,Nyb+2,FACE_VAR,blockCount])

    call MPI_WIN_SHARED_QUERY(facey_win, shared_id+1,facey_size,disp_unit,facey_ptr,ierr)
    call C_F_POINTER(facey_ptr,eastFACEY,[Nxb+2,Nyb+2,FACE_VAR,blockCount])
    end if
    call MPI_BARRIER(shared_comm,ierr)

    if(x_id > 0 .and. shared_part(myid+1-1) /= MPI_UNDEFINED) then
    call MPI_WIN_SHARED_QUERY(center_win, shared_id-1 ,center_size,disp_unit,center_ptr,ierr)
    call C_F_POINTER(center_ptr,westCENTER,[Nxb+2,Nyb+2,CENT_VAR,blockCount])

    call MPI_WIN_SHARED_QUERY(facex_win, shared_id-1,facex_size,disp_unit,facex_ptr,ierr)
    call C_F_POINTER(facex_ptr,westFACEX,[Nxb+2,Nyb+2,FACE_VAR,blockCount])

    call MPI_WIN_SHARED_QUERY(facey_win, shared_id-1,facey_size,disp_unit,facey_ptr,ierr)
    call C_F_POINTER(facey_ptr,westFACEY,[Nxb+2,Nyb+2,FACE_VAR,blockCount])
    end if
    call MPI_BARRIER(shared_comm,ierr)

    if(y_id < y_procs-1 .and. shared_part(myid+1+x_procs) /= MPI_UNDEFINED) then
    call MPI_WIN_SHARED_QUERY(center_win, shared_id+x_procs ,center_size,disp_unit,center_ptr,ierr)
    call C_F_POINTER(center_ptr,northCENTER,[Nxb+2,Nyb+2,CENT_VAR,blockCount])

    call MPI_WIN_SHARED_QUERY(facex_win, shared_id+x_procs,facex_size,disp_unit,facex_ptr,ierr)
    call C_F_POINTER(facex_ptr,northFACEX,[Nxb+2,Nyb+2,FACE_VAR,blockCount])

    call MPI_WIN_SHARED_QUERY(facey_win, shared_id+x_procs,facey_size,disp_unit,facey_ptr,ierr)
    call C_F_POINTER(facey_ptr,northFACEY,[Nxb+2,Nyb+2,FACE_VAR,blockCount])
    end if
    call MPI_BARRIER(shared_comm,ierr)

    if(y_id > 0 .and. shared_part(myid+1-x_procs) /= MPI_UNDEFINED) then
    call MPI_WIN_SHARED_QUERY(center_win, shared_id-x_procs ,center_size,disp_unit,center_ptr,ierr)
    call C_F_POINTER(center_ptr,southCENTER,[Nxb+2,Nyb+2,CENT_VAR,blockCount])

    call MPI_WIN_SHARED_QUERY(facex_win, shared_id-x_procs,facex_size,disp_unit,facex_ptr,ierr)
    call C_F_POINTER(facex_ptr,southFACEX,[Nxb+2,Nyb+2,FACE_VAR,blockCount])

    call MPI_WIN_SHARED_QUERY(facey_win, shared_id-x_procs,facey_size,disp_unit,facey_ptr,ierr)
    call C_F_POINTER(facey_ptr,southFACEY,[Nxb+2,Nyb+2,FACE_VAR,blockCount])
    end if
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

    allocate(localCENTER(Nxb+2,Nyb+2,CENT_VAR,blockCount))
    allocate(localFACEX(Nxb+2,Nyb+2,FACE_VAR,blockCount))
    allocate(localFACEY(Nxb+2,Nyb+2,FACE_VAR,blockCount))

    allocate(eastORIGIN(Nyb+2))
    allocate(westORIGIN(Nyb+2))
    allocate(northORIGIN(Nxb+2))
    allocate(southORIGIN(Nxb+2))

    RMA_size   = (Nyb+2+Nxb+2+Nyb+2+Nxb+2)*sizeof(A)

    disp_unit  = sizeof(A)

#ifdef MPI_RMA_ACTIVE
    call MPI_WIN_ALLOCATE(RMA_size,disp_unit,mpi_info_key,solver_comm,RMA_ptr,RMA_win,ierr)
    call C_F_POINTER(RMA_ptr,dataTARGET,[Nyb+2+Nyb+2+Nxb+2+Nxb+2])
#endif

#ifdef MPI_RMA_PASSIVE
    call MPI_WIN_ALLOCATE(RMA_size,disp_unit,MPI_INFO_NULL,solver_comm,RMA_ptr,RMA_win,ierr)
    call C_F_POINTER(RMA_ptr,dataTARGET,[Nyb+2+Nyb+2+Nxb+2+Nxb+2])
#endif

#endif

    !call cpu_time(start)
    !start = omp_get_wtime()
    start = MPI_Wtime()

end subroutine MPIsolver_init 
