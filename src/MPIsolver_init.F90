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

    !_______________code for AMR - still in debug_________!

    blockCount = ((nblockx*nblocky)/procs)
    checkSumMPI = blockCount*procs

    if(procs /= nblockx*nblocky) then

        call MPI_FINALIZE(ierr)

        if (myid == 0) then
        print *,"RUNTIME ERROR: The number of MPI processes must be equal to the total number of blocks."
        end if

        call exit(status)

    end if


    !if (checkSumMPI /= nblockx*nblocky) then
       
    !    call MPI_FINALIZE(ierr)

    !    if (myid == 0) then 
    !    print *,"RUNTIME ERROR: The number blocks should be greater than and exactly divisible by total number of MPI processes."
    !    end if

    !    call exit(status)

    !end if

    blockID(1:blockCount) = (/(I,I=1,blockCount)/)
 
    call morton_sort(blockCount,myid,procs,blockLC)

    !_________Define Communication Based On Grid______________!
    call MPI_COMM_SPLIT(solver_comm,myid/nblockx,myid,x_comm,ierr)
    call MPI_COMM_SPLIT(solver_comm,mod(myid,nblockx),myid,y_comm,ierr)

    call MPI_COMM_RANK(x_comm,x_id,ierr)
    call MPI_COMM_SIZE(x_comm,x_procs,ierr)

    call MPI_COMM_RANK(y_comm,y_id,ierr)
    call MPI_COMM_size(y_comm,y_procs,ierr)

#ifdef MPI_DIST
    allocate(solnData(CENT_VAR,Nxb+2,Nyb+2))
    allocate(facexData(FACE_VAR,Nxb+2,Nyb+2))
    allocate(faceyData(FACE_VAR,Nxb+2,Nyb+2))
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
    center_size = CENT_VAR*(Nxb+2)*(Nyb+2)*sizeof(A)
    facex_size  = FACE_VAR*(Nxb+2)*(Nyb+2)*sizeof(A)
    facey_size  = FACE_VAR*(Nxb+2)*(Nyb+2)*sizeof(A) 

    disp_unit = sizeof(A)

    call MPI_WIN_ALLOCATE_SHARED(center_size,disp_unit,mpi_info_key,shared_comm,center_ptr,center_win,ierr)
    call MPI_WIN_ALLOCATE_SHARED(facex_size,disp_unit,mpi_info_key,shared_comm,facex_ptr,facex_win,ierr)
    call MPI_WIN_ALLOCATE_SHARED(facey_size,disp_unit,mpi_info_key,shared_comm,facey_ptr,facey_win,ierr)

    !__________________Point to local chunk of the shared data_______________________________!
    call MPI_WIN_SHARED_QUERY(center_win, shared_id, center_size, disp_unit, center_ptr,ierr)
    call MPI_BARRIER(shared_comm,ierr)
    call C_F_POINTER(center_ptr, solnData,[CENT_VAR,Nxb+2,Nyb+2])

    call MPI_WIN_SHARED_QUERY(facex_win, shared_id, facex_size, disp_unit, facex_ptr,ierr)
    call MPI_BARRIER(shared_comm,ierr)
    call C_F_POINTER(facex_ptr, facexData, [FACE_VAR,Nxb+2,Nyb+2])

    call MPI_WIN_SHARED_QUERY(facey_win, shared_id, facey_size, disp_unit, facey_ptr,ierr)
    call MPI_BARRIER(shared_comm,ierr)
    call C_F_POINTER(facey_ptr, faceyData, [FACE_VAR,Nxb+2,Nyb+2])

    !_____________________Point to the neighbour's data________________________________!
    if(x_id < x_procs-1 .and. shared_part(myid+1+1) /= MPI_UNDEFINED) then
    call MPI_WIN_SHARED_QUERY(center_win, shared_id+1 ,center_size,disp_unit,center_ptr,ierr)
    call C_F_POINTER(center_ptr,eastCENTER,[CENT_VAR,Nxb+2,Nyb+2])

    call MPI_WIN_SHARED_QUERY(facex_win, shared_id+1,facex_size,disp_unit,facex_ptr,ierr)
    call C_F_POINTER(facex_ptr,eastFACEX,[FACE_VAR,Nxb+2,Nyb+2])

    call MPI_WIN_SHARED_QUERY(facey_win, shared_id+1,facey_size,disp_unit,facey_ptr,ierr)
    call C_F_POINTER(facey_ptr,eastFACEY,[FACE_VAR,Nxb+2,Nyb+2])
    end if
    call MPI_BARRIER(shared_comm,ierr)

    if(x_id > 0 .and. shared_part(myid+1-1) /= MPI_UNDEFINED) then
    call MPI_WIN_SHARED_QUERY(center_win, shared_id-1 ,center_size,disp_unit,center_ptr,ierr)
    call C_F_POINTER(center_ptr,westCENTER,[CENT_VAR,Nxb+2,Nyb+2])

    call MPI_WIN_SHARED_QUERY(facex_win, shared_id-1,facex_size,disp_unit,facex_ptr,ierr)
    call C_F_POINTER(facex_ptr,westFACEX,[FACE_VAR,Nxb+2,Nyb+2])

    call MPI_WIN_SHARED_QUERY(facey_win, shared_id-1,facey_size,disp_unit,facey_ptr,ierr)
    call C_F_POINTER(facey_ptr,westFACEY,[FACE_VAR,Nxb+2,Nyb+2])
    end if
    call MPI_BARRIER(shared_comm,ierr)

    if(y_id < y_procs-1 .and. shared_part(myid+1+x_procs) /= MPI_UNDEFINED) then
    call MPI_WIN_SHARED_QUERY(center_win, shared_id+x_procs ,center_size,disp_unit,center_ptr,ierr)
    call C_F_POINTER(center_ptr,northCENTER,[CENT_VAR,Nxb+2,Nyb+2])

    call MPI_WIN_SHARED_QUERY(facex_win, shared_id+x_procs,facex_size,disp_unit,facex_ptr,ierr)
    call C_F_POINTER(facex_ptr,northFACEX,[FACE_VAR,Nxb+2,Nyb+2])

    call MPI_WIN_SHARED_QUERY(facey_win, shared_id+x_procs,facey_size,disp_unit,facey_ptr,ierr)
    call C_F_POINTER(facey_ptr,northFACEY,[FACE_VAR,Nxb+2,Nyb+2])
    end if
    call MPI_BARRIER(shared_comm,ierr)

    if(y_id > 0 .and. shared_part(myid+1-x_procs) /= MPI_UNDEFINED) then
    call MPI_WIN_SHARED_QUERY(center_win, shared_id-x_procs ,center_size,disp_unit,center_ptr,ierr)
    call C_F_POINTER(center_ptr,southCENTER,[CENT_VAR,Nxb+2,Nyb+2])

    call MPI_WIN_SHARED_QUERY(facex_win, shared_id-x_procs,facex_size,disp_unit,facex_ptr,ierr)
    call C_F_POINTER(facex_ptr,southFACEX,[FACE_VAR,Nxb+2,Nyb+2])

    call MPI_WIN_SHARED_QUERY(facey_win, shared_id-x_procs,facey_size,disp_unit,facey_ptr,ierr)
    call C_F_POINTER(facey_ptr,southFACEY,[FACE_VAR,Nxb+2,Nyb+2])
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

    allocate(solnData(CENT_VAR,Nxb+2,Nyb+2))
    allocate(facexData(FACE_VAR,Nxb+2,Nyb+2))
    allocate(faceyData(FACE_VAR,Nxb+2,Nyb+2))

    allocate(dataTARGET(Nyb+2+Nyb+2+Nxb+2+Nxb+2))

    allocate(eastORIGIN(Nyb+2))
    allocate(westORIGIN(Nyb+2))
    allocate(northORIGIN(Nxb+2))
    allocate(southORIGIN(Nxb+2))

    RMA_size   = (Nyb+2+Nxb+2+Nyb+2+Nxb+2)*sizeof(A)

    disp_unit  = sizeof(A)

#ifdef MPI_RMA_ACTIVE
    call MPI_WIN_CREATE(dataTARGET,RMA_size,disp_unit,mpi_info_key,solver_comm,RMA_win,ierr)
#endif

#ifdef MPI_RMA_PASSIVE
    call MPI_WIN_CREATE(dataTARGET,RMA_size,disp_unit,MPI_INFO_NULL,solver_comm,RMA_win,ierr)
#endif

#endif

    !call cpu_time(start)
    !start = omp_get_wtime()
    start = MPI_Wtime()

end subroutine MPIsolver_init 
