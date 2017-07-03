subroutine MPI_applyBC_ORIG(u_ex)

#include "Solver.h"

      use MPI_data

      implicit none

      real, dimension(:,:), intent(inout) :: u_ex
      integer :: status(MPI_STATUS_SIZE)

      if (x_procs > 1) then
        
       if(mod(x_id,2) == 0) then
           
             if(x_id == 0) then

                 call MPI_SENDRECV(u_ex(Nxb+1,:), Nyb+2, MPI_REAL, mod(x_id+1,x_procs), 1,&
                                   u_ex(Nxb+2,:), Nyb+2, MPI_REAL, mod(x_id+1,x_procs), 2,x_comm, status, ierr) 

             else if(x_id == nblockx-1) then
            
                call MPI_SENDRECV(u_ex(2,:), Nyb+2, MPI_REAL, mod(x_id-1+x_procs,x_procs), 3,&
                                  u_ex(1,:), Nyb+2, MPI_REAL, mod(x_id-1+x_procs,x_procs), 4,x_comm,status, ierr)

             else
                call MPI_SENDRECV(u_ex(Nxb+1,:), Nyb+2, MPI_REAL, mod(x_id+1,x_procs), 1,&
                                  u_ex(Nxb+2,:), Nyb+2, MPI_REAL, mod(x_id+1,x_procs), 2, x_comm, status, ierr) 
                  
                call MPI_SENDRECV(u_ex(2,:), Nyb+2, MPI_REAL, mod(x_id-1+x_procs,x_procs), 3,&
                                  u_ex(1,:), Nyb+2, MPI_REAL, mod(x_id-1+x_procs,x_procs), 4,x_comm,status, ierr)
                                   

             end if

       else if (mod(x_id,2) == 1) then

             if(x_id == nblockx-1) then
           
               call MPI_SENDRECV(u_ex(2,:), Nyb+2,MPI_REAL, mod(x_id-1+x_procs,x_procs), 2,&
                                 u_ex(1,:), Nyb+2,MPI_REAL, mod(x_id-1+x_procs,x_procs), 1,x_comm, status, ierr)

             else
             
               call MPI_SENDRECV(u_ex(2,:), Nyb+2, MPI_REAL, mod(x_id-1+x_procs,x_procs), 2,&
                                 u_ex(1,:), Nyb+2, MPI_REAL, mod(x_id-1+x_procs,x_procs), 1,x_comm,status, ierr)

               call MPI_SENDRECV(u_ex(Nxb+1,:), Nyb+2, MPI_REAL, mod(x_id+1,x_procs), 4,&
                                 u_ex(Nxb+2,:), Nyb+2, MPI_REAL, mod(x_id+1,x_procs), 3,x_comm,status,ierr)
                  

             end if

       end if

     end if

     !! Second dimension !!

     if (y_procs > 1) then

       if(mod(y_id,2) == 0) then

             if(y_id == 0) then

                  call MPI_SENDRECV(u_ex(:,Nyb+1), Nxb+2, MPI_REAL, mod(y_id+1,y_procs), 5,&
                                    u_ex(:,Nyb+2), Nxb+2, MPI_REAL, mod(y_id+1,y_procs), 6,y_comm,status, ierr)

             else if(y_id == nblocky-1) then

                  call MPI_SENDRECV(u_ex(:,2), Nxb+2, MPI_REAL, mod(y_id-1+y_procs,y_procs), 7,&
                                    u_ex(:,1), Nxb+2, MPI_REAL, mod(y_id-1+y_procs,y_procs), 8,y_comm,status,ierr)

             else 
                  call MPI_SENDRECV(u_ex(:,Nyb+1), Nxb+2, MPI_REAL, mod(y_id+1,y_procs), 5,&
                                    u_ex(:,Nyb+2), Nxb+2, MPI_REAL, mod(y_id+1,y_procs), 6,y_comm,status,ierr)

                  call MPI_SENDRECV(u_ex(:,2), Nxb+2, MPI_REAL, mod(y_id-1+y_procs,y_procs), 7,&
                                    u_ex(:,1), Nxb+2, MPI_REAL, mod(y_id-1+y_procs,y_procs), 8,y_comm,status,ierr)
                  

             end if

       else if (mod(y_id,2) == 1) then

             if(y_id == nblocky-1) then

                  call MPI_SENDRECV(u_ex(:,2), Nxb+2, MPI_REAL, mod(y_id-1+y_procs,y_procs), 6,&
                                    u_ex(:,1), Nxb+2, MPI_REAL, mod(y_id-1+y_procs,y_procs), 5,y_comm,status, ierr)

             else

                  call MPI_SENDRECV(u_ex(:,2), Nxb+2, MPI_REAL, mod(y_id-1+y_procs,y_procs), 6,&
                                    u_ex(:,1), Nxb+2, MPI_REAL, mod(y_id-1+y_procs,y_procs), 5,y_comm,status, ierr)

 
                  call MPI_SENDRECV(u_ex(:,Nyb+1), Nxb+2, MPI_REAL, mod(y_id+1,y_procs), 8,&
                                    u_ex(:,Nyb+2), Nxb+2, MPI_REAL, mod(y_id+1,y_procs), 7,y_comm,status, ierr)
                  

             end if

       end if            

    end if
  
end subroutine MPI_applyBC_ORIG
