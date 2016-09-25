subroutine morton_sort(blockCount,myid,procs,blockLC)

#include "Solver.h"

   implicit none

   integer, dimension(:,:), intent(inout) :: blockLC
   integer,intent(in) :: blockCount 
   integer,intent(in) :: myid
   integer,intent(in) :: procs

   integer, dimension(nblockx*nblocky) :: morton_blks
   integer, dimension(nblockx*nblocky) :: temp_blks
   integer, dimension(nblockx*nblocky,2) :: temp_id, morton_id
   integer :: i,k
   logical, dimension(nblockx*nblocky) :: mask

   mask = .true.

   temp_blks = (/(I,I=0,nblockx*nblocky-1)/) 

   temp_id(:,1) = mod(temp_blks,nblockx)
   temp_id(:,2) = temp_blks/nblockx

   k = 1

   do i=1,nblockx*nblocky

      if (mask(i)) then

         morton_blks(k) = temp_blks(i)         
         mask(i) = .false.
         k = k+1

         if (temp_id(i,2) == temp_id(i+1,2) .and. mask(i+1)) then

         morton_blks(k) = temp_blks(i+1)
         mask(i+1) = .false.
         k = k+1

         end if

         if(i+nblockx <= nblockx*nblocky) then        

            if(temp_id(i,1) == temp_id(i+nblockx,1) .and. mask(i+nblockx)) then

            morton_blks(k) = temp_blks(i+nblockx)
            mask(i+nblockx) = .false.
            k = k+1

            end if

         end if

         if(i+1+nblockx <= nblockx*nblocky) then

            if(temp_id(i,2) == temp_id(i+1+nblockx,2)-1 .and. temp_id(i,1) == temp_id(i+1+nblockx,1)-1 .and. mask(i+1+nblockx)) then

            morton_blks(k) = temp_blks(i+1+nblockx)
            mask(i+1+nblockx) = .false.
            k = k+1

            end if
          
         end if
   
      else 
           
          continue  

      end if

   end do

   blockLC(:,1) = morton_blks

   do i=0,procs-1

      blockLC(i*blockCount+1:i*blockCount+1+blockCount,2) = i

   end do

   
end subroutine morton_sort
