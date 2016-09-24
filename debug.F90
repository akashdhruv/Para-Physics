program debug

#include "Solver.h"

   implicit none

   integer, dimension(nblockx*nblocky) :: morton_blks

   integer, dimension(nblockx*nblocky) :: temp_blks
   integer, dimension(nblockx,nblocky,2) :: temp_id, morton_id
   integer :: i,j,k
   logical, dimension(nblockx*nblocky) :: mask

   mask = .true.

   temp_blks = (/(I,I=0,nblockx*nblocky-1)/) 

   k = 1

   do j=1,nblocky
    do i=1,nblockx

      temp_id(i,j,1) = mod(temp_blks(k),nblockx)
      temp_id(i,j,2) = temp_blks(k)/nblockx

      k = k+1

    end do
  end do
  
  morton_id(1,1,1) = temp_id(1,1,1)
  morton_id(1,1,2) = temp_id(1,1,2)

  mask = .false.

  do j=1,nblocky
     do i=

end program debug
