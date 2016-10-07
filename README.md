## Fortran based Hybrid MPI-OpenMP, Block Structured Computational Multi-Physics Solver 

### Important Information

  1. This is the latest version of ParaSolve, for previous versions look at other ParaSolve repositories.
  2. Single phase Heat Equation, Navier Stokes and Poisson Solver modules are working.
  3. Immersed boundary module is working.
  4. Conjugate Heat Transfer is in beta.
  5. AMR module not yet implemented.
  6. Multiphase module is still in beta.

### Execution instructions

  1. Download the source code 
  2. Make sure you have make utility and the latest version of GNU and MPI installed
  3. Edit the Makefile to include your MPI path.

     ~~~terminal 
        make
        mpirun -n [number_of_procs] ./Solver 
     ~~~

  5. Note that the total number of MPI processes must be equal to the total number of blocks (nblockx X nblocky) defined in Solver.h

  4. To plot results, edit the python file (plot.py) to match your grid size and simply type (make sure k = nblockx and d = nblocky)

     ~~~terminal
        python plot.py
     ~~~ 

### Author - Akash V. Dhruv  
### License - Refer LICENSE.md
