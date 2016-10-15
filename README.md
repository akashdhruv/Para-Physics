## Fortran based Hybrid MPI-OpenMP, Block Structured Computational Multi-Physics Solver 

### Important Information

  1. This is the latest version of ParaSolve, for previous versions look at other ParaSolve repositories.
  2. Single phase Heat Equation, Navier Stokes and Poisson Solver modules are working.
  3. Immersed boundary module is working.
  4. Conjugate Heat Transfer is in beta.
  5. AMR module not yet implemented.
  6. Multiphase module is in beta.

### Execution instructions

  1. Download the source code 
  2. Make sure you have make utility and the latest version of GNU and MPI installed
  3. Edit the Makefile to include your MPI path.

     ~~~terminal 
        make
        mpirun -n [number_of_procs] ./Solver 
     ~~~

  4. Note that the total number of MPI processes must be equal to the total number of blocks (nblockx X nblocky) defined in Solver.h

  5. You can also optimize Poisson solver by multi-threading using OpenMP on each MPI process. The number of threads are specified 
     in Solver.h

  6. To plot results, edit the python file (plot.py) to match your grid size and simply type (make sure k = nblockx and d = nblocky)

     ~~~terminal
        python plot.py
     ~~~ 

### Software outline

  1. The header file Solver.h handles simulation parameters like number of grid points per block, total number of blocks, module on/off, etc.

  2. The main program file is Solver.F90, which calls the functions Solver_init, Solver_evolve and Solver_finalize, which in turn calls
     module specific subroutines (module_init, module_solver, module_finalize).

  3. The best way to understand the code is to start with Solver.F90 and follow function calls one by one.

  4. All Grid data is stored in separate multi-dimensional arrays for cell-centers and faces located in the module, physicaldata.F90

  5. Module specific data is stored in files named module_data.F90. The interface for a given module, module_interface.F90 contains function
     definitions.

### Solver capabilities

<p align="center">
  <img src="./images/Image_3.png" width="900"/>
  <img src="./images/Image_4.png" width="1000"/>
</p>
<p align="center">
  Fig 1. Velocity Streamlines of the air flow around a solid cylinder and temperature distribution demonstarting Conjugate Heat Transfer
  (using immersed boundary and level set methods)
</p>

<p align="center">
  <img src="./images/Image_5.png" width="900"/>
  <img src="./images/Image_6.png" width="1000"/>
  <img src="./images/Image_7.png" width="1000"/>
</p>
<p align="center">
  Fig 2. Velocity Streamlines and Pressure and Temperature Contours for Lid Driven Cavity flow (Re = 1000)
</p>
 
### Author - Akash V. Dhruv  
### License - Refer LICENSE.md
