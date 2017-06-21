## Fortran based Hybrid MPI-OpenMP, Block Structured Computational Multi-Physics Solver 

### Important Information

  1. This is version 1.0 of **Para-Physics**, a computational multi-physics solver.
  2. Single phase heat equation, Navier-Stokes and Poisson solver modules are working.
  3. Immersed boundary module is working.
  4. Conjugate Heat Transfer is in beta.
  5. AMR module not yet implemented.
  6. Multiphase module is in beta.

### Execution instructions

  1. Download the source code 
  2. Make sure you have make utility and the latest version of GNU and MPI installed
  3. Edit the Makefile in ./src to include your MPI path.

     ~~~terminal 
        cd src
        make
        cp Solver ../working/.
        cd ../working
        mpirun -n [number_of_procs] ./Solver 
     ~~~

  4. Note that the total number of MPI processes must be equal to the total number of blocks (nblockx X nblocky) defined in Solver.h

  5. You can also optimize Poisson solver and heat equation by multi-threading using OpenMP on each MPI process. The number of threads are specified 
     in Solver.h

  6. To plot results, edit the python file (plot.py) to match your grid size and simply type (make sure k = nblockx and d = nblocky)

     ~~~terminal
        python plot.py
     ~~~ 

  7. Use make clean in src and working folder to remove rebuildables.

### Software outline

  1. The header file Solver.h handles simulation parameters like number of grid points per block, total number of blocks, module on/off, etc.

  2. The main program file is Solver.F90, which calls the functions Solver_init, Solver_evolve and Solver_finalize, which in turn calls
     module specific subroutines (module_init, module_solver, module_finalize).

  3. The best way to understand the code is to start with Solver.F90 and follow function calls one by one.

  4. All Grid data is stored in separate multi-dimensional arrays for cell-centers and faces located in the module, physicaldata.F90

  5. Module specific data is stored in files named module_data.F90. The interface for a given module, module_interface.F90 contains function
     definitions.

### Solver capabilities

#### CASE 1 - Air flow ovear a solid cylinder with Conjugate Heat Transfer (CHT), Re = 500
##### (Grid - 800 x 400) (4 x 4 MPI processes) (2 OpenMP threads per process)

<p align="center">
  <img src="./images/Image_3.png" width="700"/>
  <img src="./images/Image_4.png" width="700"/>
</p>

<p align="center">
  Figure 1. Velocity Streamlines and Temperature Contours, Re = 500, t = 30.0 s. 
</p>

<p align="center">
  <img src="./images/Image_8.png" width="700"/>
  <img src="./images/Image_9.png" width="700"/>
</p>

<p align="center">
  Figure 2. Vorticity Contours, Re = 500, t=60.0 s
</p>

#### CASE 2 - Lid Driven Cavity Flow, Re = 1000
##### (Grid - 100 x 80) (2 x 2 MPI processes) (4 OpenMP threads per process)

<p align="center">
  <img src="./images/Image_5.png" width="700"/>
  <img src="./images/Image_6.png" width="700"/>
  <img src="./images/Image_7.png" width="700"/>
</p>

<p align="center">
  Figure 3. Velocity Streamlines and Pressure and Temperature Contours for Lid Driven Cavity flow, Re = 1000.
</p>

#### CASE 3 - Conjugate Heat Transfer between fluid and multiple immersed boundaries, Re = 500
##### (Grid - 800 x 400) (4 x 4 MPI processes) (2 OpenMP threads per process)

<p align="center">
  <img src="./images/Image_10.png" width="700"/>
</p>

<p align="center">
  Figure 4. Conjugate Heat Transfer between air and three solid cylinders, Re = 500, t=25 s.
</p>

### Author - Akash V. Dhruv  
### License - Refer LICENSE.md

### Cite as - Akash Dhruv. (2017, June 21). akidhruv/Para-Physics: Para-Physics. Zenodo. http://doi.org/10.5281/zenodo.815018
