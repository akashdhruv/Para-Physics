## PARA-PHYSICS, A parallel block structured computational multi-physics solver 

### Important Information

  1. This is version 2.0 of the software, previous versions can be found in the old releases  
  2. The software can solve Incompressible Navier-Stokes and Heat advection diffusion equations for both singlephase and multiphase problems
  3. Immersed boundary method is implemented to compute flow over solid bodies
  4. Progress is underway to add new physics modules

### Parallelization Options

  1. Standard MPI distributed memory model
  2. MPI Remote Memory Access (RMA) model
  3. MPI Shared Memory (SHM) model
  4. Hyper-threading using OpenMP

### Cache Optimization Features

  1. Loop tiling for linear algebra solvers
  2. Array padding for spatial data
  3. Blocked data storage for stencil computations

### Compilation and Execution Instructions

  1. Download the source code.
  2. Make sure you have make utility and the latest version of GNU and MPI-3 or higher installed
  3. Edit the Makefile in ./src to include your MPI path

     ~~~terminal 
        cd src
        make
        cp Solver ../working/.
        cd ../working
        mpirun -n [number_of_procs] ./Solver 
     ~~~
  4. You can also optimize Poisson solver and heat equation by hyper-threading using OpenMP on each core. The number of threads are specified 
     in Solver.h

  5. To plot results, edit the python file (plot.py) to match your grid size and simply type (make sure k = nblockx and d = nblocky)

     ~~~terminal
        python plot.py
     ~~~ 

## Examples

### 1. Wake suppression through conjugate heat transfer 

**<p align="center">
  <img src="./images/Vort_WakeSup.png" width="700"/>
  Figure 1. Vorticity contours at t = 100 units for Re = 2500 (a) Without heat transfer (b) With heat transfer (Block size - 20 x 20) (Num blocks - 8 x 4)
</p>**

**<p align="center">
  <img src="./images/Temp_Dens.png" width="400"/>  
  Figure 2. (a) Non-dimensional temperature distribution (b) Density ratio for conjugate heat transfer at Re = 2500
</p>**

### 2. Lid Driven Cavity

**<p align="center">
  <img src="./images/Ghia_Comparison.png" width="700"/>
  Figure 3. (a) Numerical solution, (Block size - 20 x 20) (Num blocks - 6 x 6) (b) Reference solution by [Ghia .et .al]
</p>**

### 3. Conjugate heat transfer between fluid and multiple immersed boundaries

**<p align="center">
  <img src="./images/Image_10.png" width="700"/>
  Figure 4. Conjugate heat transfer between air and three solid cylinders, Re = 500, t = 25 units, (Block size - 20 x 20) (Num blocks - 40 x 40)
</p>**

### 4. MPI shared memory vs distributed memory model

**<p align="center">
  <img src="./images/XeonPhi_Sp.png" width="700"/>
  Figure 5. Speed up comparison on a 64 core Xeon Phi KNL processor  
</p>**

MPI SHM gives almost twice as much speed up in comparison to traditional distributed memory model.

##__________________________________________
### Author - Akash V. Dhruv  
### License - Refer LICENSE.md

### Cite as - Akash Dhruv. (2017, June 21). akidhruv/Para-Physics: Para-Physics. Zenodo. http://doi.org/10.5281/zenodo.815018

[Ghia .et .al]: https://pdfs.semanticscholar.org/211b/45b6a06336a72ca064a6e59b14ebc520211c.pdf
