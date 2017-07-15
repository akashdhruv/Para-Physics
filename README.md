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

 1. Download the source code
 2. Make sure you have make utility and the latest version of GNU and MPI-3 or higher installed
 3. Edit the Makefile in ./src to include your MPI path

   ~~~terminal
	src
	make
	cp Solver ../working/.
	cd ../working
	mpirun -n [number_of_procs] ./Solver
   ~~~

 4. You can also optimize Poisson solver and heat equation by hyper-threading using OpenMP on each core. The number of threads are specified in Solver.h

 5. To plot results, edit the python file (plot.py) to match your grid size and simply type (make sure k = nblockx and d = nblocky)

	~~~terminal
	python plot.py
	~~~

### Author - Akash V. Dhruv
### License - Refer LICENSE.md
### Cite as - Akash Dhruv. (2017, June 21). akidhruv/Para-Physics: Para-Physics. Zenodo. http://doi.org/10.5281/zenodo.815018
