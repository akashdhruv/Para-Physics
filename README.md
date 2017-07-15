<h2> PARA-PHYSICS, A parallel block structured computational multi-physics solver </h2>

<h3> Important Information </h3>
<p align="justify"> <font size="4"> 
<ol>
	<li>This is version 2.0 of the software, previous versions can be found in the old releases</li>
	<li>The software can solve Incompressible Navier-Stokes and Heat advection diffusion equations for both singlephase and multiphase problems</li>
	<li>Immersed boundary method is implemented to compute flow over solid bodies</li>
	<li>Progress is underway to add new physics modules</li>
</ol>
</font> </p>

<h3> Parallelization Options </h3>
<p align="justify"> <font size="4">
<ol>
	<li>Standard MPI distributed memory model</li>
	<li>MPI Remote Memory Access (RMA) model</li>
	<li>MPI Shared Memory (SHM) model</li>
	<li>Hyper-threading using OpenMP</li>
</ol>
</font> </p>

<h3> Cache Optimization Features </h3>
<p align="justify"> <font size="4">
<ol>
	<li>Loop tiling for linear algebra solvers</li>
	<li>Array padding for spatial data</li>
	<li>Blocked data storage for stencil computations</li>
</ol>
</font> </p>

<h3> Compilation and Execution Instructions </h3>
<p align="justify"> <font size="4">
<ol>
	<li>Download the source code</li>
	<li>Make sure you have make utility and the latest version of GNU and MPI-3 or higher installed</li>
	<li>Edit the Makefile in ./src to include your MPI path
<pre class="bash">
src
make
cp Solver ../working/.
cd ../working
mpirun -n [number_of_procs] ./Solver
</pre>
</li>
	<li>You can also optimize Poisson solver and heat equation by hyper-threading using OpenMP on each core. The number of threads are specified 
     in Solver.h</li>
   <li>To plot results, edit the python file (plot.py) to match your grid size and simply type (make sure k = nblockx and d = nblocky)

<pre class="bash">
python plot.py
</pre>
</li>

</ol>
</font> </p>

<h3> Author - Akash V. Dhruv </h3>
<h3> License - Refer LICENSE.md </h3>
<h3> Cite as - Akash Dhruv. (2017, June 21). akidhruv/Para-Physics: Para-Physics. Zenodo. http://doi.org/10.5281/zenodo.815018 </h3>
