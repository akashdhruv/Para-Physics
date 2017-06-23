#MPI_PATH = /usr/local
#MPI_PATH = /usr
#MPI_PATH = /opt/mpich2-1.4.1p1/
MPI_PATH = /c1/apps/mpich/3.1.4
#HYPRE_PATH    = /usr/local/hypre
#PGI_PATH = /opt/pgi/osx86-64/2016
#MPI_PATH = /usr/lib64/mpich 

MPIFF = $(MPI_PATH)/bin/mpif90
#MPIFF = mpiifort
#MPIFF = mpif90
#FF = ifort
#MPIFF = $(PGI_PATH)/bin/pgf90

#LIB_HYPRE = -L${HYPRE_PATH}/lib -lHYPRE

#FFLAGS = -c
FFLAGS = -c -O3 -ftree-vectorize
#FFLAGS = -fopenmp -O3 -ftree-vectorize -c
#FFLAGS = -qopenmp -c 
#FFLAGS = -Mmpi=mpich -c

#EXEFLAGS = -o
EXEFLAGS = -O3 -ftree-vectorize -o
#EXEFLAGS = -fopenmp -O3 -ftree-vectorize -o
#EXEFLAGS = -qopenmp -o
#EXEFLAGS = -Mmpi=mpich -o

EXE = Solver

#%.o %.mod: %.F90
#	$(MPIFF) $(FFLAGS) $<


Modules += Grid_interface.mod \
           Grid_data.mod \
           IO_interface.mod \
           IncompNS_interface.mod \
           IncompNS_data.mod \
           Poisson_interface.mod \
           MPI_interface.mod \
           MPI_data.mod \
           Solver_interface.mod \
           Driver_interface.mod \
           physicaldata.mod \
	   Driver_data.mod \
           HeatAD_interface.mod \
           HeatAD_data.mod \
           Multiphase_interface.mod \
           Multiphase_data.mod \
           morton_interface.mod \
           IBM_interface.mod IBM_data.mod

Objects += physicaldata.o \
           Grid_data.o \
           IncompNS_data.o \
           HeatAD_data.o \
           Driver_data.o \
           MPI_data.o \
           Multiphase_data.o \
           IBM_data.o \
           MPIsolver_init.o \
           Grid_init.o \
           IncompNS_init.o \
           Driver_init.o \
           HeatAD_init.o \
           IBM_init.o \
           Multiphase_init.o \
           Solver_init.o \
           MPIsolver_finalize.o \
           Grid_finalize.o \
           Solver_finalize.o\
           ins_rescaleVel.o \
           IO_write.o \
           IO_display.o \
           IO_display_v2.o \
           morton_sort.o \
           MPI_periodicBC.o \
           MPI_CollectResiduals.o \
           MPI_applyBC_v2.o\
           MPI_applyBC.o \
           MPI_applyBC_shared.o \
           MPI_physicalBC_vel.o \
           MPI_physicalBC_pres.o \
           MPI_physicalBC_temp.o \
           MPI_physicalBC_dfun.o \
           MPI_physicalBC_vort.o \
           mph_FillVars.o \
           mph_FillVars_ibm.o \
           mph_PressureJumps.o \
           mph_getInterfaceVelocity.o \
           IBM_ApplyForcing.o \
           heat_tempSolver.o \
           heat_tempSolver_ibm.o \
           heat_tempSolver_mph.o \
           HeatAD_solver.o \
           ibm_evolve.o \
           Poisson_solver.o \
           ins_momentum.o \
           ins_vorticity.o \
           Poisson_solver_VC.o \
           ins_momentum_VD.o \
           IncompNS_solver.o \
           Multiphase_solver.o \
           IBM_solver.o \
           Solver_evolve.o \
           Solver.o

ALL_OBJS = $(Modules) $(Objects)

LINKER_OBJS = $(filter-out $(Modules),$(ALL_OBJS))


$(EXE): $(ALL_OBJS) 
	$(MPIFF) $(EXEFLAGS) $(EXE) $(LINKER_OBJS)
	
#Poisson_solver.o: Poisson_solver.F90
#	$(MPIFF) -mmic -c Poisson_solver.F90

%.o %.mod: %.F90
	$(MPIFF) $(FFLAGS) $<

clean:
	rm -f *.dat *.mod *.o $(EXE)

shm: MPI_data.mod MPI_applyBC_shared.o MPI_interface.mod MPI_SHMtest.o
	$(MPIFF) $(EXEFLAGS) shm MPI_applyBC_shared.o MPI_data.o MPI_SHMtest.o

cleanshm:
	rm -f *.mod *.o shm	
