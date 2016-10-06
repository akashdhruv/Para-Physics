#if 0
Defining Poisson Solver Parameters
#endif

#define POIS_SOLVER 2

#if POIS_SOLVER == 1
#define POISSON_SOLVER_JACOBI
#endif

#if POIS_SOLVER == 2
#define POISSON_SOLVER_GS
#endif

#if POIS_SOLVER == 3
#define POISSON_SOLVER_GSOR
#define omega 1.1
#endif

#if 0
Defining Temperature Solver
#endif

#define TEMP_SOLVER 2

#if TEMP_SOLVER == 1
#define TEMP_SOLVER_CENTRAL
#endif

#if TEMP_SOLVER == 2
#define TEMP_SOLVER_UPWIND
#endif


#if 0
Defining Simulation Parameters - Block Size, Domain Length, etc
#endif

#define MAX_BLOCKS 80


#define Nxb 40
#define Nyb 20


#define D_xmin -1.50
#define D_ymin -0.50


#define D_xmax 1.50
#define D_ymax 0.50

#define TIME_END 20.0

#define nblockx 2
#define nblocky 4

#define MaxIt 1500


#define PRES_VAR 1
#define TEMP_VAR 2
#define DFUN_VAR 3
#define PFUN_VAR 4
#define THCO_VAR 5
#define CPRS_VAR 6
#define NRMX_VAR 7
#define NRMY_VAR 8

#define CENT_VAR 8

#define VELC_VAR 1
#define FACE_VAR 1


#if 0
Defining Flow Type
#endif

#define FLOW 2

#if FLOW == 1
#define LID_DRIVEN_FLOW
#endif

#if FLOW == 2
#define CHANNEL_FLOW
#endif

#if FLOW == 3
#define MPH_FLOW
#endif

#if 0
Multiphase On/Off
#endif

#define MPH 0

#if MPH == 1
#define MULTIPHASE
#endif

#if MPH == 0
#define SINGLEPHASE
#endif

#if 0
Navier Stokes On/Off
#endif

#define NS 1

#if NS == 1
#define INS
#endif

#if 0
Energy On/Off
#endif

#define HEAT 1

#if HEAT == 1
#define ENERGY
#endif

#if 0
Only Poisson
#endif

#define PS 0

#if PS == 1
#define ONLY_POISSON
#endif

#if 0
OpenMP threads per MPI rank
#endif

#define NTHREADS 1
