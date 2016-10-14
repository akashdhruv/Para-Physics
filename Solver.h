#if 0
_____________Defining Poisson Solver Parameters_________________
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
_________________Defining Temperature Solver____________________
#endif

#define TEMP_SOLVER 2

#if TEMP_SOLVER == 1
#define TEMP_SOLVER_CENTRAL
#endif

#if TEMP_SOLVER == 2
#define TEMP_SOLVER_UPWIND
#endif


#if 0
_______________________Grid Parameters___________________________
#endif

#define MAX_BLOCKS 80


#define Nxb 50
#define Nyb 40


#define D_xmin -0.50
#define D_ymin -0.20


#define D_xmax  0.50
#define D_ymax  0.20

#define TIME_END 20.0

#define nblockx 5
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
#define VORT_VAR 9

#define CENT_VAR 9

#define VELC_VAR 1
#define IBMF_VAR 2
#define NMXF_VAR 3
#define NMYF_VAR 4

#define FACE_VAR 4


#if 0
_______________________Defining Flow Type_____________________
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
_______________________Multiphase On/Off_____________________
#endif

#define MPH 0

#if MPH == 1
#define MULTIPHASE
#endif

#if MPH == 0
#define SINGLEPHASE
#endif

#if 0
______________________Navier Stokes On/Off___________________
#endif

#define NS 1

#if NS == 1
#define INS
#endif

#if 0
_________________________Energy On/Off______________________
#endif

#define HEAT 1

#if HEAT == 1
#define ENERGY
#endif

#if 0
_____________________Only Poisson On/Off_____________________
#endif

#define PS 0

#if PS == 1
#define ONLY_POISSON
#endif

#if 0
___________________Immersed Boundary On/Off________________
#endif

#define IB 1

#if IB == 1
#define IBM
#endif

#define NBOD 1

#if 0
_________________OpenMP threads per MPI rank_______________
#endif

#define NTHREADS 1

