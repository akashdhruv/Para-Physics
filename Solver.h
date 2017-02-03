#if 0
_____________Defining Poisson Solver Parameters_________________
#endif

#define POIS_SOLVER 3

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


#define Nxb 20
#define Nyb 20


#define D_xmin -0.02
#define D_ymin -0.02


#define D_xmax  0.02
#define D_ymax  0.02

#define TIME_END 0.025

#define nblockx 2
#define nblocky 2

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
#define VISC_VAR 10
#define TOLD_VAR 11
#define CURV_VAR 12
#define SIGP_VAR 13
#define MDOT_VAR 14
#define SMRH_VAR 15
#define SMHV_VAR 16

#define CENT_VAR 16

#define VELC_VAR 1
#define IBMF_VAR 2
#define RH1F_VAR 3
#define RH2F_VAR 4
#define AL1F_VAR 5
#define AL2F_VAR 6
#define SIGM_VAR 7
#define VELI_VAR 8

#define FACE_VAR 8


#if 0
_______________________Defining Flow Type_____________________
#endif

#define FLOW 3

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

#define MPH 1

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

#define HEAT 0

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

#define IB 0

#if IB == 1
#define IBM
#endif

#define NBOD 1

#if 0
_________________OpenMP threads per MPI rank_______________
#endif

#define NTHREADS 2

#if 0
__________________Data Operation Parameters________________
#endif

#define SUM_DATA 1
#define MIN_DATA 2
#define MAX_DATA 3


