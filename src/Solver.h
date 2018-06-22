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
#define POISSON_SOLVER_GS_SKEW
#endif

#if POIS_SOLVER == 4
#define POISSON_SOLVER_GS_SOR
#define omega 1.1
#endif

#if 0
_________________Defining Temperature Solver____________________
#endif

#define TEMP_SOLVER 1

#if TEMP_SOLVER == 1
#define TEMP_SOLVER_CENTRAL
#endif

#if TEMP_SOLVER == 2
#define TEMP_SOLVER_UPWIND
#endif


#if 0
_______________________Grid Parameters___________________________
#endif

#define MAX_BLOCKS 64


#define Nxb 16
#define Nyb 16
#define Nzb 16

#define D_xmin  0.0
#define D_ymin  0.0


#define D_xmax  1.0
#define D_ymax  1.0

#define TIME_END 200.00
#define NEND_END 1

#define nblockx 8
#define nblocky 8
#define nblockz 8

#define MaxIt 5000

#define PRES_VAR 1
#define TEMP_VAR 2
#define DFUN_VAR 3
#define PRHS_VAR 4
#define VORO_VAR 5
#define NRMX_VAR 6
#define TOLD_VAR 7
#define NRMY_VAR 8
#define VISC_VAR 9
#define SIGP_VAR 10
#define SMHV_VAR 11
#define CURV_VAR 12
#define THCO_VAR 13
#define SMRH_VAR 14
#define CPRS_VAR 15
#define PFUN_VAR 16
#define MDOT_VAR 17
#define VORT_VAR 18
#define EXCT_VAR 19
#define EROR_VAR 20

#define CENT_VAR 20

#define VELC_VAR 1
#define IBMF_VAR 2
#define UOLD_VAR 3
#define GOLD_VAR 4
#define VELI_VAR 5
#define SIGM_VAR 6
#define USTR_VAR 7
#define RH1F_VAR 8
#define AL1F_VAR 9
#define RH2F_VAR 10
#define AL2F_VAR 11

#define FACE_VAR 11

#define CENTER 1
#define FACEX 2
#define FACEY 3

#if 0
_______________________Boundary Conditions_____________________
#endif

#define NOSLIP  1
#define SLIP    2
#define NEUMANN 3
#define INFLOW  4
#define OUTFLOW 5
#define MOVLID  6

#if 0
_______________________Number of immersed boundaries_____________________
#endif

#define NBOD 1

#if 0
_________________OpenMP threads per MPI rank_______________
#endif

#define NTHREADS 1

#if 0
_________________MPI with shared memory on/off_______________
#endif

#define SHARE_TYPE 2

#if SHARE_TYPE == 1
#define MPI_DIS
#endif

#if SHARE_TYPE == 2
#define MPI_SHM
#endif

#if SHARE_TYPE == 3
#define MPI_RMA
#define MPI_RMA_ACTIVE
#endif

#if 0
__________________Data Operation Parameters________________
#endif

#define SUM_DATA 1
#define MIN_DATA 2
#define MAX_DATA 3

#if 0
__________________Test problem forward/backward facing step________________
#endif

#define STEP 1

#if STEP == 1
#define BACKWARD_FACING_STEP
#endif

#if STEP == 2
#define FORWARD_FACING_STEP
#endif
