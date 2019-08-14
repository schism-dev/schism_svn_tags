#ifdef USE_SINGLE
# define MyREAL(xinp) REAL(xinp)
# define MySQRT(xinp) SQRT(xinp)
# define MyATAN2(xinp,yinp) ATAN2(xinp,yinp)
# define MyABS(xinp) ABS(xinp)
# define MyEXP(xinp) EXP(xinp)
# define MySINH(xinp) DSINH(xinp)
# define MyCOSH(xinp) COSH(xinp)
#else
# define MyREAL(xinp) DBLE(xinp)
# define MySQRT(xinp) DSQRT(xinp)
# define MyATAN2(xinp,yinp) DATAN2(xinp,yinp)
# define MyABS(xinp) DABS(xinp)
# define MyEXP(xinp) DEXP(xinp)
# define MySINH(xinp) SINH(xinp)
# define MyCOSH(xinp) DCOSH(xinp)
#endif
#if defined ST41 || defined ST42
# define ST_DEF
#endif
#ifdef PGMCL_COUPLING
# include "cppdefs.h"
#endif
#ifdef SELFE
# define MPI_PARALL_GRID
#endif
#ifdef WWM_MPI
# define MPI_PARALL_GRID
#endif
